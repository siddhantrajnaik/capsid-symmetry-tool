"""
gate_db.py -- local database of capsid gate residues.

One SQLite file (gate_residues.sqlite) holding, for every capsid structure analysed:
  structures       one row per structure (virus, family, class, chains, how the symmetry was obtained)
  gate_residues    one row per residue near a 5-, 3- or 2-fold axis (chain, number, amino acid, distance)
  dual_positions   residue numbers that sit near a 5-fold axis in one subunit and a 3-fold axis in another
  preferences      pooled amino-acid enrichment per capsid class (Cochran-Mantel-Haenszel odds ratio, FDR)
  meta             how and when the database was built

Build it from the output folders of gate_residue_analysis.py, or add single runs to it as you go:
    python gate_db.py build  --db gate_residues.sqlite --runs db_T3 db_T4 db_pT3 --catalog representatives.csv
    python gate_db.py add    --db gate_residues.sqlite --run results/run_20261028_130000
    python gate_db.py stats  --db gate_residues.sqlite
    python gate_db.py tidy-sources --db gate_residues.sqlite    (keep only folder names in structures.source)
"""
import argparse, os, sqlite3, glob, datetime, json
import pandas as pd
import numpy as np

SCHEMA = """
CREATE TABLE IF NOT EXISTS structures (
    pdb TEXT PRIMARY KEY, virus TEXT, family TEXT, capsid_class TEXT, resolution_A REAL,
    n_chains INTEGER, n_residues INTEGER, symmetry_from TEXT, axis_fit_rms_deg REAL,
    source TEXT, added TEXT);
CREATE TABLE IF NOT EXISTS gate_residues (
    pdb TEXT, virus TEXT, family TEXT, capsid_class TEXT, gate TEXT, chain TEXT,
    auth_seq INTEGER, aa TEXT, resname TEXT, distance_A REAL, radius_A REAL, dual_5_3 INTEGER,
    cutoff_A REAL, atoms TEXT);
CREATE TABLE IF NOT EXISTS dual_positions (
    pdb TEXT, virus TEXT, family TEXT, capsid_class TEXT, auth_seq INTEGER, aa TEXT, resname TEXT,
    chains_near_5fold TEXT, min_d5 REAL, chains_near_3fold TEXT, min_d3 REAL, cutoff_A REAL);
CREATE TABLE IF NOT EXISTS preferences (
    capsid_class TEXT, cutoff_A REAL, zone TEXT, level TEXT, key TEXT, n_structures INTEGER,
    n_with_key INTEGER, n_total INTEGER, freq_gate REAL, freq_background REAL,
    odds_ratio REAL, fdr REAL, families_enriched INTEGER, families_depleted INTEGER);
CREATE TABLE IF NOT EXISTS meta (k TEXT PRIMARY KEY, v TEXT);
CREATE INDEX IF NOT EXISTS ix_gr_pdb   ON gate_residues(pdb);
CREATE INDEX IF NOT EXISTS ix_gr_aa    ON gate_residues(aa);
CREATE INDEX IF NOT EXISTS ix_gr_gate  ON gate_residues(gate);
CREATE INDEX IF NOT EXISTS ix_dp_pdb   ON dual_positions(pdb);
CREATE INDEX IF NOT EXISTS ix_st_virus ON structures(virus);
"""


def connect(path):
    con = sqlite3.connect(path)
    con.executescript(SCHEMA)
    return con


def _catalog(path):
    """optional catalogue (representatives.csv) giving family / class / resolution per PDB id"""
    if not path or not os.path.exists(path):
        return {}
    t = pd.read_csv(path)
    out = {}
    for r in t.itertuples():
        out[str(r.pdb_id).lower()] = dict(virus=getattr(r, 'group', None), family=getattr(r, 'family', None),
                                          capsid_class=getattr(r, 'T_label', None),
                                          resolution_A=getattr(r, 'res', None))
    return out


def add_run(con, run_dir, catalog=None, cutoff=None, atoms=None, source=None):
    """import one output folder of gate_residue_analysis.py; returns the number of structures imported"""
    cat = catalog or {}
    used_f = os.path.join(run_dir, 'structures_used.csv')
    gates_f = os.path.join(run_dir, 'gate_residues_all.csv')
    dist_f = os.path.join(run_dir, 'residue_axis_distances.csv.gz')
    if not os.path.exists(used_f):
        raise FileNotFoundError(f'{run_dir}: no structures_used.csv (the run saved nothing)')
    used = pd.read_csv(used_f)
    log = os.path.join(run_dir, 'run.log')
    if cutoff is None or atoms is None:
        cutoff, atoms = 8.0, 'sc'
        if os.path.exists(log):
            for line in open(log, encoding='utf-8', errors='replace'):
                if line.startswith('arguments:'):
                    parts = line.split()
                    if '--cutoff' in parts:
                        cutoff = float(parts[parts.index('--cutoff') + 1])
                    if '--atoms' in parts:
                        atoms = parts[parts.index('--atoms') + 1]
                    break
    now = datetime.datetime.now().isoformat(timespec='seconds')
    cur = con.cursor()

    def vname(pdb, group):
        """a run started from a bare PDB id names the virus after the id; take the real name from the catalogue"""
        c = cat.get(str(pdb).lower(), {})
        return c.get('virus') or group if str(group).lower() == str(pdb).lower() else group
    for r in used.itertuples():
        c = cat.get(str(r.pdb).lower(), {})
        cur.execute('INSERT OR REPLACE INTO structures VALUES (?,?,?,?,?,?,?,?,?,?,?)',
                    (r.pdb, vname(r.pdb, r.group), c.get('family') or getattr(r, 'family', None),
                     c.get('capsid_class') or getattr(r, 'T_label', None), c.get('resolution_A'),
                     int(r.n_chains), int(r.n_residues), r.symmetry_from, float(r.axis_fit_rms_deg),
                     source or os.path.abspath(run_dir), now))
    pdbs = list(used.pdb.astype(str))
    meta = used.set_index('pdb')[['group']].to_dict('index')
    q = ','.join('?' * len(pdbs))
    cur.execute(f'DELETE FROM gate_residues WHERE pdb IN ({q})', pdbs)
    cur.execute(f'DELETE FROM dual_positions WHERE pdb IN ({q})', pdbs)

    if os.path.exists(gates_f) and os.path.getsize(gates_f) > 5:
        g = pd.read_csv(gates_f)
        rad = {}
        if os.path.exists(dist_f):                      # radial position of each residue, if available
            R = pd.read_csv(dist_f, usecols=['pdb', 'chain', 'auth_seq', 'r5_sc', 'r3_sc', 'r2_sc'])
            for n in (5, 3, 2):
                for row in R.itertuples():
                    rad[(row.pdb, row.chain, row.auth_seq, f'{n}-fold')] = getattr(row, f'r{n}_sc')
        rows = []
        for r in g.itertuples():
            cls = cat.get(str(r.pdb).lower(), {}).get('capsid_class')
            fam = cat.get(str(r.pdb).lower(), {}).get('family')
            rows.append((r.pdb, vname(r.pdb, r.group), fam, cls, r.gate, r.chain, int(r.auth_seq), r.aa, r.resname,
                         float(r.distance_A), rad.get((r.pdb, r.chain, r.auth_seq, r.gate)),
                         int(bool(r.dual_5_3)), cutoff, atoms))
        cur.executemany('INSERT INTO gate_residues VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)', rows)

    dual_f = os.path.join(run_dir, 'dual_gate_residues.csv')
    if os.path.exists(dual_f) and os.path.getsize(dual_f) > 5:
        d = pd.read_csv(dual_f)
        rows = []
        for r in d.itertuples():
            cls = cat.get(str(r.pdb).lower(), {}).get('capsid_class') or getattr(r, 'T_label', None)
            fam = cat.get(str(r.pdb).lower(), {}).get('family') or getattr(r, 'family', None)
            rows.append((r.pdb, vname(r.pdb, r.group), fam, cls, int(r.auth_seq), r.aa, r.resname,
                         r.chains_near_5fold, float(r.min_d5), r.chains_near_3fold, float(r.min_d3), cutoff))
        cur.executemany('INSERT INTO dual_positions VALUES (?,?,?,?,?,?,?,?,?,?,?,?)', rows)
    con.commit()
    return len(used)


def add_preferences(con, run_dir, capsid_class):
    """import the pooled amino-acid preferences of one capsid class"""
    cur = con.cursor()
    cur.execute('DELETE FROM preferences WHERE capsid_class = ?', (capsid_class,))
    n = 0
    for f in sorted(glob.glob(os.path.join(run_dir, 'pooled_enrichment_cutoff*.csv'))):
        cutoff = float(os.path.basename(f).split('cutoff')[1].split('.csv')[0])
        p = pd.read_csv(f)
        cf = os.path.join(run_dir, f'consistency_across_viruses_cutoff{cutoff:g}.csv')
        if os.path.exists(cf):
            c = pd.read_csv(cf)[['zone', 'level', 'key', 'families_enriched', 'families_depleted']]
            p = p.merge(c, on=['zone', 'level', 'key'], how='left')
        for r in p.itertuples():
            cur.execute('INSERT INTO preferences VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                        (capsid_class, cutoff, r.zone, r.level, r.key, int(r.n_structures),
                         int(r.fg_with_key), int(r.fg_total), float(r.freq_fg), float(r.freq_bg),
                         None if pd.isna(r.cmh_odds_ratio) else float(r.cmh_odds_ratio),
                         None if pd.isna(r.cmh_fdr) else float(r.cmh_fdr),
                         None if pd.isna(getattr(r, 'families_enriched', np.nan)) else int(r.families_enriched),
                         None if pd.isna(getattr(r, 'families_depleted', np.nan)) else int(r.families_depleted)))
            n += 1
    con.commit()
    return n


def tidy_sources(con):
    """keep only the run-folder name in structures.source (drops local paths before sharing the database)"""
    cur = con.cursor()
    rows = cur.execute('SELECT pdb, source FROM structures').fetchall()
    n = 0
    for pdb, src in rows:
        short = src.replace('\\', '/').rstrip('/').split('/')[-1] if src else src
        if short != src:
            cur.execute('UPDATE structures SET source = ? WHERE pdb = ?', (short, pdb)); n += 1
    con.commit()
    return n


def set_meta(con, **kw):
    cur = con.cursor()
    for k, v in kw.items():
        cur.execute('INSERT OR REPLACE INTO meta VALUES (?,?)', (k, json.dumps(v) if not isinstance(v, str) else v))
    con.commit()


def summary(con):
    g = lambda q: con.execute(q).fetchone()[0]
    return dict(structures=g('SELECT COUNT(*) FROM structures'),
                viruses=g('SELECT COUNT(DISTINCT virus) FROM structures'),
                families=g('SELECT COUNT(DISTINCT family) FROM structures'),
                gate_residues=g('SELECT COUNT(*) FROM gate_residues'),
                dual_positions=g('SELECT COUNT(*) FROM dual_positions'),
                preference_rows=g('SELECT COUNT(*) FROM preferences'),
                built=(con.execute("SELECT v FROM meta WHERE k='built'").fetchone() or [''])[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['build', 'add', 'stats', 'tidy-sources'])
    ap.add_argument('--db', default='gate_residues.sqlite')
    ap.add_argument('--runs', nargs='*', default=[])
    ap.add_argument('--run')
    ap.add_argument('--catalog', help='representatives.csv (adds family / class / resolution)')
    ap.add_argument('--class-of', nargs='*', default=[], help='run_dir=CLASS pairs for the preferences table')
    a = ap.parse_args()
    con = connect(a.db)
    cat = _catalog(a.catalog)
    if a.cmd == 'tidy-sources':
        print(f'{tidy_sources(con)} source paths shortened to their folder name'); return
    if a.cmd == 'stats':
        print(json.dumps(summary(con), indent=1)); return
    runs = a.runs if a.cmd == 'build' else [a.run]
    classes = dict(x.split('=', 1) for x in a.class_of)
    total = 0
    for rd in runs:
        total += add_run(con, rd, cat)
        if rd in classes:
            add_preferences(con, rd, classes[rd])
        print(f'imported {rd}')
    set_meta(con, built=datetime.datetime.now().isoformat(timespec='seconds'),
             note='Gate residues = atoms within the cutoff of a true icosahedral 5-, 3- or 2-fold axis. '
                  'Pooled preferences need >= 5 viruses per class.')
    print(f'{total} structures in {a.db}'); print(json.dumps(summary(con), indent=1))


if __name__ == '__main__':
    main()
