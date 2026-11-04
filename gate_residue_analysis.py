#!/usr/bin/env python
"""
gate_residue_analysis.py  --  which amino acids guard the 5-fold and 3-fold gates of icosahedral capsids?

For every capsid structure the program
  1. reads the PDB/mmCIF entry (RCSB, downloaded and cached) and the icosahedral operators of its biological assembly,
  2. rebuilds the icosahedral rotation group, finds the capsid centre and all 5-fold (6), 3-fold (10) and 2-fold (15) axes
     directly from those operators (so it works whatever orientation the entry is deposited in),
  3. measures, for every residue of the deposited chains, the shortest distance from its side-chain atoms
     (or all heavy atoms) to the nearest 5-, 3- and 2-fold axis,
  4. defines "gate-guarding" residues (distance <= cutoff to the axis), "dual-gate" residues (same sequence position
     near a 5-fold axis in one quasi-equivalent subunit and near a 3-fold axis in another, or the same one) and
     "gatekeepers" (the k residues closest to the axis),
  5. tests each amino acid (and physico-chemical class) for enrichment against a background of capsid residues that are
     not near any axis (optionally matched for radial position in the shell), per structure (Fisher exact + shrunk log2 enrichment), pooled over
     structures (Cochran-Mantel-Haenszel, so differences in composition between viruses do not matter) and for
     consistency across viruses (sign test / Wilcoxon on per-virus log2 enrichments, and across virus families).

Input   : --pdb 4ftb 2z2q ...   and/or   --csv list.csv  (columns: pdb_id, [group], [family], [T_label])
          and/or --model my_model.pdb ...  (local PDB/mmCIF files, e.g. your own reconstructions; need BIOMT/assembly operators or the full capsid).
          Rows sharing the same `group` are alternative structures of the same virus: the first one that can be
          analysed is used, the rest are skipped (keeps one structure per virus, avoids pseudo-replication).
Output  : <out>/residue_axis_distances.csv, dual_gate_residues.csv, gatekeepers.csv, enrichment_per_structure.csv,
          pooled_enrichment_cutoff<c>.csv, consistency_across_viruses_cutoff<c>.csv, structures_used.csv,
          figures (heatmap, forest plot) and report.txt

Examples:
  python gate_residue_analysis.py --pdb 4ftb --residue 182 --out fhv
  python gate_residue_analysis.py --csv representatives.csv --cutoff 8 --sens 6 10 --out batch

Requires: numpy pandas scipy statsmodels matplotlib gemmi requests
Caveats : the 3-fold axes used are the TRUE icosahedral 3-fold axes. Quasi-3-fold axes of T=3 capsids (A-B-C centres)
          are not analysed. Entries that deposit the whole capsid are handled by finding the symmetry from coordinates (needs the full set of 60 copies).
"""
import argparse, gzip, os, sys, io, json, traceback, datetime, warnings, faulthandler, atexit
import numpy as np, pandas as pd
import gemmi, requests
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

AA3TO1 = dict(ALA='A', CYS='C', ASP='D', GLU='E', PHE='F', GLY='G', HIS='H', ILE='I', LYS='K', LEU='L',
              MET='M', ASN='N', PRO='P', GLN='Q', ARG='R', SER='S', THR='T', VAL='V', TRP='W', TYR='Y', MSE='M')
AA_ORDER = list('AVLIMFWYCGPSTNQHKRDE')                      # grouped by chemistry for plots
CLASS = {**{a: 'aliphatic (AVLI)' for a in 'AVLI'}, 'M': 'Met', **{a: 'aromatic (FWY)' for a in 'FWY'},
         **{a: 'polar (STNQ)' for a in 'STNQ'}, **{a: 'positive (KRH)' for a in 'KRH'},
         **{a: 'negative (DE)' for a in 'DE'}, 'G': 'Gly', 'P': 'Pro', 'C': 'Cys'}
BACKBONE = {'N', 'C', 'O', 'OXT', 'CA'}
AXIS_N = {5: 6, 3: 10, 2: 15}
SHRINK = 5.0        # pseudo-observations (drawn at the background frequency) used for per-structure log2 enrichment


# ----------------------------------------------------------------------------- structure + symmetry
def fetch_cif(pid, cache, max_mb):
    os.makedirs(cache, exist_ok=True)
    p = os.path.join(cache, pid.lower() + '.cif.gz')
    if not os.path.exists(p):
        r = requests.get(f'https://files.rcsb.org/download/{pid.upper()}.cif.gz', timeout=180)
        r.raise_for_status()
        if len(r.content) > max_mb * 1e6:
            raise RuntimeError(f'file {len(r.content)/1e6:.0f} MB > --max-mb {max_mb}')
        open(p, 'wb').write(r.content)
    return p


def assembly_operators(st):
    """operators (R, t) of the first assembly with >= 30 distinct operators (the capsid)"""
    for a in st.assemblies:
        ops = []
        for g in a.generators:
            for o in g.operators:
                ops.append((np.array(o.transform.mat.tolist()), np.array(o.transform.vec.tolist())))
        uniq = []
        for R, t in ops:
            if not any(np.allclose(R, u[0], atol=1e-3) and np.allclose(t, u[1], atol=0.05) for u in uniq):
                uniq.append((R, t))
        if len(uniq) >= 30:
            return uniq
    return None


def close_group(ops, cap=70):
    G = [(np.eye(3), np.zeros(3))]
    def has(M):
        return any(np.allclose(M[0], g[0], atol=1e-2) and np.allclose(M[1], g[1], atol=0.1) for g in G)
    changed = True
    while changed and len(G) <= cap:
        changed = False
        for a in list(G):
            for g in ops:
                m = (g[0] @ a[0], g[0] @ a[1] + g[1])
                if not has(m):
                    G.append(m); changed = True
    if len(G) > 60:
        raise RuntimeError(f'operator closure gives {len(G)} elements (noisy operators)')
    return G


def _rot(axis, ang):
    k = np.asarray(axis, float); k /= np.linalg.norm(k)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


def _axes_exact(Rs, ang_tol=2e-2, dot=0.995):
    axes = {5: [], 3: [], 2: []}
    I3 = np.eye(3)
    for R in Rs:
        if np.allclose(R, I3, atol=ang_tol):
            continue
        n = next((k for k in range(2, 7) if np.allclose(np.linalg.matrix_power(R, k), I3, atol=ang_tol)), None)
        if n not in axes:
            continue
        u = np.linalg.svd(R - I3)[2][-1]
        u = u / np.linalg.norm(u)
        if u[np.argmax(np.abs(u) > 1e-6)] < 0:
            u = -u
        if not any(abs(u @ v) > dot for v in axes[n]):
            axes[n].append(u)
    return {n: np.array(v) for n, v in axes.items()}


def _ideal_axes():
    phi = (1 + 5 ** 0.5) / 2
    gens = [_rot((0, 1, phi), 2 * np.pi / 5), _rot((1, 1, 1), 2 * np.pi / 3)]
    G = [np.eye(3)]
    changed = True
    while changed:
        changed = False
        for a in list(G):
            for g in gens:
                m = g @ a
                if not any(np.allclose(m, x, atol=1e-6) for x in G):
                    G.append(m); changed = True
    assert len(G) == 60
    return _axes_exact(G, ang_tol=1e-5)


def icosahedral_axes(G, noisy=False):
    """returns (centre, {5: 6 axes, 3: 10 axes, 2: 15 axes}, axis_fit_rms_deg).
    noisy=False: axes read directly from the (exact) operators.
    noisy=True : operators come from superposing deposited chains and are only approximately symmetric; the six
                 measured 5-fold axes are used to orient an ideal icosahedron, which supplies all axes."""
    if len(G) != 60:
        raise RuntimeError(f'operator group has {len(G)} elements, expected 60 (rotational icosahedral group)')
    if any(abs(np.linalg.det(R) - 1) > (0.1 if noisy else 1e-2) for R, _ in G):
        raise RuntimeError('improper rotation among operators')
    A = np.vstack([np.eye(3) - R for R, _ in G]); b = np.concatenate([t for _, t in G])
    c = np.linalg.lstsq(A, b, rcond=None)[0]
    if not noisy:
        if max(np.linalg.norm(R @ c + t - c) for R, t in G) > 1.0:
            raise RuntimeError('operators do not share a common fixed point')
        axes = _axes_exact([R for R, _ in G])
        for n, k in AXIS_N.items():
            if len(axes[n]) != k:
                raise RuntimeError(f'found {len(axes[n])} {n}-fold axes, expected {k}')
        return c, axes, 0.0
    M5 = []
    for R, _ in G:
        th = np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))
        if min(abs(th - 72), abs(th - 144)) < 12:
            u = np.linalg.svd(R - np.eye(3))[2][-1]; u /= np.linalg.norm(u)
            if not any(abs(u @ v) > 0.94 for v in M5):
                M5.append(u)
    if len(M5) != 6:
        raise RuntimeError(f'found {len(M5)} five-fold axes among the operators, expected 6')
    M5 = np.array(M5); ideal = _ideal_axes(); D5 = ideal[5]
    def frame(u, v):
        e2 = v - (v @ u) * u; e2 /= np.linalg.norm(e2)
        return np.column_stack([u, e2, np.cross(u, e2)])
    a2 = M5[1] * np.sign(M5[0] @ M5[1]); d2 = D5[1] * np.sign(D5[0] @ D5[1])
    Q = frame(M5[0], a2) @ frame(D5[0], d2).T
    for _ in range(6):
        T = []
        for k in D5 @ Q.T:
            m = M5[np.argmax(np.abs(M5 @ k))]; T.append(m * np.sign(m @ k))
        U, _, Vt = np.linalg.svd(D5.T @ np.array(T))
        Q = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
    ang = [np.degrees(np.arccos(np.clip(np.abs(M5 @ k).max(), -1, 1))) for k in D5 @ Q.T]
    rms = float(np.sqrt(np.mean(np.square(ang))))
    return c, {n: v @ Q.T for n, v in ideal.items()}, rms


def full_capsid_model(st):
    """model with the complete biological assembly (also works when the entry already holds the whole capsid)"""
    if not st.assemblies:
        return st[0]
    return gemmi.make_assembly(st.assemblies[0], st[0], gemmi.HowToNameCopiedChain.AddNumber)


def ops_from_coordinates(model, max_err=25.0, min_gap=8.0):
    """Entries that deposit the complete capsid (or a part with < 30 operators): find the 60 icosahedral operators by
    superposing one subunit on every other copy of the same protein. A true symmetry operator maps the set of all
    subunit centroids onto itself, a quasi-equivalent (A->B) superposition does not; the 60 best operators (smallest
    worst-case centroid mismatch) must be separated from the 61st by >= min_gap A.
    Returns (operators, names of one chain per symmetry orbit = one icosahedral asymmetric unit, mismatch of #60)."""
    from scipy.spatial import cKDTree
    chains = []
    for ch in model:
        d = {}
        for res in ch:
            if res.name in AA3TO1:
                ca = res.find_atom('CA', '*')
                if ca is not None:
                    d[(res.seqid.num, res.seqid.icode)] = np.array(ca.pos.tolist())
        if len(d) >= 50:
            chains.append((ch.name, ch[0].entity_id, d))
    if len(chains) < 60:
        raise RuntimeError(f'only {len(chains)} protein chains in the assembly, not a complete capsid')
    ent = pd.Series([e for _, e, _ in chains]).value_counts().index[0]
    chains = [c for c in chains if c[1] == ent]
    C = np.array([np.mean(list(d.values()), axis=0) for _, _, d in chains])
    ctree = cKDTree(C)
    ref = chains[0][2]
    cand = []
    for name, _, d in chains:
        keys = [k for k in ref if k in d]
        if len(keys) < 30:
            continue
        X = np.array([ref[k] for k in keys]); Y = np.array([d[k] for k in keys])
        mx, my = X.mean(0), Y.mean(0)
        U, _, Vt = np.linalg.svd((X - mx).T @ (Y - my))
        R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
        t = my - R @ mx
        cand.append((ctree.query(C @ R.T + t)[0].max(), R, t))
    cand.sort(key=lambda x: x[0])
    if len(cand) < 61:
        raise RuntimeError(f'only {len(cand)} superposable copies of the reference subunit')
    e60, e61 = cand[59][0], cand[60][0]
    if e60 > max_err or e61 - e60 < min_gap:
        raise RuntimeError(f'no clear set of 60 symmetry operators (mismatch #60 = {e60:.1f} A, #61 = {e61:.1f} A)')
    ops = [(R, t) for _, R, t in cand[:60]]
    claimed, reps = set(), []
    for i in range(len(C)):
        if i in claimed:
            continue
        reps.append(chains[i][0])
        for R, t in ops:
            claimed.add(int(ctree.query(R @ C[i] + t)[1]))
    return ops, set(reps), e60


# ----------------------------------------------------------------------------- residue geometry
def residue_table(model, c, axes, pid, keep=None):
    X, ridx, issc, meta = [], [], [], []
    for chain in model:
        if keep is not None and chain.name not in keep:
            continue
        for res in chain:
            a1 = AA3TO1.get(res.name)
            if a1 is None:
                continue
            atoms = [a for a in res if a.element.name != 'H' and a.altloc in ('\0', 'A')]
            if not atoms:
                continue
            names = [a.name.strip() for a in atoms]
            sc = [(n not in BACKBONE) for n in names]
            if not any(sc):                                     # Gly or truncated: use CA as the "side chain"
                sc = [(n == 'CA') for n in names]
                if not any(sc):
                    sc = [True] * len(names)
            i = len(meta)
            ca = next((a.pos.tolist() for a in atoms if a.name.strip() == 'CA'), atoms[0].pos.tolist())
            ls = res.label_seq if res.label_seq is not None else res.seqid.num
            meta.append(dict(pdb=pid, chain=chain.name, auth_seq=res.seqid.num, icode=res.seqid.icode.strip(),
                             label_seq=int(ls), entity=res.entity_id, resname=res.name, aa=a1,
                             r_ca=float(np.linalg.norm(np.array(ca) - c))))
            for a, s in zip(atoms, sc):
                X.append(a.pos.tolist()); ridx.append(i); issc.append(s)
    X = np.array(X) - c; ridx = np.array(ridx); issc = np.array(issc)
    df = pd.DataFrame(meta)
    rr = np.linalg.norm(X, axis=1)
    for n, U in axes.items():
        s = X @ U.T
        d = np.sqrt(np.maximum(rr[:, None] ** 2 - s ** 2, 0)).min(axis=1)      # distance to nearest n-fold axis line
        for tag, mask in (('sc', issc), ('all', np.ones_like(issc))):
            dd = np.full(len(df), np.inf); rad = np.full(len(df), np.nan)
            sel = np.where(mask)[0]
            order = np.argsort(d[sel])                                          # smallest distance first
            for k in order[::-1]:                                               # last write = smallest
                dd[ridx[sel[k]]] = d[sel[k]]; rad[ridx[sel[k]]] = rr[sel[k]]
            df[f'd{n}_{tag}'] = dd
            if tag == 'sc':
                df[f'r{n}_sc'] = rad
    df['pos'] = df.entity.astype(str) + ':' + df.label_seq.astype(str)
    return df


# ----------------------------------------------------------------------------- foreground / background
def foreground_sets(df, cutoff, mode, k):
    d = lambda n: df[f'd{n}_{mode}']
    near5, near3, near2 = d(5) <= cutoff, d(3) <= cutoff, d(2) <= cutoff
    sets = {'5-fold': df[near5], '3-fold': df[near3], '2-fold': df[near2]}
    p5, p3 = set(df.pos[near5]), set(df.pos[near3])
    dual = sorted(p5 & p3)
    sets['dual 5/3-fold position'] = df[df.pos.isin(dual)].drop_duplicates('pos')
    for n in (5, 3):
        best = df.assign(dd=d(n)).sort_values('dd').drop_duplicates('pos').head(k)
        sets[f'{n}-fold gatekeepers (closest {k})'] = best
    bg = df[~(near5 | near3 | near2)]
    return sets, bg


def background_for(fg, bg, kind):
    if kind == 'all' or len(fg) < 3:
        return bg
    lo, hi = fg.r_ca.quantile(0.05), fg.r_ca.quantile(0.95)
    return bg[(bg.r_ca >= lo) & (bg.r_ca <= hi)]


def enrichment(fg, bg, meta, zone):
    rows = []
    for keyname, col in (('aa', 'aa'), ('class', 'cls')):
        f, b = fg[col].value_counts(), bg[col].value_counts()
        keys = AA_ORDER if keyname == 'aa' else list(dict.fromkeys(CLASS[a] for a in AA_ORDER))
        Nf, Nb = len(fg), len(bg)
        for kk in keys:
            a, c_ = int(f.get(kk, 0)), int(b.get(kk, 0))
            if Nf == 0 or Nb == 0:
                continue
            p = stats.fisher_exact([[a, Nf - a], [c_, Nb - c_]])[1]
            pb = (c_ + 0.5) / (Nb + 1)                      # background frequency
            l2 = np.log2(((a + SHRINK * pb) / (Nf + SHRINK)) / pb)   # shrunk toward 'no enrichment' when few residues
            rows.append(dict(**meta, zone=zone, level=keyname, key=kk, n_fg=a, N_fg=Nf, n_bg=c_, N_bg=Nb,
                             freq_fg=a / Nf, freq_bg=c_ / Nb, log2_enrichment=l2, fisher_p=p))
    return rows


# ----------------------------------------------------------------------------- pooled statistics
def _cmh(g):
    """common odds ratio over structures (Cochran-Mantel-Haenszel); Fisher exact when only one informative structure"""
    tabs = []
    for _, r in g.iterrows():
        t = np.array([[r.n_fg, r.N_fg - r.n_fg], [r.n_bg, r.N_bg - r.n_bg]], float)
        if t[0].sum() > 0 and t[1].sum() > 0 and t[:, 0].sum() > 0 and t[:, 1].sum() > 0:
            tabs.append(t)
    if len(tabs) >= 2:
        try:
            st_ = StratifiedTable(tabs)
            lo, hi = st_.oddsratio_pooled_confint()
            return len(tabs), st_.oddsratio_pooled, lo, hi, st_.test_null_odds(correction=True).pvalue
        except Exception:
            return len(tabs), np.nan, np.nan, np.nan, np.nan
    if len(tabs) == 1:
        orr, p = stats.fisher_exact(tabs[0]); return 1, orr, np.nan, np.nan, p
    return 0, np.nan, np.nan, np.nan, np.nan


def pooled(enr, loo=True):
    out = []
    for (zone, level, key), g in enr.groupby(['zone', 'level', 'key']):
        n, orr, lo, hi, p = _cmh(g)
        rec = dict(zone=zone, level=level, key=key, n_structures=len(g), n_strata_informative=n,
                   fg_total=int(g.N_fg.sum()), fg_with_key=int(g.n_fg.sum()),
                   freq_fg=g.n_fg.sum() / max(g.N_fg.sum(), 1), freq_bg=g.n_bg.sum() / max(g.N_bg.sum(), 1),
                   cmh_odds_ratio=orr, or_ci_low=lo, or_ci_high=hi, cmh_p=p)
        fams = g.family.unique()
        if loo and len(fams) >= 3 and n >= 2:                       # does one virus family drive the result?
            loo = [_cmh(g[g.family != f]) for f in fams]
            ors = [x[1] for x in loo if np.isfinite(x[1])]
            if ors:
                rec['or_leave_one_family_out_min'] = min(ors); rec['or_leave_one_family_out_max'] = max(ors)
                rec['loo_families_p_lt_0.05'] = sum(1 for x in loo if x[4] < 0.05)
                rec['n_families'] = len(fams)
        out.append(rec)
    o = pd.DataFrame(out)
    o['cmh_fdr'] = np.nan
    for (zone, level), idx in o.groupby(['zone', 'level']).groups.items():
        ok = o.loc[idx, 'cmh_p'].notna()
        if ok.sum():
            o.loc[ok[ok].index, 'cmh_fdr'] = multipletests(o.loc[ok[ok].index, 'cmh_p'], method='fdr_bh')[1]
    o['log2_or'] = np.log2(o.cmh_odds_ratio.replace(0, np.nan))
    return o


def consistency(enr):
    """Is the direction of the effect shared between viruses and between families?
    per virus : shrunk log2 enrichment (n_enriched / n_depleted, Wilcoxon on viruses)
    per family: counts of all its viruses summed first, then the sign of the family's log2 enrichment
                (exact sign test over families; this is the conservative one because related viruses are not independent)"""
    out = []
    for (zone, level, key), g in enr.groupby(['zone', 'level', 'key']):
        v = g.log2_enrichment.values
        rec = dict(zone=zone, level=level, key=key, n_structures=len(v), median_log2_enrichment=np.median(v),
                   n_enriched=int((v > 0).sum()), n_depleted=int((v < 0).sum()))
        if len(v) >= 6:
            try:
                rec['wilcoxon_p'] = stats.wilcoxon(v).pvalue
            except Exception:
                pass
        f = g.groupby('family')[['n_fg', 'N_fg', 'n_bg', 'N_bg']].sum()
        pb = (f.n_bg + 0.5) / (f.N_bg + 1)
        l2 = np.log2(((f.n_fg + SHRINK * pb) / (f.N_fg + SHRINK)) / pb)
        e_, d_ = int((l2 > 0).sum()), int((l2 < 0).sum())
        rec.update(n_families=len(f), families_enriched=e_, families_depleted=d_, family_median_log2=float(l2.median()))
        if e_ + d_ >= 6:
            rec['family_sign_p'] = stats.binomtest(e_, e_ + d_, 0.5).pvalue
        out.append(rec)
    o = pd.DataFrame(out)
    for col, new in (('wilcoxon_p', 'wilcoxon_fdr'), ('family_sign_p', 'family_sign_fdr')):
        o[new] = np.nan
        if col in o:
            for (zone, level), idx in o.groupby(['zone', 'level']).groups.items():
                ok = o.loc[idx, col].notna()
                if ok.sum():
                    o.loc[ok[ok].index, new] = multipletests(o.loc[ok[ok].index, col], method='fdr_bh')[1]
    return o


# ----------------------------------------------------------------------------- figures
def figures(enr, pol, zones, out, cutoff, min_fg=8):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
    e = enr[(enr.level == 'aa')]
    fam = e.groupby(['family', 'zone', 'key'])[['n_fg', 'N_fg', 'n_bg', 'N_bg']].sum().reset_index()
    nv = enr.drop_duplicates('pdb').groupby('family').size()
    pb = (fam.n_bg + 0.5) / (fam.N_bg + 1)
    fam['l2'] = np.log2(((fam.n_fg + SHRINK * pb) / (fam.N_fg + SHRINK)) / pb)
    order = sorted(nv.index, key=lambda f: (-nv[f], f)); names = [f'{f} (n={nv[f]})' for f in order]
    pbf = (fam.n_bg + 0.5) / (fam.N_bg + 1)
    if not ((fam.N_fg >= min_fg) & (fam.N_fg * pbf >= 3)).any():
        print('heatmap skipped: no family has enough gate residues to judge (use more viruses)')
        heat = False
    else:
        heat = True
    fig, axs = plt.subplots(1, len(zones), figsize=(3.6 * len(zones) + 1.5, max(3, 0.25 * len(order) + 1.4)), sharey=True)
    axs = np.atleast_1d(axs)
    for ax, z in zip(axs, zones):
        fz = fam[fam.zone == z]
        M = fz.pivot(index='family', columns='key', values='l2').reindex(index=order, columns=AA_ORDER)
        Nf = fz.pivot(index='family', columns='key', values='N_fg').reindex(index=order, columns=AA_ORDER)
        Pb = ((fz.pivot(index='family', columns='key', values='n_bg') + 0.5) /
              (fz.pivot(index='family', columns='key', values='N_bg') + 1)).reindex(index=order, columns=AA_ORDER)
        M = M.where((Nf >= min_fg) & (Nf * Pb >= 3))   # blank when too few residues or < 3 expected of this amino acid
        im = ax.imshow(M.values, cmap='RdBu_r', vmin=-3, vmax=3, aspect='auto')
        ax.set_xticks(range(20)); ax.set_xticklabels(AA_ORDER); ax.set_title(z, loc='left')
        ax.set_yticks(range(len(order))); ax.set_yticklabels(names, fontsize=7)
        ax.set_xlabel('amino acid')
    cb = fig.colorbar(im, ax=axs, shrink=0.6, pad=0.01); cb.set_label('log2 enrichment vs background, residues of all viruses of a family pooled (0 = none; blank = too few residues to judge)')
    if heat:
        fig.savefig(os.path.join(out, f'heatmap_enrichment_cutoff{cutoff}.png'), dpi=200, bbox_inches='tight')
    plt.close(fig)
    p = pol[(pol.level == 'aa') & pol.cmh_odds_ratio.notna()]
    fig, axs = plt.subplots(1, len(zones), figsize=(2.9 * len(zones) + 0.5, 4.6), sharey=True)
    axs = np.atleast_1d(axs)
    for ax, z in zip(axs, zones):
        q = p[p.zone == z].set_index('key').reindex(AA_ORDER)
        y = np.arange(20)[::-1]
        sig = (q.cmh_fdr < 0.05).values
        ax.hlines(y, np.log2(q.or_ci_low), np.log2(q.or_ci_high), color='0.6', lw=1)
        ax.scatter(q.log2_or[~sig], y[~sig], facecolor='white', edgecolor='k', zorder=3, s=22)
        ax.scatter(q.log2_or[sig], y[sig], color='#b2182b', zorder=3, s=26)
        ax.axvline(0, color='k', lw=0.6); ax.set_title(z, loc='left'); ax.set_xlabel('log2 pooled odds ratio (95% CI)')
        ax.set_yticks(y); ax.set_yticklabels(AA_ORDER)
    axs[0].set_ylabel('amino acid  (filled: FDR < 0.05)')
    fig.savefig(os.path.join(out, f'forest_pooled_cutoff{cutoff}.png'), dpi=200, bbox_inches='tight'); plt.close(fig)


# ----------------------------------------------------------------------------- main
class Tee:
    """copy everything printed to the screen into <out>/run.log (so a failed run always leaves a trace)"""
    def __init__(self, path):
        self.f = open(path, 'a', encoding='utf-8', buffering=1); self.o = sys.__stdout__
    def write(self, s):
        try:
            if self.o is not None:
                self.o.write(s)
        except Exception:
            pass
        self.f.write(s)
    def flush(self):
        try:
            if self.o is not None:
                self.o.flush()
        except Exception:
            pass
        self.f.flush()


def write_gate_outputs(R, used, skipped, out, cutoff, mode):
    """results that are useful for ONE structure: which residues guard each gate, which positions guard both.
    Called after every structure, so a stopped or crashed run still leaves everything finished so far."""
    R.round(3).to_csv(os.path.join(out, 'residue_axis_distances.csv.gz'), index=False)
    if used:
        pd.DataFrame(used).to_csv(os.path.join(out, 'structures_used.csv'), index=False)
    pd.DataFrame(skipped, columns=['pdb', 'group', 'reason']).to_csv(os.path.join(out, 'structures_skipped.csv'), index=False)
    rows, lines = [], [f'Gate residues: side-chain atoms within {cutoff:g} A of a true icosahedral axis'
                       if mode == 'sc' else f'Gate residues: any heavy atom within {cutoff:g} A of a true icosahedral axis',
                       'Format: residue[chain] distance-to-axis in A, closest first', '']
    for pid, df in R.groupby('pdb', sort=False):
        D = {n: df[f'd{n}_{mode}'] for n in (5, 3, 2)}
        dual = set(df.pos[D[5] <= cutoff]) & set(df.pos[D[3] <= cutoff])
        lines.append(f'=== {pid.upper()}  ({df.group.iloc[0]}) ===')
        for n in (5, 3, 2):
            sel = df[D[n] <= cutoff].assign(dist=D[n][D[n] <= cutoff]).sort_values('dist')
            for r in sel.itertuples():
                rows.append(dict(pdb=pid, group=df.group.iloc[0], gate=f'{n}-fold', chain=r.chain, auth_seq=r.auth_seq,
                                 aa=r.aa, resname=r.resname, distance_A=round(r.dist, 2), dual_5_3=(r.pos in dual)))
            lines.append(f'{n}-fold gate ({len(sel)} residues): ' +
                         (', '.join(f'{r.aa}{r.auth_seq}[{r.chain}] {r.dist:.1f}' for r in sel.itertuples()) or 'none within cutoff'))
        det = []
        for pos in sorted(dual, key=lambda p: int(df.auth_seq[df.pos == p].iloc[0])):
            g = df[df.pos == pos]; n5 = g[D[5][g.index] <= cutoff]; n3 = g[D[3][g.index] <= cutoff]
            det.append(f'{g.aa.iloc[0]}{g.auth_seq.iloc[0]} (5-fold: {",".join(n5.chain)} {D[5][n5.index].min():.1f} A; '
                       f'3-fold: {",".join(n3.chain)} {D[3][n3.index].min():.1f} A)')
        lines.append('Same residue number near BOTH a 5-fold and a 3-fold axis: ' + ('; '.join(det) or 'none'))
        lines.append('')
    pd.DataFrame(rows).to_csv(os.path.join(out, 'gate_residues_all.csv'), index=False)
    open(os.path.join(out, 'gate_summary.txt'), 'w', encoding='utf-8').write('\n'.join(lines))
    dual_rows = []                      # also written here, so runs with fewer than 5 structures get the file too
    for pid, df in R.groupby('pdb', sort=False):
        D5, D3 = df[f'd5_{mode}'], df[f'd3_{mode}']
        for pos in sorted(set(df.pos[D5 <= cutoff]) & set(df.pos[D3 <= cutoff])):
            g = df[df.pos == pos]; n5 = g[D5[g.index] <= cutoff]; n3 = g[D3[g.index] <= cutoff]
            dual_rows.append(dict(pdb=pid, group=df.group.iloc[0], family=df.family.iloc[0], T_label=df.T_label.iloc[0],
                                  position=pos, auth_seq=int(g.auth_seq.iloc[0]), aa=g.aa.iloc[0], resname=g.resname.iloc[0],
                                  chains_near_5fold=','.join(n5.chain), min_d5=round(float(D5[n5.index].min()), 2),
                                  chains_near_3fold=','.join(n3.chain), min_d3=round(float(D3[n3.index].min()), 2)))
    pd.DataFrame(dual_rows, columns=['pdb', 'group', 'family', 'T_label', 'position', 'auth_seq', 'aa', 'resname',
                                     'chains_near_5fold', 'min_d5', 'chains_near_3fold', 'min_d3']
                 ).to_csv(os.path.join(out, 'dual_gate_residues.csv'), index=False)


def local_name(path):
    n = os.path.basename(path)
    for ext in ('.gz', '.cif', '.mmcif', '.pdb', '.ent'):
        if n.lower().endswith(ext):
            n = n[:-len(ext)]
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pdb', nargs='*', default=[]); ap.add_argument('--csv')
    ap.add_argument('--model', nargs='*', default=[], help='local PDB/mmCIF model files (e.g. your own reconstructions)')
    ap.add_argument('--resume', action='store_true', help='skip structures already finished in --out/per_structure')
    ap.add_argument('--out', default='gate_out'); ap.add_argument('--cache', default='pdbcache')
    ap.add_argument('--cutoff', type=float, default=8.0, help='max axis distance (A) of a gate-guarding residue')
    ap.add_argument('--sens', type=float, nargs='*', default=[6.0, 10.0], help='extra cutoffs for sensitivity')
    ap.add_argument('--atoms', choices=['sc', 'all'], default='sc', help='side-chain atoms (sc) or all heavy atoms')
    ap.add_argument('--background', choices=['matched', 'all'], default='matched')
    ap.add_argument('--gatekeepers', type=int, default=5)
    ap.add_argument('--max-mb', type=float, default=25)
    ap.add_argument('--residue', type=int, help='print axis distances of this residue number in every chain')
    ap.add_argument('--no-plots', action='store_true')
    ap.add_argument('--reuse', help='residue_axis_distances.csv.gz from an earlier run: redo only the statistics (new cutoff/background/atoms)')
    ap.add_argument('--min-fg', type=int, default=15, help='heatmap: blank cells of families with fewer foreground residues')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(os.path.join(a.out, 'per_structure'), exist_ok=True)
    # a native crash (access violation, missing DLL entry point ...) kills the process without a Python exception;
    # faulthandler writes the Python traceback of the moment of death into the run folder
    fh_path = os.path.join(a.out, 'faulthandler.log')
    fh_file = open(fh_path, 'a', encoding='utf-8')
    faulthandler.enable(file=fh_file, all_threads=True)
    def _drop_empty_faultlog():
        try:
            fh_file.flush()
            if os.path.getsize(fh_path) == 0:
                faulthandler.disable(); fh_file.close(); os.remove(fh_path)
        except Exception:
            pass
    atexit.register(_drop_empty_faultlog)
    sys.stdout = sys.stderr = Tee(os.path.join(a.out, 'run.log'))
    print(f'gate_residue_analysis  {datetime.datetime.now():%Y-%m-%d %H:%M:%S}  python {sys.version.split()[0]}')
    print('arguments:', ' '.join(sys.argv[1:]))

    jobs = [dict(pdb_id=p, group=p.lower(), family='NA', T_label='NA') for p in a.pdb]
    if a.csv:
        t = pd.read_csv(a.csv)
        for r in t.itertuples():
            jobs.append(dict(pdb_id=r.pdb_id, group=getattr(r, 'group', r.pdb_id), family=getattr(r, 'family', 'NA'),
                             T_label=getattr(r, 'T_label', 'NA')))
    for m in a.model:
        jobs.append(dict(pdb_id=local_name(m).lower(), group=local_name(m), family='NA', T_label='NA', path=m))
    if not jobs and not a.reuse:
        sys.exit('nothing to do: give --pdb, --csv or --model')
    done_groups, frames, used, skipped = set(), [], [], []
    n_jobs = len(jobs)
    for ji, j in enumerate([] if a.reuse else jobs, 1):
        if j['group'] in done_groups:
            continue
        pid = j['pdb_id'].lower()
        pfile = os.path.join(a.out, 'per_structure', f'{pid}.csv.gz'); mfile = os.path.join(a.out, 'per_structure', f'{pid}.json')
        if a.resume and os.path.exists(pfile) and os.path.exists(mfile):
            frames.append(pd.read_csv(pfile, dtype={'pdb': str})); used.append(json.load(open(mfile))); done_groups.add(j['group'])
            print(f'[{ji}/{n_jobs}] {pid}: already finished, reused'); continue
        print(f'[{ji}/{n_jobs}] {pid} ...', flush=True)
        try:
            src = j['path'] if j.get('path') else fetch_cif(pid, a.cache, a.max_mb)
            st = gemmi.read_structure(src); st.setup_entities()
            ops, keep, model = assembly_operators(st), None, st[0]
            done, fit_rms, mism = False, 0.0, 0.0
            if ops is not None:
                try:
                    c, axes, fit_rms = icosahedral_axes(close_group(ops)); done = True
                except RuntimeError:
                    pass
            if not done:                         # complete capsid (or noisy operators): use coordinates
                model = full_capsid_model(st)
                ops, keep, mism = ops_from_coordinates(model)
                c, axes, fit_rms = icosahedral_axes(ops, noisy=True)
                if fit_rms > 5.0:
                    raise RuntimeError(f'icosahedral axes fit poorly (rms {fit_rms:.1f} deg)')
            df = residue_table(model, c, axes, pid, keep)
            if len(df) < 50:
                raise RuntimeError('too few protein residues')
        except Exception as ex:
            skipped.append((pid, j['group'], str(ex)[:160])); print(f'[skip] {pid} ({j["group"]}): {ex}', file=sys.stderr)
            if frames:
                write_gate_outputs(pd.concat(frames, ignore_index=True), used, skipped, a.out, a.cutoff, a.atoms)
            continue
        df['cls'] = df.aa.map(CLASS); df['group'] = j['group']; df['family'] = j['family']; df['T_label'] = j['T_label']
        frames.append(df); done_groups.add(j['group'])
        used.append(dict(pdb=pid, group=j['group'], family=j['family'], T_label=j['T_label'], n_chains=df.chain.nunique(),
                         n_residues=len(df), n_ops=len(ops), symmetry_from='coordinates' if keep is not None else 'assembly operators', axis_fit_rms_deg=round(fit_rms, 2), centroid_mismatch_A=round(mism, 1), center_x=c[0], center_y=c[1], center_z=c[2]))
        print(f'[ok]   {pid} ({j["group"]}): {len(df)} residues, {df.chain.nunique()} chains')
        df.round(3).to_csv(pfile, index=False); json.dump(used[-1], open(mfile, 'w'))
        write_gate_outputs(pd.concat(frames, ignore_index=True), used, skipped, a.out, a.cutoff, a.atoms)   # save as we go
    if a.reuse:
        R = pd.read_csv(a.reuse, dtype={'pdb': str})
        write_gate_outputs(R, [], [], a.out, a.cutoff, a.atoms)
    else:
        if not frames:
            why = '; '.join(f'{p}: {r}' for p, _, r in skipped[:5])
            pd.DataFrame(skipped, columns=['pdb', 'group', 'reason']).to_csv(os.path.join(a.out, 'structures_skipped.csv'), index=False)
            sys.exit('no structure could be analysed (nothing was saved except run.log). Reasons: ' + (why or 'none recorded'))
        R = pd.concat(frames, ignore_index=True)
        write_gate_outputs(R, used, skipped, a.out, a.cutoff, a.atoms)
    print(f'\nSaved gate_summary.txt and gate_residues_all.csv for {R.pdb.nunique()} structure(s) in {a.out}', flush=True)

    if a.residue is not None:
        q = R[R.auth_seq == a.residue][['pdb', 'chain', 'resname', 'auth_seq', 'label_seq', f'd5_{a.atoms}', f'd3_{a.atoms}', f'd2_{a.atoms}']]
        print(f'\nResidue {a.residue}: distance (A) of the {"side chain" if a.atoms == "sc" else "nearest heavy atom"} to the nearest axis')
        print(q.round(1).to_string(index=False))
        q.round(2).to_csv(os.path.join(a.out, f'residue_{a.residue}_query.csv'), index=False)

    report = []
    small = R.pdb.nunique() < 5
    for ci, cutoff in enumerate([] if small else [a.cutoff] + [s for s in a.sens if s != a.cutoff]):
        rows, dual_rows, gk_rows = [], [], []
        for pid, df in R.groupby('pdb', sort=False):
            meta = dict(pdb=pid, group=df.group.iloc[0], family=df.family.iloc[0], T_label=df.T_label.iloc[0])
            sets, bg = foreground_sets(df, cutoff, a.atoms, a.gatekeepers)
            for z, fg in sets.items():
                rows += enrichment(fg, background_for(fg, bg, a.background), meta, z)
            if ci == 0:
                d = sets['dual 5/3-fold position']
                for r in d.itertuples():
                    n5 = df[(df.pos == r.pos) & (df[f'd5_{a.atoms}'] <= cutoff)]
                    n3 = df[(df.pos == r.pos) & (df[f'd3_{a.atoms}'] <= cutoff)]
                    dual_rows.append(dict(**meta, position=r.pos, auth_seq=r.auth_seq, aa=r.aa, resname=r.resname,
                                          chains_near_5fold=','.join(n5.chain), min_d5=n5[f'd5_{a.atoms}'].min(),
                                          chains_near_3fold=','.join(n3.chain), min_d3=n3[f'd3_{a.atoms}'].min()))
                for z in [k for k in sets if 'gatekeepers' in k]:
                    n = int(z[0])
                    for r in sets[z].itertuples():
                        gk_rows.append(dict(**meta, zone=z, chain=r.chain, auth_seq=r.auth_seq, aa=r.aa,
                                            distance_to_axis=getattr(r, f'd{n}_{a.atoms}'), r_from_center=getattr(r, f'r{n}_sc')))
        enr = pd.DataFrame(rows)
        pol, con = pooled(enr, loo=(ci == 0)), consistency(enr)
        if ci == 0:
            enr.round(4).to_csv(os.path.join(a.out, 'enrichment_per_structure.csv'), index=False)
            pd.DataFrame(dual_rows).round(2).to_csv(os.path.join(a.out, 'dual_gate_residues.csv'), index=False)
            pd.DataFrame(gk_rows).round(2).to_csv(os.path.join(a.out, 'gatekeepers.csv'), index=False)
            if not a.no_plots and R.pdb.nunique() >= 5:
                figures(enr, pol, ['5-fold', '3-fold', 'dual 5/3-fold position'], a.out, cutoff, a.min_fg)
        tag = f'cutoff{cutoff:g}'
        pol.round(5).to_csv(os.path.join(a.out, f'pooled_enrichment_{tag}.csv'), index=False)
        con.round(5).to_csv(os.path.join(a.out, f'consistency_across_viruses_{tag}.csv'), index=False)
        report.append(f'\n=== cutoff {cutoff:g} A, atoms={a.atoms}, background={a.background}: pooled CMH, FDR<0.05, amino acids ===')
        s = pol[(pol.level == 'aa') & (pol['cmh_fdr'] < 0.05)].sort_values(['zone', 'log2_or'])
        report.append(s[['zone', 'key', 'n_structures', 'fg_with_key', 'fg_total', 'freq_fg', 'freq_bg', 'cmh_odds_ratio', 'cmh_fdr']]
                      .round(3).to_string(index=False) if len(s) else 'none')
    if small:
        report.append(f'Only {R.pdb.nunique()} structure(s) analysed. Pooled enrichment statistics and figures need at least 5 different viruses '
                      '(with fewer, p-values would be meaningless), so they were not computed.\n'
                      'The result for these structures is in gate_summary.txt and gate_residues_all.csv: the residues at each gate, '
                      'with distances, and the positions near both a 5-fold and a 3-fold axis.')
    txt = '\n'.join(report)
    open(os.path.join(a.out, 'report.txt'), 'w', encoding='utf-8').write(txt)
    print(txt)
    must = ['gate_summary.txt', 'gate_residues_all.csv', 'residue_axis_distances.csv.gz', 'report.txt']
    missing = [f for f in must if not os.path.exists(os.path.join(a.out, f))]
    if missing:
        sys.exit('ERROR: expected output files are missing: ' + ', '.join(missing))
    print('\nDONE. Files written:', ', '.join(sorted(f for f in os.listdir(a.out) if os.path.isfile(os.path.join(a.out, f)) and not (f == 'faulthandler.log' and os.path.getsize(os.path.join(a.out, f)) == 0))))


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    try:
        main()
    except SystemExit as e:
        if e.code not in (None, 0):
            print('\nSTOPPED:', e.code, file=sys.stderr)
        raise
    except Exception:
        traceback.print_exc()
        print('\nSTOPPED because of the error above.', file=sys.stderr)
        sys.exit(1)
