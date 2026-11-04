r"""
gate_compare.py -- compare two groups of viruses (capsid types or families) at the 5-fold gate, the 3-fold gate or
at positions that guard both gates: which amino acids does each group favour there, and is the difference real?

Unit of comparison is the VIRUS, not the residue. For each virus (its best-resolution structure) we take the residues
at the gate and compute each amino acid's SHARE of them ("8% of this virus's 5-fold gate residues are Asn"). Shares,
not raw counts, because capsid types differ in how many residues sit near an axis: with counts, the group with bigger
gates would look richer in every amino acid. The difference of the two groups' average shares is tested with a
permutation test (10,000 label shuffles, fixed seed) per amino acid, Benjamini-Hochberg FDR over the 20 amino acids. Also reported: in how many viruses of each group the
amino acid is present at the gate at all. Viruses with no residue at that gate are left out of the share test.
Reads gate_residues.sqlite read-only; computes nothing that is stored.

  python gate_compare.py --a "family=Nodaviridae" --b "family=Caliciviridae" --gate 5-fold
  python gate_compare.py --a "class=T=3" --b "class=T=4" --gate both --out t3_vs_t4.png
  python gate_compare.py --a "family=Picornaviridae" --b rest --gate 3-fold
"""
import argparse, os, sqlite3
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "gate_residues.sqlite")
AA_ORDER = "AVLIMFWYCGPSTNQHKRDE"
AA_NAME = dict(zip(AA_ORDER, ("Alanine", "Valine", "Leucine", "Isoleucine", "Methionine", "Phenylalanine", "Tryptophan",
                              "Tyrosine", "Cysteine", "Glycine", "Proline", "Serine", "Threonine", "Asparagine",
                              "Glutamine", "Histidine", "Lysine", "Arginine", "Aspartate", "Glutamate")))
GATES = {"5-fold": "the 5-fold gate", "3-fold": "the 3-fold gate", "both": "BOTH gates (same residue number)"}
REST = "rest"
N_PERM = 10000
COL_A, COL_B = "#1d5fa8", "#b2512f"


def connect(path=DB):
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    return con


def group_choices(con):
    """-> list of (label, spec) for the drop-downs; spec is 'class=<x>' or 'family=<x>'. Counts are the viruses the
    comparison will actually use (one structure per virus, so a virus with structures of two types counts once)."""
    out = []
    for (c,) in con.execute("SELECT DISTINCT capsid_class FROM structures WHERE capsid_class IS NOT NULL ORDER BY 1"):
        out.append((f"capsid type {c}  ({len(_members(con, 'class=' + c))} viruses)", f"class={c}"))
    fams = []
    for (f,) in con.execute("SELECT DISTINCT family FROM structures WHERE family IS NOT NULL AND family != 'Unclassified'"):
        n = len(_members(con, "family=" + f))
        if n >= 2:
            fams.append((n, f))
    for n, f in sorted(fams, key=lambda t: (-t[0], t[1])):
        out.append((f"family {f}  ({n} viruses)", f"family={f}"))
    return out


def label(spec):
    if spec == REST:
        return "all other viruses"
    kind, val = spec.split("=", 1)
    return ("capsid type " if kind == "class" else "") + val


def _members(con, spec, exclude=()):
    """{virus: pdb}, one structure per virus (best resolution)"""
    rows = con.execute("SELECT virus, pdb, capsid_class, family, resolution_A FROM structures ORDER BY virus, "
                       "CASE WHEN resolution_A IS NULL THEN 1 ELSE 0 END, resolution_A, pdb").fetchall()
    pick = {}
    for r in rows:
        pick.setdefault(r["virus"], r)
    if spec == REST:
        return {v: r["pdb"] for v, r in pick.items() if v not in exclude}
    kind, val = spec.split("=", 1)
    col = "capsid_class" if kind == "class" else "family"
    return {v: r["pdb"] for v, r in pick.items() if r[col] == val}


def _counts(con, pdbs, gate, dmax):
    """{pdb: {aa: n}} of residues at the gate within dmax A (every pdb present, possibly empty)"""
    out = {p: {} for p in pdbs}
    if not pdbs:
        return out
    q = ",".join("?" * len(pdbs))
    if gate == "both":
        rows = con.execute(f"SELECT pdb, aa FROM dual_positions WHERE pdb IN ({q}) AND min_d5 <= ? AND min_d3 <= ?",
                           (*pdbs, dmax, dmax))
    else:
        rows = con.execute(f"SELECT pdb, aa FROM gate_residues WHERE pdb IN ({q}) AND gate = ? AND distance_A <= ?",
                           (*pdbs, gate, dmax))
    for r in rows:
        if r["aa"] in AA_NAME:
            out[r["pdb"]][r["aa"]] = out[r["pdb"]].get(r["aa"], 0) + 1
    return out


def compare(con, spec_a, spec_b, gate="5-fold", dmax=8.0):
    """-> dict(a, b, gate, dmax, na, nb, ua, ub, rows=[per amino acid, sorted by share difference A - B], ...)"""
    from statsmodels.stats.multitest import multipletests
    ma = _members(con, spec_a)
    mb = _members(con, spec_b, exclude=set(ma)) if spec_b == REST else _members(con, spec_b)
    overlap = set(ma) & set(mb)
    ca, cb = _counts(con, list(ma.values()), gate, dmax), _counts(con, list(mb.values()), gate, dmax)
    use_a = [c for c in ca.values() if c]                   # viruses with at least one residue at this gate
    use_b = [c for c in cb.values() if c]

    # viruses x amino acids matrix of shares; permutation test on the difference of the group means
    S = np.array([[c.get(a, 0) / sum(c.values()) for a in AA_ORDER] for c in use_a + use_b]).reshape(-1, len(AA_ORDER))
    n1 = len(use_a)
    obs = S[:n1].mean(0) - S[n1:].mean(0) if n1 and len(use_b) else np.zeros(len(AA_ORDER))
    pvals = np.ones(len(AA_ORDER))
    if n1 >= 2 and len(use_b) >= 2:
        rng = np.random.default_rng(20261031)               # fixed seed: the same comparison always gives the same p
        hits = np.zeros(len(AA_ORDER))
        for _ in range(N_PERM):
            idx = rng.permutation(len(S))
            d = S[idx[:n1]].mean(0) - S[idx[n1:]].mean(0)
            hits += np.abs(d) >= np.abs(obs) - 1e-12
        pvals = (hits + 1) / (N_PERM + 1)
    rows = []
    for j, a in enumerate(AA_ORDER):
        ka = sum(1 for c in ca.values() if c.get(a))
        kb = sum(1 for c in cb.values() if c.get(a))
        rows.append(dict(aa=a, name=AA_NAME[a], sa=float(S[:n1, j].mean()) if n1 else 0.0,
                         sb=float(S[n1:, j].mean()) if len(use_b) else 0.0, ka=ka, kb=kb,
                         fa=ka / len(ma) if ma else 0.0, fb=kb / len(mb) if mb else 0.0, p=float(pvals[j])))
    fdr = multipletests([r["p"] for r in rows], method="fdr_bh")[1]
    for r, f in zip(rows, fdr):
        r["fdr"] = float(f)
        r["diff"] = r["sa"] - r["sb"]
    rows.sort(key=lambda r: (-r["diff"], r["aa"]))

    def nfam(m):
        if not m:
            return 0
        q = ",".join("?" * len(m))
        return len({r[0] for r in con.execute(f"SELECT DISTINCT family FROM structures WHERE pdb IN ({q})", list(m.values()))})

    size = lambda use: float(np.mean([sum(c.values()) for c in use])) if use else 0.0
    return dict(a=spec_a, b=spec_b, gate=gate, dmax=dmax, na=len(ma), nb=len(mb), ua=len(use_a), ub=len(use_b),
                rows=rows, overlap=len(overlap), fam_a=nfam(ma), fam_b=nfam(mb), size_a=size(use_a), size_b=size(use_b))


def summary(res):
    """the answer in plain words"""
    A, B, g = label(res["a"]), label(res["b"]), GATES[res["gate"]]
    na, nb, ua, ub = res["na"], res["nb"], res["ua"], res["ub"]
    head = f"{A} ({na} viruses) vs {B} ({nb} viruses), at {g}, within {res['dmax']:g} A of the axis."
    if not na or not nb:
        return head + "  One of the groups has no viruses."
    empty = [(G, n) for G, u, n in ((A, ua, na), (B, ub, nb)) if u == 0]
    if empty:
        why = ("  (Pseudo T=3 capsids are built from three different proteins, so the same residue number cannot "
               "sit at both gates.)" if res["gate"] == "both" and any("pseudo" in G for G, _ in empty) else "")
        return (head + "  " + " and ".join(f"{G} has no residue at this gate in any of its {n} viruses" for G, n in empty)
                + ", so there is nothing to compare." + why)

    def pair(r, w, o):
        return (f"{r['name']} {100 * r['s' + w]:.0f}% vs {100 * r['s' + o]:.0f}% "
                f"(in {r['k' + w]}/{res['n' + w]} vs {r['k' + o]}/{res['n' + o]} viruses)")

    sig = [r for r in res["rows"] if r["fdr"] < 0.05]
    parts = [head, f"A virus has on average {res['size_a']:.0f} residues at this gate in {A} and {res['size_b']:.0f} in {B};"
                   f" the % below is each amino acid's share of them."]
    for grp, rows, w, o in ((A, [r for r in sig if r["diff"] > 0], "a", "b"), (B, [r for r in sig if r["diff"] < 0], "b", "a")):
        if rows:
            parts.append(f"Favoured more in {grp}: " + "; ".join(pair(r, w, o) for r in rows) + ".")
    if not sig:
        top = [r for r in res["rows"] if abs(r["diff"]) >= 0.03][:4]
        parts.append("No difference is statistically clear (FDR < 0.05)." + (
            " Largest gaps, not yet significant: " + "; ".join(
                f"{r['name']} {100 * r['sa']:.0f}% vs {100 * r['sb']:.0f}%" for r in top) + "." if top else ""))
    shared = [r["name"] for r in sorted(res["rows"], key=lambda r: -min(r["fa"], r["fb"])) if min(r["fa"], r["fb"]) >= 0.5]
    if shared:
        parts.append("Found at this gate in at least half the viruses of BOTH groups: " + ", ".join(shared) + ".")
    if min(ua, ub) < 5:
        parts.append(f"Caution: only {min(ua, ub)} viruses with residues at this gate in the smaller group, "
                     f"too few for firm conclusions.")
    if (res["fam_a"] == 1 and na > 1) or (res["fam_b"] == 1 and nb > 1):
        parts.append("Viruses of one family share ancestry, so a family difference can be inherited rather than selected.")
    if res["overlap"]:
        parts.append(f"Note: {res['overlap']} viruses are in both groups.")
    return "  ".join(parts)


def draw(fig, res):
    """paired bars: average share of the gate residues, group A vs group B; returns the summary text"""
    fig.clear()
    fig.set_layout_engine("constrained")
    ax = fig.add_subplot(111)
    rows = res["rows"]
    y = np.arange(len(rows))[::-1]
    h = 0.38
    ax.barh(y + h / 2, [100 * r["sa"] for r in rows], h, color=COL_A, label=f"{label(res['a'])}  ({res['ua']} of {res['na']} viruses have residues at this gate)")
    ax.barh(y - h / 2, [100 * r["sb"] for r in rows], h, color=COL_B, label=f"{label(res['b'])}  ({res['ub']} of {res['nb']} viruses have residues at this gate)")
    top = max([100 * max(r["sa"], r["sb"]) for r in rows] + [1])
    for yi, r in zip(y, rows):
        if r["fdr"] < 0.05:
            ax.axhspan(yi - 0.5, yi + 0.5, color="#fff1c9", zorder=0)
            ax.text(100 * max(r["sa"], r["sb"]) + top * 0.015, yi, f"*  FDR {r['fdr']:.3f}", va="center", fontsize=8)
    ax.set_yticks(y, [f"{r['aa']}  {r['name']}" for r in rows], fontsize=8)
    ax.set_xlim(0, top * 1.3)
    ax.set_xlabel("average share of the gate residues (%)\n"
                  "top = more typical of the first group, bottom = of the second, yellow * = clear difference")
    ax.grid(axis="x", color="#eee")
    ax.set_axisbelow(True)
    ax.legend(loc="center right", fontsize=8, frameon=True, framealpha=0.95)
    ax.set_title(f"{label(res['a'])}  vs  {label(res['b'])}:  {GATES[res['gate']]}, within {res['dmax']:g} A",
                 loc="left", fontsize=10, weight="bold")
    return summary(res)


def write_csv(res, path):
    import csv
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["aa", "name", f"mean_share_pct {label(res['a'])}", f"mean_share_pct {label(res['b'])}",
                    "viruses_with_it_A", "n_A", "viruses_with_it_B", "n_B", "p_permutation", "fdr"])
        for r in res["rows"]:
            w.writerow([r["aa"], r["name"], round(100 * r["sa"], 2), round(100 * r["sb"], 2), r["ka"], res["na"],
                        r["kb"], res["nb"], f"{r['p']:.4g}", f"{r['fdr']:.4g}"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=DB)
    ap.add_argument("--a", required=True, help='group A: "class=T=3" or "family=Nodaviridae"')
    ap.add_argument("--b", required=True, help='group B: same form, or "rest" for all other viruses')
    ap.add_argument("--gate", default="5-fold", choices=list(GATES))
    ap.add_argument("--dmax", type=float, default=8.0, help="distance limit in A (at most the stored 8 A cutoff)")
    ap.add_argument("--out", help="PNG/PDF for the chart; the table is written next to it as .csv")
    a = ap.parse_args()
    con = connect(a.db)
    res = compare(con, a.a, a.b, a.gate, a.dmax)
    print(summary(res))
    print(f"\n{'aa':<4}{'share A':>9}{'share B':>9}{'present A':>12}{'present B':>12}{'FDR':>8}")
    for r in res["rows"]:
        print(f"{r['aa']:<4}{100 * r['sa']:>8.1f}%{100 * r['sb']:>8.1f}%{r['ka']:>7}/{res['na']:<4}"
              f"{r['kb']:>7}/{res['nb']:<4}{r['fdr']:>8.3f}")
    if a.out:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(12, 9))
        draw(fig, res)
        fig.savefig(a.out, dpi=150)
        write_csv(res, os.path.splitext(a.out)[0] + ".csv")
        print("written:", a.out)


if __name__ == "__main__":
    main()
