r"""
gate_stemleaf.py -- back-to-back stem-and-leaf plot of every virus in gate_residues.sqlite.

One leaf per virus: the distance of its residue closest to the 5-fold axis (left side) and to the 3-fold axis (right
side). Stem = whole angstroms, leaf = the tenths digit, so "2 | 8" is 2.8 A. Each leaf is coloured by the chemistry
of the amino acid that sits closest (Met is its own colour). Pick one amino acid to ask "how close does Asn get in
each virus?" instead of "what is closest". One structure per virus (the best-resolution one); viruses with nothing
within 8 A are counted underneath. Reads the database read-only.

  python gate_stemleaf.py                         -> text plot in the console
  python gate_stemleaf.py --class "T=3" --aa M --out met_T3.png
"""
import argparse, os, sqlite3
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "gate_residues.sqlite")
ALL = "all classes"
MAXD = 8.0                                               # gate residues are stored out to 8 A
CLASS_OF = {**{a: "aliphatic (A V L I)" for a in "AVLI"}, "M": "Met", **{a: "aromatic (F W Y)" for a in "FWY"},
            **{a: "polar (S T N Q)" for a in "STNQ"}, **{a: "positive (K R H)" for a in "KRH"},
            **{a: "negative (D E)" for a in "DE"}, "G": "Gly / Pro / Cys", "P": "Gly / Pro / Cys", "C": "Gly / Pro / Cys"}
CLASS_COLOR = {"aliphatic (A V L I)": "#7a8794", "Met": "#d4a017", "aromatic (F W Y)": "#7b4fa8",
               "polar (S T N Q)": "#1f7a5a", "positive (K R H)": "#1d5fa8", "negative (D E)": "#b2512f",
               "Gly / Pro / Cys": "#a0a0a0"}
AA_NAME = dict(zip("AVLIMFWYCGPSTNQHKRDE", ("Alanine", "Valine", "Leucine", "Isoleucine", "Methionine", "Phenylalanine",
                                            "Tryptophan", "Tyrosine", "Cysteine", "Glycine", "Proline", "Serine",
                                            "Threonine", "Asparagine", "Glutamine", "Histidine", "Lysine", "Arginine",
                                            "Aspartate", "Glutamate")))


def connect(path=DB):
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    return con


def _viruses(con, cls):
    """{virus: (pdb, family)} using the best-resolution structure of each virus"""
    rows = con.execute("SELECT virus, pdb, family, capsid_class FROM structures ORDER BY virus, "
                       "CASE WHEN resolution_A IS NULL THEN 1 ELSE 0 END, resolution_A, pdb").fetchall()
    pick = {}
    for r in rows:
        pick.setdefault(r["virus"], r)
    return {v: (r["pdb"], r["family"]) for v, r in pick.items() if cls == ALL or r["capsid_class"] == cls}


def data(con, cls=ALL, aa=None):
    """-> dict(gate -> list of dict(d, virus, pdb, family, aa, resnum), n=viruses, missing=gate -> [virus])"""
    vir = _viruses(con, cls)
    by_pdb = {p: (v, f) for v, (p, f) in vir.items()}
    out, missing = {}, {}
    for gate in ("5-fold", "3-fold"):
        q = ("SELECT pdb, aa, auth_seq, distance_A FROM gate_residues WHERE gate = ?" + (" AND aa = ?" if aa else "") +
             " ORDER BY distance_A")
        best = {}
        for r in con.execute(q, (gate, aa) if aa else (gate,)):
            if r["pdb"] in by_pdb and r["pdb"] not in best:
                best[r["pdb"]] = r
        out[gate] = sorted((dict(d=float(r["distance_A"]), virus=by_pdb[p][0], pdb=p, family=by_pdb[p][1],
                                 aa=r["aa"], resnum=r["auth_seq"]) for p, r in best.items()), key=lambda x: x["d"])
        missing[gate] = sorted(by_pdb[p][0] for p in by_pdb if p not in best)
    return dict(gates=out, n=len(vir), missing=missing, cls=cls, aa=aa)


def _split(d):
    """2.86 -> (2, 8): stem = whole A, leaf = tenths digit (truncated, as in a classic stem-and-leaf)"""
    t = int(np.floor(round(d, 2) * 10 + 1e-9))
    return t // 10, t % 10


def text(res):
    """classic back-to-back stem-and-leaf as plain text"""
    g5, g3 = res["gates"]["5-fold"], res["gates"]["3-fold"]
    rows = {s: ([], []) for s in range(int(MAXD) + 1)}
    for side, g in ((0, g5), (1, g3)):
        for x in g:
            s, l = _split(x["d"])
            rows[min(s, int(MAXD))][side].append(str(l))
    w = max([len(v[0]) for v in rows.values()] + [8])
    what = f"closest {AA_NAME[res['aa']]}" if res["aa"] else "closest residue"
    lines = [f"Distance of the {what} to the symmetry axis, one leaf per virus ({res['cls']}, {res['n']} viruses)",
             "stem = whole A, leaf = tenths:  '2 | 8' = 2.8 A", "",
             f"{'5-fold axis':>{w}}   |   3-fold axis", f"{'-' * w}---+---{'-' * w}"]
    for s in sorted(rows):
        left = "".join(sorted(rows[s][0], reverse=True))
        right = "".join(sorted(rows[s][1]))
        lines.append(f"{left:>{w}} {s:>2}| {right}")
    lines.append("")
    for gate, g in (("5-fold", g5), ("3-fold", g3)):
        lines.append(f"{gate}: {len(g)} viruses plotted, median {np.median([x['d'] for x in g]):.1f} A"
                     if g else f"{gate}: no virus has it within {MAXD:g} A")
        if res["missing"][gate]:
            lines.append(f"   not within {MAXD:g} A in {len(res['missing'][gate])} viruses")
    return "\n".join(lines)


def draw(fig, res, show_letters=False):
    """back-to-back stem-and-leaf figure; leaves coloured by amino-acid chemistry. Returns a one-line description."""
    from matplotlib.patches import Patch
    fig.clear()
    fig.set_layout_engine(None)
    ax = fig.add_axes([0.02, 0.10, 0.96, 0.76]); ax.axis("off")
    g5, g3 = res["gates"]["5-fold"], res["gates"]["3-fold"]
    nstem = int(MAXD) + 1
    rows = {s: ([], []) for s in range(nstem)}
    for side, g in ((0, g5), (1, g3)):
        for x in g:
            s, l = _split(x["d"])
            rows[min(s, nstem - 1)][side].append((l, x))
    for s in rows:
        rows[s][0].sort(key=lambda t: (t[1]["d"], t[1]["virus"]))
        rows[s][1].sort(key=lambda t: (t[1]["d"], t[1]["virus"]))
    longest = max([len(v[0]) for v in rows.values()] + [len(v[1]) for v in rows.values()] + [10])
    step = 0.47 / (longest + 1)                            # horizontal room per leaf on each side
    spacing_pt = step * 0.96 * fig.get_figwidth() * 72        # room per leaf, in points
    fs = max(5.0, min(15.0, 1.15 * spacing_pt))
    for i, s in enumerate(sorted(rows)):
        y = 1 - (i + 0.5) / nstem
        ax.text(0.5, y, f"{s}", ha="center", va="center", fontsize=fs + 2, weight="bold", family="monospace")
        for side, sign in ((0, -1), (1, 1)):
            for k, (leaf, x) in enumerate(rows[s][side]):
                xx = 0.5 + sign * (0.025 + (k + 0.5) * step)
                col = CLASS_COLOR[CLASS_OF.get(x["aa"], "Gly / Pro / Cys")]
                ax.text(xx, y, x["aa"] if show_letters else str(leaf), ha="center", va="center", fontsize=fs,
                        family="monospace", color=col, weight="bold")
    ax.plot([0.5 - 0.015, 0.5 - 0.015], [0, 1], color="#333", lw=1.2)
    ax.plot([0.5 + 0.015, 0.5 + 0.015], [0, 1], color="#333", lw=1.2)
    ax.text(0.25, 1.03, f"5-fold axis   ({len(g5)} viruses)", ha="center", fontsize=11, weight="bold", color="#1d5fa8")
    ax.text(0.75, 1.03, f"3-fold axis   ({len(g3)} viruses)", ha="center", fontsize=11, weight="bold", color="#b2512f")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.06)
    what = f"closest {AA_NAME[res['aa']]}" if res["aa"] else "closest residue"
    fig.suptitle(f"Stem-and-leaf: distance of each virus's {what} to the symmetry axis  -  {res['cls']}, "
                 f"{res['n']} viruses, one leaf per virus", x=0.02, ha="left", fontsize=12, weight="bold")
    leaf_txt = ("leaf = the amino acid letter (stem 2 = 2.0-2.9 A)" if show_letters
                else "leaf = tenths digit: stem 2 with leaf 8 = 2.8 A")
    fig.text(0.02, 0.925, f"stem (middle) = whole angstroms, {leaf_txt}.  Leaves run outwards from the stem in order of "
                          f"distance; smaller stems = closer to the axis (in the pore).", fontsize=9, color="#444")
    handles = [Patch(color=c, label=k) for k, c in CLASS_COLOR.items()]
    if not res["aa"]:
        fig.legend(handles=handles, loc="lower center", ncol=len(handles), fontsize=8, frameon=False,
                   title="colour = chemistry of the closest amino acid", title_fontsize=8)
    miss = [f"{g}: {len(res['missing'][g])} viruses have no {'such residue' if res['aa'] else 'residue'} within "
            f"{MAXD:g} A" for g in ("5-fold", "3-fold") if res["missing"][g]]
    if miss:
        fig.text(0.02, 0.895, "Not plotted - " + ";  ".join(miss), ha="left", fontsize=8, color="#666")
    med = lambda g: f"{np.median([x['d'] for x in g]):.1f} A" if g else "-"
    return (f"{res['cls']}  |  {what}  |  {res['n']} viruses  |  median 5-fold {med(g5)}, 3-fold {med(g3)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=DB)
    ap.add_argument("--class", dest="cls", default=ALL)
    ap.add_argument("--aa", help="one-letter amino acid (default: whatever is closest)")
    ap.add_argument("--letters", action="store_true", help="figure: write the amino acid letter as the leaf")
    ap.add_argument("--out", help="PNG/PDF for the figure (the text version is always printed)")
    a = ap.parse_args()
    con = connect(a.db)
    res = data(con, a.cls, a.aa.upper() if a.aa else None)
    print(text(res))
    if a.out:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(16, 9))
        draw(fig, res, a.letters)
        fig.savefig(a.out, dpi=150)
        print("written:", a.out)


if __name__ == "__main__":
    main()
