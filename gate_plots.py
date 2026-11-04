r"""
gate_plots.py -- figures from gate_residues.sqlite: which amino acids recur at the 5-fold and 3-fold gates,
which sit at positions that guard BOTH gates, and which are favoured there compared with the rest of the shell.

Reads only what gate_db.py stored (gate_residues, dual_positions, preferences); computes no new statistics.

  A  favoured at both gates?   log2 odds ratio at the 5-fold gate (x) against the 3-fold gate (y), one point per
                               amino acid; upper-right = enriched at both. Ringed = FDR < 0.05 at both gates.
  B  enrichment heatmap        log2 odds ratio for 5-fold, 3-fold and dual 5/3 positions; * = FDR < 0.05.
  C  how often it recurs       % of viruses that have the amino acid at the 5-fold gate, at the 3-fold gate,
                               and at a dual position (same residue number near both axes).
  D  dual positions by family  the dual-position count of each amino acid split by virus family, so a pattern
                               carried by many families can be told apart from one family of related viruses.

  python gate_plots.py                      -> all capsid classes, 8 A cutoff, plots\gate_plots_<time>.png
  python gate_plots.py --class "T=3" --cutoff 6 --out my.png
"""
import argparse, datetime, os, sqlite3
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "gate_residues.sqlite")
AA_ORDER = list("AVLIMFWYCGPSTNQHKRDE")
ZONES = [("5-fold", "5-fold"), ("3-fold", "3-fold"), ("dual 5/3-fold position", "dual 5/3")]
CLASS_COLOR = ["#1d5fa8", "#b2512f", "#1f7a5a", "#7b4fa8", "#8a7a1f"]
GATE_COLOR = {"5-fold": "#1d5fa8", "3-fold": "#b2512f", "dual": "#e0a030"}
ALL = "all classes"


def classes(con):
    return [r[0] for r in con.execute("SELECT DISTINCT capsid_class FROM structures WHERE capsid_class IS NOT NULL ORDER BY 1")]


def cutoffs(con):
    return [r[0] for r in con.execute("SELECT DISTINCT cutoff_A FROM preferences ORDER BY 1")] or [8.0]


def _log2(o):
    return float(np.clip(np.log2(o), -3, 3)) if o and o > 0 else -3.0


def build(fig, con, cls=ALL, cutoff=8.0):
    """draw the four panels into a matplotlib Figure; returns a one-line description of what was plotted"""
    import matplotlib
    from matplotlib.patches import Rectangle
    fig.clear()
    cls_list = classes(con) if cls == ALL else [cls]
    where = "" if cls == ALL else " AND capsid_class=?"
    arg = () if cls == ALL else (cls,)
    fig.set_layout_engine("constrained")                  # keeps titles and labels apart at any window size
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.15])
    axA, axB, axC, axD = (fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
                          fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]))

    # preferences: (class, zone, aa) -> (odds ratio, fdr, n_structures)
    pref = {}
    for c, z, k, o, f, n in con.execute("SELECT capsid_class, zone, key, odds_ratio, fdr, n_structures FROM preferences "
                                        "WHERE level='aa' AND cutoff_A=?" + where, (cutoff,) + arg):
        pref[(c, z, k)] = (o, f, n)

    # ---- A: 5-fold vs 3-fold enrichment
    axA.axhline(0, color="#999", lw=0.8); axA.axvline(0, color="#999", lw=0.8)
    axA.add_patch(Rectangle((0, 0), 3.3, 3.3, color="#fbe9df", zorder=0))
    axA.text(3.2, 3.2, "favoured at\nBOTH gates", ha="right", va="top", fontsize=9, color="#b2512f", weight="bold")
    for ci, c in enumerate(cls_list):
        col = CLASS_COLOR[ci % len(CLASS_COLOR)]
        for a in AA_ORDER:
            p5, p3 = pref.get((c, "5-fold", a)), pref.get((c, "3-fold", a))
            if not (p5 and p3):
                continue
            x, y = _log2(p5[0]), _log2(p3[0])
            sig5, sig3 = p5[1] is not None and p5[1] < 0.05, p3[1] is not None and p3[1] < 0.05
            axA.scatter(x, y, s=70 if (sig5 or sig3) else 32, color=col, alpha=0.85, zorder=3,
                        edgecolor="black" if (sig5 and sig3) else "white", linewidth=1.6 if (sig5 and sig3) else 0.6,
                        label=c if a == AA_ORDER[0] else None)
            axA.annotate(a, (x, y), xytext=(4, 3), textcoords="offset points", fontsize=8, color=col, zorder=4)
    axA.set_xlim(-3.3, 3.3); axA.set_ylim(-3.3, 3.3)
    axA.set_xlabel("log2 odds ratio at the 5-fold gate"); axA.set_ylabel("log2 odds ratio at the 3-fold gate")
    axA.set_title("A  Favoured at both gates?", loc="left", fontsize=11, weight="bold")
    if len(cls_list) > 1:
        axA.legend(fontsize=8, loc="upper left", frameon=True, framealpha=0.9)
    axA.text(-3.2, -3.2, "big = FDR<0.05 at one gate; ringed = both", fontsize=7.5, color="#555", va="bottom")

    # ---- B: heatmap of log2 OR (rows amino acids, columns class x zone)
    cols = [(c, z, lab) for c in cls_list for z, lab in ZONES]
    M = np.full((len(AA_ORDER), len(cols)), np.nan)
    for j, (c, z, _) in enumerate(cols):
        for i, a in enumerate(AA_ORDER):
            p = pref.get((c, z, a))
            if p:
                M[i, j] = _log2(p[0])
                if p[1] is not None and p[1] < 0.05:
                    axB.text(j, i, "*", ha="center", va="center", fontsize=12, color="black", weight="bold")
    cmap = matplotlib.colormaps["RdBu_r"].copy(); cmap.set_bad("#e6e6e6")
    im = axB.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=-2.5, vmax=2.5, aspect="auto")
    axB.set_yticks(range(len(AA_ORDER)), AA_ORDER, fontsize=8)
    axB.set_xticks(range(len(cols)), [lab if len(cls_list) == 1 else f"{lab}\n{c}" for c, _, lab in cols],
                   fontsize=7.5, rotation=0 if len(cls_list) == 1 else 60, ha="center" if len(cls_list) == 1 else "right")
    for j in range(len(ZONES), len(cols), len(ZONES)):
        axB.axvline(j - 0.5, color="white", lw=3)
    fig.colorbar(im, ax=axB, fraction=0.04, pad=0.02).set_label("log2 odds ratio (red = favoured)", fontsize=8)
    axB.set_title("B  Enrichment vs rest of shell  (* FDR<0.05, grey = none)",
                  loc="left", fontsize=11, weight="bold")

    # ---- C: % of viruses having the amino acid at each gate / at a dual position
    # gate residues are stored out to 8 A, so C and D follow the chosen cutoff up to 8 A
    dlim = min(cutoff, 8.0)
    dual_w = " AND min_d5 <= ? AND min_d3 <= ?" + where
    dual_a = (dlim, dlim) + arg
    nv = con.execute("SELECT COUNT(DISTINCT virus) FROM structures WHERE 1=1" + where, arg).fetchone()[0] or 1
    cnt = {g: dict(con.execute("SELECT aa, COUNT(DISTINCT virus) FROM gate_residues WHERE gate=? AND distance_A <= ?" +
                               where + " GROUP BY aa", (g, dlim) + arg).fetchall()) for g in ("5-fold", "3-fold")}
    cnt["dual"] = dict(con.execute("SELECT aa, COUNT(DISTINCT virus) FROM dual_positions WHERE 1=1" + dual_w +
                                   " GROUP BY aa", dual_a).fetchall())
    xs = np.arange(len(AA_ORDER)); w = 0.27
    for k, (g, lab) in enumerate((("5-fold", "at the 5-fold gate"), ("3-fold", "at the 3-fold gate"),
                                  ("dual", "at a dual 5/3 position"))):
        axC.bar(xs + (k - 1) * w, [100 * cnt[g].get(a, 0) / nv for a in AA_ORDER], w, color=GATE_COLOR[g], label=lab)
    axC.set_xticks(xs, AA_ORDER); axC.set_ylabel(f"% of viruses (n = {nv})")
    axC.legend(fontsize=8, frameon=False)
    axC.set_title("C  How often each amino acid recurs", loc="left", fontsize=11, weight="bold")
    axC.grid(axis="y", color="#eee"); axC.set_axisbelow(True)

    # ---- D: dual positions per amino acid, split by family
    rows = con.execute("SELECT aa, IFNULL(family,'unknown'), COUNT(DISTINCT virus) FROM dual_positions WHERE 1=1" + dual_w +
                       " GROUP BY 1, 2", dual_a).fetchall()
    fam_tot = {}
    for a, f, n in rows:
        fam_tot[f] = fam_tot.get(f, 0) + n
    top = [f for f, _ in sorted(fam_tot.items(), key=lambda t: -t[1])][:7]
    pal = matplotlib.colormaps["tab10"]
    bottom = np.zeros(len(AA_ORDER))
    for fi, f in enumerate(top + ["other families"]):
        v = np.array([sum(n for a2, f2, n in rows if a2 == a and (f2 == f if f in top else f2 not in top))
                      for a in AA_ORDER], float)
        if v.sum():
            axD.bar(xs, v, 0.7, bottom=bottom, color=pal(fi) if f in top else "#bbbbbb", label=f)
            bottom += v
    nfam = {a: len({f for a2, f, n in rows if a2 == a}) for a in AA_ORDER}
    for i, a in enumerate(AA_ORDER):
        if bottom[i]:
            axD.text(i, bottom[i] + 0.3, f"{nfam[a]}f", ha="center", fontsize=7, color="#555")
    axD.set_xticks(xs, AA_ORDER); axD.set_ylabel("viruses with a dual 5/3 position")
    axD.set_ylim(0, max(bottom.max(), 1) * 1.45)           # headroom for the family legend
    axD.legend(fontsize=7, frameon=False, ncol=4, loc="upper center")
    axD.set_title("D  Residues at BOTH gates, by family  (Nf = families)",
                  loc="left", fontsize=11, weight="bold")
    axD.grid(axis="y", color="#eee"); axD.set_axisbelow(True)

    ndual = con.execute("SELECT COUNT(*), COUNT(DISTINCT virus) FROM dual_positions WHERE 1=1" + dual_w, dual_a).fetchone()
    desc = (f"{cls}  |  cutoff {cutoff:g} A  |  {nv} viruses  |  {ndual[0]} dual 5/3 positions in {ndual[1]} viruses"
            + ("  (C, D: within 8 A, the largest stored)" if cutoff > 8 else ""))
    fig.suptitle("Gate residues of icosahedral capsids  -  " + desc, fontsize=12, x=0.01, ha="left")
    return desc


AA_NAME = dict(A="Alanine", V="Valine", L="Leucine", I="Isoleucine", M="Methionine", F="Phenylalanine",
               W="Tryptophan", Y="Tyrosine", C="Cysteine", G="Glycine", P="Proline", S="Serine", T="Threonine",
               N="Asparagine", Q="Glutamine", H="Histidine", K="Lysine", R="Arginine", D="Aspartate", E="Glutamate")


def _answer(rows):
    """plain-English sentence from [(aa, or5, fdr5, or3, fdr3)] already sorted best-first"""
    sig = lambda f: f is not None and f < 0.05
    both_sig = [a for a, o5, f5, o3, f3 in rows if o5 > 1 and o3 > 1 and sig(f5) and sig(f3)]
    both_fav = [a for a, o5, f5, o3, f3 in rows if o5 > 1 and o3 > 1 and a not in both_sig]
    avoid = [r for r in rows if r[1] < 1 and r[3] < 1]
    avoid_sig = [r[0] for r in reversed(avoid) if sig(r[2]) and sig(r[4])]
    avoid_other = [r[0] for r in reversed(avoid) if r[0] not in avoid_sig]
    j = lambda xs: ", ".join(xs) if xs else "none"
    s = f"Clearly favoured at BOTH gates (significant at both): {j(both_sig)}.   "
    s += f"Favoured at both, not yet significant: {j(both_fav)}.   "
    if avoid_sig:
        s += f"Avoided at both: {j(avoid_sig)} (significant)" + (f", also {j(avoid_other)}" if avoid_other else "") + "."
    else:
        s += f"Avoided at both (none significant yet): {j(avoid_other)}."
    return s


def build_simple(fig, con, cls=ALL, cutoff=8.0):
    """one easy question: which amino acids are favoured at BOTH the 5-fold and the 3-fold gate?
    Paired bars of log2 odds ratio per amino acid, sorted so 'favoured at both' is on top."""
    import textwrap
    fig.clear()
    fig.set_layout_engine("constrained")
    cls_list = classes(con) if cls == ALL else [cls]
    pref = {}
    for c, z, k, o, f, n in con.execute("SELECT capsid_class, zone, key, odds_ratio, fdr, n_structures FROM preferences "
                                        "WHERE level='aa' AND cutoff_A=? AND zone IN ('5-fold', '3-fold')", (cutoff,)):
        pref[(c, z, k)] = (o, f, n)
    axes = fig.subplots(1, len(cls_list), squeeze=False)[0]
    small = fig.get_figwidth() * fig.dpi < 1100                 # embedded in a window rather than saved big
    fs = 7.5 if small else 9
    afs = fs + 1                                                # answer text size
    wrap = max(60, int(fig.get_figwidth() * 72 / (afs * 0.6)))  # characters per line across the whole figure
    descs, answers = [], []
    for ax, c in zip(axes, cls_list):
        rows = []
        for a in AA_ORDER:
            p5, p3 = pref.get((c, "5-fold", a)), pref.get((c, "3-fold", a))
            if p5 and p3 and p5[0] is not None and p3[0] is not None:
                rows.append((a, p5[0], p5[1], p3[0], p3[1]))
        if not rows:
            ax.set_axis_off()
            ax.text(0.5, 0.5, "no pooled statistics for this selection\n(needs >= 5 viruses)",
                    ha="center", va="center", fontsize=fs + 1, color="#555", transform=ax.transAxes)
            ax.set_title(f"{c}", loc="left", fontsize=fs + 2, weight="bold")
            descs.append(f"{c}: no data")
            continue
        rows.sort(key=lambda r: -min(_log2(r[1]), _log2(r[3])))
        nvir = next(pref[(c, "5-fold", r[0])][2] for r in rows)
        y = np.arange(len(rows))[::-1]                          # first row at the top
        h = 0.38
        for yi, (a, o5, f5, o3, f3) in zip(y, rows):
            if o5 > 1 and o3 > 1:
                ax.axhspan(yi - 0.5, yi + 0.5, color="#fbe9df", zorder=0, lw=0)
            for dy, o, f, col in ((h / 2, o5, f5, GATE_COLOR["5-fold"]), (-h / 2, o3, f3, GATE_COLOR["3-fold"])):
                v = _log2(o); sg = f is not None and f < 0.05
                ax.barh(yi + dy, v, h, color=col, edgecolor="black" if sg else col, linewidth=1.2 if sg else 0, zorder=2)
                if sg:
                    ax.text(v + (0.08 if v >= 0 else -0.08), yi + dy, "*", ha="left" if v >= 0 else "right",
                            va="center", fontsize=fs + 2, weight="bold", zorder=3)
        ax.axvline(0, color="#333", lw=0.9, zorder=1)
        ax.set_yticks(y, [f"{r[0]}  {AA_NAME[r[0]]}" for r in rows], fontsize=fs)
        ax.set_ylim(-0.6, len(rows) - 0.4)
        ax.set_xlim(-3.6, 3.6)
        ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])
        ax.tick_params(axis="x", labelsize=fs)
        ax.grid(axis="x", color="#eee"); ax.set_axisbelow(True)
        ax.set_xlabel("← avoided at the gate      favoured at the gate →\n"
                      "(log2 odds ratio vs the rest of the shell)", fontsize=fs)
        title = f"{c}  ({nvir} viruses)"
        ax.set_title(title, loc="left", fontsize=fs + 1.5, weight="bold")
        answers.append(textwrap.fill(f"{title}:   " + _answer(rows), wrap, subsequent_indent="        "))
        descs.append(f"{c}: {nvir} viruses")
    # shared legend
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=GATE_COLOR["5-fold"], label="at the 5-fold gate"),
                        Patch(color=GATE_COLOR["3-fold"], label="at the 3-fold gate"),
                        Patch(facecolor="white", edgecolor="black", label="* = significant (FDR < 0.05)"),
                        Patch(color="#fbe9df", label="shaded row = favoured at both gates")],
               loc="outside lower center", ncol=4, fontsize=fs, frameon=False)
    # headline question + the plain-English answer, one paragraph per class (suptitle so constrained layout
    # reserves room for it at any window size)
    head = f"WHICH AMINO ACIDS ARE FAVOURED AT BOTH THE 5-FOLD AND THE 3-FOLD GATE?   (cutoff {cutoff:g} A)"
    note = "Top rows = best at both gates.  A bar at -3 means the amino acid was never seen at that gate."
    fig.suptitle("\n".join([head] + answers + [note]), fontsize=afs, x=0.01, ha="left", linespacing=1.5)
    return f"{cls}  |  cutoff {cutoff:g} A  |  " + "; ".join(descs)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=DB)
    ap.add_argument("--class", dest="cls", default=ALL, help='capsid class, e.g. "T=3" (default: all classes)')
    ap.add_argument("--cutoff", type=float, default=8.0)
    ap.add_argument("--out", help="PNG or PDF to write (default plots\\gate_plots_<time>.png)")
    ap.add_argument("--simple", action="store_true",
                    help="one easy figure: which amino acids are favoured at BOTH the 5-fold and 3-fold gate")
    a = ap.parse_args()
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    out = a.out or os.path.join(HERE, "plots", f"gate_plots_{datetime.datetime.now():%Y%m%d_%H%M%S}.png")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig = plt.figure(figsize=(16, 11))
    print((build_simple if a.simple else build)(fig, con, a.cls, a.cutoff))
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("written:", out)


if __name__ == "__main__":
    main()
