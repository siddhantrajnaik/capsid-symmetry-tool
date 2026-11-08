r"""
gate_chemistry.py -- pore chemistry: the charge, hydropathy and make-up of the residues that line the 5-fold and the
3-fold pores of each virus, compared across capsid types or families.

Unit of comparison is the VIRUS (its best-resolution structure). The database holds the residues of the deposited
asymmetric unit (ASU) that lie within 8 A of an axis. In the assembled capsid each of these residues occurs `fold`
times around the pore (5 times at the 5-fold axis, 3 times at the 3-fold axis), so counts and the net charge are
multiplied by the fold to describe the whole pore; averages (hydropathy, fractions, volume) are taken over the ASU
residues, which gives the same value as averaging over all copies. Net charge counts Lys and Arg as +1 and Asp and
Glu as -1; His is counted as 0 (mostly uncharged at pH 7), and terminal charges and local pKa shifts are ignored.
Viruses with no residue within the distance limit are kept in the table (n_residues 0, averages NaN) but left out of
the plot and the statistics. Reads gate_residues.sqlite read-only; computes nothing that is stored.

  python gate_chemistry.py --gate 5-fold
  python gate_chemistry.py --gate 3-fold --by family --metric hydropathy --out pore3_hydropathy.png
  python gate_chemistry.py --gate 5-fold --dmax 5 --class "T=3" --metric frac_aromatic
"""
import argparse, os, sqlite3
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "gate_residues.sqlite")
ALL = "all classes"
MAXD = 8.0                                               # gate residues are stored out to 8 A
DMAX_CHOICES = (3.0, 4.0, 5.0, 6.0, 7.0, 8.0)
FOLD = {"5-fold": 5, "3-fold": 3}
GATES = {"5-fold": "the 5-fold axis", "3-fold": "the 3-fold axis"}
BY = {"capsid_class": "capsid type", "family": "family"}
MIN_FAMILY = 3                                           # families with fewer viruses are pooled
OTHER = "other / unclassified families"
NO_CLASS = "class not assigned"

# Kyte J, Doolittle RF (1982) J Mol Biol 157:105-132, hydropathy index
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
# Zamyatnin AA (1972) Prog Biophys Mol Biol 24:107-123, mean volume of the residue in solution (A^3)
VOLUME = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4, G=60.1, H=153.2, I=166.7,
              L=166.7, K=168.6, M=162.9, F=189.9, P=112.7, S=89.0, T=116.1, W=227.8, Y=193.6, V=140.0)
POSITIVE, NEGATIVE = set("KR"), set("DE")
AROMATIC, POLAR, CHARGED = set("FWY"), set("STNQ"), set("KRDE")

METRICS = {
    "net_charge": "net charge of the pore (sum over all copies)",
    "positives": "positive residues lining the pore (Lys + Arg, all copies)",
    "negatives": "negative residues lining the pore (Asp + Glu, all copies)",
    "n_lining": "residues lining the pore (all copies)",
    "n_residues": "residues near the axis in one asymmetric unit",
    "hydropathy": "mean hydropathy (Kyte-Doolittle; > 0 = water-repelling)",
    "frac_aromatic": "share of aromatic residues (F W Y), %",
    "frac_polar": "share of polar residues (S T N Q), %",
    "frac_charged": "share of charged residues (K R D E), %",
    "volume": "mean residue volume (A^3; larger = bulkier)",
}
SHORT = {"net_charge": "have a median net charge of", "positives": "have a median of", "negatives": "have a median of",
         "n_lining": "are lined by a median of", "n_residues": "have a median of",
         "hydropathy": "have a median Kyte-Doolittle hydropathy of", "frac_aromatic": "have a median aromatic share of",
         "frac_polar": "have a median polar share of", "frac_charged": "have a median charged share of",
         "volume": "have a median residue volume of"}
LABEL_SHORT = {"net_charge": "net charge", "positives": "number of positive residues",
               "negatives": "number of negative residues", "n_lining": "number of lining residues",
               "n_residues": "number of residues per asymmetric unit", "hydropathy": "hydropathy",
               "frac_aromatic": "aromatic share", "frac_polar": "polar share", "frac_charged": "charged share",
               "volume": "residue volume"}
PERCENT = {"frac_aromatic", "frac_polar", "frac_charged"}    # stored as 0-1, shown as %
COLUMNS = ["virus", "pdb", "family", "capsid_class", "gate", "fold", "n_residues", "n_lining", "net_charge",
           "positives", "negatives", "hydropathy", "frac_aromatic", "frac_polar", "frac_charged", "volume", "closest"]
COL_POS, COL_NEG, COL_ZERO, COL_DOT = "#1d5fa8", "#b2512f", "#7a8794", "#1f7a5a"


def connect(path=DB):
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    return con


def _viruses(con, cls):
    """{virus: row} using the best-resolution structure of each virus (same choice as gate_compare)"""
    rows = con.execute("SELECT virus, pdb, family, capsid_class, resolution_A FROM structures ORDER BY virus, "
                       "CASE WHEN resolution_A IS NULL THEN 1 ELSE 0 END, resolution_A, pdb").fetchall()
    pick = {}
    for r in rows:
        pick.setdefault(r["virus"], r)
    return {v: r for v, r in pick.items() if cls == ALL or r["capsid_class"] == cls}


def data(con, gate="5-fold", dmax=8.0, cls=ALL):
    """-> DataFrame, one row per virus, with the chemistry of its pore at `gate` within `dmax` A of the axis"""
    if gate not in FOLD:
        raise ValueError(f"gate must be one of {list(FOLD)}")
    dmax = float(min(dmax, MAXD))
    fold = FOLD[gate]
    vir = _viruses(con, cls)
    res = {r["pdb"]: [] for r in vir.values()}
    if res:
        q = ",".join("?" * len(res))
        for r in con.execute(f"SELECT pdb, chain, auth_seq, aa, distance_A FROM gate_residues WHERE pdb IN ({q}) "
                             f"AND gate = ? AND distance_A <= ? ORDER BY distance_A, chain, auth_seq",
                             (*res, gate, dmax)):
            if r["aa"] in KD:                            # standard amino acids only
                res[r["pdb"]].append(r)
    out = []
    for v, s in vir.items():
        rr = res[s["pdb"]]
        aa = [r["aa"] for r in rr]
        n = len(aa)
        pos = sum(a in POSITIVE for a in aa)
        neg = sum(a in NEGATIVE for a in aa)
        frac = lambda grp: sum(a in grp for a in aa) / n if n else np.nan
        c = rr[0] if rr else None
        out.append(dict(
            virus=v, pdb=s["pdb"], family=s["family"] or "Unclassified", capsid_class=s["capsid_class"] or NO_CLASS,
            gate=gate, fold=fold, n_residues=n, n_lining=fold * n, net_charge=fold * (pos - neg),
            positives=fold * pos, negatives=fold * neg,
            hydropathy=float(np.mean([KD[a] for a in aa])) if n else np.nan,
            frac_aromatic=frac(AROMATIC), frac_polar=frac(POLAR), frac_charged=frac(CHARGED),
            volume=float(np.mean([VOLUME[a] for a in aa])) if n else np.nan,
            closest=(f"{c['aa']}{c['auth_seq']} {c['distance_A']:.1f} A" if c else f"no residue within {dmax:g} A")))
    df = pd.DataFrame(out, columns=COLUMNS).sort_values(["capsid_class", "family", "virus"]).reset_index(drop=True)
    df.attrs.update(gate=gate, dmax=dmax, cls=cls)
    return df


def _groups(df, by, metric=None):
    """-> list of (group label, sub-DataFrame) for viruses with residues at the gate. Capsid types in name order;
    families from the highest to the lowest median of `metric`, pooled small families last"""
    use = df[df["n_residues"] > 0]
    if by == "family":
        size = use.groupby("family")["virus"].count()
        keep = {f for f, n in size.items() if n >= MIN_FAMILY and f != "Unclassified"}
        key = use["family"].where(use["family"].isin(keep), OTHER)
    else:
        key = use["capsid_class"]
    out = [(g, use[key == g]) for g in sorted(key.unique(), key=lambda g: (g == OTHER, g == NO_CLASS, g))]
    if by == "family" and metric:
        out.sort(key=lambda t: (t[0] == OTHER, -t[1][metric].median(), t[0]))
    return out


def _fmt(metric, x):
    """a metric value as written in the text and on the chart"""
    if x is None or not np.isfinite(x):
        return "-"
    if metric in PERCENT:
        return f"{100 * x:.0f}%"
    if metric == "hydropathy":
        return f"{x:+.2f}"
    if metric == "volume":
        return f"{x:.0f} A^3"
    if metric == "net_charge":
        return "0" if x == 0 else f"{x:+g}"
    return f"{x:g}"


def _kruskal(groups, metric):
    """Kruskal-Wallis across the named groups with >= 5 viruses; -> (p, [labels]) or (None, labels)"""
    from scipy.stats import kruskal
    test = [(g, d[metric].dropna().to_numpy()) for g, d in groups if g not in (OTHER, NO_CLASS)]
    test = [(g, x) for g, x in test if len(x) >= 5]
    if len(test) < 2:
        return None, [g for g, _ in test]
    allx = np.concatenate([x for _, x in test])
    if np.all(allx == allx[0]):                          # all values identical: no difference to test
        return 1.0, [g for g, _ in test]
    return float(kruskal(*[x for _, x in test]).pvalue), [g for g, _ in test]


def _p_words(p):
    if p < 0.001:
        return "p < 0.001, a clear difference"
    if p < 0.05:
        return f"p = {p:.3f}, a difference unlikely to be chance"
    return f"p = {p:.2f}, no clear difference"


def _one(g, d, metric, brief=False):
    x = d[metric]
    sign = (f"positive in {int((x > 0).sum())}, negative in {int((x < 0).sum())} of {len(x)} viruses"
            if metric == "net_charge" else f"{len(x)} viruses")
    if brief:
        return f"{g} ({_fmt(metric, x.median())}; {sign})"
    return f"{g} pores {SHORT[metric]} {_fmt(metric, x.median())} ({sign})"


def summary(df, by="capsid_class", metric="net_charge"):
    """the answer in plain words, 2-4 sentences"""
    gate, dmax = df.attrs.get("gate", "5-fold"), df.attrs.get("dmax", MAXD)
    groups = _groups(df, by, metric)
    empty = int((df["n_residues"] == 0).sum())
    head = f"At {GATES[gate]} (within {dmax:g} A), "
    if not groups:
        return head + f"none of the {len(df)} viruses has a residue near the axis."
    named = [(g, d) for g, d in groups if g not in (OTHER, NO_CLASS)] or groups
    if by == "family" and len(named) > 4:            # name only the two highest and the two lowest families
        out = [head + f"the families with the highest median {LABEL_SHORT[metric]} are "
               + " and ".join(_one(g, d, metric, True) for g, d in named[:2]) + ", the lowest "
               + " and ".join(_one(g, d, metric, True) for g, d in named[-2:][::-1]) + "."]
    else:
        out = [head + "; ".join(_one(g, d, metric) for g, d in named) + "."]
    p, tested = _kruskal(groups, metric)
    if p is not None:
        what = f"{len(tested)} families with at least 5 viruses" if by == "family" else ", ".join(tested)
        out.append(f"A Kruskal-Wallis test across {what} gives {_p_words(p)}"
                   + ("; viruses of one family share ancestry, so a family difference can be inherited rather "
                      "than selected." if by == "family" else "."))
    else:
        out.append("Too few groups with at least 5 viruses for a test across groups."
                   + (" Viruses of one family share ancestry." if by == "family" else ""))
    note = []
    if metric == "net_charge":
        note.append(f"charges are summed over all {FOLD[gate]} copies around the pore (Lys, Arg +1; Asp, Glu -1; "
                    f"His counted as 0 at pH 7)")
    if empty:
        note.append(f"{empty} of {len(df)} viruses have no residue within {dmax:g} A of this axis and are not included")
    if note:
        t = "; ".join(note)
        out.append(t[0].upper() + t[1:] + ".")
    return "  ".join(out)


def draw(fig, df, by="capsid_class", metric="net_charge"):
    """box plot with one dot per virus for the chosen metric, one row per group; returns a one-line description"""
    fig.clear()
    fig.set_layout_engine("constrained")
    ax = fig.add_subplot(111)
    gate, dmax = df.attrs.get("gate", "5-fold"), df.attrs.get("dmax", MAXD)
    groups = _groups(df, by, metric)
    scale = 100.0 if metric in PERCENT else 1.0
    if not groups:
        ax.axis("off")
        ax.text(0.5, 0.5, f"No virus has a residue within {dmax:g} A of {GATES[gate]}.", ha="center", va="center")
        return f"{GATES[gate]} | no residues within {dmax:g} A"
    groups = groups[::-1]                                # first group at the top
    vals = [d[metric].to_numpy(dtype=float) * scale for _, d in groups]
    y = np.arange(len(groups))
    import matplotlib
    horiz = (dict(orientation="horizontal") if tuple(int(v) for v in matplotlib.__version__.split(".")[:2]) >= (3, 10)
             else dict(vert=False))                  # older matplotlib (e.g. in the browser) lacks `orientation`
    ax.boxplot(vals, positions=y, widths=0.6, **horiz, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="#f2f2f2", edgecolor="#999"), whiskerprops=dict(color="#999"),
               capprops=dict(color="#999"), medianprops=dict(color="#222", lw=2.2), zorder=1)
    rng = np.random.default_rng(20261108)                # fixed jitter: the same data always gives the same picture
    small = len(groups) > 8
    for yi, x in zip(y, vals):
        jit = rng.uniform(-0.22, 0.22, len(x))
        if metric == "net_charge":
            col = np.where(x > 0, COL_POS, np.where(x < 0, COL_NEG, COL_ZERO))
        else:
            col = COL_DOT
        ax.scatter(x, yi + jit, s=16 if small else 26, c=col, alpha=0.75, edgecolors="white", linewidths=0.4, zorder=3)
    if metric in ("net_charge", "hydropathy"):
        ax.axvline(0, color="#555", lw=1, ls="--", zorder=0)
    lo = min(float(np.nanmin(v)) for v in vals)
    hi = max(float(np.nanmax(v)) for v in vals)
    pad = 0.05 * (hi - lo or 1)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(-0.6, len(groups) - 0.4)
    ax.set_yticks(y, [f"{g}  (n={len(d)})" for g, d in groups], fontsize=8 if small else 10)
    # median written to the right of each row, outside the axes so it never sits on a dot
    tr = ax.get_yaxis_transform()
    for yi, x in zip(y, vals):
        m = float(np.median(x))
        txt = _fmt(metric, m / scale)
        ax.text(1.01, yi, f"median {txt}", transform=tr, va="center", ha="left", fontsize=8 if small else 9,
                color="#222", clip_on=False)
    narrow = fig.get_figwidth() < 10                     # embedded small: break the explanation into lines
    sep = "\n" if narrow else ";  "
    key = ("blue = positive pore, red = negative, grey = neutral" + sep if metric == "net_charge" else "")
    ax.set_xlabel(METRICS[metric] + "\n" + key + "box = middle half of the viruses, black bar = median",
                  fontsize=9 if narrow else 10)
    ax.grid(axis="x", color="#eee")
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    cls = df.attrs.get("cls", ALL)
    ax.set_title(f"Pore chemistry at {GATES[gate]}, within {dmax:g} A, by {BY[by]}" + ("\n" if narrow else "  ")
                 + f"({cls}, {sum(len(v) for v in vals)} viruses, one dot each)", loc="left",
                 fontsize=10 if narrow else 11, weight="bold")
    p, _ = _kruskal(groups, metric)
    empty = int((df["n_residues"] == 0).sum())
    extra = [f"Kruskal-Wallis p = {p:.2g}" if p is not None else None,
             f"{empty} viruses with no residue within {dmax:g} A not shown" if empty else None]
    return (f"{GATES[gate]} | within {dmax:g} A | {METRICS[metric]} by {BY[by]} | "
            + " | ".join(e for e in extra if e)).rstrip(" |")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=DB)
    ap.add_argument("--gate", default="5-fold", choices=list(FOLD))
    ap.add_argument("--dmax", type=float, default=8.0, choices=DMAX_CHOICES, help="distance limit in A")
    ap.add_argument("--class", dest="cls", default=ALL)
    ap.add_argument("--by", default="capsid_class", choices=list(BY))
    ap.add_argument("--metric", default="net_charge", choices=list(METRICS))
    ap.add_argument("--out", help="PNG/PDF for the chart; the table is written next to it as .csv")
    a = ap.parse_args()
    con = connect(a.db)
    df = data(con, a.gate, a.dmax, a.cls)
    print(summary(df, a.by, a.metric))
    print()
    with pd.option_context("display.width", 200, "display.max_columns", 30, "display.precision", 2):
        print(df.drop(columns=["gate"]).head(15).to_string(index=False))
    if a.out:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(13, 7))
        print(draw(fig, df, a.by, a.metric))
        fig.savefig(a.out, dpi=100)
        df.to_csv(os.path.splitext(a.out)[0] + ".csv", index=False, float_format="%.4g")
        print("written:", a.out)


if __name__ == "__main__":
    main()
