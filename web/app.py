"""
Capsid Symmetry Tool - web version.

Runs the same read-only views as the desktop program over gate_residues.sqlite: Browse, Find, Plots, Compare,
Stem-and-leaf. Works both as a normal Streamlit app (streamlit run web/app.py) and fully in the browser on GitHub
Pages through stlite (Pyodide). Analysing NEW structures needs gemmi and stays in the desktop program.

Author: Siddhant Raj Naik, PhD scholar, Kusuma School of Biological Sciences (KSBS), Indian Institute of Technology Delhi
"""
import os, sys, sqlite3
import pandas as pd
import streamlit as st

HERE = os.path.dirname(os.path.abspath(__file__))
for p in (HERE, os.path.dirname(HERE)):                 # modules sit next to app.py online, one level up in the repo
    if p not in sys.path:
        sys.path.insert(0, p)
DB = next(p for p in (os.path.join(HERE, "gate_residues.sqlite"), os.path.join(os.path.dirname(HERE), "gate_residues.sqlite"))
          if os.path.exists(p))

import gate_plots, gate_compare, gate_stemleaf           # noqa: E402

CREDIT = ("**Siddhant Raj Naik**, PhD scholar, Kusuma School of Biological Sciences (KSBS), "
          "Indian Institute of Technology Delhi")
AA_ORDER = "AVLIMFWYCGPSTNQHKRDE"
AA_LABEL = {a: f"{a}  {gate_stemleaf.AA_NAME[a]}" for a in AA_ORDER}
GATE_COLOR = {"5-fold": "#1d5fa8", "3-fold": "#b2512f", "2-fold": "#5b7189"}

st.set_page_config(page_title="Capsid Symmetry Tool", page_icon="🦠", layout="wide")


def con():
    c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True, check_same_thread=False)
    c.row_factory = sqlite3.Row
    return c


@st.cache_data
def q(sql, args=()):
    c = con()
    df = pd.read_sql_query(sql, c, params=args)
    c.close()
    return df


def figure(w=13, h=8):
    from matplotlib.figure import Figure
    return Figure(figsize=(w, h), dpi=100)


# ---------------------------------------------------------------- header
s = q("SELECT (SELECT COUNT(DISTINCT virus) FROM structures) v, (SELECT COUNT(DISTINCT family) FROM structures) f, "
      "(SELECT COUNT(*) FROM gate_residues) g, (SELECT COUNT(*) FROM dual_positions) d").iloc[0]
st.title("Capsid Symmetry Tool")
st.caption("Residues guarding the 5-fold and 3-fold symmetry axes (\"gates\") of icosahedral virus capsids  ·  "
           f"{s.v} viruses · {s.f} families · {s.g:,} gate residues · {s.d} positions at both gates")

page = st.sidebar.radio("Page", ["Home", "Browse", "Find", "Plots", "Stem-and-leaf", "Compare", "Help"])
st.sidebar.markdown("---")
st.sidebar.markdown("Created by " + CREDIT)
st.sidebar.caption("Read-only web version. To analyse new PDB entries or your own models, use the desktop program "
                   "from the GitHub repository.")
classes = gate_plots.classes(con())


# ---------------------------------------------------------------- pages
if page == "Home":
    st.markdown(f"""
A virus shell (**capsid**) is built from many copies of the same protein arranged with icosahedral symmetry. It has
three kinds of symmetry axis: **5-fold** (6 axes), **3-fold** (10) and **2-fold** (15). Where an axis passes through the
shell there is a pore or a junction - a **gate**.

This tool measures, for every virus in the database, which amino acids sit at those gates, which residue numbers sit at
**both** the 5-fold and the 3-fold gate, and which amino acids are **favoured** there more than chance would give.

| Page | Question it answers |
|---|---|
| **Browse** | One virus: what sits at its gates? |
| **Find** | One amino acid: in which viruses is it at a gate? |
| **Plots** | Which amino acids are favoured at both gates, per capsid type? |
| **Stem-and-leaf** | Every virus on one chart, at the distance of its closest residue. |
| **Compare** | Two capsid types or families side by side. |

Created by {CREDIT}.
""")

elif page == "Browse":
    vir = q("""SELECT s.virus, s.pdb, s.family, s.capsid_class AS class, s.resolution_A AS resolution,
                 (SELECT COUNT(*) FROM gate_residues g WHERE g.pdb=s.pdb AND g.gate='5-fold') AS "5-fold",
                 (SELECT COUNT(*) FROM gate_residues g WHERE g.pdb=s.pdb AND g.gate='3-fold') AS "3-fold",
                 (SELECT COUNT(DISTINCT auth_seq) FROM dual_positions d WHERE d.pdb=s.pdb) AS "both gates"
               FROM structures s ORDER BY s.virus""")
    c1, c2 = st.columns([2, 1])
    txt = c1.text_input("Search virus, PDB ID or family")
    cls = c2.selectbox("Capsid class", ["all"] + classes)
    v = vir
    if txt:
        m = v.virus.str.contains(txt, case=False) | v.pdb.str.contains(txt, case=False) | \
            v.family.fillna("").str.contains(txt, case=False)
        v = v[m]
    if cls != "all":
        v = v[v["class"] == cls]
    st.dataframe(v.assign(pdb=v.pdb.str.upper()), hide_index=True, use_container_width=True, height=260)
    labels = [f"{r.virus}  ({r.pdb.upper()})" for r in v.itertuples()]
    if labels:
        pick = st.selectbox("Show the gates of", labels, index=next((i for i, l in enumerate(labels) if "4FTB" in l), 0))
        pdb = v.iloc[labels.index(pick)].pdb
        res = q("SELECT gate, chain, auth_seq, aa, resname, distance_A, dual_5_3 FROM gate_residues WHERE pdb=? "
                "ORDER BY gate, distance_A", (pdb,))
        dual = q("SELECT auth_seq, aa, chains_near_5fold, min_d5, chains_near_3fold, min_d3 FROM dual_positions "
                 "WHERE pdb=? ORDER BY auth_seq", (pdb,))
        st.subheader(pick)
        if len(dual):
            st.markdown("**Same residue number at both gates:** " + "; ".join(
                f"{r.aa}{r.auth_seq} (5-fold {r.min_d5:.1f} Å in {r.chains_near_5fold}; 3-fold {r.min_d3:.1f} Å in "
                f"{r.chains_near_3fold})" for r in dual.itertuples()))
        else:
            st.markdown("No residue number guards both gates within the cutoff.")
        fig = figure(13, 3.6)
        ax = fig.add_subplot(111)
        for i, g in enumerate(("5-fold", "3-fold", "2-fold")):
            r = res[res.gate == g]
            y = 2 - i
            ax.scatter(r.distance_A, [y] * len(r), s=45, zorder=3,
                       color=["#e0a030" if d else GATE_COLOR[g] for d in r.dual_5_3])
            last = []
            for j, row in enumerate(r.itertuples()):
                off = 0.22 if j % 2 == 0 else -0.3
                ax.annotate(f"{row.aa}{row.auth_seq}", (row.distance_A, y), xytext=(0, 14 if off > 0 else -16),
                            textcoords="offset points", ha="center", fontsize=8)
        ax.set_yticks([2, 1, 0], ["5-fold", "3-fold", "2-fold"]); ax.set_ylim(-0.7, 2.7); ax.set_xlim(-0.3, 8.3)
        ax.set_xlabel("distance of the side chain from the symmetry axis (Å)  -  orange = residue number at both gates")
        ax.grid(axis="x", color="#eee")
        fig.tight_layout()
        st.pyplot(fig)
        st.dataframe(res.rename(columns={"auth_seq": "residue", "distance_A": "distance Å", "dual_5_3": "both gates"}),
                     hide_index=True, use_container_width=True)

elif page == "Find":
    st.markdown("Pick an amino acid and a gate. Distance = how far the residue's side chain is from the symmetry axis: "
                "0-3 Å is right on the axis (in the pore), 8 Å is the outer edge of the gate.")
    c = st.columns(5)
    aa = c[0].selectbox("Amino acid", ["any"] + [AA_LABEL[a] for a in AA_ORDER], index=1 + AA_ORDER.index("N"))
    gate = c[1].selectbox("Gate", ["any", "5-fold", "3-fold", "2-fold"])
    cls = c[2].selectbox("Capsid class", ["any"] + classes)
    dmax = c[3].select_slider("Max distance (Å)", options=[1, 2, 3, 4, 5, 6, 7, 8], value=8)
    both = c[4].checkbox("only residues at BOTH gates")
    sql, args = "SELECT * FROM gate_residues WHERE distance_A <= ?", [float(dmax)]
    if aa != "any":
        sql += " AND aa = ?"; args.append(aa[:1])
    if gate != "any":
        sql += " AND gate = ?"; args.append(gate)
    if cls != "any":
        sql += " AND capsid_class = ?"; args.append(cls)
    if both:
        sql += " AND dual_5_3 = 1"
    r = q(sql + " ORDER BY distance_A", tuple(args))
    nv = q("SELECT COUNT(DISTINCT virus) n FROM structures" + (" WHERE capsid_class=?" if cls != "any" else ""),
           (cls,) if cls != "any" else ()).n[0]
    if len(r):
        fams = r.groupby(r.family.fillna("unknown")).virus.nunique().sort_values(ascending=False)
        nboth = r[r.dual_5_3 == 1].virus.nunique()
        st.info(f"Found in **{r.virus.nunique()} of {nv} viruses** ({100 * r.virus.nunique() / max(nv, 1):.0f}%), from "
                f"{len(fams)} families ({', '.join(f'{f} {n}' for f, n in fams.head(4).items())}"
                f"{', ...' if len(fams) > 4 else ''})." +
                ("" if both else f"  In {nboth} of these viruses the same residue number sits at BOTH gates."))
    else:
        st.warning("Nothing found. Try a larger distance or 'any' gate.")
    st.dataframe(r[["virus", "pdb", "family", "capsid_class", "gate", "resname", "auth_seq", "chain", "distance_A", "dual_5_3"]]
                 .rename(columns={"capsid_class": "capsid type", "auth_seq": "residue", "distance_A": "distance Å",
                                  "dual_5_3": "both gates"}), hide_index=True, use_container_width=True, height=520)
    st.download_button("Download table (CSV)", r.to_csv(index=False), "gate_search.csv", "text/csv")

elif page == "Plots":
    c = st.columns(3)
    mode = c[0].selectbox("Show", ["Simple answer", "Full figure (4 panels)"])
    cls = c[1].selectbox("Capsid class", [gate_plots.ALL] + classes, index=1 + classes.index("T=3") if "T=3" in classes else 0)
    cut = c[2].selectbox("Cutoff (Å)", [3.0, 6.0, 8.0, 10.0], index=2)
    fig = figure(16, 10 if mode.startswith("Simple") else 11)
    with st.spinner("drawing ..."):
        desc = (gate_plots.build_simple if mode.startswith("Simple") else gate_plots.build)(fig, con(), cls, cut)
    st.caption(desc)
    st.pyplot(fig)
    st.markdown("Top rows (simple answer) = favoured at both gates. Blue = 5-fold gate, orange = 3-fold gate, "
                "\\* = significant (FDR < 0.05). Odds ratios compare each amino acid at the gate with the rest of the "
                "capsid shell at the same radius, pooled over viruses (Cochran-Mantel-Haenszel).")

elif page == "Stem-and-leaf":
    c = st.columns(3)
    cls = c[0].selectbox("Capsid class", [gate_stemleaf.ALL] + classes)
    aa = c[1].selectbox("Measure", ["closest residue"] + [AA_LABEL[a] for a in AA_ORDER])
    letters = c[2].radio("Leaf shows", ["tenths digit", "amino acid letter"], horizontal=True) == "amino acid letter"
    res = gate_stemleaf.data(gate_stemleaf.connect(DB), cls, None if aa == "closest residue" else aa[:1])
    fig = figure(18, 9)
    st.caption(gate_stemleaf.draw(fig, res, letters))
    st.pyplot(fig)
    with st.expander("Text version"):
        st.code(gate_stemleaf.text(res))

elif page == "Compare":
    st.markdown("Pick two capsid types or families. For every virus the residues at the chosen gate are taken and each "
                "amino acid's share of them is computed; the groups' average shares are compared (permutation test, "
                "FDR over 20 amino acids). Yellow rows with * are differences unlikely to be chance.")
    cc = gate_compare.connect(DB)
    ch = gate_compare.group_choices(cc)
    labs = [l for l, _ in ch]; spec = dict(ch)
    rest = "all other viruses (everything not in the first group)"
    c = st.columns(4)
    a = c[0].selectbox("First group (blue)", labs, index=next((i for i, (l, s_) in enumerate(ch) if s_ == "class=T=3"), 0))
    b = c[1].selectbox("Second group (orange)", labs + [rest],
                       index=next((i for i, (l, s_) in enumerate(ch) if s_ == "class=pseudo T=3 (P=3)"), 1))
    g = c[2].selectbox("Gate", ["5-fold gate", "3-fold gate", "BOTH gates (dual positions)"])
    d = c[3].selectbox("Max distance (Å)", [8.0, 6.0, 5.0, 3.0])
    gk = {"5-fold gate": "5-fold", "3-fold gate": "3-fold"}.get(g, "both")
    with st.spinner("comparing ..."):
        res = gate_compare.compare(cc, spec[a], gate_compare.REST if b == rest else spec[b], gk, d)
    fig = figure(12, 8.5)
    st.info(gate_compare.draw(fig, res))
    st.pyplot(fig)
    st.dataframe(pd.DataFrame([dict(amino_acid=f"{r['aa']}  {r['name']}", first_pct=round(100 * r["sa"], 1),
                                    second_pct=round(100 * r["sb"], 1), first_has_it=f"{r['ka']}/{res['na']}",
                                    second_has_it=f"{r['kb']}/{res['nb']}", FDR=round(r["fdr"], 3)) for r in res["rows"]]),
                 hide_index=True, use_container_width=True)

else:
    st.markdown(f"""
### What is measured
For every subunit of the capsid the program finds the true icosahedral 5-fold, 3-fold and 2-fold axes from the
structure's own symmetry operators, then measures the shortest distance from each residue's side chain to the nearest
axis of each kind. A residue within the cutoff (8 Å by default) is a **gate residue**. A residue number that is near a
5-fold axis in one subunit and near a 3-fold axis in another is a **dual position**: the same position in the sequence
guards both gates, as Met182 does in Flock House virus (PDB 4FTB).

### Statistics
Each amino acid is compared with capsid residues of the same structure that are not near any axis and lie at the same
radius in the shell. Structures are combined with a Cochran-Mantel-Haenszel odds ratio, so one virus with an unusual
composition cannot drive the result, and p-values are corrected for testing 20 amino acids (Benjamini-Hochberg FDR).

### Limits
Distance to an axis is not pore radius. Only atomic models are read, not density maps. Quasi-3-fold axes of T=3
capsids are not analysed. Viruses of one family share ancestry, so a family pattern can be inherited rather than
selected. A conserved gate residue says where to look, not what it does.

### Credit
Created by {CREDIT}.
""")
