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

EMAIL = "siddhantrajnaik@gmail.com"
CREDIT = "Siddhant Raj Naik, PhD scholar, Kusuma School of Biological Sciences, IIT Delhi"
AA_ORDER = "AVLIMFWYCGPSTNQHKRDE"
AA_LABEL = {a: f"{a}  {gate_stemleaf.AA_NAME[a]}" for a in AA_ORDER}
GATE_COLOR = {"5-fold": "#1d5fa8", "3-fold": "#b2512f", "2-fold": "#5b7189"}

st.set_page_config(page_title="Capsid Symmetry Tool", layout="wide")


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
st.caption(f"{s.v} viruses, {s.f} families, {s.g:,} residues near a symmetry axis")

page = st.sidebar.radio("Page", ["Home", "Browse", "Find", "Plots", "Stem-and-leaf", "Compare", "Help"])
st.sidebar.markdown("---")
st.sidebar.markdown(f"Made by **Siddhant Raj Naik**  \nPhD scholar, Kusuma School of Biological Sciences, IIT Delhi  \n"
                    f"[{EMAIL}](mailto:{EMAIL})")
st.sidebar.caption("You can look things up here. To run new structures, download the desktop version from "
                   "[GitHub](https://github.com/siddhantrajnaik/capsid-symmetry-tool).")
classes = gate_plots.classes(con())


# ---------------------------------------------------------------- pages
if page == "Home":
    st.markdown(f"""
I'm Siddhant, a PhD scholar at the Kusuma School of Biological Sciences, IIT Delhi.

Icosahedral virus capsids have openings where their symmetry axes pass through the shell: five subunits meet at
each 5-fold axis and three at each 3-fold axis. I wanted to know which amino acids line these openings, and whether
the same residue ever sits at both. In Flock House virus it does. Met182 is 2.8 Å from the 5-fold axis in one
subunit and 2.5 Å from the 3-fold axis in the other two. I wasn't sure if that was a one-off, so I checked
{s.v} viruses.

For each structure from the PDB I worked out where the symmetry axes are, using the structure's own assembly
operators, and measured how far every side chain is from the nearest one. Anything within 8 Å counts as sitting
at that axis. The statistics compare those residues with the rest of the capsid shell at the same radius.

The pages on the left:

- **Browse**: pick a virus and see which residues sit near each axis.
- **Find**: pick an amino acid and see which viruses have it there.
- **Plots**: which amino acids turn up near the axes more than you'd expect.
- **Stem-and-leaf**: every virus on one chart, by how close its nearest residue gets.
- **Compare**: two capsid types or families side by side.

It all runs inside your browser, which is why the first load is slow.

If something looks wrong, or you'd like a structure added, write to me at [{EMAIL}](mailto:{EMAIL}).
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
    st.markdown("Pick an amino acid and an axis. The distance is from the side chain to the axis, so anything "
                "under about 3 Å is practically on it.")
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
    st.markdown("Blue bars are the 5-fold axis and orange bars the 3-fold axis. Bars to the right mean the amino "
                "acid shows up there more often than in the rest of the shell, and a star means the difference holds "
                "after correcting for testing all 20 amino acids. The rows at the top are the ones that are high at both.")

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
    st.markdown("Pick two groups. For each virus I take the residues near the chosen axis and work out what fraction "
                "each amino acid makes up, then compare the two groups' averages. Yellow rows are the differences "
                "that survive a permutation test.")
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
### How the numbers are made

Each PDB entry comes with the operators that build the full capsid. From those I rebuild the icosahedral group, find
the centre of the capsid and the 6 five-fold, 10 three-fold and 15 two-fold axes, and measure the shortest distance
from each residue's side-chain atoms to each kind of axis. A residue within 8 Å of an axis is counted as being at
that axis. If the same residue number is near a 5-fold axis in one subunit and near a 3-fold axis in another, I
call it a position at both axes.

To decide whether an amino acid is enriched, I compare it with residues of the same structure that are at the same
radius in the shell but not near any axis. The viruses are pooled with a Cochran-Mantel-Haenszel test, so one virus
with an odd composition can't carry the result. The p-values are corrected for testing 20 amino acids.

### Things to keep in mind

Being close to an axis doesn't mean a residue blocks the pore. Distance from the axis line isn't the same thing as
pore radius. I only use atomic models, not density maps, and I don't look at the quasi-3-fold axes of T=3 capsids.
Viruses in the same family are related, so a pattern inside one family could just be inherited.

### Contact

{CREDIT}. Email [{EMAIL}](mailto:{EMAIL}). The code is on
[GitHub](https://github.com/siddhantrajnaik/capsid-symmetry-tool).
""")
