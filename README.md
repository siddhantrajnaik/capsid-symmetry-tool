<p align="center"><img src="gate_icon.png" width="110" alt="Capsid Symmetry Tool icon"></p>

<h1 align="center">Capsid Symmetry Tool</h1>

<p align="center">
Which amino acids guard the <b>5-fold</b> and <b>3-fold</b> symmetry axes ("gates") of icosahedral virus capsids,
which residue positions sit at <b>both</b> gates, and which amino acids are <b>favoured</b> there.
</p>

<p align="center">
<a href="https://siddhantrajnaik.github.io/capsid-symmetry-tool/"><b>Open the web version</b></a> ·
<a href="#desktop-program">Desktop program</a> ·
<a href="#how-it-works">How it works</a> ·
<a href="#citation">Cite</a>
</p>

<p align="center">
<a href="https://siddhantrajnaik.github.io/capsid-symmetry-tool/"><img src="https://img.shields.io/badge/web%20version-online-1d5fa8.svg" alt="web version online"></a>
<a href="LICENSE"><img src="https://img.shields.io/badge/licence-MIT-blue.svg" alt="MIT licence"></a>
</p>

---

**Author:** Siddhant Raj Naik, PhD scholar, Kusuma School of Biological Sciences (KSBS), Indian Institute of Technology Delhi

---

## What it does

A virus shell (capsid) is built from many copies of the same protein arranged with icosahedral symmetry. It has
5-fold (6), 3-fold (10) and 2-fold (15) symmetry axes; where an axis passes through the shell there is a pore or a
junction - a **gate**. For every capsid structure the tool

1. reads the PDB/mmCIF entry and the icosahedral operators of its biological assembly,
2. rebuilds the icosahedral group and finds the capsid centre and all 5-, 3- and 2-fold axes from the operators,
3. measures the shortest distance from every residue's side chain to the nearest axis of each kind,
4. calls residues within the cutoff (8 Å by default; 3, 6 and 10 Å also stored) **gate residues**, and residue
   numbers near a 5-fold axis in one subunit and a 3-fold axis in another **dual positions** - for example Met182
   of Flock House virus (PDB 4FTB: 2.8 Å from the 5-fold axis, 2.5 Å from the 3-fold axis),
5. tests each amino acid for enrichment at the gates against capsid residues at the same radius that are not near
   any axis, pooled over viruses (Cochran-Mantel-Haenszel odds ratio, Benjamini-Hochberg FDR).

The bundled database `gate_residues.sqlite` holds **181 structures of 175 viruses from 23 families** (T=3, T=4 and
pseudo-T=3 capsids), 6,724 gate residues and 170 dual positions.

## Web version

**<https://siddhantrajnaik.github.io/capsid-symmetry-tool/>** - runs entirely in your browser (Python via
Pyodide/stlite), nothing to install. The first visit takes about a minute to load.

| Page | Question it answers |
|---|---|
| Browse | One virus: what sits at its gates? |
| Find | One amino acid: in which viruses is it at a gate? |
| Plots | Which amino acids are favoured at both gates, per capsid type (3, 6, 8, 10 Å)? |
| Stem-and-leaf | Every virus on one chart, at the distance of its closest residue. |
| Compare | Two capsid types or families side by side. |

The web version reads the database only. Analysing new PDB entries or your own models needs the desktop program.

![Web version, Plots page](docs/web_plots.jpg)

## Desktop program

**Windows:** download the repository (green *Code* button → *Download ZIP*, or `git clone`), then double-click
**`Capsid Symmetry Tool.bat`**. The first start creates a private Python environment in `.venv` (needs Python
3.11-3.13 from python.org and internet, a few minutes); later starts are instant.

**Any system:**

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt          # Windows: .venv\Scripts\pip install -r requirements.txt
.venv/bin/python gate_explorer.py                  # the window
```

The window has the same pages as the web version plus **Analyse**, which runs the analysis on new PDB IDs or local
model files and adds them to the database.

![Desktop program, Compare page](docs/desktop_compare.png)

### Command line

```bash
python gate_residue_analysis.py --pdb 4ftb --residue 182 --out results/fhv          # one structure
python gate_residue_analysis.py --csv representatives.csv --sens 3 6 10 --out results/all
python gate_plots.py --simple --class "T=3" --cutoff 3 --out plots/t3_3A.png       # favoured at both gates
python gate_compare.py --a "family=Nodaviridae" --b "family=Caliciviridae" --gate 5-fold
python gate_stemleaf.py --class "T=3" --aa M                                        # stem-and-leaf in the console
python gate_db.py add --db gate_residues.sqlite --run results/fhv                   # add a run to the database
```

## How it works

| File | Role |
|---|---|
| `gate_residue_analysis.py` | symmetry detection, axis distances, gate/dual residues, enrichment statistics, figures |
| `gate_db.py` | builds and updates `gate_residues.sqlite` from analysis output folders |
| `gate_explorer.py` | desktop window (tkinter): Browse, Find, Preferences, Plots, Compare, Analyse, Help |
| `gate_plots.py`, `gate_compare.py`, `gate_stemleaf.py` | figures and comparisons shared by the window, the web version and the command line |
| `gate_env.py` | picks the interpreter and checks every native library before a run starts |
| `web/` | the browser version (`app.py` for Streamlit, `index.html` for GitHub Pages) |
| `representatives.csv` | one row per virus with candidate PDB IDs |

**Comparisons between groups** (Compare page) use each virus's share of its gate residues per amino acid, so
capsid types with larger gates do not look richer in every amino acid, and a permutation test on the difference of
group means.

**Limits.** Distance to an axis is not pore radius. Only atomic models are read, not density maps. Quasi-3-fold
axes of T=3 capsids are not analysed. Viruses of one family share ancestry, so a family pattern can be inherited
rather than selected. A conserved gate residue tells you where to look, not what it does.

Structures come from the [RCSB Protein Data Bank](https://www.rcsb.org/); please cite the original depositions of
any structure you discuss.

## Citation

If you use this tool, please cite it (GitHub shows a *Cite this repository* button from `CITATION.cff`):

> Naik, S. R. *Capsid Symmetry Tool: residues at the 5-fold and 3-fold axes of icosahedral virus capsids.*
> Kusuma School of Biological Sciences, Indian Institute of Technology Delhi (2026).
> https://github.com/siddhantrajnaik/capsid-symmetry-tool

## Licence

MIT © 2026 Siddhant Raj Naik - see [LICENSE](LICENSE).
