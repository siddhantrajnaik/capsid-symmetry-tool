"""
Capsid Symmetry Tool -- desktop application for the 5-fold / 3-fold gate residues of icosahedral capsids.

Four views over one SQLite database (gate_residues.sqlite):
  Browse      every virus in the database, its gate residues and the positions that guard both gates
  Find        search gate residues across all viruses by amino acid, gate, distance or class
  Preferences pooled amino-acid enrichment per capsid class
  Analyse     run gate_residue_analysis.py on new PDB entries or your own model files and add them to the database
"""
import os, sys, sqlite3, csv, subprocess, threading, queue, datetime, glob, webbrowser, traceback
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "gate_residues.sqlite")
SCRIPT = os.path.join(HERE, "gate_residue_analysis.py")
RESULTS = os.path.join(HERE, "results")
import gate_env
PY = gate_env.find_python()                                 # .venv\Scripts\python.exe when present
NOWIN = 0x08000000 if os.name == "nt" else 0

C = dict(bg="#eef2f7", card="#ffffff", ink="#1b2a3a", muted="#5b7189", line="#d3dde8",
         accent="#1d5fa8", accent_soft="#e6eefa", warm="#b2512f", good="#1f7a5a", rowalt="#f7fAfc")
FONT = ("Segoe UI", 10)
FONT_B = ("Segoe UI Semibold", 10)
FONT_H = ("Segoe UI Semibold", 15)
FONT_S = ("Segoe UI", 9)
MONO = ("Consolas", 10)
GATE_COLOR = {"5-fold": "#1d5fa8", "3-fold": "#b2512f", "2-fold": "#5b7189"}
WEB_URL = "https://siddhantrajnaik.github.io/capsid-symmetry-tool/"
AA_ORDER = "AVLIMFWYCGPSTNQHKRDE"
AA_LABEL = {a: f"{a}  {n}" for a, n in zip(AA_ORDER, (
    "Alanine", "Valine", "Leucine", "Isoleucine", "Methionine", "Phenylalanine", "Tryptophan", "Tyrosine", "Cysteine",
    "Glycine", "Proline", "Serine", "Threonine", "Asparagine", "Glutamine", "Histidine", "Lysine", "Arginine",
    "Aspartate", "Glutamate"))}


# --------------------------------------------------------------------------- database access
class Store:
    def __init__(self, path):
        self.path = path
        self.con = None
        self.open()

    def open(self):
        try:
            if self.con:
                self.con.close()
        except Exception:
            pass
        self.con = sqlite3.connect(self.path, check_same_thread=False)
        self.con.row_factory = sqlite3.Row
        self.con.executescript(SCHEMA_MIN)

    def q(self, sql, args=()):
        try:
            return self.con.execute(sql, args).fetchall()
        except sqlite3.Error:
            return []

    def summary(self):
        one = lambda s: (self.q(s) or [{0: 0}])[0][0]
        return dict(structures=one("SELECT COUNT(*) FROM structures"),
                    viruses=one("SELECT COUNT(DISTINCT virus) FROM structures"),
                    families=one("SELECT COUNT(DISTINCT family) FROM structures"),
                    residues=one("SELECT COUNT(*) FROM gate_residues"),
                    dual=one("SELECT COUNT(DISTINCT pdb || auth_seq) FROM dual_positions"))

    def classes(self):
        return [r[0] for r in self.q("SELECT DISTINCT capsid_class FROM structures WHERE capsid_class IS NOT NULL ORDER BY 1")]


SCHEMA_MIN = """
CREATE TABLE IF NOT EXISTS structures (pdb TEXT PRIMARY KEY, virus TEXT, family TEXT, capsid_class TEXT,
    resolution_A REAL, n_chains INTEGER, n_residues INTEGER, symmetry_from TEXT, axis_fit_rms_deg REAL,
    source TEXT, added TEXT);
CREATE TABLE IF NOT EXISTS gate_residues (pdb TEXT, virus TEXT, family TEXT, capsid_class TEXT, gate TEXT,
    chain TEXT, auth_seq INTEGER, aa TEXT, resname TEXT, distance_A REAL, radius_A REAL, dual_5_3 INTEGER,
    cutoff_A REAL, atoms TEXT);
CREATE TABLE IF NOT EXISTS dual_positions (pdb TEXT, virus TEXT, family TEXT, capsid_class TEXT, auth_seq INTEGER,
    aa TEXT, resname TEXT, chains_near_5fold TEXT, min_d5 REAL, chains_near_3fold TEXT, min_d3 REAL, cutoff_A REAL);
CREATE TABLE IF NOT EXISTS preferences (capsid_class TEXT, cutoff_A REAL, zone TEXT, level TEXT, key TEXT,
    n_structures INTEGER, n_with_key INTEGER, n_total INTEGER, freq_gate REAL, freq_background REAL,
    odds_ratio REAL, fdr REAL, families_enriched INTEGER, families_depleted INTEGER);
CREATE TABLE IF NOT EXISTS meta (k TEXT PRIMARY KEY, v TEXT);
"""


# --------------------------------------------------------------------------- small widgets
def card(parent, title=None):
    outer = tk.Frame(parent, bg=C["line"], padx=1, pady=1)
    inner = tk.Frame(outer, bg=C["card"])
    inner.pack(fill="both", expand=True)
    if title:
        tk.Label(inner, text=title, font=FONT_B, bg=C["card"], fg=C["ink"]).pack(anchor="w", padx=12, pady=(10, 2))
    return outer, inner


def make_tree(parent, columns, widths, anchors=None, height=12):
    wrap = tk.Frame(parent, bg=C["card"])
    tv = ttk.Treeview(wrap, columns=columns, show="headings", height=height, style="Gate.Treeview")
    k = min(parent.winfo_fpixels("1i") / 96.0, 1.25)       # column widths are pixels, fonts scale with Windows display scaling
    for col, w in zip(columns, (int(w * k) for w in widths)):
        tv.heading(col, text=col, anchor="w")
        tv.column(col, width=w, anchor=(anchors or {}).get(col, "w"), stretch=(w > 120 * k))
    vs = ttk.Scrollbar(wrap, orient="vertical", command=tv.yview)
    tv.configure(yscrollcommand=vs.set)
    tv.pack(side="left", fill="both", expand=True)
    vs.pack(side="right", fill="y")
    tv.tag_configure("odd", background=C["rowalt"])
    tv.tag_configure("dual", foreground=C["warm"], font=FONT_B)
    return wrap, tv


def fill_tree(tv, rows):
    tv.delete(*tv.get_children())
    for i, r in enumerate(rows):
        tags = ["odd"] if i % 2 else []
        if isinstance(r, (list, tuple)) and r and r[-1] == "DUAL":
            r = r[:-1]; tags.append("dual")
        tv.insert("", "end", values=r, tags=tags)


def export_tree(tv, name):
    cols = tv["columns"]
    rows = [tv.item(i, "values") for i in tv.get_children()]
    if not rows:
        messagebox.showinfo("Nothing to export", "The table is empty."); return
    p = filedialog.asksaveasfilename(defaultextension=".csv", initialfile=name, filetypes=[("CSV", "*.csv")])
    if not p:
        return
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(cols); w.writerows(rows)
    messagebox.showinfo("Saved", f"{len(rows)} rows written to\n{p}")


class DistanceChart(tk.Canvas):
    """residues drawn on a 0 -> cutoff distance scale, one lane per gate"""
    def __init__(self, parent, height=190):
        super().__init__(parent, height=height, bg=C["card"], highlightthickness=0)
        self.rows = []
        self.bind("<Configure>", lambda e: self.redraw())

    def show(self, rows, cutoff):
        self.rows, self.cutoff = rows, max(cutoff or 8, 1)
        self.redraw()

    def redraw(self):
        self.delete("all")
        W = self.winfo_width() or 600
        left, right, top, lane = 74, 24, 24, 80
        axis_y = top + 3 * lane - 6                         # tick labels sit here, legend one row below
        self.configure(height=axis_y + 40)
        if not self.rows:
            self.create_text(W / 2, 60, text="select a virus to see its gate residues",
                             fill=C["muted"], font=FONT_S)
            return
        span = max(W - left - right, 80)
        x = lambda d: left + span * min(d / self.cutoff, 1.0)
        for t in (0, self.cutoff / 2, self.cutoff):
            self.create_line(x(t), top - 6, x(t), axis_y - 10, fill=C["line"])
            self.create_text(x(t), axis_y, text=f"{t:g}", fill=C["muted"], font=FONT_S)
        self.create_text(left + span / 2, 10, text="distance from the symmetry axis (A)", fill=C["muted"], font=FONT_S)
        tiers = (-13, 13, -25, 25, -37, 37)                         # label rows above / below the lane line
        for i, gate in enumerate(("5-fold", "3-fold", "2-fold")):
            y = top + 26 + i * lane
            self.create_text(8, y, text=gate, anchor="w", fill=GATE_COLOR[gate], font=FONT_B)
            sel = sorted([r for r in self.rows if r["gate"] == gate], key=lambda r: r["distance_A"])
            if not sel:
                self.create_text(left, y, text="no residue within the cutoff", anchor="w", fill=C["muted"], font=FONT_S)
                continue
            self.create_line(left, y, left + span, y, fill=C["line"])
            last = [-1e9] * len(tiers)                      # right edge of the last label in each tier
            for r in sel:
                cx = x(r["distance_A"])
                dual = bool(r["dual_5_3"])
                self.create_oval(cx - 5, y - 5, cx + 5, y + 5, fill=(C["warm"] if dual else GATE_COLOR[gate]),
                                 outline=C["card"], width=1)
                lab = f'{r["aa"]}{r["auth_seq"]}'
                half = 4.3 * len(lab) + 4
                k = next((j for j in range(len(tiers)) if cx - half > last[j]), None)
                if k is None:                               # all rows busy: take the one that frees up first
                    k = min(range(len(tiers)), key=lambda j: last[j])
                last[k] = cx + half
                self.create_text(cx, y + tiers[k], text=lab, fill=C["ink"], font=FONT_S)
        self.create_text(8, axis_y + 22, anchor="w", font=FONT_S, fill=C["warm"],
                         text="orange = same residue at 5-fold and 3-fold gates")


class BarChart(tk.Canvas):
    """log2 odds ratio per amino acid, filled bars = significant"""
    def __init__(self, parent, height=260):
        super().__init__(parent, height=height, bg=C["card"], highlightthickness=0)
        self.rows = []
        self.bind("<Configure>", lambda e: self.redraw())

    def show(self, rows):
        self.rows = rows
        self.redraw()

    def redraw(self):
        import math
        self.delete("all")
        W, H = self.winfo_width() or 700, int(self["height"])
        if not self.rows:
            self.create_text(W / 2, H / 2, text="no pooled statistics for this selection", fill=C["muted"], font=FONT_S)
            return
        vals = []
        for r in self.rows:
            o = r["odds_ratio"]
            if o is None or o <= 0:
                vals.append((r["key"], None, r["fdr"]))
            else:
                vals.append((r["key"], math.log2(o), r["fdr"]))
        m = max([abs(v) for _, v, _ in vals if v is not None] + [1])
        pad, top, bot = 30, 24, 34
        mid = pad + (W - 2 * pad) / 2
        half = (W - 2 * pad) / 2
        n = len(vals)
        bh = max(min((H - top - bot) / max(n, 1), 20), 8)
        self.create_line(mid, top - 6, mid, top + n * bh + 4, fill=C["muted"])
        self.create_text(mid, 10, text="enriched at the gate  ->", anchor="w", fill=C["muted"], font=FONT_S)
        self.create_text(mid, 10, text="<-  depleted", anchor="e", fill=C["muted"], font=FONT_S)
        for i, (k, v, fdr) in enumerate(vals):
            y = top + i * bh + bh / 2
            self.create_text(mid - half - 6, y, text=k, anchor="e", fill=C["ink"], font=FONT_S)
            if v is None:
                self.create_text(mid + 6, y, text="not estimable", anchor="w", fill=C["muted"], font=FONT_S)
                continue
            w = half * min(abs(v) / m, 1.0)
            sig = fdr is not None and fdr < 0.05
            col = C["accent"] if v > 0 else C["muted"]
            x0, x1 = (mid, mid + w) if v > 0 else (mid - w, mid)
            if sig:
                self.create_rectangle(x0, y - bh * 0.32, x1, y + bh * 0.32, fill=col, outline="")
            else:
                self.create_rectangle(x0, y - bh * 0.32, x1, y + bh * 0.32, fill=C["card"], outline=col)
            self.create_text(x1 + 4 if v > 0 else x0 - 4, y, anchor="w" if v > 0 else "e",
                             text=f'{2 ** v:.1f}x', fill=C["muted"], font=FONT_S)
        self.create_text(W / 2, H - 12, font=FONT_S, fill=C["muted"],
                         text="bar length = pooled odds ratio vs non-gate residues; filled = false-discovery rate < 0.05")


# --------------------------------------------------------------------------- application
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Capsid Symmetry Tool")
        self.geometry("1180x780")
        self.minsize(1000, 680)
        try:
            self.state("zoomed")                            # open maximised: the tables need the height
        except tk.TclError:
            pass
        self.configure(bg=C["bg"])
        try:
            self.iconbitmap(os.path.join(HERE, "gate_icon.ico"))
        except tk.TclError:
            pass
        self.store = Store(DB)
        self.q = queue.Queue()
        self.proc = None
        self.last_out = None
        self._style()
        self._header()
        body = tk.Frame(self, bg=C["bg"]); body.pack(fill="both", expand=True)
        self._nav(body)
        self.stack = tk.Frame(body, bg=C["bg"]); self.stack.pack(side="left", fill="both", expand=True, padx=(0, 14), pady=10)
        self.views = {}
        self.c_fig = None
        self.k_fig = None
        for name, builder in (("Browse", self.view_browse), ("Find", self.view_find),
                              ("Preferences", self.view_pref), ("Plots", self.view_plots), ("Pores", self.view_pores), ("Compare", self.view_compare), ("Analyse", self.view_run), ("Help", self.view_help)):
            fr = tk.Frame(self.stack, bg=C["bg"]); self.views[name] = fr
            builder(fr)
        self.status = tk.StringVar(value="")
        tk.Label(self, textvariable=self.status, bg=C["bg"], fg=C["muted"], font=FONT_S, anchor="w").pack(fill="x", padx=18, pady=(0, 8))
        self.show("Browse")
        self.refresh_all()
        self.protocol("WM_DELETE_WINDOW", self.close)
        self.after(150, self.pump)

    # ---------- chrome
    def _style(self):
        s = ttk.Style(self)
        try:
            s.theme_use("clam")
        except tk.TclError:
            pass
        s.configure("Gate.Treeview", background=C["card"], fieldbackground=C["card"], foreground=C["ink"],
                    rowheight=int(25 * self.winfo_fpixels("1i") / 96.0), font=FONT, borderwidth=0)
        s.configure("Gate.Treeview.Heading", background=C["accent_soft"], foreground=C["ink"], font=FONT_B,
                    relief="flat", padding=(6, 5))
        s.map("Gate.Treeview.Heading", background=[("active", C["accent_soft"])])
        s.map("Gate.Treeview", background=[("selected", C["accent"])], foreground=[("selected", "#ffffff")])
        s.configure("TCombobox", fieldbackground=C["card"], background=C["card"])
        s.configure("Gate.TCheckbutton", background=C["card"], font=FONT)
        s.configure("Vertical.TScrollbar", background=C["bg"], troughcolor=C["bg"], borderwidth=0, arrowsize=12)

    def _header(self):
        h = tk.Frame(self, bg=C["accent"], height=66); h.pack(fill="x"); h.pack_propagate(False)
        tk.Label(h, text="Capsid Symmetry Tool", font=("Segoe UI Semibold", 17), bg=C["accent"], fg="white").pack(side="left", padx=18)
        tk.Label(h, text="residues guarding the 5-fold and 3-fold axes of icosahedral capsids",
                 font=FONT_S, bg=C["accent"], fg="#cfe0f5").pack(side="left", pady=(6, 0))
        self.chips = tk.Frame(h, bg=C["accent"]); self.chips.pack(side="right", padx=16)

    def _nav(self, body):
        nav = tk.Frame(body, bg=C["bg"], width=176); nav.pack(side="left", fill="y", padx=(14, 14), pady=10); nav.pack_propagate(False)
        self.nav_btns = {}
        items = [("Browse", "Viruses and their gates"), ("Find", "Search gate residues"),
                 ("Preferences", "Amino-acid statistics"),
                 ("Plots", "Residues at both gates"), ("Pores", "Charge and chemistry of each opening"), ("Compare", "Capsid types or families"), ("Analyse", "Add new structures"), ("Help", "What this measures")]
        for name, sub in items:
            b = tk.Frame(nav, bg=C["bg"], cursor="hand2")
            b.pack(fill="x", pady=2)
            lab = tk.Label(b, text=name, font=FONT_B, bg=C["bg"], fg=C["ink"], anchor="w")
            lab.pack(fill="x", padx=12, pady=(7, 0))
            sl = tk.Label(b, text=sub, font=FONT_S, bg=C["bg"], fg=C["muted"], anchor="w", wraplength=150, justify="left")
            sl.pack(fill="x", padx=12, pady=(0, 7))
            for w in (b, lab, sl):
                w.bind("<Button-1>", lambda e, n=name: self.show(n))
            self.nav_btns[name] = (b, lab, sl)

    def show(self, name):
        for n, fr in self.views.items():
            fr.pack_forget()
        self.views[name].pack(fill="both", expand=True)
        if name == "Plots":
            self.after(50, self.draw_plots)
        if name == "Pores" and self.k_fig is None:
            self.after(50, self.draw_pores)
        if name == "Compare" and self.c_fig is None:
            self.after(50, self.do_compare)
        for n, (b, lab, sl) in self.nav_btns.items():
            on = (n == name)
            for w in (b, lab, sl):
                w.configure(bg=C["card"] if on else C["bg"])
            lab.configure(fg=C["accent"] if on else C["ink"])

    def chip(self, text):
        tk.Label(self.chips, text=text, font=FONT_S, bg="#17508e", fg="white", padx=9, pady=3).pack(side="left", padx=3)

    def refresh_all(self):
        for w in self.chips.winfo_children():
            w.destroy()
        s = self.store.summary()
        if s["structures"]:
            self.chip(f'{s["viruses"]} viruses'); self.chip(f'{s["families"]} families')
            self.chip(f'{s["residues"]:,} gate residues'); self.chip(f'{s["dual"]} dual positions')
        else:
            self.chip("database empty - use Analyse")
        self.reload_browse()
        self.reload_pref_classes()
        self.status.set(f'database: {DB}' + ("" if s["structures"] else "   (no structures yet)"))

    # ---------- view 1: browse
    def view_browse(self, root):
        top = tk.Frame(root, bg=C["bg"]); top.pack(fill="x")
        tk.Label(top, text="Search", bg=C["bg"], fg=C["muted"], font=FONT_S).pack(side="left")
        self.b_search = tk.StringVar()
        e = tk.Entry(top, textvariable=self.b_search, font=FONT, bg=C["card"], relief="flat", width=26)
        e.pack(side="left", padx=(6, 14), ipady=4)
        e.bind("<KeyRelease>", lambda ev: self.reload_browse())
        tk.Label(top, text="Capsid class", bg=C["bg"], fg=C["muted"], font=FONT_S).pack(side="left")
        self.b_class = ttk.Combobox(top, state="readonly", width=18, font=FONT, values=["all"])
        self.b_class.set("all"); self.b_class.pack(side="left", padx=6)
        self.b_class.bind("<<ComboboxSelected>>", lambda ev: self.reload_browse())
        self.b_dual = tk.BooleanVar(value=False)
        ttk.Checkbutton(top, text="only viruses with a dual 5/3 position", variable=self.b_dual,
                        command=self.reload_browse).pack(side="left", padx=12)
        tk.Button(top, text="Export list", command=lambda: export_tree(self.b_tree, "viruses.csv"),
                  font=FONT_S, bg=C["card"], fg=C["ink"], relief="flat", padx=10, pady=3).pack(side="right")

        pan = tk.PanedWindow(root, orient="horizontal", bg=C["bg"], sashwidth=10, bd=0)
        pan.pack(fill="both", expand=True, pady=10)
        o1, i1 = card(pan, "Viruses in the database")
        w, self.b_tree = make_tree(i1, ("virus", "PDB", "family", "class", "5-fold", "3-fold", "dual"),
                                   (230, 60, 150, 110, 60, 60, 50),
                                   {"5-fold": "center", "3-fold": "center", "dual": "center", "PDB": "center"}, height=20)
        w.pack(fill="both", expand=True, padx=10, pady=(0, 10))
        self.b_tree.bind("<<TreeviewSelect>>", self.on_pick_virus)
        pan.add(o1, minsize=420)

        # chart above, residue table below; drag the divider between them to give either one more room
        right = tk.PanedWindow(pan, orient="vertical", bg=C["bg"], sashwidth=10, bd=0)
        o2, i2 = card(right, "Gate residues")
        self.b_title = tk.Label(i2, text="select a virus", font=FONT, bg=C["card"], fg=C["muted"], anchor="w", justify="left")
        self.b_title.pack(fill="x", padx=12, pady=(0, 4))
        self.b_pdb = None
        tk.Button(i2, text="View in 3D (opens in your browser)", command=self.open_3d, font=FONT_S, bg=C["accent_soft"],
                  fg=C["accent"], relief="flat", padx=8, cursor="hand2").pack(anchor="e", padx=12, pady=(0, 4))
        self.b_chart = DistanceChart(i2); self.b_chart.pack(fill="x", padx=6, pady=(0, 6))
        right.add(o2, minsize=260, height=455, stretch="never")
        o3, i3 = card(right, "Residues within the cutoff  (drag the bar above to resize)")
        w2, self.b_res = make_tree(i3, ("gate", "residue", "chain", "distance A", "dual 5/3"),
                                   (85, 130, 70, 110, 125), {"distance A": "center", "dual 5/3": "center"}, height=11)
        w2.pack(fill="both", expand=True, padx=10, pady=(0, 6))
        tk.Button(i3, text="Export residues", command=lambda: export_tree(self.b_res, "gate_residues.csv"),
                  font=FONT_S, bg=C["card"], fg=C["accent"], relief="flat").pack(anchor="e", padx=12, pady=(0, 8))
        right.add(o3, minsize=230, stretch="always")
        pan.add(right, minsize=430)

    def reload_browse(self):
        if not hasattr(self, "b_tree"):
            return
        txt = f'%{self.b_search.get().strip()}%' if self.b_search.get().strip() else "%"
        cls = self.b_class.get()
        sql = """SELECT s.virus, s.pdb, s.family, s.capsid_class,
                   (SELECT COUNT(*) FROM gate_residues g WHERE g.pdb=s.pdb AND g.gate='5-fold') n5,
                   (SELECT COUNT(*) FROM gate_residues g WHERE g.pdb=s.pdb AND g.gate='3-fold') n3,
                   (SELECT COUNT(DISTINCT auth_seq) FROM dual_positions d WHERE d.pdb=s.pdb) nd
                 FROM structures s
                 WHERE (s.virus LIKE ? OR s.pdb LIKE ? OR IFNULL(s.family,'') LIKE ?)"""
        args = [txt, txt, txt]
        if cls and cls != "all":
            sql += " AND s.capsid_class = ?"; args.append(cls)
        sql += " ORDER BY s.virus"
        rows = self.store.q(sql, tuple(args))
        if self.b_dual.get():
            rows = [r for r in rows if r["nd"]]
        out = [(r["virus"], r["pdb"].upper(), r["family"] or "", r["capsid_class"] or "", r["n5"], r["n3"], r["nd"] or "")
               for r in rows]
        fill_tree(self.b_tree, out)
        self.b_class["values"] = ["all"] + self.store.classes()
        self.status.set(f"{len(out)} viruses listed")

    def open_3d(self):
        if not self.b_pdb:
            messagebox.showinfo("Pick a virus first", "Select a virus in the list, then press View in 3D."); return
        webbrowser.open(f"{WEB_URL}?pdb={self.b_pdb}")

    def on_pick_virus(self, _ev=None):
        sel = self.b_tree.selection()
        if not sel:
            return
        pdb = self.b_tree.item(sel[0], "values")[1].lower()
        self.b_pdb = pdb
        st = self.store.q("SELECT * FROM structures WHERE pdb=?", (pdb,))
        res = self.store.q("SELECT * FROM gate_residues WHERE pdb=? ORDER BY gate, distance_A", (pdb,))
        dual = self.store.q("SELECT * FROM dual_positions WHERE pdb=? ORDER BY auth_seq", (pdb,))
        if st:
            s = st[0]
            cutoff = res[0]["cutoff_A"] if res else 8
            bits = [f'{s["virus"]}   |   PDB {s["pdb"].upper()}   |   {s["family"] or "family not recorded"}',
                    f'{s["capsid_class"] or "class not recorded"}' +
                    (f'   |   {s["resolution_A"]:.2f} A' if s["resolution_A"] else "") +
                    f'   |   {s["n_chains"]} deposited chains, {s["n_residues"]} residues   |   symmetry from {s["symmetry_from"]}']
            if dual:
                bits.append("guards both gates: " + ", ".join(
                    f'{d["aa"]}{d["auth_seq"]} (5-fold {d["min_d5"]:.1f} A in {d["chains_near_5fold"]}; '
                    f'3-fold {d["min_d3"]:.1f} A in {d["chains_near_3fold"]})' for d in dual))
            else:
                bits.append("no residue number guards both gates within the cutoff")
            self.b_title.configure(text="\n".join(bits), fg=C["ink"])
            self.b_chart.show([dict(gate=r["gate"], aa=r["aa"], auth_seq=r["auth_seq"],
                                    distance_A=r["distance_A"], dual_5_3=r["dual_5_3"]) for r in res], cutoff)
        rows = [(r["gate"], f'{r["aa"]}{r["auth_seq"]} ({r["resname"]})', r["chain"], f'{r["distance_A"]:.2f}',
                 "yes" if r["dual_5_3"] else "", "DUAL" if r["dual_5_3"] else None) for r in res]
        fill_tree(self.b_res, [tuple(x for x in r if x is not None) + (("DUAL",) if r[-1] else ()) for r in
                               [(a, b, c, d, e, f) for a, b, c, d, e, f in rows]])

    # ---------- view 2: find
    def view_find(self, root):
        o, i = card(root, "Find: where does an amino acid sit at the gates, in every virus?")
        tk.Label(i, text="Pick an amino acid and a gate, then Search. Distance = how far the residue's side chain is from the "
                         "symmetry axis: 0-3 A is right on the axis (in the pore), 8 A is the outer edge of the gate.",
                 bg=C["card"], fg=C["muted"], font=FONT_S, anchor="w", justify="left", wraplength=1100).pack(fill="x", padx=12)
        bar = tk.Frame(i, bg=C["card"]); bar.pack(fill="x", padx=12, pady=(6, 4))
        def lab(c, t):
            tk.Label(bar, text=t, bg=C["card"], fg=C["muted"], font=FONT_S).grid(row=0, column=c, sticky="w")
        lab(0, "amino acid")
        self.f_aa = ttk.Combobox(bar, width=18, state="readonly", font=FONT, values=["any"] + [AA_LABEL[a] for a in AA_ORDER])
        self.f_aa.set("any"); self.f_aa.grid(row=1, column=0, padx=(0, 12))
        lab(1, "gate (symmetry axis)")
        self.f_gate = ttk.Combobox(bar, width=12, state="readonly", font=FONT, values=["any", "5-fold", "3-fold", "2-fold"])
        self.f_gate.set("any"); self.f_gate.grid(row=1, column=1, padx=(0, 12))
        lab(2, "capsid type")
        self.f_class = ttk.Combobox(bar, width=18, state="readonly", font=FONT, values=["any"])
        self.f_class.set("any"); self.f_class.grid(row=1, column=2, padx=(0, 12))
        lab(3, "max distance (A)")
        self.f_dist = ttk.Combobox(bar, width=5, font=FONT, values=["3", "5", "8"])
        self.f_dist.set("8"); self.f_dist.grid(row=1, column=3, padx=(0, 12))
        lab(4, "virus or family contains")
        self.f_text = tk.StringVar()
        en = tk.Entry(bar, textvariable=self.f_text, width=16, font=FONT, bg=C["bg"], relief="flat")
        en.grid(row=1, column=4, padx=(0, 12), ipady=3); en.bind("<Return>", lambda ev: self.do_find())
        self.f_dual = tk.BooleanVar(value=False)
        ttk.Checkbutton(bar, text="only residues at BOTH gates", variable=self.f_dual,
                        style="Gate.TCheckbutton", command=self.do_find).grid(row=2, column=0, columnspan=3, sticky="w", pady=(6, 0))
        tk.Button(bar, text="Search", command=self.do_find, font=FONT_B, bg=C["accent"], fg="white",
                  relief="flat", padx=16, pady=4).grid(row=1, column=5, padx=6)
        for wdg in (self.f_aa, self.f_gate, self.f_class, self.f_dist):
            wdg.bind("<<ComboboxSelected>>", lambda ev: self.do_find())
        quick = tk.Frame(i, bg=C["card"]); quick.pack(fill="x", padx=12, pady=(0, 6))
        tk.Label(quick, text="quick questions:", bg=C["card"], fg=C["muted"], font=FONT_S).pack(side="left")
        for txt, preset in (("all residues at BOTH gates", dict(dual=True)),
                            ("on the 5-fold axis (< 3 A)", dict(gate="5-fold", dist="3")),
                            ("on the 3-fold axis (< 3 A)", dict(gate="3-fold", dist="3")),
                            ("Asn at both gates", dict(aa="N", dual=True))):
            tk.Button(quick, text=txt, command=lambda p=preset: self.find_preset(**p), font=FONT_S, bg=C["accent_soft"],
                      fg=C["accent"], relief="flat", padx=8, pady=2, cursor="hand2").pack(side="left", padx=4)
        tk.Button(quick, text="Export table", command=lambda: export_tree(self.f_tree, "gate_search.csv"), font=FONT_S,
                  bg=C["card"], fg=C["accent"], relief="flat").pack(side="right")
        self.f_answer = tk.Label(i, text="The answer appears here in words after you search.", bg=C["accent_soft"],
                                 fg=C["ink"], font=FONT, anchor="w", justify="left", wraplength=1100, padx=10, pady=6)
        self.f_answer.pack(fill="x", padx=12, pady=(0, 6))
        self.f_answer.bind("<Configure>", lambda e: self.f_answer.configure(wraplength=max(e.width - 24, 300)))
        cols = ("virus", "PDB", "family", "capsid type", "gate", "residue", "subunit", "distance A", "at both gates?")
        w, self.f_tree = make_tree(i, cols, (190, 50, 105, 120, 55, 70, 60, 75, 100),
                                   {"distance A": "center", "at both gates?": "center", "PDB": "center",
                                    "subunit": "center", "gate": "center"}, height=18)
        w.pack(fill="both", expand=True, padx=10, pady=(0, 4))
        self.f_tree.bind("<Double-1>", self.find_open_virus)
        self.f_count = tk.Label(i, text="Orange rows: the same residue number also sits at the other gate (5-fold and 3-fold).   "
                                        "Double-click a row to open that virus in Browse.",
                                bg=C["card"], fg=C["muted"], font=FONT_S, anchor="w")
        self.f_count.pack(fill="x", padx=12, pady=(0, 8))
        o.pack(fill="both", expand=True)

    def find_preset(self, aa="any", gate="any", dist="8", dual=False):
        self.f_aa.set(AA_LABEL.get(aa, "any")); self.f_gate.set(gate); self.f_dist.set(dist)
        self.f_dual.set(dual); self.f_class.set("any"); self.f_text.set("")
        self.do_find()

    def find_open_virus(self, _ev=None):
        sel = self.f_tree.selection()
        if not sel:
            return
        pdb = str(self.f_tree.item(sel[0], "values")[1])
        self.b_search.set(""); self.b_class.set("all"); self.b_dual.set(False); self.reload_browse()
        for iid in self.b_tree.get_children():
            if str(self.b_tree.item(iid, "values")[1]) == pdb:
                self.show("Browse"); self.b_tree.selection_set(iid); self.b_tree.focus(iid); self.b_tree.see(iid)
                return

    def do_find(self):
        try:
            dmax = float(self.f_dist.get())
        except ValueError:
            messagebox.showwarning("Check the distance", "Max distance must be a number, for example 8"); return
        aa = self.f_aa.get()[:1] if self.f_aa.get() != "any" else None
        cls = self.f_class.get()
        sql = "SELECT * FROM gate_residues WHERE distance_A <= ?"
        args = [dmax]
        if aa:
            sql += " AND aa = ?"; args.append(aa)
        if self.f_gate.get() != "any":
            sql += " AND gate = ?"; args.append(self.f_gate.get())
        if cls != "any":
            sql += " AND capsid_class = ?"; args.append(cls)
        if self.f_dual.get():
            sql += " AND dual_5_3 = 1"
        if self.f_text.get().strip():
            sql += " AND (virus LIKE ? OR IFNULL(family,'') LIKE ?)"
            args += [f'%{self.f_text.get().strip()}%'] * 2
        sql += " ORDER BY distance_A LIMIT 4000"
        rows = self.store.q(sql, tuple(args))
        fill_tree(self.f_tree, [(r["virus"], r["pdb"].upper(), r["family"] or "", r["capsid_class"] or "", r["gate"],
                                 f'{(r["resname"] or r["aa"]).title()}{r["auth_seq"]}', r["chain"],
                                 f'{r["distance_A"]:.2f}', "yes" if r["dual_5_3"] else "")
                                + (("DUAL",) if r["dual_5_3"] else ()) for r in rows])
        self.f_class["values"] = ["any"] + self.store.classes()
        # the answer, in words
        what = AA_LABEL[aa].split("  ", 1)[1] if aa else "Any amino acid"
        where = {"any": "at any gate", "5-fold": "at the 5-fold gate", "3-fold": "at the 3-fold gate",
                 "2-fold": "at the 2-fold axis"}[self.f_gate.get()]
        cond = f"{what} {where} (within {dmax:g} A of the axis)" + (", at BOTH gates" if self.f_dual.get() else "")
        if cls != "any":
            cond += f", {cls} capsids"
        if self.f_text.get().strip():
            cond += f', "{self.f_text.get().strip()}"'
        if not rows:
            self.f_answer.configure(text=f"{cond}:  nothing found. Try a larger distance or 'any' gate.", bg="#fbe9df")
            return
        q = "SELECT COUNT(DISTINCT virus) c FROM structures" + (" WHERE capsid_class=?" if cls != "any" else "")
        nv_all = self.store.q(q, (cls,) if cls != "any" else ())[0]["c"]
        vir = {r["virus"] for r in rows}
        fams = {}
        for r in rows:
            fams.setdefault(r["family"] or "unknown", set()).add(r["virus"])
        both = {r["virus"] for r in rows if r["dual_5_3"]}
        topf = ", ".join(f"{f} {len(v)}" for f, v in sorted(fams.items(), key=lambda t: -len(t[1]))[:4])
        txt = (f"{cond}:  found in {len(vir)} of {nv_all} viruses ({100 * len(vir) / max(nv_all, 1):.0f}%), "
               f"from {len(fams)} famil{'y' if len(fams) == 1 else 'ies'} ({topf}{', ...' if len(fams) > 4 else ''}).")
        if not self.f_dual.get():
            txt += f"  In {len(both)} of these viruses the same residue number sits at BOTH the 5-fold and 3-fold gate."
        if len(fams) == 1 and len(vir) > 1:
            txt += "  All from one family: probably inherited from a common ancestor, not a general rule."
        txt += f"  ({len(rows)} rows{', first 4000 shown' if len(rows) == 4000 else ''})"
        self.f_answer.configure(text=txt, bg=C["accent_soft"])

    # ---------- view 3: preferences
    def view_pref(self, root):
        o, i = card(root, "Which amino acids are over-represented at the gates")
        bar = tk.Frame(i, bg=C["card"]); bar.pack(fill="x", padx=12, pady=(2, 8))
        self.p_class = ttk.Combobox(bar, width=16, state="readonly", font=FONT, values=[])
        self.p_zone = ttk.Combobox(bar, width=24, state="readonly", font=FONT,
                                   values=["5-fold", "3-fold", "2-fold", "dual 5/3-fold position"])
        self.p_zone.set("5-fold")
        self.p_cut = ttk.Combobox(bar, width=6, state="readonly", font=FONT, values=["3", "6", "8", "10"])
        self.p_cut.set("8")
        for lab, wdg in (("capsid class", self.p_class), ("zone", self.p_zone), ("cutoff (A)", self.p_cut)):
            tk.Label(bar, text=lab, bg=C["card"], fg=C["muted"], font=FONT_S).pack(side="left", padx=(0, 4))
            wdg.pack(side="left", padx=(0, 14))
            wdg.bind("<<ComboboxSelected>>", lambda ev: self.reload_pref())
        tk.Button(bar, text="Export", command=lambda: export_tree(self.p_tree, "preferences.csv"), font=FONT_S,
                  bg=C["card"], fg=C["accent"], relief="flat").pack(side="right")
        split = tk.Frame(i, bg=C["card"]); split.pack(fill="both", expand=True, padx=10, pady=(0, 10))
        self.p_chart = BarChart(split); self.p_chart.pack(side="left", fill="both", expand=True, padx=(0, 8))
        w, self.p_tree = make_tree(split, ("aa", "at gate", "background", "odds ratio", "FDR", "families +/-"),
                                   (40, 80, 90, 85, 75, 90),
                                   {c: "center" for c in ("aa", "at gate", "background", "odds ratio", "FDR", "families +/-")},
                                   height=20)
        w.pack(side="left", fill="both")
        self.p_note = tk.Label(i, text="", bg=C["card"], fg=C["muted"], font=FONT_S, anchor="w", justify="left", wraplength=1000)
        self.p_note.pack(fill="x", padx=12, pady=(0, 8))
        o.pack(fill="both", expand=True)

    def reload_pref_classes(self):
        if not hasattr(self, "p_class"):
            return
        cl = [r[0] for r in self.store.q("SELECT DISTINCT capsid_class FROM preferences ORDER BY 1")]
        self.p_class["values"] = cl
        if cl and self.p_class.get() not in cl:
            self.p_class.set(cl[0])
        self.reload_pref()

    def reload_pref(self):
        if not hasattr(self, "p_tree"):
            return
        rows = self.store.q("""SELECT * FROM preferences WHERE capsid_class=? AND zone=? AND cutoff_A=? AND level='aa'
                               ORDER BY CASE WHEN odds_ratio IS NULL THEN 1 ELSE 0 END, odds_ratio DESC""",
                            (self.p_class.get(), self.p_zone.get(), float(self.p_cut.get() or 8)))
        self.p_chart.show(rows)
        fill_tree(self.p_tree, [(r["key"], f'{100*r["freq_gate"]:.1f}%', f'{100*r["freq_background"]:.1f}%',
                                 "-" if r["odds_ratio"] is None else f'{r["odds_ratio"]:.2f}',
                                 "-" if r["fdr"] is None else f'{r["fdr"]:.3f}',
                                 "" if r["families_enriched"] is None else f'{r["families_enriched"]}/{r["families_depleted"]}')
                                for r in rows])
        n = rows[0]["n_structures"] if rows else 0
        tot = rows[0]["n_total"] if rows else 0
        if rows:
            msg = (f'{n} viruses, {tot} residues in this zone. Odds ratio above 1 = more common at the gate than in the rest of '
                   f'the capsid shell at the same radius. "families +/-" counts virus families in which the amino acid is enriched '
                   f'or depleted: a high odds ratio with a split family count is a pattern, not a rule.')
        else:
            nv = (self.store.q('SELECT COUNT(DISTINCT pdb) FROM structures WHERE capsid_class=?', (self.p_class.get(),))
                  or [[0]])[0][0]
            nd = (self.store.q('SELECT COUNT(*) FROM dual_positions WHERE capsid_class=?', (self.p_class.get(),))
                  or [[0]])[0][0]
            if self.p_zone.get().startswith('dual') and nv and not nd:
                msg = (f'No residue number guards both gates in this capsid class: none was found in its {nv} viruses, '
                       f'so there is nothing to test. In pseudo-T=3 capsids the three subunits are different proteins, '
                       f'so a shared residue number is not expected.')
            else:
                msg = 'No pooled statistics stored for this selection. They are computed only for runs with 5 or more viruses.'
        self.p_note.configure(text=msg)

    # ---------- view: plots (figure code lives in gate_plots.py, shared with the command line)
    def view_plots(self, root):
        o, i = card(root, "Which residues recur at the 5-fold and 3-fold gates, and which are favoured there")
        bar = tk.Frame(i, bg=C["card"]); bar.pack(fill="x", padx=12, pady=(2, 6))
        self.g_class = ttk.Combobox(bar, width=18, state="readonly", font=FONT)
        self.g_cut = ttk.Combobox(bar, width=6, state="readonly", font=FONT, values=["3", "6", "8", "10"]); self.g_cut.set("8")
        self.g_mode = ttk.Combobox(bar, width=22, state="readonly", font=FONT,
                                   values=["Simple answer", "Full figure (4 panels)", "Stem-and-leaf (leaf = tenths)",
                                           "Stem-and-leaf (leaf = amino acid)"]); self.g_mode.set("Simple answer")
        self.g_aa = ttk.Combobox(bar, width=16, state="readonly", font=FONT, values=["closest residue"] + [AA_LABEL[a] for a in AA_ORDER])
        self.g_aa.set("closest residue")
        for lab, wdg in (("show", self.g_mode), ("capsid class", self.g_class), ("cutoff (A)", self.g_cut),
                         ("stem-and-leaf amino acid", self.g_aa)):
            tk.Label(bar, text=lab, bg=C["card"], fg=C["muted"], font=FONT_S).pack(side="left", padx=(0, 4))
            wdg.pack(side="left", padx=(0, 14))
            wdg.bind("<<ComboboxSelected>>", lambda ev: self.draw_plots())
        tk.Button(bar, text="Save figure (PNG / PDF / SVG)...", command=self.save_plots, font=FONT_S,
                  bg=C["card"], fg=C["accent"], relief="flat").pack(side="right")
        self.g_host = tk.Frame(i, bg=C["card"]); self.g_host.pack(fill="both", expand=True, padx=8, pady=(0, 8))
        self.g_fig = None
        o.pack(fill="both", expand=True)

    def draw_plots(self):
        import sqlite3
        try:
            import gate_plots
            if self.g_fig is None:
                from matplotlib.figure import Figure
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
                self.g_fig = Figure(figsize=(16, 11), dpi=72)
                self.g_canvas = FigureCanvasTkAgg(self.g_fig, master=self.g_host)
                NavigationToolbar2Tk(self.g_canvas, self.g_host).pack(side="bottom", fill="x")
                self.g_canvas.get_tk_widget().pack(fill="both", expand=True)
            con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
            vals = [gate_plots.ALL] + gate_plots.classes(con)
            self.g_class["values"] = vals
            if self.g_class.get() not in vals:
                self.g_class.set(gate_plots.ALL)
            mode = self.g_mode.get()
            if mode.startswith("Stem"):
                import gate_stemleaf
                aa = None if self.g_aa.get() == "closest residue" else self.g_aa.get()[:1]
                con2 = gate_stemleaf.connect(DB)
                res = gate_stemleaf.data(con2, self.g_class.get(), aa); con2.close()
                desc = gate_stemleaf.draw(self.g_fig, res, show_letters=mode.endswith("amino acid)"))
            else:
                simple = mode.startswith("Simple") and hasattr(gate_plots, "build_simple")
                draw = gate_plots.build_simple if simple else gate_plots.build
                desc = draw(self.g_fig, con, self.g_class.get(), float(self.g_cut.get() or 8))
            con.close()
            self.g_canvas.draw_idle()
            self.status.set("Plots: " + desc)
        except Exception as e:
            self.status.set(f"could not draw the plots: {e}")
            traceback.print_exc()

    def save_plots(self):
        if self.g_fig is None:
            return
        os.makedirs(os.path.join(HERE, "plots"), exist_ok=True)
        cls = (self.g_class.get() or "all").replace("=", "").replace(" ", "_").replace("(", "").replace(")", "")
        p = filedialog.asksaveasfilename(initialdir=os.path.join(HERE, "plots"), defaultextension=".png",
                                         initialfile=f"gate_plots_{cls}_cutoff{self.g_cut.get()}.png",
                                         filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
        if p:
            self.g_fig.savefig(p, dpi=200, bbox_inches="tight")
            messagebox.showinfo("Saved", p)

    # ---------- view: compare two groups (logic in gate_compare.py, shared with the command line)
    def view_pores(self, root):
        import gate_chemistry
        o, i = card(root, "Pores: what lines each opening (charge, hydropathy, bulk)")
        tk.Label(i, text="Each residue near the axis counts once for every copy around the pore: five times at a 5-fold "
                         "axis, three times at a 3-fold axis. His is counted as neutral.",
                 bg=C["card"], fg=C["muted"], font=FONT_S, anchor="w", justify="left", wraplength=1100).pack(fill="x", padx=12)
        bar = tk.Frame(i, bg=C["card"]); bar.pack(fill="x", padx=12, pady=(6, 4))
        self.k_gate = ttk.Combobox(bar, width=8, state="readonly", font=FONT, values=["5-fold", "3-fold"]); self.k_gate.set("5-fold")
        self.k_metric = ttk.Combobox(bar, width=44, state="readonly", font=FONT, values=list(gate_chemistry.METRICS.values()))
        self.k_metric.set(gate_chemistry.METRICS["net_charge"])
        self.k_by = ttk.Combobox(bar, width=12, state="readonly", font=FONT, values=["capsid class", "family"]); self.k_by.set("capsid class")
        self.k_dist = ttk.Combobox(bar, width=4, state="readonly", font=FONT, values=["3", "4", "5", "6", "7", "8"]); self.k_dist.set("8")
        for lab, wdg in (("axis", self.k_gate), ("measure", self.k_metric), ("group by", self.k_by), ("within (A)", self.k_dist)):
            tk.Label(bar, text=lab, bg=C["card"], fg=C["muted"], font=FONT_S).pack(side="left", padx=(0, 4))
            wdg.pack(side="left", padx=(0, 12))
            wdg.bind("<<ComboboxSelected>>", lambda ev: self.draw_pores())
        tk.Button(bar, text="Export table", command=self.export_pores, font=FONT_S, bg=C["card"], fg=C["accent"],
                  relief="flat").pack(side="right")
        self.k_answer = tk.Label(i, text="", bg=C["accent_soft"], fg=C["ink"], font=FONT, anchor="w", justify="left",
                                 wraplength=1100, padx=10, pady=6)
        self.k_answer.pack(fill="x", padx=12, pady=(0, 6))
        self.k_answer.bind("<Configure>", lambda e: self.k_answer.configure(wraplength=max(e.width - 24, 300)))
        self.k_host = tk.Frame(i, bg=C["card"]); self.k_host.pack(fill="both", expand=True, padx=8, pady=(0, 8))
        self.k_df = None
        o.pack(fill="both", expand=True)

    def draw_pores(self):
        import gate_chemistry
        try:
            if self.k_fig is None:
                from matplotlib.figure import Figure
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
                self.k_fig = Figure(figsize=(13, 7), dpi=72)
                self.k_canvas = FigureCanvasTkAgg(self.k_fig, master=self.k_host)
                self.k_canvas.get_tk_widget().pack(fill="both", expand=True)
            metric = next(k for k, v in gate_chemistry.METRICS.items() if v == self.k_metric.get())
            by = "family" if self.k_by.get() == "family" else "capsid_class"
            con = gate_chemistry.connect(DB)
            self.k_df = gate_chemistry.data(con, self.k_gate.get(), float(self.k_dist.get()))
            con.close()
            self.k_answer.configure(text=gate_chemistry.summary(self.k_df, by, metric))
            self.status.set("Pores: " + gate_chemistry.draw(self.k_fig, self.k_df, by, metric))
            self.k_canvas.draw_idle()
        except Exception as e:
            self.k_answer.configure(text=f"could not draw: {e}")
            traceback.print_exc()

    def export_pores(self):
        if self.k_df is None:
            return
        p = filedialog.asksaveasfilename(defaultextension=".csv", initialfile=f"pore_chemistry_{self.k_gate.get()}.csv",
                                         filetypes=[("CSV", "*.csv")])
        if p:
            self.k_df.to_csv(p, index=False)
            messagebox.showinfo("Saved", p)

    def view_compare(self, root):
        o, i = card(root, "Compare: which amino acids does one group of viruses favour at the gates, and the other?")
        tk.Label(i, text="Pick two capsid types or families. For every virus the tool takes the residues at the chosen gate "
                         "and works out what share each amino acid makes up; the groups' average shares are then compared. "
                         "Yellow rows with * are differences unlikely to be chance (FDR < 0.05).",
                 bg=C["card"], fg=C["muted"], font=FONT_S, anchor="w", justify="left", wraplength=1100).pack(fill="x", padx=12)
        bar = tk.Frame(i, bg=C["card"]); bar.pack(fill="x", padx=12, pady=(6, 4))
        def lab(c, t):
            tk.Label(bar, text=t, bg=C["card"], fg=C["muted"], font=FONT_S).grid(row=0, column=c, sticky="w")
        lab(0, "first group (blue)"); lab(2, "second group (orange)"); lab(3, "gate"); lab(4, "max distance (A)")
        self.c_a = ttk.Combobox(bar, width=34, state="readonly", font=FONT)
        self.c_a.grid(row=1, column=0, padx=(0, 6))
        tk.Label(bar, text="vs", bg=C["card"], fg=C["ink"], font=FONT_B).grid(row=1, column=1, padx=4)
        self.c_b = ttk.Combobox(bar, width=34, state="readonly", font=FONT)
        self.c_b.grid(row=1, column=2, padx=(6, 12))
        self.c_gate = ttk.Combobox(bar, width=22, state="readonly", font=FONT,
                                   values=["5-fold gate", "3-fold gate", "BOTH gates (dual positions)"])
        self.c_gate.set("5-fold gate"); self.c_gate.grid(row=1, column=3, padx=(0, 12))
        self.c_dist = ttk.Combobox(bar, width=5, state="readonly", font=FONT, values=["8", "6", "5", "3"])
        self.c_dist.set("8"); self.c_dist.grid(row=1, column=4, padx=(0, 12))
        tk.Button(bar, text="Compare", command=self.do_compare, font=FONT_B, bg=C["accent"], fg="white",
                  relief="flat", padx=16, pady=4).grid(row=1, column=5, padx=6)
        tk.Button(bar, text="swap", command=self.swap_compare, font=FONT_S, bg=C["accent_soft"], fg=C["accent"],
                  relief="flat", padx=8).grid(row=1, column=6, padx=4)
        for wdg in (self.c_a, self.c_b, self.c_gate, self.c_dist):
            wdg.bind("<<ComboboxSelected>>", lambda ev: self.do_compare())
        quick = tk.Frame(i, bg=C["card"]); quick.pack(fill="x", padx=12, pady=(0, 6))
        tk.Label(quick, text="quick comparisons:", bg=C["card"], fg=C["muted"], font=FONT_S).pack(side="left")
        for txt, a, b, g in (("T=3 vs pseudo T=3", "class=T=3", "class=pseudo T=3 (P=3)", "5-fold gate"),
                             ("T=3 vs T=4", "class=T=3", "class=T=4", "5-fold gate"),
                             ("Nodaviridae vs Caliciviridae", "family=Nodaviridae", "family=Caliciviridae", "5-fold gate"),
                             ("Picornaviridae vs all others", "family=Picornaviridae", "rest", "3-fold gate")):
            tk.Button(quick, text=txt, command=lambda a=a, b=b, g=g: self.compare_preset(a, b, g), font=FONT_S,
                      bg=C["accent_soft"], fg=C["accent"], relief="flat", padx=8, pady=2, cursor="hand2").pack(side="left", padx=4)
        tk.Button(quick, text="Save figure...", command=self.save_compare, font=FONT_S, bg=C["card"], fg=C["accent"],
                  relief="flat").pack(side="right")
        tk.Button(quick, text="Export table", command=lambda: export_tree(self.c_tree, "compare.csv"), font=FONT_S,
                  bg=C["card"], fg=C["accent"], relief="flat").pack(side="right", padx=6)
        self.c_answer = tk.Label(i, text="", bg=C["accent_soft"], fg=C["ink"], font=FONT, anchor="w", justify="left",
                                 wraplength=1100, padx=10, pady=6)
        self.c_answer.pack(fill="x", padx=12, pady=(0, 6))
        self.c_answer.bind("<Configure>", lambda e: self.c_answer.configure(wraplength=max(e.width - 24, 300)))
        split = tk.PanedWindow(i, orient="horizontal", bg=C["card"], sashwidth=8, bd=0)
        split.pack(fill="both", expand=True, padx=8, pady=(0, 8))
        self.c_host = tk.Frame(split, bg=C["card"]); split.add(self.c_host, minsize=420, stretch="always")
        tw = tk.Frame(split, bg=C["card"])
        cols = ("amino acid", "1st %", "2nd %", "1st has it", "2nd has it", "FDR")
        w, self.c_tree = make_tree(tw, cols, (115, 60, 60, 80, 80, 55),
                                   {c: "center" for c in cols[1:]}, height=20)
        w.pack(fill="both", expand=True)
        split.add(tw, minsize=300, stretch="never")
        self.c_fig = None; self.c_spec = {}; self.c_res = None
        o.pack(fill="both", expand=True)

    def _compare_choices(self):
        import gate_compare
        con = gate_compare.connect(DB)
        ch = gate_compare.group_choices(con); con.close()
        self.c_spec = {lab: spec for lab, spec in ch}
        self.c_spec["all other viruses (everything not in the first group)"] = gate_compare.REST
        labs = [lab for lab, _ in ch]
        self.c_a["values"] = labs
        self.c_b["values"] = labs + ["all other viruses (everything not in the first group)"]
        by_spec = {v: k for k, v in self.c_spec.items()}
        if self.c_a.get() not in self.c_spec:
            self.c_a.set(by_spec.get("class=T=3", labs[0] if labs else ""))
        if self.c_b.get() not in self.c_spec:
            self.c_b.set(by_spec.get("class=pseudo T=3 (P=3)", labs[1] if len(labs) > 1 else ""))
        return by_spec

    def compare_preset(self, a, b, gate):
        by_spec = self._compare_choices()
        if a in by_spec and b in by_spec:
            self.c_a.set(by_spec[a]); self.c_b.set(by_spec[b]); self.c_gate.set(gate)
            self.do_compare()

    def swap_compare(self):
        a, b = self.c_a.get(), self.c_b.get()
        if self.c_spec.get(b) == "rest":
            return
        self.c_a.set(b); self.c_b.set(a); self.do_compare()

    def do_compare(self):
        import gate_compare
        try:
            self._compare_choices()
            if self.c_fig is None:
                from matplotlib.figure import Figure
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
                self.c_fig = Figure(figsize=(9, 8), dpi=72)
                self.c_canvas = FigureCanvasTkAgg(self.c_fig, master=self.c_host)
                self.c_canvas.get_tk_widget().pack(fill="both", expand=True)
            gate = {"5-fold gate": "5-fold", "3-fold gate": "3-fold"}.get(self.c_gate.get(), "both")
            self.c_answer.configure(text="comparing ...", bg=C["accent_soft"]); self.update_idletasks()
            con = gate_compare.connect(DB)
            res = gate_compare.compare(con, self.c_spec[self.c_a.get()], self.c_spec[self.c_b.get()], gate,
                                       float(self.c_dist.get() or 8))
            con.close()
            self.c_res = res
            self.c_answer.configure(text=gate_compare.draw(self.c_fig, res))
            self.c_canvas.draw_idle()
            fill_tree(self.c_tree, [(f'{r["aa"]}  {r["name"]}', f'{100 * r["sa"]:.1f}%', f'{100 * r["sb"]:.1f}%',
                                     f'{r["ka"]}/{res["na"]}', f'{r["kb"]}/{res["nb"]}', f'{r["fdr"]:.3f}')
                                    + (("DUAL",) if r["fdr"] < 0.05 else ()) for r in res["rows"]])
            self.status.set(f"Compare: {gate_compare.label(res['a'])} ({res['na']} viruses) vs "
                            f"{gate_compare.label(res['b'])} ({res['nb']} viruses)")
        except Exception as e:
            self.c_answer.configure(text=f"could not compare: {e}", bg="#fbe9df")
            traceback.print_exc()

    def save_compare(self):
        if self.c_fig is None or self.c_res is None:
            return
        import gate_compare
        os.makedirs(os.path.join(HERE, "plots"), exist_ok=True)
        clean = lambda t: "".join(ch for ch in t if ch.isalnum() or ch in "-_") or "group"
        name = f"compare_{clean(gate_compare.label(self.c_res['a']))}_vs_{clean(gate_compare.label(self.c_res['b']))}_{self.c_res['gate']}"
        p = filedialog.asksaveasfilename(initialdir=os.path.join(HERE, "plots"), defaultextension=".png",
                                         initialfile=name + ".png",
                                         filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
        if p:
            self.c_fig.savefig(p, dpi=200, bbox_inches="tight")
            gate_compare.write_csv(self.c_res, os.path.splitext(p)[0] + ".csv")
            messagebox.showinfo("Saved", f"{p}\n\nand the table as {os.path.splitext(p)[0] + '.csv'}")

    # ---------- view 4: analyse
    def view_run(self, root):
        o, i = card(root, "Analyse new structures and add them to the database")
        g = tk.Frame(i, bg=C["card"]); g.pack(fill="x", padx=12, pady=(2, 6)); g.columnconfigure(1, weight=1)
        self.r_pdb = tk.StringVar(value="4ftb")
        self.r_models = tk.StringVar()
        self.r_csv = tk.StringVar()
        self.r_res = tk.StringVar(value="")
        self.r_cut = tk.StringVar(value="8")
        rows = [("PDB IDs (space separated)", self.r_pdb, None),
                ("Your own model files (.pdb / .cif)", self.r_models, self.pick_models),
                ("or a CSV list of entries", self.r_csv, self.pick_csv)]
        for n, (lab, var, cb) in enumerate(rows):
            tk.Label(g, text=lab, bg=C["card"], fg=C["muted"], font=FONT_S).grid(row=n, column=0, sticky="w", pady=3)
            tk.Entry(g, textvariable=var, font=FONT, bg=C["bg"], relief="flat").grid(row=n, column=1, sticky="ew", padx=8, ipady=4)
            if cb:
                tk.Button(g, text="Browse", command=cb, font=FONT_S, bg=C["card"], fg=C["accent"], relief="flat").grid(row=n, column=2)
        opt = tk.Frame(i, bg=C["card"]); opt.pack(fill="x", padx=12, pady=4)
        tk.Label(opt, text="gate cutoff (A)", bg=C["card"], fg=C["muted"], font=FONT_S).pack(side="left")
        tk.Entry(opt, textvariable=self.r_cut, width=5, font=FONT, bg=C["bg"], relief="flat").pack(side="left", padx=(6, 16), ipady=3)
        tk.Label(opt, text="residue number to locate", bg=C["card"], fg=C["muted"], font=FONT_S).pack(side="left")
        tk.Entry(opt, textvariable=self.r_res, width=8, font=FONT, bg=C["bg"], relief="flat").pack(side="left", padx=(6, 16), ipady=3)
        self.r_add = tk.BooleanVar(value=True)
        ttk.Checkbutton(opt, text="add the results to the database when it finishes", variable=self.r_add,
                        style="Gate.TCheckbutton").pack(side="left")
        bt = tk.Frame(i, bg=C["card"]); bt.pack(fill="x", padx=12, pady=6)
        self.r_go = tk.Button(bt, text="Run analysis", command=self.run, font=FONT_B, bg=C["accent"], fg="white",
                              relief="flat", padx=18, pady=5); self.r_go.pack(side="left")
        self.r_stop = tk.Button(bt, text="Stop", command=self.stop, font=FONT, bg=C["card"], fg=C["warm"],
                                relief="flat", padx=12, state="disabled"); self.r_stop.pack(side="left", padx=8)
        for txt, cmd in (("Open results folder", self.open_out), ("Check setup", self.check),
                         ("Use bundled virus list", self.bundled)):
            tk.Button(bt, text=txt, command=cmd, font=FONT_S, bg=C["card"], fg=C["ink"], relief="flat", padx=10).pack(side="left", padx=4)
        self.r_status = tk.Label(i, text="idle", bg=C["card"], fg=C["muted"], font=FONT_S, anchor="w")
        self.r_status.pack(fill="x", padx=12)
        self.r_log = scrolledtext.ScrolledText(i, height=16, font=MONO, wrap="none", bg="#0f1a26", fg="#d8e2ec",
                                               insertbackground="white", relief="flat")
        self.r_log.pack(fill="both", expand=True, padx=10, pady=(6, 10))
        self.log("Enter a PDB ID or pick your own model file, then press Run analysis.\n"
                 "Every run writes run.log into its own results folder, so a failure always leaves a trace.\n")
        o.pack(fill="both", expand=True)

    def log(self, t):
        self.r_log.insert("end", t); self.r_log.see("end")
        if self.last_out:
            try:
                with open(os.path.join(self.last_out, "gui_session.log"), "a", encoding="utf-8") as fh:
                    fh.write(t)
            except Exception:
                pass

    def pick_models(self):
        ps = filedialog.askopenfilenames(filetypes=[("Models", "*.pdb *.cif *.mmcif *.ent *.gz"), ("All", "*.*")])
        if ps:
            self.r_models.set(";".join(ps)); self.r_pdb.set("")

    def pick_csv(self):
        p = filedialog.askopenfilename(filetypes=[("CSV", "*.csv")])
        if p:
            self.r_csv.set(p); self.r_pdb.set("")

    def bundled(self):
        p = os.path.join(HERE, "representatives.csv")
        if not os.path.exists(p):
            messagebox.showinfo("Not found", "representatives.csv is not next to the program."); return
        self.r_csv.set(p); self.r_pdb.set("")
        self.log("Bundled list selected: one structure per virus, 176 viruses. A full run takes about an hour; "
                 "results are saved structure by structure.\n")

    def open_out(self):
        p = self.last_out or RESULTS
        os.makedirs(p, exist_ok=True)
        try:
            os.startfile(p)
        except Exception:
            webbrowser.open("file://" + p)

    def run(self):
        if not os.path.exists(SCRIPT):
            messagebox.showerror("Missing file", f"gate_residue_analysis.py was not found next to this program:\n{SCRIPT}")
            return
        ids = [x for x in self.r_pdb.get().replace(",", " ").split() if x]
        models = [m for m in self.r_models.get().split(";") if m.strip()]
        csvf = self.r_csv.get().strip()
        if not (ids or models or csvf):
            messagebox.showwarning("Input needed", "Enter PDB IDs, choose model files, or pick a CSV list."); return
        for m in models:
            if not os.path.exists(m):
                messagebox.showwarning("File not found", m); return
        try:
            cut = float(self.r_cut.get())
            if not 1 <= cut <= 40:
                raise ValueError
        except ValueError:
            messagebox.showwarning("Check the cutoff", "The cutoff must be a number between 1 and 40 A."); return
        if self.r_res.get().strip() and not self.r_res.get().strip().lstrip("-").isdigit():
            messagebox.showwarning("Check the residue number", "The residue number must be a whole number."); return
        out = os.path.join(RESULTS, "run_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        try:
            os.makedirs(out, exist_ok=True)
        except Exception as e:
            messagebox.showerror("Cannot create the results folder", f"{out}\n\n{e}"); return
        self.last_out = out
        cmd = [PY, "-u", SCRIPT, "--cache", os.path.join(HERE, "pdbcache"), "--out", out,
               "--cutoff", str(cut), "--atoms", "sc", "--background", "matched"]
        if ids:
            cmd += ["--pdb"] + ids
        if models:
            cmd += ["--model"] + models
        if csvf:
            cmd += ["--csv", csvf]
        if self.r_res.get().strip():
            cmd += ["--residue", self.r_res.get().strip()]
        self.log("\n$ " + " ".join(cmd) + "\n")
        self.r_status.configure(text="running ...", fg=C["accent"])
        self.r_go.configure(state="disabled"); self.r_stop.configure(state="normal")
        threading.Thread(target=self.worker, args=(cmd, out, bool(self.r_add.get())), daemon=True).start()

    def worker(self, cmd, out, add_to_db):
        tail = []
        self.q.put("[environment check] importing numpy, scipy, pandas, statsmodels, matplotlib, gemmi, requests ...\n")
        good, rep, bad = gate_env.probe(PY)
        if not good:
            msg = (f"Environment check failed for {PY}\nFailing module: {bad}\n{rep}"
                   f"Fix: delete the .venv folder and start 'Capsid Symmetry Tool.bat' again to rebuild it.\n")
            try:
                with open(os.path.join(out, "run.log"), "a", encoding="utf-8") as fh:
                    fh.write(f"gate_explorer  {datetime.datetime.now():%Y-%m-%d %H:%M:%S}\n" + msg)
            except Exception:
                pass
            self.q.put(msg)
            self.q.put(("done", -2, False, out, msg, ""))
            return
        self.q.put(rep + "[environment check] ok\n")
        try:
            self.proc = subprocess.Popen(cmd, cwd=HERE, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                                         creationflags=NOWIN, env=gate_env.child_env(PY))
            for line in self.proc.stdout:
                self.q.put(line); tail = (tail + [line])[-25:]
            rc = self.proc.wait()
        except Exception as e:
            rc = -1; tail.append(f"could not start the analysis: {e}\n"); self.q.put(tail[-1])
        ok = rc == 0 and os.path.exists(os.path.join(out, "gate_summary.txt"))
        added = ""
        if ok and add_to_db:
            try:
                import gate_db
                con = gate_db.connect(DB)
                cat = gate_db._catalog(os.path.join(HERE, "representatives.csv"))
                n = gate_db.add_run(con, out, cat)
                gate_db.set_meta(con, built=datetime.datetime.now().isoformat(timespec="seconds"))
                con.close()
                added = f"{n} structure(s) added to the database"
                self.q.put(f"\n[database] {added}: {DB}\n")
            except Exception as e:
                added = f"could not add to the database: {e}"
                self.q.put(f"\n[database] {added}\n{traceback.format_exc()}\n")
        self.q.put(("done", rc, ok, out, "".join(tail[-12:]), added))

    def stop(self):
        if self.proc and self.proc.poll() is None:
            self.proc.kill()
            self.log("\n[stopped - structures finished so far are saved in the results folder]\n")

    def finished(self, rc, ok, out, tail, added):
        self.r_go.configure(state="normal"); self.r_stop.configure(state="disabled")
        n = len([f for f in os.listdir(out) if os.path.isfile(os.path.join(out, f))]) if os.path.isdir(out) else 0
        if ok:
            self.r_status.configure(text=f"finished - {n} files in {out}. {added}", fg=C["good"])
            self.log(f"\n[finished] {n} files in {out}\n")
            self.store.open(); self.refresh_all()
            if messagebox.askyesno("Finished", f"Results saved in:\n{out}\n\n{added}\n\nOpen gate_summary.txt?"):
                try:
                    os.startfile(os.path.join(out, "gate_summary.txt"))
                except Exception:
                    pass
        else:
            self.r_status.configure(text=f"failed (exit code {rc}) - see run.log in {out}", fg=C["warm"])
            self.log(f"\n[failed, exit code {rc}] see {os.path.join(out, 'run.log')}\n")
            messagebox.showerror("Python environment is broken - analysis not started" if rc == -2 else "The analysis did not finish",
                                 f"Exit code {rc}. Files saved so far: {n}.\n\nLog: {os.path.join(out, 'run.log')}\n\nLast lines:\n{tail}")

    def check(self):
        self.log("\n--- setup check ---\n")
        threading.Thread(target=self.check_worker, args=(RESULTS,), daemon=True).start()

    def check_worker(self, base):
        P = self.q.put
        P(f"python: {PY}\nprogram: {SCRIPT} ({'found' if os.path.exists(SCRIPT) else 'MISSING'})\n"
          f"database: {DB} ({'found' if os.path.exists(DB) else 'will be created'})\n")
        P(gate_env.probe(PY)[1])
        try:
            import urllib.request
            urllib.request.urlopen(urllib.request.Request("https://files.rcsb.org/download/4FTB.cif.gz", method="HEAD"), timeout=20)
            P("  ok   files.rcsb.org reachable\n")
        except Exception as e:
            P(f"  FAIL files.rcsb.org not reachable ({e}); PDB IDs cannot be downloaded, local model files still work\n")
        try:
            os.makedirs(base, exist_ok=True)
            t = os.path.join(base, "_write_test.txt"); open(t, "w").write("x"); os.remove(t)
            P(f"  ok   can write to {base}\n")
        except Exception as e:
            P(f"  FAIL cannot write to {base}: {e}\n")
        P("--- end of check ---\n")

    # ---------- view 5: help
    def view_help(self, root):
        o, i = card(root, "Help: what the Capsid Symmetry Tool measures")
        t = tk.Text(i, wrap="word", font=FONT, bg=C["card"], fg=C["ink"], relief="flat", padx=14, pady=8, height=30)
        sb = ttk.Scrollbar(i, orient="vertical", command=t.yview); t.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y"); t.pack(fill="both", expand=True)
        t.tag_configure("h", font=FONT_B, foreground=C["accent"], spacing1=8, spacing3=3)
        t.tag_configure("b", font=FONT_B)
        def put(tag, s):
            t.insert("end", s + "\n", tag)
        put("h", "In one paragraph")
        put("", "A virus shell (capsid) is built from many copies of the same protein arranged with icosahedral symmetry. "
                "It has three kinds of symmetry axis: 5-fold (6 axes), 3-fold (10) and 2-fold (15). Where an axis passes "
                "through the shell there is a pore or a junction - a gate. This tool measures which amino acids sit at "
                "those gates in every virus in the database, which residue numbers sit at BOTH the 5-fold and the 3-fold "
                "gate, and which amino acids are favoured there more than chance would give.")
        put("h", "Which page answers which question")
        put("", "Browse - one virus: what sits at its gates? The chart puts each residue at its distance from the axis "
                "(left = on the axis). Orange = the same residue number is at both gates.\n"
                "Find - one amino acid: in which viruses is it at a gate? The blue box under the buttons gives the answer "
                "in words; double-click a row to open that virus.\n"
                "Preferences - the statistics table: odds ratio above 1 = more common at the gate than in the rest of the "
                "shell at the same radius; FDR below 0.05 = unlikely to be chance.\n"
                "Pores - charge, hydropathy and bulk of what lines each opening, by capsid type or family.\n"
                "Browse also has a View in 3D button that opens the whole capsid in your browser.\n"
                "Compare - two capsid types or families side by side: which amino acids each group puts at a gate, "
                "as a share of its gate residues, and whether the difference is real (yellow * = FDR below 0.05).\n"
                "Plots - the overall answer. 'Simple answer' ranks amino acids from favoured at both gates (top) to "
                "avoided (bottom); 'Full figure' adds how often each one recurs and which virus families it comes from; 'Stem-and-leaf' puts every virus on one chart as a leaf at the distance of its closest residue (or of one amino acid you pick).\n"
                "Analyse - add a new PDB entry or your own model to the database.")
        put("h", "Words used on the screens")
        put("", "gate = the region within the cutoff (8 A by default) of a symmetry axis.   distance A = shortest distance "
                "from the residue's side chain to the axis, in angstrom.   subunit / chain = which protein copy.   "
                "dual position / at both gates = the same residue number is near a 5-fold axis in one copy and near a "
                "3-fold axis in another.   capsid type T=3, T=4, pseudo T=3 = how many protein copies make up one "
                "triangular face (3 or 4 copies; pseudo T=3 = three different proteins).")
        put("h", "The measurement")
        put("", "For every subunit of the capsid the program finds the true icosahedral 5-fold, 3-fold and 2-fold axes, then "
                "measures the shortest distance from each residue's side chain to the nearest axis of each kind. A residue "
                "within the cutoff (8 A by default) is called a gate residue. A residue number that is near a 5-fold axis in "
                "one subunit and near a 3-fold axis in another is a dual position: the same position in the sequence guards "
                "both gates, as Met182 does in Flock House virus.")
        put("h", "What you can do with it")
        put("", "Browse: see the gate residues of any virus in the database and which position guards both gates.\n"
                "Plots: see which amino acids are favoured at both gates, as one ranked chart or a four-panel figure.\n"
                "Find: ask where one amino acid sits across every virus, or list all dual positions in a family.\n"
                "Preferences: see which amino acids are over-represented at each gate, pooled over viruses.\n"
                "Analyse: run a new PDB entry or your own model and add it to the database, so your own structures can be "
                "compared with the published ones under identical rules.")
        put("h", "How the statistics work")
        put("", "Each amino acid is compared against capsid residues of the same structure that are not near any axis and lie "
                "at the same radius in the shell. Structures are combined with a Cochran-Mantel-Haenszel odds ratio, so a "
                "virus with an unusual composition cannot drive the result, and p-values are corrected for testing 20 amino "
                "acids. Pooled statistics are only computed when a run contains 5 or more viruses.")
        put("h", "Limits you should keep in mind")
        put("", "Distance to an axis is not pore radius: a residue near the axis is not proven to block anything.\n"
                "Only atomic models are read, never density maps.\n"
                "Quasi-3-fold axes of T=3 capsids are not analysed, only true icosahedral axes.\n"
                "A conserved gate residue tells you where to look, not what it does; that still needs an experiment.")
        put("h", "Your own model files")
        put("", "A local .pdb or .cif file must either contain the icosahedral operators (REMARK 350 or the mmCIF assembly "
                "records) or the complete capsid, otherwise the symmetry axes cannot be established and the file is skipped "
                "with the reason written to run.log.")
        t.configure(state="disabled")
        o.pack(fill="both", expand=True)

    # ---------- plumbing
    def pump(self):
        try:
            while True:
                m = self.q.get_nowait()
                if isinstance(m, tuple):
                    self.finished(*m[1:])
                else:
                    self.log(m)
                    if m.startswith("["):
                        part = m.split("]")[0].strip("[")
                        if "/" in part:
                            self.r_status.configure(text=f"running: structure {part}", fg=C["accent"])
        except queue.Empty:
            pass
        self.after(150, self.pump)

    def close(self):
        if self.proc and self.proc.poll() is None:
            if not messagebox.askyesno("Run in progress", "Stop the running analysis and close?"):
                return
            self.stop()
        self.destroy()


if __name__ == "__main__":
    App().mainloop()
