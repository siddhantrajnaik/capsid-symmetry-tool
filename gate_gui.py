"""Gate Residue Analysis - desktop launcher (tkinter front end for gate_residue_analysis.py)."""
import os, sys, subprocess, threading, queue, datetime, re, glob
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "gate_residue_analysis.py")
RESULTS = os.path.join(HERE, "results")
import gate_env
PY = gate_env.find_python()                                 # .venv\Scripts\python.exe when present
NOWIN = 0x08000000 if os.name == "nt" else 0               # CREATE_NO_WINDOW


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Gate Residue Analysis - 5-fold / 3-fold gates of icosahedral capsids")
        self.geometry("900x720")
        self.q = queue.Queue(); self.proc = None; self.last_out = None; self.outdir = None
        f = ttk.Frame(self, padding=10); f.pack(fill="both", expand=True); f.columnconfigure(1, weight=1)

        r = 0
        ttk.Label(f, text="PDB IDs (space or comma separated):").grid(row=r, column=0, sticky="w")
        self.pdb = tk.StringVar(value="4ftb"); ttk.Entry(f, textvariable=self.pdb).grid(row=r, column=1, columnspan=2, sticky="ew", pady=2); r += 1
        ttk.Label(f, text="Your own model files (.pdb / .cif):").grid(row=r, column=0, sticky="w")
        self.models = tk.StringVar(); ttk.Entry(f, textvariable=self.models).grid(row=r, column=1, sticky="ew", pady=2)
        ttk.Button(f, text="Add files...", command=self.pick_models).grid(row=r, column=2, sticky="e"); r += 1
        ttk.Label(f, text="or CSV list (pdb_id, group, family, T_label):").grid(row=r, column=0, sticky="w")
        self.csv = tk.StringVar(); ttk.Entry(f, textvariable=self.csv).grid(row=r, column=1, sticky="ew", pady=2)
        bf = ttk.Frame(f); bf.grid(row=r, column=2, sticky="e")
        ttk.Button(bf, text="Browse", command=self.pick_csv).pack(side="left")
        ttk.Button(bf, text="Use bundled list", command=self.bundled).pack(side="left", padx=3); r += 1

        ttk.Label(f, text="Residue number to query (optional):").grid(row=r, column=0, sticky="w")
        self.res = tk.StringVar(value="182"); ttk.Entry(f, textvariable=self.res, width=10).grid(row=r, column=1, sticky="w", pady=2); r += 1
        ttk.Label(f, text="Gate cutoff (A from axis):").grid(row=r, column=0, sticky="w")
        self.cut = tk.StringVar(value="8"); ttk.Entry(f, textvariable=self.cut, width=10).grid(row=r, column=1, sticky="w", pady=2); r += 1
        ttk.Label(f, text="Sensitivity cutoffs (A, used with 5+ viruses):").grid(row=r, column=0, sticky="w")
        self.sens = tk.StringVar(value="6 10"); ttk.Entry(f, textvariable=self.sens, width=14).grid(row=r, column=1, sticky="w", pady=2); r += 1
        ttk.Label(f, text="Atoms used for distance:").grid(row=r, column=0, sticky="w")
        self.atoms = tk.StringVar(value="sc"); af = ttk.Frame(f); af.grid(row=r, column=1, sticky="w")
        ttk.Radiobutton(af, text="side chain", variable=self.atoms, value="sc").pack(side="left")
        ttk.Radiobutton(af, text="all heavy atoms", variable=self.atoms, value="all").pack(side="left", padx=8); r += 1
        ttk.Label(f, text="Background:").grid(row=r, column=0, sticky="w")
        self.bg = tk.StringVar(value="matched"); bf2 = ttk.Frame(f); bf2.grid(row=r, column=1, sticky="w")
        ttk.Radiobutton(bf2, text="radius-matched", variable=self.bg, value="matched").pack(side="left")
        ttk.Radiobutton(bf2, text="all other residues", variable=self.bg, value="all").pack(side="left", padx=8); r += 1

        ttk.Label(f, text="Save results in (a new time-stamped\nfolder is made inside it for every run):").grid(row=r, column=0, sticky="w")
        self.base = tk.StringVar(value=RESULTS); ttk.Entry(f, textvariable=self.base).grid(row=r, column=1, sticky="ew", pady=2)
        ttk.Button(f, text="Browse", command=self.pick_out).grid(row=r, column=2, sticky="e"); r += 1
        self.resume = tk.BooleanVar(value=False)
        ttk.Checkbutton(f, text="Resume the newest run (skip structures it already finished)", variable=self.resume).grid(row=r, column=1, sticky="w"); r += 1

        bt = ttk.Frame(f); bt.grid(row=r, column=0, columnspan=3, pady=8, sticky="w")
        self.run_btn = ttk.Button(bt, text="Run analysis", command=self.run); self.run_btn.pack(side="left")
        self.stop_btn = ttk.Button(bt, text="Stop", command=self.stop, state="disabled"); self.stop_btn.pack(side="left", padx=6)
        ttk.Button(bt, text="Open results folder", command=self.open_out).pack(side="left")
        ttk.Button(bt, text="Open summary", command=self.open_summary).pack(side="left", padx=6)
        ttk.Button(bt, text="Check setup", command=self.check).pack(side="left")
        ttk.Button(bt, text="Help", command=self.help).pack(side="left", padx=6); r += 1
        self.status = tk.StringVar(value="Ready."); ttk.Label(f, textvariable=self.status, foreground="#245").grid(row=r, column=0, columnspan=3, sticky="w"); r += 1

        self.log = scrolledtext.ScrolledText(f, height=22, font=("Consolas", 9), wrap="none")
        self.log.grid(row=r, column=0, columnspan=3, sticky="nsew"); f.rowconfigure(r, weight=1)
        self.say("Needs internet the first time each PDB entry is used (cached afterwards).\nEvery run writes run.log into its results folder.\n")
        self.protocol("WM_DELETE_WINDOW", self.close)
        self.after(150, self.pump)

    # ---- helpers
    def say(self, t):
        self.log.insert("end", t); self.log.see("end")
        if self.outdir:
            try:
                with open(os.path.join(self.outdir, "gui_session.log"), "a", encoding="utf-8") as fh: fh.write(t)
            except Exception: pass
    def pick_csv(self):
        p = filedialog.askopenfilename(filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if p: self.csv.set(p); self.pdb.set("")
    def pick_models(self):
        ps = filedialog.askopenfilenames(filetypes=[("Models", "*.pdb *.cif *.mmcif *.ent *.gz"), ("All", "*.*")])
        if ps: self.models.set(";".join(ps)); self.pdb.set("")
    def bundled(self):
        self.csv.set(os.path.join(HERE, "representatives.csv")); self.pdb.set("")
        self.say("Bundled list = one best entry per virus (T=3, T=4, pseudo-T=3). A full run takes roughly an hour; results are saved as it goes.\n")
    def pick_out(self):
        p = filedialog.askdirectory()
        if p: self.base.set(p)
    def open_out(self):
        p = self.last_out or self.base.get(); os.makedirs(p, exist_ok=True); os.startfile(p)
    def open_summary(self):
        p = os.path.join(self.last_out or "", "gate_summary.txt")
        if os.path.exists(p): os.startfile(p)
        else: messagebox.showinfo("No summary yet", "gate_summary.txt appears in the results folder after the first structure is analysed.")
    def help(self):
        messagebox.showinfo("Help", "Gate residue = side-chain atom within the cutoff of a true 5-fold or 3-fold axis.\n"
            "gate_summary.txt lists, for every structure, the residues at each gate and the residue numbers near BOTH a 5-fold and a 3-fold axis.\n"
            "Pooled statistics and figures need 5 or more different viruses.\n"
            "Your own models need the icosahedral operators (REMARK 350 / mmCIF assembly) or the full capsid in the file.\n"
            "If something fails, run.log in the results folder says why.")
    def close(self):
        if self.proc and self.proc.poll() is None and not messagebox.askyesno("Run in progress", "Stop the running analysis and close?"):
            return
        self.stop(); self.destroy()

    # ---- run
    def run(self):
        ids = [x for x in self.pdb.get().replace(",", " ").split() if x]
        models = [m for m in self.models.get().split(";") if m.strip()]
        csv = self.csv.get().strip()
        if not (ids or models or csv):
            messagebox.showwarning("Input needed", "Enter PDB IDs, add model files, or choose a CSV list."); return
        base = self.base.get().strip() or RESULTS
        if self.resume.get():
            runs = sorted(glob.glob(os.path.join(base, "run_*")))
            if not runs: messagebox.showwarning("Nothing to resume", "No earlier run folder found in " + base); return
            out = runs[-1]
        else:
            out = os.path.join(base, "run_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        try:
            os.makedirs(out, exist_ok=True)
        except Exception as e:
            messagebox.showerror("Cannot create results folder", f"{out}\n\n{e}"); return
        self.outdir = self.last_out = out
        cmd = [PY, "-u", SCRIPT, "--cache", os.path.join(HERE, "pdbcache"), "--out", out,
               "--cutoff", self.cut.get().strip() or "8", "--atoms", self.atoms.get(), "--background", self.bg.get()]
        if ids: cmd += ["--pdb"] + ids
        if models: cmd += ["--model"] + models
        if csv: cmd += ["--csv", csv]
        if self.sens.get().strip(): cmd += ["--sens"] + self.sens.get().split()
        if self.res.get().strip(): cmd += ["--residue", self.res.get().strip()]
        if self.resume.get(): cmd += ["--resume"]
        self.say("\n$ " + " ".join(cmd) + "\n")
        self.status.set("Running ...")
        self.run_btn.config(state="disabled"); self.stop_btn.config(state="normal")
        threading.Thread(target=self.worker, args=(cmd, out), daemon=True).start()

    def worker(self, cmd, out):
        tail = []
        self.q.put("[environment check] importing numpy, scipy, pandas, statsmodels, matplotlib, gemmi, requests ...\n")
        good, rep, bad = gate_env.probe(PY)
        if not good:
            msg = (f"Environment check failed for {PY}\nFailing module: {bad}\n{rep}"
                   f"Fix: delete the .venv folder and start the launcher .bat again to rebuild it.\n")
            try:
                with open(os.path.join(out, "run.log"), "a", encoding="utf-8") as fh:
                    fh.write(f"gate_gui  {datetime.datetime.now():%Y-%m-%d %H:%M:%S}\n" + msg)
            except Exception:
                pass
            self.q.put(msg)
            self.q.put(("done", -2, False, out, msg))
            return
        self.q.put(rep + "[environment check] ok\n")
        try:
            self.proc = subprocess.Popen(cmd, cwd=HERE, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                         text=True, encoding="utf-8", errors="replace", creationflags=NOWIN, env=gate_env.child_env(PY))
            for line in self.proc.stdout:
                self.q.put(line); tail.append(line); tail = tail[-25:]
            rc = self.proc.wait()
        except Exception as e:
            rc = -1; tail.append(f"could not start the analysis: {e}\n"); self.q.put(tail[-1])
        ok = rc == 0 and os.path.exists(os.path.join(out, "gate_summary.txt"))
        self.q.put(("done", rc, ok, out, "".join(tail[-12:])))

    def stop(self):
        if self.proc and self.proc.poll() is None:
            self.proc.kill(); self.say("\n[stopped by user - structures finished so far are saved in the results folder]\n")

    def finished(self, rc, ok, out, tail):
        self.run_btn.config(state="normal"); self.stop_btn.config(state="disabled")
        n = len([f for f in os.listdir(out) if os.path.isfile(os.path.join(out, f))]) if os.path.isdir(out) else 0
        if ok:
            self.status.set(f"Finished. {n} files saved in {out}")
            self.say(f"\n[finished OK - {n} files in {out}]\n")
            if messagebox.askyesno("Finished", f"Results saved in:\n{out}\n\nOpen gate_summary.txt now?"):
                os.startfile(os.path.join(out, "gate_summary.txt"))
        else:
            self.status.set(f"FAILED (exit code {rc}). See run.log in {out}")
            self.say(f"\n[FAILED, exit code {rc}] details: {os.path.join(out, 'run.log')}\n")
            messagebox.showerror("Python environment is broken - analysis not started" if rc == -2 else "Analysis did not finish", f"Exit code {rc}. Files saved so far: {n}.\n\nFull log: {os.path.join(out, 'run.log')}\n\nLast lines:\n{tail}")

    def check(self):
        self.say("\n--- setup check ---\n")
        threading.Thread(target=self.check_worker, args=(self.base.get().strip() or RESULTS,), daemon=True).start()

    def check_worker(self, base):
        P = self.q.put
        P(f"python used for the analysis: {PY}\n")
        P(gate_env.probe(PY)[1])
        try:
            import urllib.request
            urllib.request.urlopen(urllib.request.Request("https://files.rcsb.org/download/4FTB.cif.gz", method="HEAD"), timeout=20)
            P("  ok   internet: files.rcsb.org reachable\n")
        except Exception as e:
            P(f"  FAIL internet: files.rcsb.org not reachable ({e}) - PDB IDs cannot be downloaded (local model files still work)\n")
        try:
            os.makedirs(base, exist_ok=True)
            t = os.path.join(base, "_write_test.txt"); open(t, "w").write("x"); os.remove(t)
            P(f"  ok   can write to {base}\n")
        except Exception as e:
            P(f"  FAIL cannot write to the results folder: {e}\n")
        P("--- end of check ---\n")

    def pump(self):
        try:
            while True:
                m = self.q.get_nowait()
                if isinstance(m, tuple): self.finished(*m[1:])
                else:
                    self.say(m)
                    g = re.match(r"\[(\d+)/(\d+)\] (\S+)", m)
                    if g: self.status.set(f"Running: structure {g.group(1)} of {g.group(2)} ({g.group(3)})")
        except queue.Empty: pass
        self.after(150, self.pump)


if __name__ == "__main__":
    App().mainloop()
