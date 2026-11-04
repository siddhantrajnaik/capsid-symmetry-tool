"""Interpreter selection and pre-flight import probe shared by gate_explorer.py and gate_gui.py.

The probe runs in a child process of the SAME interpreter that will run the analysis. It imports every native
dependency and makes one small BLAS/LAPACK call (a bare `import numpy` can succeed while the first matrix product
kills the process, which is how the 0xC06D007F failures looked). A native crash therefore shows up here as a
non-zero exit code plus the name of the module that was being loaded, instead of as a run that dies silently.
"""
import os, sys, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
NOWIN = 0x08000000 if os.name == "nt" else 0               # CREATE_NO_WINDOW

_PROBE = r"""
import sys
def step(name, fn):
    print('TRY ' + name, flush=True)
    r = fn()
    print('OK  ' + name + ' ' + str(r or ''), flush=True)
def numpy_():
    import numpy as np
    a = np.random.default_rng(0).random((6, 6)); b = a @ a
    np.linalg.inv(a + 6 * np.eye(6)); np.linalg.svd(b)
    return np.__version__
def scipy_():
    import scipy, scipy.linalg, scipy.spatial, scipy.stats
    scipy.linalg.solve([[3., 1.], [1., 2.]], [1., 2.])
    return scipy.__version__
def pandas_():
    import pandas as pd
    pd.DataFrame({'a': [1., 2.]}).corr()
    return pd.__version__
def statsmodels_():
    import statsmodels
    from statsmodels.stats.contingency_tables import StratifiedTable
    from statsmodels.stats.multitest import multipletests
    multipletests([0.01, 0.02, 0.5], method='fdr_bh')
    return statsmodels.__version__
def matplotlib_():
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure(); plt.close('all')
    return matplotlib.__version__
def gemmi_():
    import gemmi
    return gemmi.__version__
def requests_():
    import requests
    return requests.__version__
for n, f in [('numpy', numpy_), ('scipy', scipy_), ('pandas', pandas_), ('statsmodels', statsmodels_),
             ('matplotlib', matplotlib_), ('gemmi', gemmi_), ('requests', requests_)]:
    try:
        step(n, f)
    except Exception as e:
        print('FAIL ' + n + ' ' + type(e).__name__ + ': ' + str(e), flush=True)
        sys.exit(2)
"""


def find_python():
    """the interpreter that does the work: the tool's own .venv if it exists, else the one running this program"""
    cand = os.path.join(HERE, ".venv", "Scripts", "python.exe")
    if os.path.exists(cand):
        return cand
    py = sys.executable
    if py.lower().endswith("pythonw.exe"):                  # the analysis runs under python.exe (console, piped output)
        c = py[:-len("pythonw.exe")] + "python.exe"
        py = c if os.path.exists(c) else py
    return py


def child_env(py):
    """environment for child processes: put the interpreter's own DLL folders first on PATH. A conda env that is
    started without `conda activate` otherwise cannot resolve forwarded DLL imports (libiomp5md.dll -> libomp.dll)."""
    env = dict(os.environ)
    root = os.path.dirname(os.path.abspath(py))
    extra = [d for d in (root, os.path.join(root, "Scripts"), os.path.join(root, "Library", "bin"),
                         os.path.join(root, "Library", "mingw-w64", "bin"), os.path.join(root, "Library", "usr", "bin"))
             if os.path.isdir(d)]
    env["PATH"] = os.pathsep.join(extra + [env.get("PATH", "")])
    return env


def probe(py, timeout=120):
    """-> (ok, report_text, failing_module_or_None)"""
    try:
        r = subprocess.run([py, "-u", "-c", _PROBE], capture_output=True, text=True, timeout=timeout, errors="replace",
                           stdin=subprocess.DEVNULL, creationflags=NOWIN, env=child_env(py))
    except Exception as e:
        return False, f"could not start {py}: {e}\n", "python"
    lines = [l for l in r.stdout.splitlines() if l.strip()]
    report = "".join(("  ok   " + l[4:] if l.startswith("OK  ") else "  FAIL " + l[5:] if l.startswith("FAIL ") else "")
                     + ("\n" if l.startswith(("OK  ", "FAIL ")) else "") for l in lines)
    if r.returncode == 0:
        return True, report, None
    bad = None
    for l in lines:
        if l.startswith("FAIL "):
            bad = l.split()[1]
    if bad is None:                                         # process died without a Python exception
        tried = [l[4:].strip() for l in lines if l.startswith("TRY ")]
        done = {l[4:].split()[0] for l in lines if l.startswith("OK  ")}
        pend = [t for t in tried if t not in done]
        bad = pend[-1] if pend else "python"
        report += (f"  FAIL {bad}: the interpreter died with exit code {r.returncode & 0xFFFFFFFF} (0x{r.returncode & 0xFFFFFFFF:08X}) "
                   f"while loading or first using it - a native DLL problem, not a Python error.\n")
    tail = r.stderr.strip()[-600:]
    if tail:
        report += "  stderr: " + tail.replace("\n", "\n          ") + "\n"
    return False, report, bad
