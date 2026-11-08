"""
3D view of a capsid with its axis residues coloured, using the PDBe build of Mol* (the viewer the PDB itself uses).
The structure is read from the RCSB PDB as BinaryCIF, so a whole icosahedral capsid draws in a few seconds.

Residue numbers and chains are the author ones (auth_asym_id / auth_seq_id) stored in gate_residues.sqlite. In the
whole-capsid view Mol* applies a selection to every symmetry copy of a chain, so one residue lights up at all of its
copies and the rings around each 5-fold and 3-fold opening become visible.
"""
import json

MOLSTAR = "https://cdn.jsdelivr.net/npm/pdbe-molstar@3.12.0/build"
COLOURS = {"5-fold": "#1d5fa8", "3-fold": "#e0662f", "2-fold": "#8a96a3", "both": "#e8b100"}


def html(pdb, residues, whole_capsid=True, height=620):
    """residues: iterable of dicts with chain, auth_seq, gate ('5-fold'/'3-fold'/'2-fold'/'both')"""
    sel, seen = [], set()
    for r in residues:
        key = (r["chain"], int(r["auth_seq"]))
        if key in seen:
            continue
        seen.add(key)
        sel.append({"auth_asym_id": r["chain"], "auth_seq_id": int(r["auth_seq"]), "color": COLOURS[r["gate"]],
                    "representation": "spacefill", "representationColor": COLOURS[r["gate"]], "sideChain": True})
    options = {
        "moleculeId": pdb.lower(),
        # whole-entry BinaryCIF from RCSB (has the assembly operators); the PDBe model server is the plugin default
        # but is sometimes unavailable, and this file is small and served with CORS
        "customData": {"url": f"https://models.rcsb.org/{pdb.lower()}.bcif", "format": "cif", "binary": True},
        "bgColor": "#ffffff",
        "hideControls": True,
        "sequencePanel": False,
        "pdbeLink": False,
        "landscape": True,
        "visualStyle": {"polymer": {"type": "cartoon"}},
        "hideCanvasControls": ["animation", "trajectory", "controlInfo"],
    }
    if whole_capsid:
        options["assemblyId"] = "preferred"
    return f"""
<link rel="stylesheet" href="{MOLSTAR}/pdbe-molstar-light.css">
<script src="{MOLSTAR}/pdbe-molstar-plugin.js"></script>
<div id="msg" style="font:13px 'Segoe UI',sans-serif;color:#5b7189;padding:4px 2px">
  Loading {pdb.upper()} from the PDB{' (whole capsid, can take 10-30 s)' if whole_capsid else ''} ...</div>
<div id="view" style="position:relative;width:100%;height:{height - 30}px"></div>
<script>
  const plugin = new PDBeMolstarPlugin();
  const sel = {json.dumps(sel)};
  plugin.render(document.getElementById('view'), {json.dumps(options)});
  plugin.events.loadComplete.subscribe(async (ok) => {{
    const msg = document.getElementById('msg');
    if (!ok) {{ msg.textContent = 'Could not load {pdb.upper()} from the PDB. Try again in a minute.'; return; }}
    if (sel.length) await plugin.visual.select({{ data: sel, nonSelectedColor: '#cfd6de' }});
    msg.textContent = 'Drag to rotate, scroll to zoom, right-drag to move. Click a residue to see its name.';
  }});
</script>
"""
