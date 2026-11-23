"""Execute notebooks locally with nbclient; report errors and outputs.

Usage: .venv/bin/python tests/run_notebooks.py [prefix ...]
"""
import pathlib, sys, time
import nbformat
from nbclient import NotebookClient

ROOT = pathlib.Path(__file__).resolve().parents[1]
nbs = sorted((ROOT / "notebooks").glob("*.ipynb"))
if sys.argv[1:]:
    nbs = [n for n in nbs if any(n.name.startswith(a) for a in sys.argv[1:])]
bad = 0
for p in nbs:
    nb = nbformat.read(p, 4)
    t = time.time()
    try:
        NotebookClient(nb, timeout=600, kernel_name="python3").execute()
        status = "OK"
    except Exception as e:  # noqa: BLE001
        status = "ERROR " + str(e).strip().splitlines()[-1][:300]
        bad += 1
    print(f"{p.name:<40} {status} ({time.time() - t:.1f}s)")
sys.exit(1 if bad else 0)
