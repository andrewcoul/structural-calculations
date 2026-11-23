#!/usr/bin/env python3
"""Build standalone Colab notebooks from src/*.py.

Source format ("percent" cells):

    # %% [markdown]
    # # Title
    # text...

    # %%
    span = 24.0  #@param {type:"number"}

Each notebook gets a hidden "Engine" cell (inserted after the first
markdown cell) that embeds the structcalc library as a zip, so the
notebook runs in Colab with no Drive mount and no pip install.

Usage:
    python build.py              # build all (empty outputs)
    python build.py --execute    # build and save default example outputs
    python build.py 01 05        # build sources whose name starts with 01 or 05
"""

from __future__ import annotations

import base64
import hashlib
import io
import json
import pathlib
import sys
import zipfile

ROOT = pathlib.Path(__file__).resolve().parent
SRC = ROOT / "src"
OUT = ROOT / "notebooks"
LIB = ROOT / "lib" / "structcalc"

STEEL_DATA = {"W.csv", "HSS_rect.csv", "HSS_round.csv"}


def lib_zip(include_steel_data: bool) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for p in sorted(LIB.rglob("*")):
            if p.is_dir() or "__pycache__" in p.parts:
                continue
            if p.name in STEEL_DATA and not include_steel_data:
                continue
            # fixed timestamp -> reproducible output
            info = zipfile.ZipInfo(str(p.relative_to(LIB.parent)), date_time=(2026, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, p.read_bytes())
    return buf.getvalue()


def engine_cell(source_text: str) -> str:
    blob = lib_zip(True)  # steel shape data is small (~40 KB zipped); always embed
    digest = hashlib.sha1(blob).hexdigest()[:10]
    b64 = base64.b64encode(blob).decode()
    lines = [b64[i:i + 4000] for i in range(0, len(b64), 4000)]
    zip_lit = "(\n" + "\n".join(f'    "{ln}"' for ln in lines) + "\n)"
    return f'''#@title Engine — run this cell first (embedded structcalc library) {{ display-mode: "form" }}
# Unpacks the embedded structcalc library and imports it (no pip install or Drive mount needed).
# Source: structural-calculations (lib/structcalc). Rebuild notebooks with: python build.py
import base64 as _b64, io as _io, os as _os, sys as _sys, tempfile as _tf, zipfile as _zf
_ZIP = {zip_lit}
_dir = _os.path.join(_tf.gettempdir(), "structcalc_{digest}")
if not _os.path.isdir(_os.path.join(_dir, "structcalc")):
    _zf.ZipFile(_io.BytesIO(_b64.b64decode(_ZIP))).extractall(_dir)
for _m in [m for m in _sys.modules if m == "structcalc" or m.startswith("structcalc.")]:
    del _sys.modules[_m]
if _dir not in _sys.path:
    _sys.path.insert(0, _dir)
import structcalc
print(f"structcalc {{structcalc.__version__}} ready (build {digest})")
'''


def parse_cells(text: str):
    cells = []
    cur_type, cur = None, []

    def flush():
        if cur_type is None:
            return
        body = "\n".join(cur).strip("\n")
        if cur_type == "markdown":
            md = []
            for ln in body.splitlines():
                if ln.startswith("# "):
                    md.append(ln[2:])
                elif ln.strip() == "#":
                    md.append("")
                else:
                    md.append(ln)
            body = "\n".join(md)
        if body.strip():
            cells.append((cur_type, body))

    for ln in text.splitlines():
        if ln.startswith("# %%"):
            flush()
            cur_type = "markdown" if "[markdown]" in ln else "code"
            cur = []
        else:
            if cur_type is None:
                continue  # header lines before the first cell are ignored
            cur.append(ln)
    flush()
    return cells


def to_source_lines(body: str):
    lines = body.splitlines(keepends=True)
    if lines and lines[-1].endswith("\n"):
        lines[-1] = lines[-1].rstrip("\n")
    return lines


def build_one(src: pathlib.Path) -> pathlib.Path:
    text = src.read_text(encoding="utf-8")
    cells = parse_cells(text)
    title = src.stem
    for t, body in cells:
        if t == "markdown":
            for ln in body.splitlines():
                if ln.startswith("# "):
                    title = ln[2:].strip()
                    break
            break
    nb_cells = []
    inserted = False
    for t, body in cells:
        if t == "markdown":
            nb_cells.append({"cell_type": "markdown", "metadata": {}, "source": to_source_lines(body)})
        else:
            meta = {}
            if "#@title" in body.splitlines()[0]:
                meta = {"cellView": "form"}
            nb_cells.append({"cell_type": "code", "metadata": meta, "execution_count": None,
                             "outputs": [], "source": to_source_lines(body)})
        if not inserted and t == "markdown":
            nb_cells.append({"cell_type": "code", "metadata": {"cellView": "form"},
                             "execution_count": None, "outputs": [],
                             "source": to_source_lines(engine_cell(text))})
            inserted = True
    nb = {
        "nbformat": 4,
        "nbformat_minor": 0,
        "metadata": {
            "colab": {"name": f"{src.stem}.ipynb", "provenance": [], "toc_visible": True},
            "kernelspec": {"name": "python3", "display_name": "Python 3"},
            "language_info": {"name": "python"},
        },
        "cells": nb_cells,
    }
    OUT.mkdir(exist_ok=True)
    out = OUT / f"{src.stem}.ipynb"
    out.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    return out


def execute_notebook(path: pathlib.Path) -> None:
    import nbformat
    from nbclient import NotebookClient

    nb = nbformat.read(path, as_version=4)
    NotebookClient(nb, timeout=300, kernel_name="python3").execute()
    # Strip transient execution timestamps for clean git diffs
    for cell in nb.cells:
        if cell.cell_type == "code":
            cell.metadata.pop("execution", None)
    nbformat.write(nb, path)


def main(argv):
    execute = "--execute" in argv
    prefixes = [a for a in argv if a != "--execute"]
    srcs = sorted(SRC.glob("*.py"))
    if prefixes:
        srcs = [s for s in srcs if any(s.name.startswith(a) for a in prefixes)]
    for s in srcs:
        out = build_one(s)
        if execute:
            execute_notebook(out)
        print(f"built {out.relative_to(ROOT)} ({out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main(sys.argv[1:])
