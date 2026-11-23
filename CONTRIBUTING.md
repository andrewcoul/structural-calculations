# Contributing

This repository (`structural-calculations`) builds standalone Google Colab calculation notebooks from plain Python sources and a shared `structcalc` library.

## Module layout

Each calculation sheet consists of four files:

| File | Purpose |
| --- | --- |
| `lib/structcalc/calcs/<name>.py` | Engineering module exposing `run(**inputs) -> dict` with at least `report` (`structcalc.Report`). |
| `src/NN_<name>.py` | Percent-format notebook source with Colab form cells (`#@title`, `#@param`), a **Calculate** cell, an optional plot cell, and an optional HTML export cell. |
| `tests/test_<name>.py` | `pytest` verification suite with at least one hand-checked or published example per primary limit state. |
| `notebooks/NN_<name>.ipynb` | Generated standalone notebook built by `python build.py`. Never edit `.ipynb` files by hand. |

## Conventions

- **Units:** US customary throughout (`kip`, `ft`, `in`, `ksi`/`psi`, `plf`/`psf`). State units in every form helper line and report row.
- **Code basis:** IBC 2024, ASCE 7-22, AISC 360-22, NDS 2024 + NDS Supplement 2024, SDPWS 2021, ACI 318-19.
- **Checks:** Use `report.check(name, demand, capacity, unit, ref="...")`. Ratios `<= 1.0` pass (`OK`).
- **Intermediate values:** Record key intermediate quantities with `report.result(name, value, unit, ref="...")`.
- **Code tables:** Tabulated values not in the shared data files live in the module as small dicts with table references in comments, plus form overrides where appropriate.
- **Assumptions & limits:** Document scope limits with `report.note(...)`.
- **Colab forms:**
  - Single-line string inputs cannot hold newlines; separate list items with `;` or `,`.
  - Numeric dropdowns **must** include `{type:"raw"}` (for example, `bar_bot = 8  #@param [3, 4, 5, 6, 7, 8] {type:"raw"}`) so Colab preserves integer/float types when edited.
- **Notebook copy:** The top markdown cell contains `# Title`, a bold code-basis line, a concise scope summary, the **How to use** line, and the engineering disclaimer.
- **Shared core modules** (review carefully before modifying): `report.py`, `beam.py`, `beaminput.py`, `loads.py`, `steel.py`, `concrete.py`, `concrete_member.py`, `foundation.py`, `wood.py`, `wood_conn.py`, `build.py`.

## Development commands

```bash
.venv/bin/python -m pytest -q tests          # run unit tests
.venv/bin/python build.py                    # rebuild all notebooks (or: build.py 01 05)
.venv/bin/python tests/run_notebooks.py      # execute built notebooks locally
```
