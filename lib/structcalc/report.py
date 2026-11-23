"""Calculation report: inputs, intermediate results, and code checks.

Usage
-----
    r = Report("Steel Beam", code="AISC 360-22 (LRFD)")
    r.input("Span", 24, "ft")
    r.result("Mu", 180.2, "kip-ft", ref="ASCE 7-22 2.3.1 combo 2")
    r.check("Flexure", demand=180.2, capacity=252.0, unit="kip-ft", ref="AISC F2")
    r.show()

A check passes when demand <= capacity (ratio <= 1.0). For checks where
"bigger is better" (e.g. a factor of safety), pass ``demand`` as the
required value and ``capacity`` as the provided value.
"""

from __future__ import annotations

import html
import math
from dataclasses import dataclass, field


def _fmt(v, digits=4):
    if v is None:
        return ""
    if isinstance(v, str):
        return v
    if isinstance(v, bool):
        return "Yes" if v else "No"
    if isinstance(v, int):
        return f"{v:,d}"
    try:
        x = float(v)
    except (TypeError, ValueError):
        return str(v)
    if math.isnan(x):
        return "n/a"
    if math.isinf(x):
        return "inf"
    if abs(x) < 1e-9:
        return "0"
    ax = abs(x)
    if ax >= 1e5 or ax < 1e-3:
        return f"{x:.{digits - 1}e}"
    # significant-figure style formatting
    decimals = max(0, digits - 1 - int(math.floor(math.log10(ax))))
    s = f"{x:,.{decimals}f}"
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s


@dataclass
class Item:
    name: str
    value: object
    unit: str = ""
    ref: str = ""
    note: str = ""


@dataclass
class Check:
    name: str
    demand: float
    capacity: float
    unit: str = ""
    ref: str = ""
    note: str = ""

    @property
    def ratio(self) -> float:
        if self.capacity == 0:
            return math.inf if self.demand > 0 else 0.0
        return abs(self.demand) / abs(self.capacity)

    @property
    def ok(self) -> bool:
        return self.ratio <= 1.0 + 1e-9


@dataclass
class Report:
    title: str
    code: str = ""
    project: str = ""
    inputs: list = field(default_factory=list)
    results: list = field(default_factory=list)
    checks: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    sections: list = field(default_factory=list)

    # -- builders -----------------------------------------------------
    def input(self, name, value, unit="", note=""):
        self.inputs.append(Item(name, value, unit, note=note))
        return value

    def result(self, name, value, unit="", ref="", note=""):
        self.results.append(Item(name, value, unit, ref, note))
        return value

    def check(self, name, demand, capacity, unit="", ref="", note=""):
        c = Check(name, float(demand), float(capacity), unit, ref, note)
        self.checks.append(c)
        return c

    def note(self, text):
        self.notes.append(text)

    def section(self, title, rows, columns=None):
        """Add an extra table (list of dicts or list of lists)."""
        self.sections.append((title, rows, columns))

    # -- summary --------------------------------------------------------
    @property
    def max_ratio(self) -> float:
        return max((c.ratio for c in self.checks), default=0.0)

    @property
    def ok(self) -> bool:
        return all(c.ok for c in self.checks)

    def governing(self):
        if not self.checks:
            return None
        return max(self.checks, key=lambda c: c.ratio)

    # -- output ---------------------------------------------------------
    def to_text(self) -> str:
        out = []
        bar = "=" * 78
        out += [bar, self.title.upper()]
        if self.code:
            out.append(f"Code basis: {self.code}")
        if self.project:
            out.append(f"Project: {self.project}")
        out.append(bar)
        if self.inputs:
            out.append("INPUTS")
            for i in self.inputs:
                out.append(f"  {i.name:<38} {_fmt(i.value):>14} {i.unit}")
        if self.results:
            out.append("RESULTS")
            for i in self.results:
                ref = f"  [{i.ref}]" if i.ref else ""
                out.append(f"  {i.name:<38} {_fmt(i.value):>14} {i.unit}{ref}")
        for title, rows, columns in self.sections:
            out.append(title.upper())
            rows, columns = _normalize_rows(rows, columns)
            out.append("  " + " | ".join(str(c) for c in columns))
            for r in rows:
                out.append("  " + " | ".join(_fmt(v) for v in r))
        if self.checks:
            out.append("CODE CHECKS")
            out.append(f"  {'Check':<32} {'Demand':>11} {'Capacity':>11} {'Ratio':>7}  Status")
            for c in self.checks:
                out.append(
                    f"  {c.name:<32} {_fmt(c.demand):>11} {_fmt(c.capacity):>11} "
                    f"{c.ratio:7.3f}  {'OK' if c.ok else 'NG'}  {c.unit} [{c.ref}]"
                )
            g = self.governing()
            out.append(bar)
            out.append(
                f"SUMMARY: {'PASS' if self.ok else 'FAIL'} - governing: {g.name} "
                f"(ratio {g.ratio:.3f})"
            )
        if self.notes:
            out.append("NOTES")
            for n in self.notes:
                out.append(f"  - {n}")
        out.append(bar)
        return "\n".join(out)

    def to_html(self) -> str:
        e = html.escape
        css = (
            "<style>.sc{font-family:system-ui,Segoe UI,Arial,sans-serif;font-size:13px;"
            "max-width:980px}.sc h2{margin:4px 0;font-size:18px}.sc h3{margin:14px 0 4px;"
            "font-size:14px;border-bottom:1px solid #8884}.sc table{border-collapse:collapse;"
            "width:100%}.sc td,.sc th{padding:3px 8px;border-bottom:1px solid #8882;"
            "text-align:left}.sc td.n{text-align:right;font-variant-numeric:tabular-nums}"
            ".sc .ok{color:#0a7d32;font-weight:600}.sc .ng{color:#c62828;font-weight:700}"
            ".sc .ref{color:#888;font-size:12px}.sc .banner{padding:8px 12px;border-radius:6px;"
            "margin:10px 0;font-weight:600}.sc .bok{background:#0a7d3222;border:1px solid #0a7d32}"
            ".sc .bng{background:#c6282822;border:1px solid #c62828}</style>"
        )
        h = [css, '<div class="sc">', f"<h2>{e(self.title)}</h2>"]
        if self.code:
            h.append(f'<div class="ref">Code basis: {e(self.code)}</div>')
        if self.project:
            h.append(f'<div class="ref">Project: {e(self.project)}</div>')
        if self.checks:
            g = self.governing()
            cls = "bok" if self.ok else "bng"
            word = "PASS" if self.ok else "FAIL"
            h.append(
                f'<div class="banner {cls}">{word} &mdash; governing check: {e(g.name)} '
                f"(ratio {g.ratio:.3f})</div>"
            )

        def items_table(title, items, with_ref):
            if not items:
                return
            h.append(f"<h3>{e(title)}</h3><table>")
            for i in items:
                ref = f'<td class="ref">{e(i.ref)}</td>' if with_ref else ""
                note = f' <span class="ref">{e(i.note)}</span>' if i.note else ""
                h.append(
                    f"<tr><td>{e(i.name)}{note}</td><td class='n'>{e(_fmt(i.value))}</td>"
                    f"<td>{e(i.unit)}</td>{ref}</tr>"
                )
            h.append("</table>")

        items_table("Inputs", self.inputs, False)
        items_table("Results", self.results, True)
        for title, rows, columns in self.sections:
            rows, columns = _normalize_rows(rows, columns)
            h.append(f"<h3>{e(title)}</h3><table><tr>")
            h += [f"<th>{e(str(c))}</th>" for c in columns]
            h.append("</tr>")
            for r in rows:
                h.append("<tr>" + "".join(f"<td class='n'>{e(_fmt(v))}</td>" for v in r) + "</tr>")
            h.append("</table>")
        if self.checks:
            h.append(
                "<h3>Code checks</h3><table><tr><th>Check</th><th>Demand</th>"
                "<th>Capacity</th><th>Unit</th><th>Ratio</th><th>Status</th><th>Reference</th></tr>"
            )
            for c in self.checks:
                st = '<span class="ok">OK</span>' if c.ok else '<span class="ng">NG</span>'
                note = f' <span class="ref">{e(c.note)}</span>' if c.note else ""
                h.append(
                    f"<tr><td>{e(c.name)}{note}</td><td class='n'>{e(_fmt(c.demand))}</td>"
                    f"<td class='n'>{e(_fmt(c.capacity))}</td><td>{e(c.unit)}</td>"
                    f"<td class='n'>{c.ratio:.3f}</td><td>{st}</td>"
                    f"<td class='ref'>{e(c.ref)}</td></tr>"
                )
            h.append("</table>")
        if self.notes:
            h.append("<h3>Notes</h3><ul>")
            h += [f"<li>{e(n)}</li>" for n in self.notes]
            h.append("</ul>")
        h.append(
            '<div class="ref" style="margin-top:10px">Design aid only. A licensed engineer '
            "must verify all inputs and results.</div></div>"
        )
        return "".join(h)

    def _repr_html_(self):
        return self.to_html()

    def show(self):
        try:
            from IPython import get_ipython
            from IPython.display import HTML, display

            if get_ipython() is not None:
                display(HTML(self.to_html()))
                return
        except ImportError:
            pass
        print(self.to_text())

    def save_html(self, path):
        with open(path, "w", encoding="utf-8") as f:
            f.write(f"<!doctype html><meta charset='utf-8'><title>{html.escape(self.title)}"
                    f"</title>{self.to_html()}")
        return path

    def as_dict(self):
        return {
            "title": self.title,
            "ok": self.ok,
            "max_ratio": self.max_ratio,
            "checks": {c.name: {"demand": c.demand, "capacity": c.capacity,
                                "ratio": c.ratio, "ok": c.ok} for c in self.checks},
            "results": {i.name: i.value for i in self.results},
        }


def _normalize_rows(rows, columns):
    rows = list(rows)
    if rows and isinstance(rows[0], dict):
        columns = columns or list(rows[0].keys())
        rows = [[r.get(c, "") for c in columns] for r in rows]
    return rows, (columns or [])
