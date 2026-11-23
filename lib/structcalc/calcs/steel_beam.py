"""Steel beam: AISC 360-22 flexure + shear, ASCE 7-22 combinations,
IBC 1604.3 deflection. W shapes and rectangular/round HSS, strong axis.
"""

from __future__ import annotations

import numpy as np

from .. import loads, steel
from ..beam import Beam
from ..beaminput import apply_loads, deflection_checks, make_supports, parse_loads
from ..report import Report

BRACING = ["Continuous (both flanges)", "Top flange continuous", "At supports only", "At interval"]


def _brace_points(length, supports, mode, interval):
    xs = {0.0, length} | {x for x, _ in supports}
    if mode == "At interval" and interval > 0:
        n = int(np.floor(length / interval + 1e-9))
        xs |= {i * interval for i in range(n + 1)}
    return sorted(x for x in xs if 0 <= x <= length + 1e-9)


def _segment_ratios(res, p, Fy, method, brace_pts, mode, supports, length):
    """Return list of segment dicts with flexural capacity and ratio."""
    phi = steel.factor("b", method)
    sup_x = sorted(x for x, _ in supports)
    single_fixed = len(supports) == 1
    rows = []
    Mn0 = steel.flexure(p, Fy, "x", 0.0, 1.0).Mn
    for a, b in zip(brace_pts[:-1], brace_pts[1:]):
        if b - a < 1e-6:
            continue
        m = (res.x >= a - 1e-9) & (res.x <= b + 1e-9)
        Mseg = res.M[m]
        xseg = res.x[m]
        Mpos = max(float(Mseg.max()), 0.0)
        Mneg = max(float(-Mseg.min()), 0.0)
        Lb = (b - a) * 12
        # free end segment of a cantilever / overhang -> Cb = 1.0
        free = (a < sup_x[0] - 1e-9) or (b > sup_x[-1] + 1e-9) or single_fixed
        if mode == "Continuous (both flanges)":
            Lb_pos = Lb_neg = 0.0
        elif mode == "Top flange continuous":
            Lb_pos, Lb_neg = 0.0, Lb
        else:
            Lb_pos = Lb_neg = Lb
        Mabs = np.abs(Mseg)
        if free:
            Cb = 1.0
        else:
            def at(f):
                return float(np.abs(np.interp(a + f * (b - a), xseg, Mseg)))
            Cb = steel.Cb_factor(at(0.25), at(0.5), at(0.75), float(Mabs.max()))
        fpos = steel.flexure(p, Fy, "x", Lb_pos, Cb if Lb_pos else 1.0)
        fneg = steel.flexure(p, Fy, "x", Lb_neg, Cb if Lb_neg else 1.0)
        rpos = Mpos * 12 / (phi * fpos.Mn) if Mpos > 0 else 0.0
        rneg = Mneg * 12 / (phi * fneg.Mn) if Mneg > 0 else 0.0
        if rpos >= rneg:
            rows.append(dict(a=a, b=b, Mu=Mpos, Lb=Lb_pos / 12, Cb=Cb, cap=phi * fpos.Mn / 12,
                             ratio=rpos, ls=fpos.limit_state, sign="+"))
        else:
            rows.append(dict(a=a, b=b, Mu=Mneg, Lb=Lb_neg / 12, Cb=Cb, cap=phi * fneg.Mn / 12,
                             ratio=rneg, ls=fneg.limit_state, sign="-"))
    if not rows:
        rows.append(dict(a=0, b=length, Mu=0, Lb=0, Cb=1, cap=phi * Mn0 / 12, ratio=0,
                         ls="-", sign="+"))
    return rows


def analyze(section, Fy=50.0, method="LRFD", span_type="Simple span", span=20.0,
            overhang_left=0.0, overhang_right=0.0, custom_supports="",
            w=None, extra_loads="", self_weight=True, bracing="At supports only",
            brace_interval=0.0, live_defl=360, total_defl=240, E=steel.E_STEEL):
    """Analyze one section. Returns dict(report, governing, combos, ...)."""
    w = dict(w or {})
    p = steel.shape(section)
    length, supports = make_supports(span_type, span, overhang_left, overhang_right, custom_supports)
    beam = Beam(length, E, p["Ix"], supports)
    apply_loads(beam, w, parse_loads(extra_loads))
    if self_weight:
        beam.self_weight = p["weight"] / 1000.0
    per_case = beam.solve()
    combos = loads.lrfd() if method.upper() == "LRFD" else loads.asd()
    results = beam.combos(combos)
    brace_pts = _brace_points(length, supports, bracing, brace_interval)
    phiv_Vn, Vn, vnote = steel.shear(p, Fy, "x", method)

    r = Report(f"Steel Beam - {p['name']}", code=f"AISC 360-22 ({method.upper()}), ASCE 7-22, IBC 2024")
    r.input("Section", p["name"])
    r.input("Fy", Fy, "ksi")
    r.input("E", E, "ksi")
    r.input("Span type", span_type)
    r.input("Total length", length, "ft")
    r.input("Supports", ", ".join(f"{x:g} ft {t}" for x, t in supports))
    for c, v in w.items():
        if v:
            r.input(f"Uniform load {c}", v * 1000, "plf")
    if extra_loads.strip():
        r.input("Extra loads", "; ".join(ln.strip() for ln in extra_loads.splitlines() if ln.strip()))
    r.input("Self weight included", bool(self_weight))
    r.input("Bracing", bracing + (f" @ {brace_interval:g} ft" if bracing == "At interval" else ""))

    gov = {"flex": None, "shear": None}
    table = []
    for name, f, res in results:
        segs = _segment_ratios(res, p, Fy, method, brace_pts, bracing, supports, length)
        s = max(segs, key=lambda d: d["ratio"])
        V = res.Vabs
        table.append({"Combination": name, "M+ (k-ft)": res.Mmax, "M- (k-ft)": res.Mmin,
                      "V (k)": V, "Flex ratio": s["ratio"], "Shear ratio": V / phiv_Vn})
        if gov["flex"] is None or s["ratio"] > gov["flex"][1]["ratio"]:
            gov["flex"] = (name, s, res)
        if gov["shear"] is None or V > gov["shear"][1]:
            gov["shear"] = (name, V, res)

    fname, seg, fres = gov["flex"]
    word = "phi" if method.upper() == "LRFD" else "1/Omega"
    r.result("Ix", p["Ix"], "in^4")
    r.result("Zx", p["Zx"], "in^3")
    r.result("Sx", p["Sx"], "in^3")
    Mp = Fy * p["Zx"] / 12
    r.result(f"{word}*Mp", steel.factor("b", method) * Mp, "kip-ft", "AISC F2.1")
    f0 = steel.flexure(p, Fy, "x", 0.0)
    if p["kind"] == "W":
        fx = steel.flexure_W_major(p, Fy, 0.0)
        r.result("Lp", fx.Lp / 12, "ft", "AISC F2-5")
        r.result("Lr", fx.Lr / 12, "ft", "AISC F2-6")
    r.result("Section compactness", f0.detail)
    r.result("Governing flexure combination", fname)
    r.result("Governing segment", f"{seg['a']:g} to {seg['b']:g} ft ({seg['sign']} moment)")
    r.result("Lb (governing segment)", seg["Lb"], "ft")
    r.result("Cb", seg["Cb"], "", "AISC F1-1")
    r.result("Flexure limit state", seg["ls"])
    r.check("Flexure", seg["Mu"], seg["cap"], "kip-ft", "AISC F2/F3/F7/F8",
            note=f"[{fname}]")
    sname, Vu, sres = gov["shear"]
    r.check("Shear", Vu, phiv_Vn, "kip", "AISC G2/G4/G5", note=f"[{sname}] {vnote}")

    deflection_checks(r, per_case, length, supports, [
        ("Live (L)", {"L": 1.0}, live_defl),
        ("Roof live/snow (Lr or S)", {"Lr": 1.0, "S": 1.0}, live_defl),
        ("Total (D+L+Lr+S)", {"D": 1.0, "L": 1.0, "Lr": 1.0, "S": 1.0}, total_defl),
    ])
    reac_rows = []
    for x in sorted(per_case[next(iter(per_case))].reactions):
        row = {"Support x (ft)": x}
        for c, res in per_case.items():
            row[f"R {c} (k)"] = res.reactions[x][0]
        reac_rows.append(row)
    r.section("Service reactions by load case (kip, + up)", reac_rows)
    r.section("Load combinations", table)
    r.note("Unbraced lengths: compression flange braced at the points listed in Bracing.")
    r.note("Deflection: unfactored loads. Roof live and snow are added together (conservative).")
    if p["kind"] == "W" and "NONCOMPACT" in f0.detail:
        r.note("Web is noncompact: AISC F4 not implemented. Verify by hand.")
    return {"report": r, "section": p, "governing": gov, "per_case": per_case,
            "combos": results, "beam": beam, "length": length, "supports": supports}


def run(section="W18X35", auto_select=False, depth_filter="", **kw):
    """Analyze ``section``; if auto_select, return the lightest passing W shape."""
    if not auto_select:
        return analyze(section, **kw)
    names = steel.lightest("W", lambda p: True, nominal_depth=depth_filter or None)
    tried = []
    for n in names:
        out = analyze(n, **kw)
        tried.append((n, out["report"].max_ratio))
        if out["report"].ok:
            out["report"].note(f"Auto-select: lightest passing W shape is {n} "
                               f"({len(tried)} shapes checked).")
            out["tried"] = tried
            return out
    out = analyze(names[-1], **kw)
    out["report"].note("Auto-select: no W shape passed. Heaviest shape shown.")
    out["tried"] = tried
    return out
