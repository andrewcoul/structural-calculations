"""Wood beam: NDS 2024 (ASD or LRFD) bending, shear, bearing, deflection.

Sawn lumber (NDS Supplement Tables 4A/4B/4D), DF glulam (Table 5A) and
generic LVL/PSL. Analysis with structcalc.beam (any span type, any loads),
ASCE 7-22 load combinations, CD (ASD) or KF*phi*lambda (LRFD) per
combination.
"""

from __future__ import annotations

import numpy as np

from .. import loads, wood
from ..beam import Beam
from ..beaminput import apply_loads, deflection_checks, make_supports, parse_loads
from ..report import Report

BRACING = ["Continuous (compression edge braced)", "Top edge braced (sheathing on top)",
           "At supports only", "At interval"]

SAWN_SIZES = ["2x6", "2x8", "2x10", "2x12", "(2)2x8", "(2)2x10", "(2)2x12", "(3)2x10",
              "(3)2x12", "4x6", "4x8", "4x10", "4x12", "4x14", "(4)2x12", "6x8", "6x10", "6x12",
              "6x14", "6x16", "8x12", "8x14", "8x16", "8x18"]
GLULAM_SIZES = [f"{b}x{d:g}" for b in (3.125, 5.125, 6.75, 8.75)
                for d in np.arange(6, 37.5, 1.5)]
SCL_SIZES = [f"({n})1.75x{d}" for n in (1, 2, 3, 4)
             for d in (7.25, 9.25, 9.5, 11.25, 11.875, 14, 16, 18, 20, 24)]


def candidates(material):
    if material == "Sawn lumber":
        names = SAWN_SIZES
    elif material == "Glulam":
        names = GLULAM_SIZES
    else:
        names = SCL_SIZES
    secs = []
    for n in names:
        s = wood.parse_size(n, material)
        secs.append((s.A, s.d, n))
    return [n for _, _, n in sorted(secs)]


def _brace_points(length, supports, mode, interval):
    xs = {0.0, length} | {x for x, _ in supports}
    if mode == "At interval" and interval > 0:
        n = int(np.floor(length / interval + 1e-9))
        xs |= {i * interval for i in range(n + 1)}
    return sorted(x for x in xs if 0 <= x <= length + 1e-9)


def analyze(size="2x10", material="Sawn lumber", species="DF-L", grade="No.2", method="ASD",
            span_type="Simple span", span=12.0, overhang_left=0.0, overhang_right=0.0,
            custom_supports="", w=None, extra_loads="", self_weight=True, density=35.0,
            bracing="Continuous (compression edge braced)", brace_interval=0.0,
            le_load="general", wet=False, temp="T <= 100F", incised=False, repetitive=False,
            bearing_length=3.5, live_defl=360, total_defl=240, Kcr=1.0, override=None):
    """Analyze one wood beam. Loads in kip/ft (w dict) and kip (extra loads)."""
    w = dict(w or {})
    method = method.upper()
    sec = wood.parse_size(size, material)
    m = wood.Member(material, species, grade, sec, wet=wet, temp=temp, incised=incised,
                    repetitive=repetitive, override=override or {})
    length, supports = make_supports(span_type, span, overhang_left, overhang_right, custom_supports)
    E_adj = m.E_adj()
    beam = Beam(length, E_adj / 1000.0, sec.Ix, supports)
    apply_loads(beam, w, parse_loads(extra_loads))
    if self_weight:
        beam.self_weight = density * sec.A / 144.0 / 1000.0
    per_case = beam.solve()
    combos = loads.asd() if method == "ASD" else loads.lrfd()
    results = beam.combos(combos)
    brace_pts = _brace_points(length, supports, bracing, brace_interval)
    sup_x = sorted(x for x, _ in supports)
    single_fixed = len(supports) == 1

    r = Report(f"Wood Beam - {size} {material if material != 'Sawn lumber' else species + ' ' + grade}",
               code=f"NDS 2024 ({method}), ASCE 7-22, IBC 2024")
    r.input("Section", f"{size} ({sec.b:g} x {sec.d:g} in actual)")
    r.input("Material", material if material != "Sawn lumber" else f"{species} {grade}")
    r.input("Span type", span_type)
    r.input("Total length", length, "ft")
    r.input("Supports", ", ".join(f"{x:g} ft {t}" for x, t in supports))
    for c, v in w.items():
        if v:
            r.input(f"Uniform load {c}", v * 1000, "plf")
    if extra_loads.strip():
        r.input("Extra loads", extra_loads.strip())
    r.input("Self weight", f"{density:g} pcf" if self_weight else "not included")
    r.input("Bracing", bracing + (f" @ {brace_interval:g} ft" if bracing == "At interval" else ""))
    r.input("Service", ("wet" if wet else "dry") + f", {temp}" + (", incised" if incised else "")
            + (", repetitive (Cr)" if repetitive else ""))
    wood.reference_inputs(r, m)
    r.result("A / Sx / Ix", f"{sec.A:.2f} / {sec.Sx:.2f} / {sec.Ix:.1f}", "in² / in³ / in⁴")
    r.section("Adjustment factors (constant)", m.factor_table())

    d_over_b = sec.d / sec.ply_b if sec.plies > 1 else sec.d / sec.b
    CV = wood.volume_factor(max(min(b - a for a, b in zip(sup_x[:-1], sup_x[1:])) if len(sup_x) > 1
                                else length, 1.0), sec) if material == "Glulam" else 1.0
    if material == "Glulam":
        r.result("CV (volume factor)", CV, "", "NDS 5.3.6")

    best = {"flex": None, "shear": None, "bear": None}
    table = []
    for name, f, res in results:
        Fb_star = m.Fb_star(method, f)
        Emin_p = m.Emin_adj(method)
        seg_best = None
        for a, b in zip(brace_pts[:-1], brace_pts[1:]):
            if b - a < 1e-6:
                continue
            msk = (res.x >= a - 1e-9) & (res.x <= b + 1e-9)
            Mpos = max(float(res.M[msk].max()), 0.0)
            Mneg = max(float(-res.M[msk].min()), 0.0)
            lu = (b - a) * 12
            free = (a < sup_x[0] - 1e-9) or (b > sup_x[-1] + 1e-9) or single_fixed
            if bracing.startswith("Continuous"):
                lu_pos = lu_neg = 0.0
            elif bracing.startswith("Top edge"):
                lu_pos, lu_neg = 0.0, lu
            else:
                lu_pos = lu_neg = lu
            for Mx, lux, sign in ((Mpos, lu_pos, "+"), (Mneg, lu_neg, "-")):
                if Mx <= 0:
                    continue
                if lux <= 0 or d_over_b <= 2.0:
                    CL, le, RB = 1.0, 0.0, 0.0
                else:
                    lt = "cantilever_uniform" if (free and le_load == "uniform") else le_load
                    le = wood.le_beam(lux, sec.d, lt)
                    CL, RB, _ = wood.beam_stability(Fb_star, Emin_p, le, sec.d, sec.b)
                Fb_base = Fb_star
                if material == "Glulam" and sign == "-":
                    Fb_base = Fb_star * m.ref.get("Fb_neg", m.ref["Fb"]) / m.ref["Fb"]
                Fbp = Fb_base * (min(CL, CV) if material == "Glulam" else CL)
                fb = Mx * 12000 / sec.Sx
                ratio = fb / Fbp
                if seg_best is None or ratio > seg_best["ratio"]:
                    seg_best = dict(ratio=ratio, fb=fb, Fbp=Fbp, CL=CL, le=le, RB=RB, M=Mx,
                                    a=a, b=b, sign=sign)
        Fvp = m.Fv_adj(method, f)
        V = res.Vabs
        fv = 1.5 * V * 1000 / sec.A
        Rmax = max(abs(v[0]) for v in res.reactions.values())
        table.append({"Combination": name, "CD/λ": wood.CD_for_combo(f) if method == "ASD"
                      else wood.lambda_for_combo(f), "M+ (k-ft)": round(res.Mmax, 6), "M- (k-ft)": round(res.Mmin, 6),
                      "V (k)": V, "fb/Fb'": seg_best["ratio"] if seg_best else 0.0,
                      "fv/Fv'": fv / Fvp})
        if seg_best and (best["flex"] is None or seg_best["ratio"] > best["flex"][1]["ratio"]):
            best["flex"] = (name, seg_best, res, Fb_star)
        if best["shear"] is None or fv / Fvp > best["shear"][1] / best["shear"][2]:
            best["shear"] = (name, fv, Fvp, res)
        # bearing at each support
        for x, (R, _) in res.reactions.items():
            at_end = x <= 1e-6 or x >= length - 1e-6
            Cb = wood.bearing_area_factor(bearing_length, at_end)
            Fcp = m.Fc_perp_adj(method, Cb)
            fcp = abs(R) * 1000 / (sec.b * bearing_length)
            if best["bear"] is None or fcp / Fcp > best["bear"][1] / best["bear"][2]:
                best["bear"] = (name, fcp, Fcp, x, abs(R), Cb)

    fname, seg, fres, Fb_star = best["flex"]
    r.result("Governing flexure combination", fname)
    r.result("Fb* (all factors except CL/CV)", Fb_star, "psi")
    r.result("Governing segment", f"{seg['a']:g} to {seg['b']:g} ft ({seg['sign']} moment)")
    if seg["le"]:
        r.result("Effective length le", seg["le"] / 12, "ft", "NDS Table 3.3.3")
        r.result("Slenderness RB", seg["RB"], "", "NDS 3.3.3.6 (RB <= 50)")
        if seg["RB"] > 50:
            r.check("Beam slenderness RB", seg["RB"], 50, "", "NDS 3.3.3.7")
    r.result("CL (beam stability)", seg["CL"], "", "NDS 3.3.3")
    if d_over_b <= 2.0 and not bracing.startswith("Continuous"):
        r.note("d/b <= 2: no lateral support required, CL = 1.0 (NDS 4.4.1.2).")
    r.check("Bending fb", seg["fb"], seg["Fbp"], "psi", "NDS 3.3", note=f"[{fname}] M={seg['M']:.2f} k-ft")
    sname, fv, Fvp, _ = best["shear"]
    r.check("Shear fv = 1.5V/A", fv, Fvp, "psi", "NDS 3.4.2", note=f"[{sname}]")
    bname, fcp, Fcp, bx, R, Cb = best["bear"]
    r.check(f"Bearing at support x={bx:g} ft", fcp, Fcp, "psi", "NDS 3.10.2",
            note=f"[{bname}] R={R:.2f} k, lb={bearing_length:g} in, Cb={Cb:.3f}")
    r.result("Required bearing length (governing support)", R * 1000 / (sec.b * Fcp), "in")
    deflection_checks(r, per_case, length, supports, [
        ("Live (L)", {"L": 1.0}, live_defl),
        ("Roof live/snow (Lr or S)", {"Lr": 1.0, "S": 1.0}, live_defl),
        (f"Total ({'Kcr·' if Kcr != 1 else ''}D+L+Lr+S)", {"D": Kcr, "L": 1.0, "Lr": 1.0, "S": 1.0},
         total_defl),
    ])
    r.result("E' (deflection)", E_adj, "psi")
    reac = []
    for x in sorted(per_case[next(iter(per_case))].reactions):
        row = {"Support x (ft)": x}
        for c, res in per_case.items():
            row[f"R {c} (k)"] = res.reactions[x][0]
        reac.append(row)
    r.section("Service reactions by load case (kip, + up)", reac)
    r.section("Load combinations", table)
    r.note("Shear uses the maximum V (conservative); NDS 3.4.3.1 permits loads within d of the "
           "support to be ignored.")
    if Kcr != 1.0:
        r.note(f"Long-term deflection: Kcr = {Kcr:g} applied to dead load (NDS 3.5.2).")
    if material.startswith(("LVL", "PSL")):
        r.note("SCL values are generic — verify design values against the manufacturer's evaluation report.")
    if material == "Glulam":
        r.note("Glulam: CV and CL are not combined; the lesser governs (NDS 5.3.6). "
               "Negative moment uses Fbx- (compression zone stressed in tension).")
    brace_note = bracing + (f" @ {brace_interval:g} ft" if bracing == "At interval" else "")
    show_braces = brace_pts if bracing in ("At interval", "At supports only") else None
    return {"report": r, "member": m, "section": sec, "governing": best, "per_case": per_case,
            "combos": results, "beam": beam, "length": length, "supports": supports,
            "brace_pts": show_braces, "bracing_note": f"Bracing: {brace_note}"}


def run(size="2x10", auto_select=False, **kw):
    if not auto_select:
        return analyze(size, **kw)
    material = kw.get("material", "Sawn lumber")
    tried = []
    for n in candidates(material):
        try:
            out = analyze(n, **kw)
        except ValueError:
            continue  # e.g. timber values not tabulated for this species
        tried.append((n, out["report"].max_ratio))
        if out["report"].ok:
            out["report"].note(f"Auto-select: lightest passing size is {n} ({len(tried)} sizes checked).")
            out["tried"] = tried
            return out
    if not tried:
        raise ValueError("no candidate size could be analyzed")
    out = analyze(tried[-1][0], **kw)
    out["report"].note("Auto-select: no candidate size passed. Largest candidate shown.")
    out["tried"] = tried
    return out
