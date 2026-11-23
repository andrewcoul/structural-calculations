"""Wood column / stud: NDS 2024 axial + biaxial bending (NDS 3.7, 3.9).

Axial loads by case (kip), end eccentricity (in), lateral wind pressure
on the strong axis (stud walls), weak-axis bracing by sheathing, nailed
built-up columns (NDS 15.3, Kf = 0.6), bearing on the plate (NDS 3.10).
"""

from __future__ import annotations

import math

from .. import loads, wood
from ..report import Report

WEAK_BRACING = ["Same as strong axis", "Braced by sheathing (stud wall)", "Braced at interval"]


def _cfu(sec: wood.Section, material):
    """Flat use factor for bending about the weak axis (Table 4A/4B adj.)."""
    if material != "Sawn lumber" or sec.kind != "sawn":
        return 1.0
    w, t = sec.nominal_d, sec.nominal_b
    if t <= 3:
        return 1.0 if w <= 3 else 1.1 if w <= 4 else 1.15 if w <= 8 else 1.2
    return 1.0 if w <= 4 else 1.05 if w <= 8 else 1.1


def analyze(size="2x6", material="Sawn lumber", species="DF-L", grade="Stud", method="ASD",
            height=10.0, Ke=1.0, weak_bracing="Braced by sheathing (stud wall)",
            weak_interval=0.0, P=None, ex=0.0, ey=0.0, wind_psf=0.0, trib_width=16.0,
            stud_wall=True, wet=False, temp="T <= 100F", incised=False, Kf=0.6,
            plate_Fc_perp=0.0, override=None):
    """P: {case: kip} axial (+ compression). wind_psf: strength-level wind
    pressure on the strong axis (ASCE 7-22). trib_width in inches."""
    P = {k: v for k, v in (P or {}).items() if v}
    method = method.upper()
    sec = wood.parse_size(size, material)
    m = wood.Member(material, species, grade, sec, wet=wet, temp=temp, incised=incised,
                    repetitive=stud_wall, override=override or {})
    L_in = height * 12
    le1 = Ke * L_in
    if weak_bracing.startswith("Braced by sheathing"):
        le2 = 0.0
    elif weak_bracing == "Braced at interval" and weak_interval > 0:
        le2 = Ke * min(weak_interval * 12, L_in)
    else:
        le2 = le1
    d1, d2 = sec.d, sec.b
    s1, s2 = le1 / d1, le2 / d2 if le2 else 0.0
    c = wood.column_c(material)
    Emin_p = m.Emin_adj(method)
    built_up = sec.plies > 1 and material == "Sawn lumber"

    r = Report(f"Wood Column - {size} {species + ' ' + grade if material == 'Sawn lumber' else material}",
               code=f"NDS 2024 ({method}), ASCE 7-22")
    r.input("Section", f"{size} ({sec.b:g} x {sec.d:g} in actual)")
    r.input("Material", material if material != "Sawn lumber" else f"{species} {grade}")
    r.input("Height", height, "ft")
    r.input("Ke", Ke)
    r.input("Weak-axis bracing", weak_bracing + (f" @ {weak_interval:g} ft" if weak_bracing == "Braced at interval" else ""))
    for k, v in P.items():
        r.input(f"Axial load {k}", v, "kip")
    if ex or ey:
        r.input("Eccentricity ex / ey", f"{ex:g} / {ey:g}", "in")
    if wind_psf:
        r.input("Wind pressure (strength level)", wind_psf, "psf")
        r.input("Tributary width", trib_width, "in")
    r.input("Stud wall (Cr = 1.15 on Fb)", bool(stud_wall))
    wood.reference_inputs(r, m)
    r.result("A / Sx / Sy", f"{sec.A:.2f} / {sec.Sx:.2f} / {sec.Sy:.2f}", "in² / in³")
    r.result("le1/d1 (strong axis)", s1, "", "NDS 3.7.1")
    r.result("le2/d2 (weak axis)", s2 if le2 else "braced", "", "NDS 3.7.1")
    smax = max(s1, s2)
    r.check("Slenderness le/d", smax, 50.0, "", "NDS 3.7.1.4")

    cases = set(P) | ({"W"} if wind_psf else set())
    combos = loads.asd() if method == "ASD" else loads.lrfd()
    combos = loads.active(combos, cases)
    w_plf = wind_psf * trib_width / 12.0  # lb/ft (strength level)
    Mx_W = w_plf * height**2 / 8.0 * 12.0  # lb-in, strength level
    Cfu = _cfu(sec, material)
    rows, gov = [], None
    for name, f in combos:
        Pc = sum(P.get(k, 0.0) * v for k, v in f.items()) * 1000.0  # lb
        Mx = abs(Pc) * ey + abs(f.get("W", 0.0)) * Mx_W  # both add (conservative)
        My = abs(Pc) * ex
        Fc_star = m.Fc_star(method, f)
        CP1, FcE1 = wood.column_stability(Fc_star, Emin_p, s1, c)
        CP2, FcE2 = wood.column_stability(Fc_star, Emin_p, s2, c) if le2 else (1.0, math.inf)
        if built_up and le2:
            CP2 *= Kf
        CP = min(CP1, CP2)
        Fcp = Fc_star * CP
        Fb_star = m.Fb_star(method, f)
        # strong-axis bending stability: compression edge braced in stud walls
        dob = sec.d / sec.b
        if le2 == 0.0 or dob <= 2.0 or (material == "Sawn lumber" and dob <= 4.0 and stud_wall):
            CL, FbE = 1.0, math.inf
        else:
            le_b = wood.le_beam(L_in if not le2 else le2, sec.d, "general")
            CL, _, FbE = wood.beam_stability(Fb_star, Emin_p, le_b, sec.d, sec.b)
        Fb1 = Fb_star * CL
        Fb2 = Fb_star / m.ref["Fb"] * m.ref.get("Fby", m.ref["Fb"]) * Cfu
        A = sec.A
        if Pc >= 0:
            fc = Pc / A
            fb1, fb2 = Mx / sec.Sx, My / sec.Sy
            ok_terms = fc < FcE1 and fc < FcE2 and fb1 < FbE
            if ok_terms:
                t2 = fb2 / (Fb2 * (1 - fc / FcE2 - (fb1 / FbE) ** 2)) if fb2 else 0.0
                ratio = (fc / Fcp) ** 2 + fb1 / (Fb1 * (1 - fc / FcE1)) + t2
                eq = "NDS 3.9-3"
                if fb2 and fc / FcE2 + (fb1 / FbE) ** 2 >= 1.0:
                    ratio = math.inf
            else:
                ratio = math.inf
                eq = "fc >= FcE (unstable)"
            axial = fc / Fcp
        else:
            ft = -Pc / A
            Ftp = m.Ft_adj(method, f)
            fb1, fb2 = Mx / sec.Sx, My / sec.Sy
            ratio = ft / Ftp + fb1 / Fb_star + fb2 / Fb2
            eq = "NDS 3.9-1 (tension)"
            fc, axial, CP = -ft, ft / Ftp, CP
        rows.append({"Combination": name, "P (lb)": Pc, "Mx (lb-in)": Mx, "My (lb-in)": My,
                     "CP": CP, "Fc' (psi)": Fcp, "Axial ratio": axial, "Interaction": ratio})
        if gov is None or ratio > gov["ratio"]:
            gov = dict(name=name, ratio=ratio, eq=eq, fc=fc, Fcp=Fcp, CP=CP, Fc_star=Fc_star,
                       FcE1=FcE1, FcE2=FcE2, fb1=fb1, fb2=fb2, Fb1=Fb1, Fb2=Fb2, CL=CL, P=Pc,
                       axial=axial)
    if gov is None:
        raise ValueError("no loads entered")
    r.result("Governing combination", gov["name"])
    r.result("Fc* (all factors except CP)", gov["Fc_star"], "psi")
    r.result("FcE1 / FcE2", f"{gov['FcE1']:.0f} / {gov['FcE2']:.0f}" if gov["FcE2"] != math.inf
             else f"{gov['FcE1']:.0f} / braced", "psi", "NDS 3.7.1")
    r.result("CP (column stability)", gov["CP"], "", "NDS 3.7.1" + (" + 15.3 (Kf)" if built_up and le2 else ""))
    r.result("fc = P/A", gov["fc"], "psi")
    r.result("fb1 (strong) / fb2 (weak)", f"{gov['fb1']:.0f} / {gov['fb2']:.0f}", "psi")
    r.result("Fb1' / Fb2'", f"{gov['Fb1']:.0f} / {gov['Fb2']:.0f}", "psi")
    ax = max(rows, key=lambda row: row["Axial ratio"])
    r.check("Axial only fc/Fc' (or ft/Ft')", ax["Axial ratio"], 1.0, "", "NDS 3.6.3 / 3.8.1",
            note=f"[{ax['Combination']}]")
    r.check("Combined axial + bending", gov["ratio"], 1.0, "", gov["eq"], note=f"[{gov['name']}]")

    # Bearing on plate (max compression over combinations, no load duration on Fc_perp)
    Pmax = max((row["P (lb)"] for row in rows), default=0.0)
    if Pmax > 0:
        Cb = wood.bearing_area_factor(sec.b)  # bearing length along the plate = column thickness
        Fcp_plate = (plate_Fc_perp or m.ref["Fc_perp"]) * m.CM["Fc_perp"] * m.Ct["Fc_perp"] * Cb
        if method == "LRFD":
            Fcp_plate *= wood.KF["Fc_perp"] * wood.PHI["Fc_perp"]
        r.check("Bearing on plate fc_perp", Pmax / sec.A, Fcp_plate, "psi", "NDS 3.10.2/3.10.4",
                note=f"Cb={Cb:.3f}")
    r.section("Load combinations", rows)
    r.note("End eccentricity moments P*e are applied over the full height (conservative).")
    r.note("Wind bending: simple span w*h^2/8 on the strong axis, combined with axial load in "
           "the same combination.")
    if stud_wall:
        r.note("Stud wall: repetitive member factor Cr = 1.15 on Fb; strong-axis compression edge "
               "assumed braced by sheathing when weak axis is braced.")
    if built_up:
        r.note(f"Built-up column: weak-axis CP multiplied by Kf = {Kf:g} (nailed, NDS 15.3.2). "
               "Verify fastening requirements of NDS 15.3.3.")
    return {"report": r, "member": m, "section": sec, "governing": gov, "rows": rows}


def run(**kw):
    return analyze(**kw)
