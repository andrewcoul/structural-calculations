"""Wood shear wall line (segmented, full-height segments): SDPWS 2021.

* Nominal unit shear capacity vn: SDPWS 2021 Table 4.3A (wood structural
  panels, blocked). SDPWS 2021 tabulates ONE nominal value for wind and
  seismic; ASD: vn/2.0 wind, vn/2.8 seismic (4.1.4); LRFD: 0.8*vn wind,
  0.5*vn seismic.
* Specific gravity factor 1-(0.5-G) <= 1.0 (Table 4.3A footnote 3).
* Aspect ratio: h/b <= 3.5 (blocked WSP, Table 4.3.3); for h/b > 2.0 the
  capacity is multiplied by 1.25-0.125h/b (4.3.3.2) - or by 2b/h when the
  line shear is distributed in proportion to capacity (4.3.5.5.1 Exc. 1).
* Hold-down (chord) tension and compression from overturning.
* Deflection SDPWS 2021 Eq. 4.3-1 (three-term):
      delta = 8 v h^3/(E A b) + v h/(1000 Ga) + h da/b
  v plf, h ft, E psi, A in^2, b ft, Ga kips/in, da in.
* Story drift: Cd*delta/Ie vs allowable drift ratio (ASCE 7-22 12.12).

Loads: Eh (= rho*QE) and W are STRENGTH-level line shears (kip), as from
ASCE 7-22. ASD uses 0.7E and 0.6W.
"""

from __future__ import annotations

import math

from .. import wood
from ..report import Report

# Nominal unit shear capacity vn (plf) for 6/4/3/2 in. panel edge nail spacing,
# and apparent shear stiffness Ga (kips/in, OSB) - SDPWS 2021 Table 4.3A,
# wood structural panels, blocked. vn = 1.4 x (SDPWS 2015 seismic vs).
# Ga = None: not in this library - enter Ga_override for deflection.
TABLE_43A = {
    "Structural I 5/16 in. 6d": ((560, 840, 1090, 1430), (13, 18, 23, 35)),
    "Structural I 3/8 in. 8d": ((645, 1010, 1290, 1710), (19, 24, 30, 43)),
    "Structural I 7/16 in. 8d": ((715, 1105, 1415, 1875), (16, 21, 27, 40)),
    "Structural I 15/32 in. 8d": ((785, 1205, 1540, 2045), (14, 18, 24, 37)),
    "Structural I 15/32 in. 10d": ((950, 1430, 1860, 2435), None),
    "Sheathing 5/16 in. 6d": ((505, 755, 980, 1260), (13, 18, 24, 37)),
    "Sheathing 3/8 in. 6d": ((560, 840, 1090, 1430), (11, 15, 20, 32)),
    "Sheathing 3/8 in. 8d": ((615, 895, 1150, 1485), (17, 25, 31, 45)),
    "Sheathing 7/16 in. 8d": ((670, 980, 1260, 1640), (15, 22, 28, 42)),
    "Sheathing 15/32 in. 8d": ((730, 1065, 1370, 1790), (13, 19, 25, 39)),
    "Sheathing 15/32 in. 10d": ((870, 1290, 1680, 2155), (22, 30, 37, 52)),
    "Sheathing 19/32 in. 10d": ((950, 1430, 1860, 2435), (19, 26, 33, 48)),
}
SPACINGS = (6, 4, 3, 2)
SHEATHING = list(TABLE_43A)

OMEGA = {"seismic": 2.8, "wind": 2.0}   # SDPWS 2021 4.1.4 (ASD reduction factors)
PHI = {"seismic": 0.5, "wind": 0.8}     # SDPWS 2021 4.1.4 (LRFD resistance factors)


def table_43a(sheathing, spacing):
    vns, gas = TABLE_43A[sheathing]
    i = SPACINGS.index(int(spacing))
    return vns[i], (gas[i] if gas else None)


def aspect_factor(h, b, method="4.3.3.2"):
    """Return the aspect ratio factor. method: '4.3.3.2' (1.25-0.125h/b)
    or '2b/h' (4.3.5.5.1 Exception 1)."""
    ar = h / b
    if ar <= 2.0:
        return 1.0
    return min(1.0, 2 * b / h) if method == "2b/h" else 1.25 - 0.125 * ar


def parse_lengths(text):
    vals = [float(t) for t in str(text).replace(";", ",").split(",") if t.strip()]
    if not vals:
        raise ValueError("enter at least one full-height segment length")
    return vals


def analyze(sheathing="Sheathing 7/16 in. 8d", spacing=4, sides=1, framing_G=0.50,
            height=9.0, segments="8, 6", method="ASD", Eh=4.0, W=5.0, rho=1.0, sds=1.0,
            wall_D=150.0, wall_L=0.0, aspect_method="4.3.3.2", hd_offset=3.0,
            chord="(2)2x4", chord_species="DF-L", chord_grade="No.2", hd_allow=0.0,
            hd_defl=0.12, Cd=4.0, Ie=1.0, drift_limit=0.020, vn_override=0.0, Ga_override=0.0,
            seismic_wall=True, chord_trib=1.0):
    """Shear wall line. Eh, W: strength-level line shears (kip). wall_D/L: plf on wall."""
    method = method.upper()
    bs = parse_lengths(segments)
    Lsum = sum(bs)
    vn_tab, Ga_tab = table_43a(sheathing, spacing)
    vn1 = vn_override or vn_tab
    Ga1 = Ga_override or Ga_tab
    CG = min(1.0 - (0.5 - framing_G), 1.0)
    vn = vn1 * CG * sides
    Ga = Ga1 * sides if Ga1 else None

    r = Report("Wood Shear Wall (segmented)", code=f"SDPWS 2021 ({method}), ASCE 7-22, NDS 2024")
    r.input("Sheathing", f"{sheathing}, edge nailing {spacing} in. o.c., {sides} side(s), blocked")
    r.input("Framing specific gravity G", framing_G)
    r.input("Wall height h", height, "ft")
    r.input("Full-height segments b", ", ".join(f"{b:g}" for b in bs), "ft")
    r.input("Seismic line shear Eh (strength, incl. rho)", Eh, "kip")
    r.input("Wind line shear W (strength)", W, "kip")
    r.input("Wall dead / live load", f"{wall_D:g} / {wall_L:g}", "plf")
    r.input("SDS (for Ev)", sds)
    r.result("vn (Table 4.3A, one side)", vn1, "plf", "SDPWS 2021 Table 4.3A")
    r.result("Specific gravity factor 1-(0.5-G)", CG, "", "Table 4.3A note 3")
    r.result("Ga (OSB, total)", Ga if Ga else "not tabulated", "kips/in", "SDPWS 2021 Table 4.3A")
    if sides == 2:
        r.note("Two sides of the same sheathing: vn and Ga added (SDPWS 4.3.5.4.1, similar materials).")

    # capacity per segment (distribution in proportion to capacity)
    caps, rows = [], []
    for b in bs:
        af = aspect_factor(height, b, aspect_method)
        caps.append(vn * af * b)
    Ctot = sum(caps)
    ar_max = max(height / b for b in bs)
    r.check("Aspect ratio h/b (max segment)", ar_max, 3.5, "", "SDPWS 2021 Table 4.3.3")

    # (shear factor, dead factor, live factor) per ASCE 7-22 combination
    if method == "ASD":
        fac = {"seismic": 1 / OMEGA["seismic"], "wind": 1 / OMEGA["wind"]}
        up = {"seismic": (0.7, 0.6 - 0.14 * sds), "wind": (0.6, 0.6)}
        down = {"seismic": [(0.7, 1.0 + 0.14 * sds, 0.0), (0.525, 1.0 + 0.105 * sds, 0.75)],
                "wind": [(0.6, 1.0, 0.0), (0.45, 1.0, 0.75)]}
    else:
        fac = dict(PHI)
        up = {"seismic": (1.0, 0.9 - 0.2 * sds), "wind": (1.0, 0.9)}
        down = {"seismic": [(1.0, 1.2 + 0.2 * sds, 1.0)], "wind": [(1.0, 1.2, 1.0)]}
    strength = {"seismic": Eh, "wind": W}

    chord_sec = wood.parse_size(chord)
    gov_T = gov_C = None
    for kind in ("seismic", "wind"):
        Vs = strength[kind]
        if Vs <= 0:
            continue
        V = up[kind][0] * Vs
        cap = Ctot * fac[kind] / 1000.0  # kip
        r.check(f"Line shear capacity ({kind})", V, cap, "kip",
                "SDPWS 2021 4.3.3, 4.1.4",
                note=(f"ASD: {up[kind][0]:g}x{'E' if kind == 'seismic' else 'W'}, vn/{OMEGA[kind]}"
                      if method == "ASD" else f"LRFD: phi={PHI[kind]}"))
        for b, c in zip(bs, caps):
            share = c / Ctot
            arm = b - hd_offset / 12.0             # ft
            Vi = V * share                         # kip
            v = Vi / b * 1000.0                    # plf
            RM_up = up[kind][1] * wall_D / 1000.0 * b * b / 2
            T = max((Vi * height - RM_up) / arm, 0.0)
            C = max(fv * Vs * share * height / arm
                    + (fd * wall_D + fl * wall_L) / 1000.0 * chord_trib
                    for fv, fd, fl in down[kind])
            rows.append({"Case": kind, "b (ft)": b, "h/b": height / b, "Aspect factor":
                         aspect_factor(height, b, aspect_method), "V seg (k)": Vi, "v (plf)": v,
                         "T hold-down (k)": T, "C chord (k)": C})
            if gov_T is None or T > gov_T[0]:
                gov_T = (T, kind, b)
            if gov_C is None or C > gov_C[0]:
                gov_C = (C, kind, b)
    r.section("Segment forces (design level)", rows)

    if gov_T and gov_T[0] > 0:
        T, kind, b = gov_T
        r.result("Max hold-down tension", T, "kip", "SDPWS 4.3.6.1",
                 note=f"{kind}, b={b:g} ft, arm b-{hd_offset:g} in")
        if hd_allow > 0:
            r.check("Hold-down capacity", T, hd_allow, "kip", "manufacturer rating",
                    note=f"{kind}, b={b:g} ft")
        else:
            r.note("Enter the hold-down allowable (ASD) or design (LRFD) capacity to check it.")
    # compression chord: column with weak axis braced by sheathing
    if gov_C:
        from . import wood_column
        Cmax, kind, b = gov_C
        cc = wood_column.analyze(size=chord, species=chord_species, grade=chord_grade,
                                 method=method, height=height - 0.375,
                                 weak_bracing="Braced by sheathing (stud wall)",
                                 P={"W" if kind == "wind" else "E": Cmax / (0.6 if method == "ASD" and kind == "wind"
                                                                              else 0.7 if method == "ASD" else 1.0)},
                                 stud_wall=False)
        cap_ratio = cc["governing"]["axial"]
        Pallow = Cmax / cap_ratio if cap_ratio else math.inf
        r.check(f"Compression chord {chord} {chord_species} {chord_grade}", Cmax, Pallow, "kip",
                "NDS 3.7 (CD=1.6), 3.10" if method == "ASD" else "NDS 3.7 (lambda=1.0), 3.10",
                note=f"{kind}, b={b:g} ft, CP={cc['governing']['CP']:.3f}")
        Fcp = cc["member"].ref["Fc_perp"] * 1.0
        r.check("Chord bearing on plate", Cmax * 1000 / chord_sec.A, Fcp * wood.bearing_area_factor(chord_sec.b)
                * (wood.KF["Fc_perp"] * wood.PHI["Fc_perp"] if method == "LRFD" else 1.0), "psi",
                "NDS 3.10.2")

    # deflection and drift (seismic, strength level, rho = 1.0)
    if seismic_wall and Eh > 0:
        if not Ga:
            r.note("Ga not tabulated for this sheathing in this library: enter Ga_override to check drift.")
        else:
            QE = Eh / rho
            E_chord = wood.Member("Sawn lumber", chord_species, chord_grade, chord_sec).ref["E"]
            worst = None
            for b, c in zip(bs, caps):
                v = QE * c / Ctot / b * 1000.0
                d1 = 8 * v * height**3 / (E_chord * chord_sec.A * b)
                d2 = v * height / (1000.0 * Ga)
                d3 = height * hd_defl / b
                dsw = d1 + d2 + d3
                if worst is None or dsw > worst[0]:
                    worst = (dsw, b, d1, d2, d3, v)
            dsw, b, d1, d2, d3, v = worst
            r.result("Segment v for drift (strength, rho=1)", v, "plf", "ASCE 7-22 12.3.4.1")
            r.result("delta bending / shear / anchorage", f"{d1:.3f} / {d2:.3f} / {d3:.3f}", "in",
                     "SDPWS 2021 Eq. 4.3-1")
            r.result("delta_xe (elastic)", dsw, "in")
            dx = Cd * dsw / Ie
            r.check("Story drift Cd*delta/Ie", dx, drift_limit * height * 12, "in",
                    "ASCE 7-22 12.8.6, Table 12.12-1", note=f"b={b:g} ft, limit {drift_limit:g}h")
    r.note("vn values: SDPWS 2021 Table 4.3A, blocked wood structural panels, DF-L/SP framing. "
           "Footnote 2 (3/8 and 7/16 panels may use 15/32 values at 16 in. studs) is not applied.")
    r.note("Line shear is distributed to segments in proportion to capacity (incl. aspect factor).")
    r.note("Hold-down tension uses the uplift combination with reduced dead load "
           f"({'0.6-0.14SDS' if method == 'ASD' else '0.9-0.2SDS'} for seismic). Anchorage elongation "
           "da is the total hold-down slip/stretch at the drift load.")
    r.note("Check SDPWS 4.3.7.1 framing requirements (3x members) for 2 in. edge nailing, 10d at 3 in. "
           "or less, or vn > 980 plf per side in SDC D-F.")
    return {"report": r, "rows": rows, "vn": vn, "Ga": Ga, "caps": caps}


def run(**kw):
    return analyze(**kw)
