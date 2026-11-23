"""Wind loads: ASCE 7-22 Ch. 26 + Ch. 27 directional procedure (MWFRS) for an
enclosed / partially enclosed / partially open / open rectangular building,
plus Ch. 30 C&C wall zones 4 and 5 (h <= 60 ft).

ASCE 7-22 changes used here (verified):
* qz = 0.00256 Kz Kzt Ke V^2 (Eq. 26.10-1) -> Kd is NOT in qz.
* p = q Kd G Cp - qi Kd (GCpi) (Eq. 27.3-1); C&C p = qh Kd [(GCp) - (GCpi)].
* Kz = 2.41 (z/zg)^(2/alpha) for 15 ft <= z <= zg; z < 15 ft uses z = 15 ft
  (Table 26.10-1 note). Table 26.11-1: B alpha=7.5 zg=3280 ft;
  C alpha=9.8 zg=2460 ft; D alpha=11.5 zg=1935 ft.
* Ke = exp(-0.0000362 ze) (Table 26.9-1 note), ze = ground elevation, ft.

Pressures are strength level (load case W, LRFD factor 1.0, ASD 0.6).
Units: psf, ft, kip.
"""

from __future__ import annotations

import math

import numpy as np

from ..report import Report

EXPOSURE = {"B": (7.5, 3280.0), "C": (9.8, 2460.0), "D": (11.5, 1935.0)}
GCPI = {"Enclosed": 0.18, "Partially enclosed": 0.55, "Partially open": 0.18, "Open": 0.0}

# Fig. 27.3-1 roof Cp, wind normal to ridge, theta >= 10 deg.
# windward: h/L -> [(angle, Cp_min, Cp_max)]
_WW = {
    0.25: [(10, -0.7, -0.18), (15, -0.5, 0.0), (20, -0.3, 0.2), (25, -0.2, 0.3), (30, -0.2, 0.3),
           (35, 0.0, 0.4), (45, 0.4, 0.4), (60, 0.6, 0.6)],
    0.5: [(10, -0.9, -0.18), (15, -0.7, -0.18), (20, -0.4, 0.0), (25, -0.3, 0.2), (30, -0.2, 0.2),
          (35, -0.2, 0.3), (45, 0.0, 0.4), (60, 0.6, 0.6)],
    1.0: [(10, -1.3, -0.18), (15, -1.0, -0.18), (20, -0.7, -0.18), (25, -0.5, 0.0), (30, -0.3, 0.2),
          (35, -0.2, 0.2), (45, 0.0, 0.3), (60, 0.6, 0.6)],
}
_LW = {
    0.25: [(10, -0.3), (15, -0.5), (20, -0.6)],
    0.5: [(10, -0.5), (15, -0.5), (20, -0.6)],
    1.0: [(10, -0.7), (15, -0.6), (20, -0.6)],
}


def Kz(z, exposure="C"):
    alpha, zg = EXPOSURE[exposure]
    z = min(max(z, 15.0), zg)
    return 2.41 * (z / zg) ** (2.0 / alpha)


def Ke(elev_ft):
    return math.exp(-0.0000362 * elev_ft)


def qz(z, V, exposure="C", Kzt=1.0, ke=1.0):
    """Velocity pressure (psf), ASCE 7-22 Eq. 26.10-1 (no Kd)."""
    return 0.00256 * Kz(z, exposure) * Kzt * ke * V**2


def cp_leeward_wall(L_over_B):
    """Fig. 27.3-1 leeward wall Cp: 0-1: -0.5, 2: -0.3, >=4: -0.2."""
    return float(np.interp(L_over_B, [0, 1, 2, 4], [-0.5, -0.5, -0.3, -0.2]))


def _interp_hl(table, hL, fn):
    """Evaluate fn on each h/L row, then interpolate linearly in h/L."""
    keys = sorted(table)
    hL = min(max(hL, keys[0]), keys[-1])
    vals = [fn(table[k]) for k in keys]
    return tuple(float(np.interp(hL, keys, [v[i] for v in vals])) for i in range(len(vals[0])))


def _trapz(y, x):
    f = getattr(np, "trapezoid", None) or np.trapz
    return float(f(y, x))


def cp_roof_windward(theta, h_over_L):
    """(Cp_min, Cp_max) for windward roof, wind normal to ridge, theta >= 10."""
    def at_angle(rows):
        if theta >= 60:
            return (0.01 * theta, 0.01 * theta)
        a = [r[0] for r in rows]
        return (float(np.interp(theta, a, [r[1] for r in rows])),
                float(np.interp(theta, a, [r[2] for r in rows])))
    return _interp_hl(_WW, h_over_L, at_angle)


def cp_roof_leeward(theta, h_over_L):
    def at_angle(rows):
        a = [r[0] for r in rows]
        return (float(np.interp(theta, a, [r[1] for r in rows])),)
    return _interp_hl(_LW, h_over_L, at_angle)[0]


def flat_roof_zones(h, L):
    """Fig. 27.3-1 roof Cp for theta < 10 deg or wind parallel to ridge.

    Returns [(x1, x2, Cp_min, Cp_max)] with x from the windward edge (ft).
    Linear interpolation between h/L = 0.5 and 1.0 is applied zone by zone.
    """
    hL = h / L

    def table_low(x):  # h/L <= 0.5
        if x <= h:
            return -0.9
        if x <= 2 * h:
            return -0.5
        return -0.3

    def table_high(x):  # h/L >= 1.0
        return -1.3 if x <= h / 2 else -0.7

    edges = sorted({0.0, min(h / 2, L), min(h, L), min(2 * h, L), L})
    out = []
    for x1, x2 in zip(edges[:-1], edges[1:]):
        if x2 - x1 < 1e-9:
            continue
        xm = 0.5 * (x1 + x2)
        lo, hi = table_low(xm), table_high(xm)
        if hL <= 0.5:
            cp = lo
        elif hL >= 1.0:
            cp = hi
        else:
            cp = lo + (hi - lo) * (hL - 0.5) / 0.5
        out.append((x1, x2, cp, -0.18))
    return out


def gcp_wall_cc(zone, A, theta):
    """Fig. 30.3-1 wall GCp (h <= 60 ft): returns (positive, negative)."""
    A = min(max(A, 10.0), 500.0)
    t = (math.log10(A) - 1.0) / (math.log10(500.0) - 1.0)
    pos = 1.0 + (0.7 - 1.0) * t
    neg_start = -1.1 if zone == 4 else -1.4
    neg = neg_start + (-0.8 - neg_start) * t
    if theta <= 10:
        pos *= 0.9
        neg *= 0.9
    return pos, neg


def run(V=115.0, exposure="C", Kzt=1.0, ground_elev=0.0, use_Ke=True, Kd=0.85, G=0.85,
        enclosure="Enclosed", roof_type="Gable", theta=20.0, h_eave=20.0,
        B_ridge=100.0, L_across=60.0, heights="", cc_area=20.0, risk_category="II"):
    """B_ridge = plan dimension along the ridge; L_across = plan dimension
    across the ridge (for a flat roof, any two plan dimensions). ft, deg, mph."""
    ke = Ke(ground_elev) if use_Ke else 1.0
    gcpi = GCPI[enclosure]
    if roof_type == "Flat":
        theta = 0.0
    rise = 0.0 if roof_type == "Flat" else (L_across / 2.0) * math.tan(math.radians(theta))
    h = h_eave + rise / 2.0 if theta > 10 else h_eave
    ridge = h_eave + rise
    qh = qz(h, V, exposure, Kzt, ke)
    qi = qh

    r = Report("Wind Loads - MWFRS Directional Procedure + C&C Walls",
               code="ASCE 7-22 Ch. 26, 27, 30 (strength-level W)")
    r.input("Basic wind speed V", V, "mph", note=f"Risk Category {risk_category}")
    r.input("Exposure", exposure)
    r.input("Kzt", Kzt)
    r.input("Ground elevation", ground_elev, "ft")
    r.input("Kd", Kd, "", note="Table 26.6-1")
    r.input("Gust effect factor G", G, "", note="26.11.1 rigid")
    r.input("Enclosure", enclosure)
    r.input("Roof type", roof_type)
    r.input("Roof angle", theta, "deg")
    r.input("Eave height", h_eave, "ft")
    r.input("Plan dimension along ridge B", B_ridge, "ft")
    r.input("Plan dimension across ridge L", L_across, "ft")
    alpha, zg = EXPOSURE[exposure]
    r.result("alpha / zg", f"{alpha:g} / {zg:g} ft", "", "Table 26.11-1")
    r.result("Ke", ke, "", "Table 26.9-1")
    r.result("Roof rise", rise, "ft")
    r.result("Mean roof height h", h, "ft", "26.2")
    r.result("Kh", Kz(h, exposure), "", "Table 26.10-1")
    r.result("qh", qh, "psf", "Eq. 26.10-1")
    r.result("GCpi", f"+/-{gcpi:g}", "", "Table 26.13-1")

    # velocity pressure table
    if heights.strip():
        hs = sorted({float(x) for x in heights.replace(";", ",").split(",") if x.strip()})
    else:
        hs = sorted({0, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160} & set(
            range(0, int(ridge) + 1)) | {round(h_eave, 2), round(h, 2), round(ridge, 2)})
    r.section("Velocity pressure", [
        {"z (ft)": z, "Kz": Kz(z, exposure), "qz (psf)": qz(z, V, exposure, Kzt, ke),
         "p windward wall, ext. (psf)": qz(z, V, exposure, Kzt, ke) * Kd * G * 0.8} for z in hs])

    def press(q, cp):
        ext = q * Kd * G * cp
        return ext - qi * Kd * gcpi, ext + qi * Kd * gcpi  # (+GCpi), (-GCpi)

    directions = []
    # Direction 1: normal to ridge (wind along L_across); Direction 2: parallel (along B_ridge)
    for dname, L, B in (("Dir 1: wind normal to ridge", L_across, B_ridge),
                        ("Dir 2: wind parallel to ridge", B_ridge, L_across)):
        normal = dname.startswith("Dir 1")
        rows = []
        cp_lw = cp_leeward_wall(L / B)
        rows.append(("Windward wall (at h)", 0.8, qh))
        rows.append(("Leeward wall", cp_lw, qh))
        rows.append(("Side walls", -0.7, qh))
        roof_rows = []
        if normal and theta >= 10:
            ww = cp_roof_windward(theta, h / L)
            lw = cp_roof_leeward(theta, h / L)
            roof_rows.append(("Roof windward slope (Cp min)", ww[0], qh))
            if abs(ww[1] - ww[0]) > 1e-9:
                roof_rows.append(("Roof windward slope (Cp max)", ww[1], qh))
            roof_rows.append(("Roof leeward slope", lw, qh))
        else:
            for x1, x2, cmin, cmax in flat_roof_zones(h, L):
                roof_rows.append((f"Roof {x1:.1f} to {x2:.1f} ft from windward edge (Cp min)", cmin, qh))
                roof_rows.append((f"Roof {x1:.1f} to {x2:.1f} ft (Cp = -0.18 case)", cmax, qh))
        table = []
        for nm, cp, q in rows + roof_rows:
            pp, pn = press(q, cp)
            table.append({"Surface": nm, "Cp": cp, "q (psf)": q, "p with +GCpi (psf)": pp,
                          "p with -GCpi (psf)": pn})
        r.section(f"{dname}: L = {L:g} ft, B = {B:g} ft, h/L = {h / L:.3f}, L/B = {L / B:.3f}", table)

        # base shear and overturning (external pressures; internal pressures cancel)
        zz = np.linspace(0.0, h_eave, 201)
        pww = np.array([qz(z, V, exposure, Kzt, ke) for z in zz]) * Kd * G * 0.8  # psf
        Fww = _trapz(pww, zz) * B / 1000.0
        Mww = _trapz(pww * zz, zz) * B / 1000.0
        plw = -qh * Kd * G * cp_lw
        Flw = plw * h_eave * B / 1000.0
        Mlw = Flw * h_eave / 2.0
        F, M = Fww + Flw, Mww + Mlw
        wall_area = B * h_eave
        roof_proj = 0.0
        if normal and theta >= 10:
            ww = cp_roof_windward(theta, h / L)
            lw = cp_roof_leeward(theta, h / L)
            zr = h_eave + rise / 2.0
            Fr = (qh * Kd * G * max(ww) - qh * Kd * G * lw) * rise * B / 1000.0
            F += Fr
            M += Fr * zr
            roof_proj = rise * B
        elif not normal and rise > 0:
            # gable-end triangles on windward and leeward walls
            At = L * rise / 2.0 if roof_type == "Gable" else 0.0
            ft = (qh * Kd * G * 0.8 + plw) * At / 1000.0
            F += ft
            M += ft * (h_eave + rise / 3.0)
            wall_area += At
        Fmin = (16.0 * wall_area + 8.0 * roof_proj) / 1000.0
        directions.append({"Direction": dname, "Base shear (k)": F, "Overturning (k-ft)": M,
                           "Minimum 27.1.5 (k)": Fmin, "Design base shear (k)": max(F, Fmin)})
    r.section("MWFRS base shear and overturning (strength level)", directions)
    for d in directions:
        r.result(f"Base shear, {d['Direction']}", d["Design base shear (k)"], "kip",
                 "27.3 / 27.1.5 minimum")

    # C&C walls
    least = min(B_ridge, L_across)
    a = max(min(0.1 * least, 0.4 * h), 0.04 * least, 3.0)
    r.result("C&C edge zone width a", a, "ft", "Fig. 30.3-1 note")
    cc = []
    for zone in (4, 5):
        pos, neg = gcp_wall_cc(zone, cc_area, theta)
        p_pos = max(qh * Kd * (pos + gcpi), 16.0)
        p_neg = min(qh * Kd * (neg - gcpi), -16.0)
        cc.append({"Zone": zone, "A (ft^2)": cc_area, "GCp +": pos, "GCp -": neg,
                   "p + (psf)": p_pos, "p - (psf)": p_neg})
        r.result(f"C&C wall zone {zone} pressure", f"+{p_pos:.1f} / {p_neg:.1f}", "psf",
                 "Eq. 30.3-1, 16 psf min 30.2.2")
    r.section("Components and cladding - walls (h <= 60 ft)", cc)
    if h > 60:
        r.note("h > 60 ft: the C&C wall values above (Fig. 30.3-1) do not apply. Use ASCE 7-22 30.4.")
    if roof_type == "Hip":
        r.note("Hip roof: Cp from the gable values of Fig. 27.3-1 (approximation). "
               "Hip end slopes are not modeled separately.")
    if enclosure == "Partially enclosed":
        r.note("Partially enclosed: qi = qh used (conservative); 26.13.1.1 reduction factor Ri not applied.")
    r.note("Positive pressure acts toward the surface; negative is suction.")
    r.note("Fig. 27.3-1 interpolation: linear in angle and h/L; the -1.3 value for h/L >= 1 is "
           "not reduced for area (conservative).")
    r.note("Torsional load cases (Fig. 27.3-8) and roof uplift combinations are not generated.")
    r.note("Tornado loads (Ch. 32) are required for Risk Category III/IV in tornado-prone regions.")
    return {"report": r, "qh": qh, "h": h, "directions": directions, "cc": cc}


def plot_profile(V=115.0, exposure="C", Kzt=1.0, ke=1.0, Kd=0.85, G=0.85, height=60.0):
    import matplotlib.pyplot as plt

    z = np.linspace(0, height, 200)
    q = [qz(zi, V, exposure, Kzt, ke) for zi in z]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(q, z, label="qz")
    ax.plot(np.array(q) * Kd * G * 0.8, z, label="windward wall p (ext.)")
    ax.set_xlabel("psf")
    ax.set_ylabel("z (ft)")
    ax.grid(alpha=0.3)
    ax.legend()
    ax.set_title(f"Velocity pressure profile, Exposure {exposure}, V={V:g} mph")
    fig.tight_layout()
    return fig
