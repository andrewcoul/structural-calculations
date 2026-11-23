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
    # Capture Dir 1 surface pressures for the building diagram
    cp_lw1 = cp_leeward_wall(L_across / B_ridge)
    p_ww1 = press(qh, 0.8)
    p_lw1 = press(qh, cp_lw1)
    if theta >= 10:
        ww_cp = cp_roof_windward(theta, h / L_across)
        lw_cp = cp_roof_leeward(theta, h / L_across)
        p_rww1 = (min(press(qh, ww_cp[0])[0], press(qh, ww_cp[1])[0]),
                  max(press(qh, ww_cp[0])[1], press(qh, ww_cp[1])[1]))
        p_rlw1 = press(qh, lw_cp)
    else:
        z0 = flat_roof_zones(h, L_across)[0]
        p_rww1 = press(qh, z0[2])
        p_rlw1 = press(qh, flat_roof_zones(h, L_across)[-1][2])
    geom = dict(V=V, exposure=exposure, Kzt=Kzt, ke=ke, Kd=Kd, G=G, roof_type=roof_type,
                theta=theta, h_eave=h_eave, rise=rise, ridge=ridge, h=h,
                B_ridge=B_ridge, L_across=L_across, a=a,
                p_ww=p_ww1, p_lw=p_lw1, p_rww=p_rww1, p_rlw=p_rlw1)
    return {"report": r, "qh": qh, "h": h, "directions": directions, "cc": cc, "geom": geom}


def plot(out):
    """Building MWFRS pressure elevation (Dir 1), C&C wall zone plan, and velocity pressure profile."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Rectangle

    g = out["geom"]
    fig, (axb, axc, axq) = plt.subplots(1, 3, figsize=(13.2, 4.4), gridspec_kw={"width_ratios": [1.4, 0.95, 0.95]})

    # --- Panel 0: Building cross-section with Dir 1 MWFRS net pressures ---
    L, he, rd = g["L_across"], g["h_eave"], g["ridge"]
    pts = [(0, 0), (L, 0), (L, he), (L / 2, rd), (0, he)] if g["rise"] > 0 else [(0, 0), (L, 0), (L, he), (0, he)]
    axb.add_patch(Polygon(pts, closed=True, fc="#e2e8f0", ec="#1e293b", lw=1.6, zorder=2))
    axb.plot([-0.28 * L, 1.28 * L], [0, 0], color="#475569", lw=1.5, zorder=3)

    # Windward wall arrows (positive -> pushing onto left wall)
    pww_hi = max(abs(x) for x in g["p_ww"])
    for yw in np.linspace(0.18 * he, 0.85 * he, 4):
        axb.annotate("", xy=(0, yw), xytext=(-0.20 * L, yw),
                     arrowprops=dict(arrowstyle="-|>,head_length=0.4,head_width=0.24", color="#2563eb", lw=1.6),
                     zorder=5)
    axb.text(-0.22 * L, he * 0.5, f"Windward\n+{pww_hi:.1f} psf", ha="right", va="center",
             fontsize=8.0, fontweight="bold", color="#1d4ed8")

    # Leeward wall arrows (suction -> pulling away to right)
    plw_min = min(g["p_lw"])
    for yw in np.linspace(0.18 * he, 0.85 * he, 4):
        axb.annotate("", xy=(L + 0.20 * L, yw), xytext=(L, yw),
                     arrowprops=dict(arrowstyle="-|>,head_length=0.4,head_width=0.24", color="#dc2626", lw=1.6),
                     zorder=5)
    axb.text(L + 0.22 * L, he * 0.5, f"Leeward\n{plw_min:.1f} psf", ha="left", va="center",
             fontsize=8.0, fontweight="bold", color="#b91c1c")

    # Roof windward & leeward arrows
    rww_lo, rww_hi = min(g["p_rww"]), max(g["p_rww"])
    rlw_lo = min(g["p_rlw"])
    # Windward slope midpoint
    xm1, ym1 = L * 0.25, 0.5 * (he + rd)
    axb.annotate("", xy=(xm1 - 0.08 * L, ym1 + 0.22 * he), xytext=(xm1, ym1),
                 arrowprops=dict(arrowstyle="-|>,head_length=0.4,head_width=0.24", color="#7c3aed", lw=1.6),
                 zorder=5)
    axb.text(xm1 - 0.06 * L, ym1 + 0.25 * he, f"WW roof\n{rww_lo:.1f} / {rww_hi:+.1f} psf",
             ha="center", va="bottom", fontsize=7.8, fontweight="semibold", color="#6d28d9")
    # Leeward slope midpoint
    xm2, ym2 = L * 0.75, 0.5 * (he + rd)
    axb.annotate("", xy=(xm2 + 0.08 * L, ym2 + 0.22 * he), xytext=(xm2, ym2),
                 arrowprops=dict(arrowstyle="-|>,head_length=0.4,head_width=0.24", color="#dc2626", lw=1.6),
                 zorder=5)
    axb.text(xm2 + 0.06 * L, ym2 + 0.25 * he, f"LW roof\n{rlw_lo:.1f} psf",
             ha="center", va="bottom", fontsize=7.8, fontweight="semibold", color="#b91c1c")

    axb.text(L / 2, he * 0.42, f"{g['roof_type']} ({g['theta']:g}°)\nL={L:g} ft, h={g['h']:.1f} ft",
             ha="center", va="center", fontsize=8.0, color="#334155")
    axb.set_xlim(-0.50 * L, 1.50 * L)
    axb.set_ylim(-0.14 * max(rd, 10), rd + 0.62 * he)
    axb.axis("off")
    axb.set_title("MWFRS Dir 1 net pressures (governing ±GCpi)", fontsize=9.2, fontweight="semibold")

    # --- Panel 1: Plan view of C&C wall Zones 4 & 5 ---
    B = g["B_ridge"]
    a = min(g["a"], B / 2, L / 2)
    axc.add_patch(Rectangle((0, 0), B, L, fc="#f8fafc", ec="#1e293b", lw=1.5, zorder=1))
    band = max(0.06 * min(B, L), 1.5)
    # Zone 4 (interior walls) in blue, Zone 5 (corner length a) in red-orange
    for (x0, y0, w0, h0) in [(a, 0, B - 2 * a, band), (a, L - band, B - 2 * a, band),
                             (0, a, band, L - 2 * a), (B - band, a, band, L - 2 * a)]:
        axc.add_patch(Rectangle((x0, y0), w0, h0, fc="#3b82f6", alpha=0.75, ec="none", zorder=2))
    for (cx, cy) in [(0, 0), (B - a, 0), (0, L - band), (B - a, L - band)]:
        axc.add_patch(Rectangle((cx, cy), a, band, fc="#ef4444", alpha=0.85, ec="none", zorder=3))
    for (cx, cy) in [(0, 0), (B - band, 0), (0, L - a), (B - band, L - a)]:
        axc.add_patch(Rectangle((cx, cy), band, a, fc="#ef4444", alpha=0.85, ec="none", zorder=3))
    cc4, cc5 = out["cc"][0], out["cc"][1]
    axc.text(B / 2, L * 0.55, f"Zone 4 (blue):\n+{cc4['p + (psf)']:.1f} / {cc4['p - (psf)']:.1f} psf",
             ha="center", va="center", fontsize=8.0, fontweight="semibold", color="#1d4ed8")
    axc.text(B / 2, L * 0.30, f"Zone 5 (red, a={g['a']:.1f} ft):\n+{cc5['p + (psf)']:.1f} / {cc5['p - (psf)']:.1f} psf",
             ha="center", va="center", fontsize=8.0, fontweight="semibold", color="#b91c1c")
    axc.set_aspect("equal")
    axc.set_xlim(-0.08 * B, 1.08 * B)
    axc.set_ylim(-0.08 * L, 1.08 * L)
    axc.set_xlabel(f"B = {B:g} ft")
    axc.set_ylabel(f"L = {L:g} ft")
    axc.set_title("C&C wall zones (plan)", fontsize=9.2, fontweight="semibold")

    # --- Panel 2: Velocity pressure profile ---
    height = max(g["h"] * 1.5, 30.0)
    z = np.linspace(0, height, 200)
    q = [qz(zi, g["V"], g["exposure"], g["Kzt"], g["ke"]) for zi in z]
    axq.plot(q, z, label="qz", color="tab:blue")
    axq.plot(np.array(q) * g["Kd"] * g["G"] * 0.8, z, label="WW wall p (ext.)", color="tab:orange")
    axq.axhline(g["h"], color="#64748b", ls=":", lw=1.0, label=f"h = {g['h']:.1f} ft")
    axq.set_xlabel("Pressure (psf)")
    axq.set_ylabel("z (ft)")
    axq.grid(alpha=0.3)
    axq.legend(fontsize=7.6)
    axq.set_title(f"qz profile (Exp {g['exposure']}, V={g['V']:g})", fontsize=9.2, fontweight="semibold")
    fig.tight_layout()
    return fig


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
