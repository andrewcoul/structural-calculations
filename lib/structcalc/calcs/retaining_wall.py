"""Cantilever retaining wall (per foot of wall length).

Stability (service loads, IBC 2024 1807.2.3): overturning and sliding
FS >= 1.5 (1.1 with seismic), soil bearing with the resultant location.
Concrete (LRFD, ASCE 7-22 2.3.1 with lateral earth/ground water pressure H
at 1.6 where it adds to the load effect): stem flexure and shear at several
heights, heel and toe flexure and shear, minimum and temperature steel.

Geometry: x from the toe (front edge of the footing), y up from the bottom
of the footing. Front face of the stem is vertical; the back face may be
battered (t_top <= t_bot). Earth pressure for stability acts on the
vertical "virtual back" plane through the heel end; soil above the heel is
part of the resisting block.

Load kinds:
    D   concrete, soil weights, stem axial dead load
    L   stem axial live load, vertical surcharge over the heel
    H   lateral earth pressure, ground water pressure and uplift, vertical
        component of the active thrust
    HL  lateral pressure from surcharges (uniform and strip)
    E   seismic earth pressure increment and wall inertia
"""

from __future__ import annotations

import math

import numpy as np

from .. import concrete as aci
from .. import foundation as fd
from ..report import Report

GAMMA_W = 62.4  # pcf

METHODS = ["Rankine", "Coulomb", "Equivalent fluid"]
SEISMIC = ["None", "Mononobe-Okabe", "Seed-Whitman"]

LRFD_COMBOS = [
    ("1.2D + 1.6H + 1.6L", {"D": 1.2, "L": 1.6, "H": 1.6, "HL": 1.6, "E": 0.0}),
    ("0.9D + 1.6H + 1.6L(lateral)", {"D": 0.9, "L": 0.0, "H": 1.6, "HL": 1.6, "E": 0.0}),
]
LRFD_SEISMIC = [
    ("1.2D + 1.0E + 1.6H + L", {"D": 1.2, "L": 1.0, "H": 1.6, "HL": 1.0, "E": 1.0}),
    ("0.9D + 1.0E + 1.6H", {"D": 0.9, "L": 0.0, "H": 1.6, "HL": 0.0, "E": 1.0}),
]


def _lateral(Hs, hw, Kh, gamma, gamma_sat, n=400):
    """Horizontal pressure profile on a vertical plane of soil height Hs (ft).

    hw: water height above the bottom of the plane (ft, 0 = dry).
    Kh: horizontal pressure coefficient applied to the vertical effective
    stress. Returns y (height above plane bottom), p_earth, p_water (psf).
    """
    y = np.linspace(0.0, Hs, n)
    depth = Hs - y
    dw = max(Hs - hw, 0.0)  # depth to water
    gp = gamma_sat - GAMMA_W
    sv = np.where(depth <= dw, gamma * depth, gamma * dw + gp * (depth - dw))
    pw = np.where(depth > dw, GAMMA_W * (depth - dw), 0.0) if hw > 0 else np.zeros_like(y)
    return y, Kh * sv, pw


def _trapz(y, x):
    fn = getattr(np, "trapezoid", None) or np.trapz
    return float(fn(y, x))


def _resultant(y, p):
    """Force (lb/ft) and height of its line of action above y=0 (ft)."""
    F = _trapz(p, y)
    if abs(F) < 1e-12:
        return 0.0, 0.0
    return F, _trapz(p * y, y) / F


def run(H_stem=10.0, H_retained=10.0, t_top=12.0, t_bot=12.0, toe=2.0, heel=4.0, h_ftg=12.0,
        key_depth=0.0, key_width=12.0, toe_cover=0.0,
        method="Rankine", gamma=120.0, phi=30.0, beta=0.0, delta=0.0, efp=40.0,
        gamma_sat=125.0, h_water=0.0,
        surcharge=0.0, strip_q=0.0, strip_offset=2.0, strip_width=3.0,
        seismic="None", kh=0.0, wall_inertia=True,
        P_D=0.0, P_L=0.0, ecc=0.0,
        qa=3.0, qa_increase_E=1.0, mu=0.45, passive_efp=250.0, passive_ignore=1.0,
        include_toe_soil=True, fs_ot=1.5, fs_sl=1.5, fs_seis=1.1,
        gamma_c=150.0, fc=3000.0, fy=60000.0,
        stem_bar=5, stem_s=12.0, stem_cover=2.0, horiz_bar=4, horiz_s=12.0,
        toe_bar=5, toe_s=12.0, heel_bar=5, heel_s=12.0, lam=1.0):
    """Units: heights, lengths ft; thicknesses, key, covers, spacings in;
    unit weights and fluid pressures pcf; surcharges psf; angles deg;
    P_D/P_L kip/ft; ecc in (+ toward the retained side); qa ksf; fc, fy psi."""
    if H_retained > H_stem + 1e-9:
        raise ValueError("retained height cannot exceed the stem height")
    if t_top > t_bot:
        raise ValueError("t_top must not exceed t_bot (front face vertical, batter on the back)")
    hf = h_ftg / 12
    tb, tt = t_bot / 12, t_top / 12
    B = toe + tb + heel
    xb = lambda z: toe + (tb - (tb - tt) * z / H_stem)  # back face x at height z above footing top
    x_bt = xb(H_retained)
    rise = math.tan(math.radians(beta)) * (B - x_bt)
    H_vb = hf + H_retained + rise  # virtual back height
    if h_water > H_vb:
        raise ValueError("water height exceeds the retained soil height")

    # earth pressure coefficients
    if method == "Rankine":
        K = fd.ka_rankine(phi, beta)
        cos_f, sin_f = math.cos(math.radians(beta)), math.sin(math.radians(beta))
        d_mo = 0.0
    elif method == "Coulomb":
        K = fd.ka_coulomb(phi, delta, beta)
        cos_f, sin_f = math.cos(math.radians(delta)), math.sin(math.radians(delta))
        d_mo = delta
    elif method == "Equivalent fluid":
        K = efp / gamma
        cos_f, sin_f = 1.0, 0.0
        d_mo = 0.0
    else:
        raise ValueError(f"method must be one of {METHODS}")
    Kh = K * cos_f

    r = Report("Cantilever Retaining Wall", code="IBC 2024 1807.2.3, ASCE 7-22, ACI 318-19 (per ft of wall)")
    r.input("Stem height / retained height", f"{H_stem:g} / {H_retained:g}", "ft")
    r.input("Stem thickness top / bottom", f"{t_top:g} / {t_bot:g}", "in")
    r.input("Toe / heel / footing thickness", f"{toe:g} ft / {heel:g} ft / {h_ftg:g} in")
    if key_depth > 0:
        r.input("Shear key (below stem) depth x width", f"{key_depth:g} x {key_width:g}", "in")
    r.input("Earth pressure method", method)
    r.input("Soil gamma / phi / backfill slope", f"{gamma:g} pcf / {phi:g} deg / {beta:g} deg")
    if method == "Coulomb":
        r.input("Wall friction delta", delta, "deg")
    if h_water > 0:
        r.input("Water height above footing bottom / gamma_sat", f"{h_water:g} ft / {gamma_sat:g} pcf")
    if surcharge:
        r.input("Uniform surcharge (L)", surcharge, "psf")
    if strip_q:
        r.input("Strip surcharge q, offset, width", f"{strip_q:g} psf, {strip_offset:g} ft, {strip_width:g} ft")
    if seismic != "None" and kh > 0:
        r.input("Seismic method / kh", f"{seismic} / {kh:g}")
    if P_D or P_L:
        r.input("Stem axial D / L, eccentricity", f"{P_D:g} / {P_L:g} kip/ft, {ecc:g} in")
    r.input("Allowable bearing / friction / passive", f"{qa:g} ksf / {mu:g} / {passive_efp:g} pcf")
    r.input("f'c / fy", f"{fc:g} / {fy:g}", "psi")

    r.result("Footing width B", B, "ft")
    r.result("Virtual back height", H_vb, "ft")
    r.result("Active coefficient K (method)", K, "", method)
    r.result("Horizontal component Kh", Kh)

    # ---------------------------------------------------------------- forces
    F = []  # (name, kind, V kip (+down), H kip (+ toward toe), x ft, y ft)

    def add(name, kind, V=0.0, H=0.0, x=0.0, y=0.0):
        F.append((name, kind, V, H, x, y))

    gc = gamma_c / 1000
    W_stem_r = gc * tt * H_stem
    add("Stem (rectangle)", "D", W_stem_r, 0, toe + tt / 2, hf + H_stem / 2)
    W_stem_t = gc * (tb - tt) * H_stem / 2
    if W_stem_t > 0:
        add("Stem (batter)", "D", W_stem_t, 0, toe + tt + (tb - tt) / 3, hf + H_stem / 3)
    add("Footing", "D", gc * B * hf, 0, B / 2, hf / 2)
    if key_depth > 0:
        add("Shear key", "D", gc * key_width / 12 * key_depth / 12, 0, toe + key_width / 24, -key_depth / 24)
    # soil above the heel (polygon, split at the water level)
    if H_retained > 0:
        # quadrilateral: back face (z = 0 .. H_retained), heel top, virtual back, sloped surface
        poly = [(xb(0.0), hf), (B, hf), (B, H_vb), (x_bt, hf + H_retained)]
        dry = fd.clip(poly, -h_water, 0, 1) if h_water > 0 else poly
        wet = fd.clip(poly, h_water, 0, -1) if h_water > 0 else []
        A, cx_, cy_ = fd.area_centroid(dry)
        if A > 0:
            add("Soil over heel", "D", gamma / 1000 * A, 0, cx_, cy_)
        A, cx_, cy_ = fd.area_centroid(wet)
        if A > 0:
            add("Soil over heel (saturated)", "D", gamma_sat / 1000 * A, 0, cx_, cy_)
    if include_toe_soil and toe_cover > 0:
        add("Soil over toe", "D", gamma / 1000 * toe * toe_cover, 0, toe / 2, hf + toe_cover / 2)
    x_ax = toe + tt / 2 + ecc / 12
    if P_D:
        add("Stem axial D", "D", P_D, 0, x_ax, hf + H_stem)
    if P_L:
        add("Stem axial L", "L", P_L, 0, x_ax, hf + H_stem)
    if surcharge:
        add("Surcharge over heel", "L", surcharge / 1000 * (B - x_bt), 0, (B + x_bt) / 2, H_vb)

    # lateral on the virtual back
    y, pe, pw = _lateral(H_vb, h_water, Kh, gamma, gamma_sat)
    Fe, ye = _resultant(y, pe)
    add("Active earth (horizontal)", "H", 0, Fe / 1000, B, ye)
    if sin_f:
        Fev = Fe / cos_f * sin_f
        add("Active earth (vertical comp.)", "H", Fev / 1000, 0, B, ye)
    if h_water > 0:
        Fw, yw = _resultant(y, pw)
        add("Water pressure", "H", 0, Fw / 1000, B, yw)
        U = 0.5 * GAMMA_W * h_water * B / 1000
        add("Water uplift", "H", -U, 0, 2 * B / 3, 0)
    if surcharge:
        add("Uniform surcharge (lateral)", "HL", 0, Kh * surcharge * H_vb / 1000, B, H_vb / 2)
    if strip_q:
        a_vb = max(strip_offset - (B - x_bt), 0.0)
        ps = fd.strip_lateral(strip_q, a_vb, strip_width, H_vb - y)
        Fs, ys = _resultant(y, ps)
        add("Strip surcharge (lateral)", "HL", 0, Fs / 1000, B, ys)

    seis = seismic != "None" and kh > 0
    dKAE = 0.0
    if seis:
        if seismic == "Mononobe-Okabe":
            kae = fd.kae_mononobe_okabe(phi, d_mo, beta, kh)
            ka_s = fd.ka_coulomb(phi, d_mo, beta)
            dKAE = kae - ka_s
            dP = 0.5 * gamma * H_vb**2 * dKAE * math.cos(math.radians(d_mo))
            r.result("KAE / KA (Mononobe-Okabe)", f"{kae:.4f} / {ka_s:.4f}", "", "AASHTO A11.3.1")
        else:
            dP = 3 / 8 * kh * gamma * H_vb**2
            dKAE = 0.75 * kh
            r.result("Seed-Whitman dKAE = 3/4 kh", dKAE)
        add("Seismic earth increment", "E", 0, dP / 1000, B, 0.6 * H_vb)
        if wall_inertia:
            for name, kind, V, H, x, yy in list(F):
                if kind == "D" and name.startswith(("Stem (", "Footing", "Shear key", "Soil over heel")):
                    add(f"Inertia: {name}", "E", 0, kh * V, x, yy)

    r.section("Service forces per ft (V + down, H + toward toe; arms from toe / footing bottom)",
              [{"Item": n, "Kind": k, "V (k)": V, "H (k)": H, "x (ft)": x, "y (ft)": yy}
               for n, k, V, H, x, yy in F])

    def totals(fac, resist_only_D=False):
        N = Hs = Mr = Mo = 0.0
        for name, kind, V, H, x, yy in F:
            f = fac.get(kind, 0.0)
            if f == 0:
                continue
            N += f * V
            Hs += f * H
            if f * V > 0:
                Mr += f * V * x
            else:
                Mo += -f * V * x
            Mo += f * H * yy
        return N, Hs, Mr, Mo

    # ------------------------------------------------------------ stability
    Dp = toe_cover + hf + key_depth / 12
    Pp = 0.5 * passive_efp * max(Dp**2 - passive_ignore**2, 0.0) / 1000 if Dp > passive_ignore else 0.0
    r.result("Passive resistance (toe + key)", Pp, "kip/ft", note=f"depth {Dp:.2f} ft, top {passive_ignore:g} ft ignored")
    cases = [("Static", {"D": 1, "H": 1, "HL": 1}, fs_ot, fs_sl, 1.0)]
    if seis:
        cases.append(("Seismic", {"D": 1, "H": 1, "HL": 1, "E": 1}, fs_seis, fs_seis, qa_increase_E))
    stab_rows = []
    for nm, fac, fso, fss, qinc in cases:
        N, Hs, Mr, Mo = totals(fac)
        FSo = Mr / Mo if Mo > 0 else math.inf
        FSs = (mu * N + Pp) / Hs if Hs > 0 else math.inf
        r.check(f"Overturning FS ({nm})", fso, FSo, "", "IBC 1807.2.3",
                note=f"Mr={Mr:.2f}, Mo={Mo:.2f} k-ft/ft")
        r.check(f"Sliding FS ({nm})", fss, FSs, "", "IBC 1807.2.3",
                note=f"mu*N={mu * N:.2f} + Pp={Pp:.2f} vs H={Hs:.2f} k/ft")
        # bearing: include live verticals
        facb = dict(fac, L=1.0)
        Nb, Hb, Mrb, Mob = totals(facb)
        xbar = (Mrb - Mob) / Nb if Nb > 0 else -1
        e_r = xbar - B / 2  # + toward heel
        p = fd.pressure_1d(Nb, e_r, B)
        if p is None:
            r.check(f"Soil bearing ({nm})", 1.0, 0.0, "ksf", "resultant outside footing")
            qtoe = qheel = float("inf")
        else:
            qtoe, qheel, cont = p
            r.result(f"Resultant from toe ({nm})", xbar, "ft",
                     note="inside middle third" if abs(e_r) <= B / 6 else "outside middle third (partial contact)")
            r.result(f"Bearing toe / heel ({nm})", f"{qtoe:.3f} / {qheel:.3f}", "ksf")
            r.check(f"Soil bearing ({nm})", max(qtoe, qheel), qa * qinc, "ksf", "IBC 1806")
        stab_rows.append({"Case": nm, "N (k)": N, "H (k)": Hs, "FS OT": FSo, "FS slide": FSs,
                          "x resultant (ft)": xbar, "q toe (ksf)": qtoe, "q heel (ksf)": qheel})
    r.section("Stability summary (service)", stab_rows)

    # ---------------------------------------------------------- concrete LRFD
    combos = LRFD_COMBOS + (LRFD_SEISMIC if seis else [])
    db_s, Ab_s = aci.bar(stem_bar)
    As_stem = Ab_s * 12 / stem_s
    # stem lateral profile on the stem back face (soil height H_retained)
    hw_s = max(h_water - hf, 0.0)
    zs, pes, pws = _lateral(H_retained, hw_s, Kh, gamma, gamma_sat, n=600) if H_retained > 0 else (
        np.array([0.0, 1e-6]), np.zeros(2), np.zeros(2))
    p_sur = np.full_like(zs, Kh * surcharge)
    if strip_q:
        p_sur = p_sur + fd.strip_lateral(strip_q, strip_offset, strip_width, H_retained - zs)
    p_seis = np.zeros_like(zs)
    if seis and H_retained > 0:
        dPs = 0.5 * gamma * H_retained**2 * dKAE * (math.cos(math.radians(d_mo))
                                                     if seismic == "Mononobe-Okabe" else 1.0)
        p0 = dPs / (2.5 * H_retained)
        p_seis = p0 * (1 + 3 * zs / H_retained)   # inverted trapezoid, resultant at 0.6 H
    z_full = np.linspace(0, H_stem, 600)
    t_z = lambda z: tb - (tb - tt) * z / H_stem
    levels = sorted({0.0} | {f * H_retained for f in (0.25, 0.5, 0.75)})
    stem_rows = []
    stem_gov = {"M": (0.0, "", 0.0), "V": (0.0, "", 0.0)}
    for cname, fac in combos:
        p = fac["H"] * (pes + pws) + fac["HL"] * p_sur + fac["E"] * p_seis   # psf on the stem
        for z0 in levels:
            m = zs >= z0
            Vz = _trapz(p[m], zs[m]) / 1000 if m.sum() > 1 else 0.0
            Mz = _trapz(p[m] * (zs[m] - z0), zs[m]) / 1000 if m.sum() > 1 else 0.0
            if fac["E"] and seis and wall_inertia:
                mi = z_full >= z0
                w_in = kh * gc * t_z(z_full[mi])
                Vz += fac["E"] * _trapz(w_in, z_full[mi])
                Mz += fac["E"] * _trapz(w_in * (z_full[mi] - z0), z_full[mi])
            Paxial = fac["D"] * P_D + fac["L"] * P_L
            if ecc < 0:
                Mz += Paxial * (-ecc) / 12
            d = t_z(z0) * 12 - stem_cover - db_s / 2
            fl = aci.flexure_rect(As_stem, 12, d, fc, fy)
            As_req = aci.As_required(Mz * 12, 12, d, fc, fy)
            stem_rows.append({"Combination": cname, "z above ftg (ft)": z0, "t (in)": t_z(z0) * 12,
                              "d (in)": d, "Mu (k-ft/ft)": Mz, "phiMn (k-ft/ft)": fl["phiMn"] / 12,
                              "As req (in2/ft)": As_req, "Vu (k/ft)": Vz})
            if z0 == 0.0 and Mz > stem_gov["M"][0]:
                stem_gov["M"] = (Mz, cname, fl["phiMn"] / 12)
        # shear at d above the base
        d0 = tb * 12 - stem_cover - db_s / 2
        zd = min(d0 / 12, H_retained)
        m = zs >= zd
        Vd = _trapz(p[m], zs[m]) / 1000 if m.sum() > 1 else 0.0
        if fac["E"] and seis and wall_inertia:
            mi = z_full >= zd
            Vd += _trapz(kh * gc * t_z(z_full[mi]), z_full[mi])
        if Vd > stem_gov["V"][0]:
            stem_gov["V"] = (Vd, cname, d0)
    d0 = tb * 12 - stem_cover - db_s / 2
    Mu, cn, phiMn = stem_gov["M"]
    r.result("Stem d at base", d0, "in")
    r.check("Stem flexure at base", Mu, phiMn, "kip-ft/ft", "ACI 22.2 / 13.3.6", note=f"[{cn}]")
    Vu, cn, _ = stem_gov["V"]
    phiVc = aci.PHI_SHEAR * aci.Vc_oneway(12, d0, fc, As=As_stem, lam=lam)
    r.check("Stem shear at d above base", Vu, phiVc, "kip/ft", "ACI 22.5.5.1 (no stirrups)", note=f"[{cn}]")
    # minimum wall reinforcement
    rho_v = 0.0012 if (stem_bar <= 5 and fy >= 60000) else 0.0015
    rho_h = 0.0020 if (horiz_bar <= 5 and fy >= 60000) else 0.0025
    r.check("Stem vertical steel minimum", rho_v * 12 * tb * 12, As_stem, "in²/ft", "ACI 11.6.1",
            note=f"rho_l >= {rho_v}, earth face bars only")
    db_h, Ab_h = aci.bar(horiz_bar)
    As_h = 2 * Ab_h * 12 / horiz_s
    t_avg = (t_top + t_bot) / 2
    r.check("Stem horizontal steel (both faces)", rho_h * 12 * t_avg, As_h, "in²/ft", "ACI 11.6.1",
            note=f"rho_t >= {rho_h}, average thickness")
    for nm_, s_ in (("stem vertical", stem_s), ("stem horizontal", horiz_s)):
        r.check(f"Bar spacing ({nm_})", s_, min(3 * t_bot, 18.0), "in", "ACI 11.7.2.1 / 11.7.3.1")

    # heel and toe
    db_t, Ab_t = aci.bar(toe_bar)
    db_hl, Ab_hl = aci.bar(heel_bar)
    As_toe = Ab_t * 12 / toe_s
    As_heel = Ab_hl * 12 / heel_s
    d_toe = h_ftg - 3.0 - db_t / 2
    d_heel = h_ftg - 2.0 - db_hl / 2
    x_heel0 = xb(0.0)
    ftg_rows = []
    G = {"Mt": (0.0, ""), "Vt": (0.0, ""), "Mh": (0.0, ""), "Vh": (0.0, "")}
    for cname, fac in combos:
        N, Hs, Mr, Mo = totals(fac)
        if N <= 0:
            continue
        xbar = (Mr - Mo) / N
        p = fd.pressure_1d(N, xbar - B / 2, B)
        if p is None:
            r.note(f"[{cname}] factored resultant outside the footing: toe/heel not designed for it.")
            continue
        ql, qr, cont = p
        xs = np.linspace(0, B, 801)
        qu = fd.pressure_1d_at(xs, ql, qr, cont, B)   # ksf (per ft of wall)
        # toe: upward pressure minus toe concrete and soil (D factor)
        wt = fac["D"] * (gc * hf + (gamma / 1000 * toe_cover if include_toe_soil else 0.0))
        mt = xs <= toe
        up = qu + (fac["H"] * GAMMA_W * h_water * xs / B / 1000 if h_water > 0 else 0.0)
        net = up[mt] - wt
        Mt = _trapz(net * (toe - xs[mt]), xs[mt])
        xv = toe - d_toe / 12
        mv = xs <= xv
        Vt = _trapz(up[mv] - wt, xs[mv]) if xv > 0 and mv.sum() > 1 else 0.0
        # heel: downward soil/concrete/surcharge/uplift minus upward pressure
        mh = xs >= x_heel0
        xh = xs[mh]
        down = fac["D"] * gc * hf * np.ones_like(xh)
        # soil column height over the heel at x (dry/wet split)
        # the heel (x >= back face at the base) has no stem above it
        hs = H_retained + np.maximum(xh - x_bt, 0) * math.tan(math.radians(beta))
        wet_h = np.clip(h_water - hf, 0.0, hs)
        down = down + fac["D"] * ((hs - wet_h) * gamma + wet_h * gamma_sat) / 1000
        down = down + fac["L"] * surcharge / 1000
        if h_water > 0:
            down = down - fac["H"] * GAMMA_W * h_water * xh / B / 1000  # uplift (linear, 0 at toe)
        netd = down - qu[mh]
        Mh = _trapz(netd * (xh - x_heel0), xh)
        # vertical component of active thrust at the heel end
        Fev = sum(V for n, k, V, H, x, yy in F if n == "Active earth (vertical comp.)")
        Mh += fac["H"] * Fev * (B - x_heel0)
        Vh = _trapz(netd, xh) + fac["H"] * Fev
        ftg_rows.append({"Combination": cname, "q toe (ksf)": ql, "q heel (ksf)": qr,
                         "Mu toe (k-ft/ft)": Mt, "Vu toe (k/ft)": Vt, "Mu heel (k-ft/ft)": Mh,
                         "Vu heel (k/ft)": Vh})
        for k, v in (("Mt", Mt), ("Vt", Vt), ("Mh", Mh), ("Vh", Vh)):
            if v > G[k][0]:
                G[k] = (v, cname)
    As_min_f = aci.As_min_slab(12, h_ftg, fy)
    for lab, Mk, Vk, As, d, sp in (("Toe", "Mt", "Vt", As_toe, d_toe, toe_s),
                                   ("Heel", "Mh", "Vh", As_heel, d_heel, heel_s)):
        Mu, cn = G[Mk]
        fl = aci.flexure_rect(As, 12, d, fc, fy)
        As_req = aci.As_required(Mu * 12, 12, d, fc, fy)
        r.result(f"{lab}: d, As req / provided", f"{d:.2f} in, {As_req:.3f} / {As:.3f} in²/ft")
        r.check(f"{lab} flexure", Mu, fl["phiMn"] / 12, "kip-ft/ft", "ACI 22.2",
                note=f"[{cn}] {'bottom' if lab == 'Toe' else 'top'} bars")
        Vu, cn = G[Vk]
        phiVc = aci.PHI_SHEAR * aci.Vc_oneway(12, d, fc, As=As, lam=lam)
        r.check(f"{lab} shear", Vu, phiVc, "kip/ft", "ACI 22.5.5.1",
                note=f"[{cn}] {'at d from stem face' if lab == 'Toe' else 'at stem face'}")
        r.check(f"{lab} minimum steel", As_min_f, As, "in²/ft", "ACI 7.6.1.1 (each face, conservative)")
        r.check(f"Bar spacing ({lab.lower()})", sp, min(3 * h_ftg, 18.0), "in", "ACI 7.7.2.3")

    r.section("Stem design by height (LRFD)", stem_rows)
    r.section("Footing design (LRFD)", ftg_rows)
    r.note("ASCE 7-22 2.3.1: lateral earth and ground water pressure H uses 1.6 where it adds to the "
           "load effect; the vertical component of the active thrust also uses the H factor (simplification).")
    r.note("Stability uses service loads; live verticals (axial L, surcharge over the heel) are not "
           "counted as resisting overturning or sliding, but are included for bearing.")
    r.note("Earth pressure for stability acts on the vertical plane through the heel; stem design uses "
           "the pressure on the stem back face from the retained height.")
    r.note("Seismic increment acts at 0.6H (inverted trapezoid on the stem); M-O uses the moist unit "
           "weight even below water. Wall inertia includes stem, footing, key, and soil over the heel.")
    r.note("Water: uplift varies linearly from gamma_w*h_w at the heel to 0 at the toe (no water in front).")
    r.note("Strip surcharge: elastic rigid-wall solution; its vertical load over the heel is ignored.")
    r.note("Stem axial load eccentricity toward the front (negative) adds P*e to the stem moment; toward "
           "the back it is ignored (conservative). Bar cut-off, lap splices, and dowel development are "
           "not checked; use the stem-by-height table for cut-off decisions.")
    geom = dict(B=B, toe=toe, tt=tt, tb=tb, hf=hf, H_stem=H_stem, H_retained=H_retained, x_bt=x_bt,
                rise=rise, H_vb=H_vb, key_depth=key_depth / 12, key_width=key_width / 12,
                h_water=h_water, toe_cover=toe_cover, y=y, pe=pe, pw=pw,
                surcharge=surcharge, strip_q=strip_q, strip_offset=strip_offset, strip_width=strip_width,
                Dp=Dp, passive_ignore=passive_ignore, Pp=Pp, passive_efp=passive_efp,
                P_D=P_D, P_L=P_L, ecc=ecc)
    return {"report": r, "forces": F, "geom": geom, "stem_rows": stem_rows, "ftg_rows": ftg_rows,
            "stability": stab_rows}


def plot(out):
    """Wall section with service lateral, surcharge, passive, and bearing pressures."""
    import matplotlib.pyplot as plt

    g = out["geom"]
    fig, ax = plt.subplots(figsize=(8.6, 6.8))
    B, toe, tt, tb, hf = g["B"], g["toe"], g["tt"], g["tb"], g["hf"]
    # concrete
    ftg = [(0, 0), (B, 0), (B, hf), (toe + tb, hf), (toe + tt, hf + g["H_stem"]),
           (toe, hf + g["H_stem"]), (toe, hf), (0, hf)]
    ax.fill(*zip(*ftg), color="0.75", ec="k", zorder=3)
    if g["key_depth"] > 0:
        ax.fill([toe, toe + g["key_width"], toe + g["key_width"], toe],
                [0, 0, -g["key_depth"], -g["key_depth"]], color="0.75", ec="k", zorder=3)
    # soil surface
    slope = g["rise"] / max(B - g["x_bt"], 1e-9)
    ax.plot([g["x_bt"], B, B + 2.5], [hf + g["H_retained"], g["H_vb"], g["H_vb"] + 2.5 * slope],
            color="saddlebrown", lw=2, zorder=4)
    if g["toe_cover"] > 0:
        ax.plot([-1.8, toe], [hf + g["toe_cover"]] * 2, color="saddlebrown", lw=2, zorder=4)
    if g["h_water"] > 0:
        ax.plot([B - 0.5, B + 3], [g["h_water"]] * 2, color="tab:blue", ls="--", zorder=4)
        ax.text(B + 2.8, g["h_water"] + 0.15, "GWT", fontsize=8, color="tab:blue", ha="right")

    # Uniform surcharge arrows on backfill surface
    if g.get("surcharge", 0) > 0:
        xs_s = np.linspace(g["x_bt"] + 0.2, B + 1.8, 6)
        for xs in xs_s:
            ys = hf + g["H_retained"] + (xs - g["x_bt"]) * slope
            ax.annotate("", xy=(xs, ys + 0.05), xytext=(xs, ys + 0.95),
                        arrowprops=dict(arrowstyle="-|>,head_length=0.3,head_width=0.18", color="#2563eb", lw=1.2),
                        zorder=5)
        ax.plot([xs_s[0], xs_s[-1]],
                [hf + g["H_retained"] + (xs_s[0] - g["x_bt"]) * slope + 0.95,
                 hf + g["H_retained"] + (xs_s[-1] - g["x_bt"]) * slope + 0.95],
                color="#2563eb", lw=1.4, label=f"surcharge ({g['surcharge']:g} psf)", zorder=5)

    # Strip surcharge arrows
    if g.get("strip_q", 0) > 0:
        x1 = g["x_bt"] + g["strip_offset"]
        x2 = x1 + g["strip_width"]
        for xs in np.linspace(x1, x2, 4):
            ys = hf + g["H_retained"] + (xs - g["x_bt"]) * slope
            ax.annotate("", xy=(xs, ys + 0.05), xytext=(xs, ys + 1.35),
                        arrowprops=dict(arrowstyle="-|>,head_length=0.32,head_width=0.2", color="#7c3aed", lw=1.5),
                        zorder=6)
        ax.text(0.5 * (x1 + x2), hf + g["H_retained"] + (0.5 * (x1 + x2) - g["x_bt"]) * slope + 1.45,
                f"strip {g['strip_q']:g} psf", ha="center", va="bottom", fontsize=8, color="#7c3aed", fontweight="semibold")

    # Stem axial load
    if g.get("P_D", 0) or g.get("P_L", 0):
        x_ax = toe + tt / 2 + g.get("ecc", 0) / 12
        y_top = hf + g["H_stem"]
        ax.annotate("", xy=(x_ax, y_top + 0.05), xytext=(x_ax, y_top + 1.3),
                    arrowprops=dict(arrowstyle="-|>,head_length=0.45,head_width=0.26", color="#dc2626", lw=2.0),
                    zorder=7)
        ax.text(x_ax, y_top + 1.38, f"P = {g['P_D'] + g['P_L']:g} k/ft", ha="center", va="bottom",
                fontsize=8, color="#dc2626", fontweight="bold")

    # Passive pressure wedge at toe
    if g.get("Pp", 0) > 0:
        y_top_p = hf + g["toe_cover"] - g["passive_ignore"]
        y_bot_p = -g["key_depth"]
        if y_top_p > y_bot_p:
            p_top_v = g["passive_efp"] * g["passive_ignore"]
            p_bot_v = g["passive_efp"] * g["Dp"]
            sc_p = 1.4 / max(p_bot_v, 1e-9)
            ax.fill([0, -p_top_v * sc_p, -p_bot_v * sc_p, 0],
                    [y_top_p, y_top_p, y_bot_p, y_bot_p],
                    color="#0d9488", alpha=0.35, ec="#0f766e", lw=1.2,
                    label=f"passive Pp={g['Pp']:.2f} k/ft", zorder=4)

    # pressure diagram on the virtual back (scaled)
    pe, pw, y = g["pe"], g["pw"], g["y"]
    ptot = pe + pw
    sc = 2.5 / max(float(ptot.max()), 1e-9)
    ax.fill_betweenx(y, B, B + pe * sc, color="tab:orange", alpha=0.5, label="earth (psf)")
    if pw.max() > 0:
        ax.fill_betweenx(y, B + pe * sc, B + ptot * sc, color="tab:blue", alpha=0.4, label="water (psf)")
    ax.text(B + ptot[0] * sc + 0.1, 0.1, f"{ptot[0]:.0f} psf", fontsize=8)
    # bearing pressure (static)
    st = out["stability"][0]
    if math.isfinite(st["q toe (ksf)"]):
        qmax = max(st["q toe (ksf)"], st["q heel (ksf)"], 1e-9)
        s2 = 1.5 / qmax
        xs = np.linspace(0, B, 100)
        p = fd.pressure_1d(1.0, st["x resultant (ft)"] - B / 2, B)
        if p is not None:
            ql, qr, cont = p
            qq = fd.pressure_1d_at(xs, ql, qr, cont, B)
            qq = qq / max(ql, qr) * qmax
            ax.fill_between(xs, -qq * s2 - 0.3, -0.3, color="tab:green", alpha=0.4, label="bearing (ksf)")
            ax.text(0, -st["q toe (ksf)"] * s2 - 0.7, f"{st['q toe (ksf)']:.2f}", fontsize=8)
            ax.text(B - 0.6, -st["q heel (ksf)"] * s2 - 0.7, f"{st['q heel (ksf)']:.2f}", fontsize=8)
    ax.set_aspect("equal")
    ax.set_xlabel("x from toe (ft)")
    ax.set_ylabel("y from footing bottom (ft)")
    ax.set_title("Retaining wall section & service pressures (static)")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    return fig
