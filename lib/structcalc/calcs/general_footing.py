"""General (isolated rectangular spread) footing.

Soil: ASCE 7-22 2.4 ASD combinations - bearing with biaxial eccentricity
(exact partial contact), overturning, sliding, uplift.
Concrete: ASCE 7-22 2.3 LRFD combinations - ACI 318-19 flexure at the
column face, one-way shear (Table 22.5.5.1, no stirrups, size effect),
two-way shear (22.6) with unbalanced moment (8.4.4.2), bearing (22.8),
minimum steel (7.6.1.1 / 24.4.3.2), spacing (7.7.2.3), short-direction band
(13.3.3.3), development length (25.4.2.3).

Coordinates: origin at the footing center, x along Lx, y along Ly.
Sign convention: +My and +Vx move the soil resultant toward +x;
+Mx and +Vy move it toward +y. Column is concentric with the footing.
"""

from __future__ import annotations

import math

import numpy as np

from .. import concrete as aci
from .. import foundation as fd
from .. import loads
from ..report import Report

CASE_ORDER = ("D", "L", "Lr", "S", "R", "W", "E")


def parse_load_string(s):
    """'P, Mx, My, Vx, Vy' (kip, kip-ft) -> tuple of 5 floats. Missing = 0."""
    vals = [float(v) for v in str(s or "").replace(";", ",").split(",") if v.strip()]
    if len(vals) > 5:
        raise ValueError(f"load '{s}': at most 5 values (P, Mx, My, Vx, Vy)")
    return tuple(vals + [0.0] * (5 - len(vals)))


def _combine(col_loads, f):
    out = np.zeros(5)
    for c, fac in f.items():
        if c in col_loads:
            out += fac * np.array(col_loads[c])
    return out


def run(Lx=8.0, Ly=8.0, h=24.0, col_x=18.0, col_y=18.0, ped_height=0.0,
        soil_cover=1.0, gamma_c=150.0, gamma_s=110.0,
        qa=3.0, qa_increase_WE=1.0, qa_is_net=False,
        mu=0.35, passive_pcf=0.0, fs_overturn=1.0, fs_sliding=1.0,
        fc=3000.0, fy=60000.0, bar_x=5, s_x=12.0, bar_y=5, s_y=12.0, cover=3.0,
        col_loads=None, lam=1.0):
    """Analyze and design the footing.

    Units: Lx, Ly, soil_cover ft; h, col_x, col_y, ped_height, s_x, s_y, cover in;
    unit weights pcf; qa ksf; passive pcf (equivalent fluid); fc, fy psi.
    col_loads: {case: (P, Mx, My, Vx, Vy)} service loads at the top of the
    pedestal (kip, kip-ft), P + down.
    """
    col_loads = {c: tuple(v) for c, v in (col_loads or {}).items() if any(abs(x) > 0 for x in v)}
    present = [c for c in CASE_ORDER if c in col_loads]
    A = Lx * Ly
    hf = h / 12.0
    hp = ped_height / 12.0
    cx_ft, cy_ft = col_x / 12.0, col_y / 12.0
    W_ftg = gamma_c * A * hf / 1000.0
    W_ped = gamma_c * cx_ft * cy_ft * hp / 1000.0
    W_soil = gamma_s * soil_cover * (A - (cx_ft * cy_ft if hp > 0 else 0.0)) / 1000.0
    W = W_ftg + W_ped + W_soil
    h_v = hp + hf  # lever arm of pedestal-top shears to the footing base (ft)

    r = Report("General Footing", code="ACI 318-19, ASCE 7-22 (ASD soil / LRFD concrete), IBC 2024")
    r.input("Footing Lx x Ly x h", f"{Lx:g} ft x {Ly:g} ft x {h:g} in")
    r.input("Column / pedestal", f"{col_x:g} x {col_y:g} in, pedestal height {ped_height:g} in")
    r.input("Soil cover over footing", soil_cover, "ft")
    r.input("Unit weights concrete / soil", f"{gamma_c:g} / {gamma_s:g}", "pcf")
    r.input("Allowable soil pressure qa", qa, "ksf", note="net" if qa_is_net else "gross")
    r.input("qa increase for W/E combinations", qa_increase_WE)
    r.input("Friction coefficient / passive", f"{mu:g} / {passive_pcf:g} pcf")
    r.input("f'c / fy", f"{fc:g} / {fy:g}", "psi")
    r.input("Bars along x (bottom layer)", f"#{bar_x} @ {s_x:g} in")
    r.input("Bars along y (second layer)", f"#{bar_y} @ {s_y:g} in")
    r.input("Bottom cover", cover, "in")
    for c in present:
        P, Mx, My, Vx, Vy = col_loads[c]
        r.input(f"Load {c}: P, Mx, My, Vx, Vy", f"{P:g}, {Mx:g}, {My:g}, {Vx:g}, {Vy:g}", "kip, kip-ft")

    r.result("Footing weight", W_ftg, "kip")
    r.result("Pedestal weight", W_ped, "kip")
    r.result("Soil weight over footing", W_soil, "kip")

    # ------------------------------------------------------------------ ASD
    asd_rows = []
    gov = {"bear": None, "ot": None, "sl": None, "up": None}
    qa_net_add = W / A if qa_is_net else 0.0
    Pp_x = passive_pcf / 2 * ((soil_cover + hf) ** 2 - soil_cover**2) * Ly / 1000.0
    Pp_y = passive_pcf / 2 * ((soil_cover + hf) ** 2 - soil_cover**2) * Lx / 1000.0
    pressures = {}
    for name, f in loads.active(loads.asd(), present):
        P, Mx, My, Vx, Vy = (float(v) for v in _combine(col_loads, f))
        fD = f.get("D", 0.0)
        Ptot = P + fD * W
        Mxb = Mx + Vy * h_v
        Myb = My + Vx * h_v
        wind_eq = any(f.get(k, 0) != 0 for k in ("W", "E"))
        qall = qa * (qa_increase_WE if wind_eq else 1.0) + qa_net_add
        row = {"Combination": name, "P total (k)": Ptot, "ex (ft)": float("nan"),
               "ey (ft)": float("nan"), "q max (ksf)": float("nan"), "Contact %": 0.0,
               "FS overturn": float("nan"), "FS slide": float("nan")}
        if P < 0:
            ratio = -P / (fD * W) if fD * W > 0 else math.inf
            if gov["up"] is None or ratio > gov["up"][0]:
                gov["up"] = (ratio, name, -P, fD * W)
        if Ptot <= 0:
            row["q max (ksf)"] = float("inf")
            asd_rows.append(row)
            if gov["bear"] is None or gov["bear"][0] < math.inf:
                gov["bear"] = (math.inf, name, math.inf, qall, None)
            continue
        ex, ey = Myb / Ptot, Mxb / Ptot
        row["ex (ft)"], row["ey (ft)"] = ex, ey
        pr = fd.bearing_pressure(Ptot, ex, ey, Lx, Ly)
        fs_ot = min(Ptot * Lx / 2 / abs(Myb) if abs(Myb) > 1e-9 else math.inf,
                    Ptot * Ly / 2 / abs(Mxb) if abs(Mxb) > 1e-9 else math.inf)
        row["FS overturn"] = fs_ot
        if gov["ot"] is None or fs_ot < gov["ot"][0]:
            gov["ot"] = (fs_ot, name)
        fs_sl = min((mu * Ptot + Pp_x) / abs(Vx) if abs(Vx) > 1e-9 else math.inf,
                    (mu * Ptot + Pp_y) / abs(Vy) if abs(Vy) > 1e-9 else math.inf)
        row["FS slide"] = fs_sl
        if gov["sl"] is None or fs_sl < gov["sl"][0]:
            gov["sl"] = (fs_sl, name)
        if pr is None:
            row["q max (ksf)"] = float("inf")
            ratio = math.inf
            gov["bear"] = (math.inf, name, math.inf, qall, None)
        else:
            row["q max (ksf)"] = pr.qmax
            row["Contact %"] = 100 * pr.contact_ratio
            ratio = pr.qmax / qall
            pressures[name] = pr
            if gov["bear"] is None or ratio > gov["bear"][0]:
                gov["bear"] = (ratio, name, pr.qmax, qall, pr)
        asd_rows.append(row)

    if gov["bear"] is not None:
        ratio, name, qmax, qall, pr = gov["bear"]
        r.result("Governing bearing combination", name)
        if pr is not None:
            r.result("Eccentricity ex, ey", f"{pr.ex:.3f}, {pr.ey:.3f}", "ft")
            r.result("Resultant inside kern (full contact)", bool(fd.in_kern(pr.ex, pr.ey, Lx, Ly)), "",
                     "|ex|/Lx + |ey|/Ly <= 1/6")
            r.result("Contact area", 100 * pr.contact_ratio, "%")
            r.result("q min", pr.qmin, "ksf")
        r.check("Soil bearing", qmax, qall, "ksf", "IBC 1806 / ASCE 7-22 2.4", note=f"[{name}]")
    if gov["ot"] is not None and math.isfinite(gov["ot"][0]):
        r.check("Overturning (FS required vs provided)", fs_overturn, gov["ot"][0], "",
                "ASD combinations", note=f"[{gov['ot'][1]}]")
    if gov["sl"] is not None and math.isfinite(gov["sl"][0]):
        r.check("Sliding (FS required vs provided)", fs_sliding, gov["sl"][0], "",
                "ASD combinations", note=f"[{gov['sl'][1]}]")
    if gov["up"] is not None:
        _, name, up, res = gov["up"]
        r.check("Uplift (upward column load vs footing+soil weight)", up, res, "kip",
                "ASCE 7-22 2.4 (0.6D)", note=f"[{name}]")

    # ------------------------------------------------------------------ LRFD
    dbx, Abx = aci.bar(bar_x)
    dby, Aby = aci.bar(bar_y)
    dx = h - cover - dbx / 2          # in
    dy = h - cover - dbx - dby / 2    # in
    davg = (dx + dy) / 2
    r.result("d (x bars, bottom) / d (y bars)", f"{dx:.2f} / {dy:.2f}", "in")
    As_x = Abx * (Ly * 12) / s_x       # bars along x, spread over width Ly
    As_y = Aby * (Lx * 12) / s_y
    x0, y0 = cx_ft / 2, cy_ft / 2
    lrfd_rows = []
    G = {k: (0.0, "") for k in ("Mx+", "Mx-", "My+", "My-", "Vx", "Vy", "vu", "Pu")}
    punch_detail = None
    for name, f in loads.active(loads.lrfd(), present):
        P, Mx, My, Vx, Vy = (float(v) for v in _combine(col_loads, f))
        fD = f.get("D", 0.0)
        if P > G["Pu"][0]:
            G["Pu"] = (P, name)
        Ptot = P + fD * W
        if P <= 0 or Ptot <= 0:
            continue
        Mxb, Myb = Mx + Vy * h_v, My + Vx * h_v
        pr = fd.bearing_pressure(Ptot, Myb / Ptot, Mxb / Ptot, Lx, Ly)
        if pr is None:
            continue
        wu = fD * W / A  # uniform downward (ksf)

        # flexure + one-way shear, x direction (bars along x)
        out = {}
        for sgn, key in ((1, "x+"), (-1, "x-")):
            Q, Qx, _ = pr.integrals([(-x0, sgn, 0.0)])
            arm = sgn * Qx - x0 * Q
            proj = Lx / 2 - x0
            out["M" + key] = arm - wu * Ly * proj**2 / 2
            xs = x0 + dx / 12
            if xs < Lx / 2:
                Q2, _, _ = pr.integrals([(-xs, sgn, 0.0)])
                out["V" + key] = Q2 - wu * Ly * (Lx / 2 - xs)
            else:
                out["V" + key] = 0.0
        for sgn, key in ((1, "y+"), (-1, "y-")):
            Q, _, Qy = pr.integrals([(-y0, 0.0, sgn)])
            arm = sgn * Qy - y0 * Q
            proj = Ly / 2 - y0
            out["M" + key] = arm - wu * Lx * proj**2 / 2
            ys = y0 + dy / 12
            if ys < Ly / 2:
                Q2, _, _ = pr.integrals([(-ys, 0.0, sgn)])
                out["V" + key] = Q2 - wu * Lx * (Ly / 2 - ys)
            else:
                out["V" + key] = 0.0
        Mux = max(out["Mx+"], out["Mx-"])
        Muy = max(out["My+"], out["My-"])
        Vux = max(out["Vx+"], out["Vx-"])
        Vuy = max(out["Vy+"], out["Vy-"])
        # two-way shear (inches for the perimeter geometry)
        ps = fd.punching_section(0.0, 0.0, col_x, col_y, davg, (-Lx * 6, Lx * 6), (-Ly * 6, Ly * 6))
        psy = fd.punching_section(0.0, 0.0, col_y, col_x, davg, (-Ly * 6, Ly * 6), (-Lx * 6, Lx * 6))
        hx1, hx2, hy1, hy2 = ps["x1"] / 12, ps["x2"] / 12, ps["y1"] / 12, ps["y2"] / 12
        Qin, _, _ = pr.integrals([(-hx1, 1, 0), (hx2, -1, 0), (-hy1, 0, 1), (hy2, 0, -1)])
        A_in = (hx2 - hx1) * (hy2 - hy1)
        Vu2 = P - Qin + wu * A_in
        My_top = abs(My + Vx * hp) * 12  # kip-in (moment in x-z plane)
        Mx_top = abs(Mx + Vy * hp) * 12
        gvx = fd.gamma_v(ps["b1"], ps["b2"])
        gvy = fd.gamma_v(psy["b1"], psy["b2"])
        vu = (Vu2 * 1000 / (ps["bo"] * davg)
              + gvx * My_top * 1000 * ps["c"] / ps["Jc"]
              + gvy * Mx_top * 1000 * psy["c"] / psy["Jc"])
        lrfd_rows.append({"Combination": name, "Pu (k)": P, "Mu x-bars (k-ft)": Mux,
                          "Mu y-bars (k-ft)": Muy, "Vu x (k)": Vux, "Vu y (k)": Vuy,
                          "Vu 2-way (k)": Vu2, "vu (psi)": vu})
        for k, v in (("Mx+", Mux), ("My+", Muy), ("Vx", Vux), ("Vy", Vuy), ("vu", vu)):
            if v > G[k][0]:
                G[k] = (v, name)
                if k == "vu":
                    punch_detail = (ps, Vu2, My_top, Mx_top, gvx, gvy)

    if lrfd_rows:
        # flexure
        for key, As, d, b_ft, s, db, lab in (("Mx+", As_x, dx, Ly, s_x, dbx, "x"),
                                             ("My+", As_y, dy, Lx, s_y, dby, "y")):
            Mu, cname = G[key]
            b = b_ft * 12
            fl = aci.flexure_rect(As, b, d, fc, fy)
            As_req = aci.As_required(Mu * 12, b, d, fc, fy)
            As_min = aci.As_min_slab(b, h, fy)
            r.result(f"Mu at column face, bars along {lab}", Mu, "kip-ft", note=f"[{cname}]")
            r.result(f"As required / As,min / As provided ({lab})",
                     f"{As_req:.2f} / {As_min:.2f} / {As:.2f}", "in²")
            r.result(f"phi ({lab} bars), eps_t", f"{fl['phi']:.3f}, {fl['eps_t']:.4f}", "", "ACI 21.2.2")
            r.check(f"Flexure, bars along {lab}", Mu, fl["phiMn"] / 12, "kip-ft", "ACI 22.2 / 13.2.7.1",
                    note=f"[{cname}]")
            r.check(f"Minimum steel ({lab})", As_min, As, "in²", "ACI 7.6.1.1 / 24.4.3.2")
            r.check(f"Bar spacing ({lab})", s, min(3 * h, 18.0), "in", "ACI 7.7.2.3")
            avail = ((Lx if lab == "x" else Ly) * 12 - (col_x if lab == "x" else col_y)) / 2 - 3.0
            ld = aci.ld_simplified(db, fc, fy, lam=lam)
            r.check(f"Development length ({lab} bars)", ld, avail, "in", "ACI 25.4.2.3 / 13.2.8.3",
                    note="straight bars, clear spacing >= 2db, cover >= db")
        # short-direction band
        if abs(Lx - Ly) > 1e-9:
            beta = max(Lx, Ly) / min(Lx, Ly)
            if Lx > Ly:   # short-direction bars are the y bars, band width = Ly along x
                Mu, _ = G["My+"]
                b = Lx * 12
                As_tot = max(aci.As_required(Mu * 12, b, dy, fc, fy), aci.As_min_slab(b, h, fy))
                As_band_prov = Aby * Ly * 12 / s_y
            else:
                Mu, _ = G["Mx+"]
                b = Ly * 12
                As_tot = max(aci.As_required(Mu * 12, b, dx, fc, fy), aci.As_min_slab(b, h, fy))
                As_band_prov = Abx * Lx * 12 / s_x
            As_band = 2 / (beta + 1) * As_tot
            r.check("Short-direction band steel", As_band, As_band_prov, "in²", "ACI 13.3.3.3",
                    note=f"beta={beta:.2f}, band width = short side, uniform spacing assumed")
        # one-way shear
        for key, As, d, b_ft, lab in (("Vx", As_x, dx, Ly, "x"), ("Vy", As_y, dy, Lx, "y")):
            Vu, cname = G[key]
            b = b_ft * 12
            Vc = aci.Vc_oneway(b, d, fc, As=As, lam=lam)
            r.check(f"One-way shear, {lab} direction", Vu, aci.PHI_SHEAR * Vc, "kip",
                    "ACI 22.5.5.1 (Av < Av,min, size effect)", note=f"[{cname}] at d from face")
        # two-way shear
        vu, cname = G["vu"]
        if punch_detail is not None:
            ps, Vu2, Myt, Mxt, gvx, gvy = punch_detail
            beta_c = max(col_x, col_y) / min(col_x, col_y)
            vc = aci.vc_twoway(davg, ps["bo"], beta_c, ps["alpha_s"], fc, lam)
            r.result("Punching: bo / d avg", f"{ps['bo']:.1f} / {davg:.2f}", "in", "ACI 22.6.4.1")
            r.result("Punching: Vu, gamma_v (x, y)", f"{Vu2:.1f} k, {gvx:.3f}, {gvy:.3f}", "",
                     "ACI 8.4.4.2")
            r.result("lambda_s", aci.lambda_s(davg), "", "ACI 22.5.5.1.3")
            r.check("Two-way (punching) shear", vu, aci.PHI_SHEAR * vc, "psi", "ACI 22.6.5.2 / 8.4.4.2",
                    note=f"[{cname}]")
        # bearing
        Pu, cname = G["Pu"]
        A1 = col_x * col_y
        k = fd.bearing_A2_factor(col_x, col_y, h, Lx * 12, Ly * 12)
        phiBn = aci.PHI_BEARING * 0.85 * fc * A1 * k / 1000
        r.check("Bearing on footing concrete", Pu, phiBn, "kip", "ACI 22.8.3.2",
                note=f"[{cname}] sqrt(A2/A1)={k:.2f}")
    else:
        r.note("No LRFD combination puts the footing in compression: concrete checks not done.")

    r.section("ASD combinations (soil)", asd_rows)
    if lrfd_rows:
        r.section("LRFD combinations (concrete)", lrfd_rows)
    r.note("Soil pressure: rigid footing, linear pressure, no tension; partial contact solved exactly.")
    r.note("Column moments and shears are applied at the top of the pedestal; shears add V x (pedestal + "
           "footing height) at the base for soil checks and V x pedestal height at the footing top.")
    r.note("Overturning/sliding FS use ASD combinations (0.6D already reduces the resisting load); "
           "set FS required = 1.5 for traditional service-load practice.")
    r.note("Sliding is checked separately in x and y with full friction in each direction.")
    r.note("Two-way shear: perimeter at d/2 (average d); unbalanced moment = column moment at the "
           "footing top; the soil pressure inside the perimeter is subtracted from Vu only.")
    r.note("Bars along x are assumed in the bottom layer. Hooks, dowels, and the pedestal "
           "itself are not designed.")
    gb = gov["bear"]
    return {"report": r, "pressure": gb[4] if gb else None, "pressures": pressures,
            "Lx": Lx, "Ly": Ly, "col": (col_x, col_y), "W": W}


def plot(out, n=81):
    """Plan view of the governing soil pressure (ksf) with the zero-pressure line."""
    import matplotlib.pyplot as plt

    pr = out["pressure"]
    Lx, Ly = out["Lx"], out["Ly"]
    fig, ax = plt.subplots(figsize=(6.5, 6.5 * min(max(Ly / Lx, 0.4), 1.6)))
    ax.add_patch(plt.Rectangle((-Lx / 2, -Ly / 2), Lx, Ly, fill=False, lw=1.5))
    cx, cy = out["col"][0] / 12, out["col"][1] / 12
    ax.add_patch(plt.Rectangle((-cx / 2, -cy / 2), cx, cy, color="0.5", alpha=0.6))
    if pr is not None:
        xs = np.linspace(-Lx / 2, Lx / 2, n)
        ys = np.linspace(-Ly / 2, Ly / 2, n)
        X, Y = np.meshgrid(xs, ys)
        Z = pr.q(X, Y)
        cs = ax.contourf(X, Y, Z, levels=15, cmap="YlOrRd", alpha=0.85)
        fig.colorbar(cs, ax=ax, label="soil pressure (ksf)")
        if not pr.full_contact:
            ax.contour(X, Y, pr.a + pr.b * X + pr.c * Y, levels=[0], colors="k", linestyles="--")
        ax.plot([pr.ex], [pr.ey], "k+", ms=12, mew=2)
        ax.set_title(f"Governing service soil pressure\nqmax = {pr.qmax:.2f} ksf, "
                     f"{100 * pr.contact_ratio:.0f}% contact (+ = resultant, -- = zero line)", fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlabel("x (ft)")
    ax.set_ylabel("y (ft)")
    ax.set_xlim(-Lx / 2 - 0.3, Lx / 2 + 0.3)
    ax.set_ylim(-Ly / 2 - 0.3, Ly / 2 + 0.3)
    fig.tight_layout()
    return fig
