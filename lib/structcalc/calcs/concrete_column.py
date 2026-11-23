"""Concrete column: ACI 318-19 P-M interaction (strain compatibility),
biaxial bending (Bresler reciprocal load / linear contour), nonsway
slenderness by moment magnification (6.6.4), reinforcement limits and
tie/spiral detailing. ASCE 7-22 LRFD combinations.

Units at the interface: b, h, D, cover in in; lu in ft; fc/fy in psi;
loads in kip and kip-ft.

Axes: b is the column dimension along x, h along y. "Mx" bends about the
x-axis (compression on the +y face, depth = h). End moments with the same
sign put the column in single curvature.
"""

from __future__ import annotations

import math

import numpy as np

from .. import concrete as C
from .. import concrete_member as cm
from .. import loads
from ..report import Report

SHAPES = ["Rectangular tied", "Rectangular spiral", "Circular spiral", "Circular tied"]
EI_METHODS = ["(0.2 Ec Ig + Es Ise)/(1+beta_dns)", "0.4 Ec Ig/(1+beta_dns)"]


def parse_case(text):
    """'P, Mx_top, Mx_bot, My_top, My_bot' -> tuple of 5 floats (blank = 0)."""
    vals = [v.strip() for v in (text or "").replace(";", ",").split(",")]
    nums = [float(v) if v else 0.0 for v in vals if v != ""]
    nums += [0.0] * (5 - len(nums))
    if len(nums) > 5:
        raise ValueError(f"load '{text}': give at most 5 values P, Mxt, Mxb, Myt, Myb")
    return tuple(nums)


def bar_layout(shape, b, h, D, cover, tie_db, db, nx, ny, n_circ):
    """Bar centers (x, y) in inches from the column center."""
    if shape.startswith("Rectangular"):
        if nx < 2 or ny < 2:
            raise ValueError("bars per face must be at least 2 (corners included)")
        xb = b / 2 - cover - tie_db - db / 2
        yb = h / 2 - cover - tie_db - db / 2
        pts = []
        for i in range(nx):
            x = -xb + 2 * xb * i / (nx - 1)
            pts += [(x, yb), (x, -yb)]
        for j in range(1, ny - 1):
            y = -yb + 2 * yb * j / (ny - 1)
            pts += [(xb, y), (-xb, y)]
        return pts
    Rb = D / 2 - cover - tie_db - db / 2
    return [(Rb * math.cos(math.pi / 2 + math.pi / n_circ + 2 * math.pi * i / n_circ),
             Rb * math.sin(math.pi / 2 + math.pi / n_circ + 2 * math.pi * i / n_circ))
            for i in range(n_circ)]


def _magnify(Pu, Mt, Mb, klu_in, r, h_dir, EI, beta_dns, axis):
    """ACI 6.6.4 nonsway magnification about one axis. Returns dict."""
    a, c = abs(Mt), abs(Mb)
    M2 = max(a, c)
    M1 = min(a, c)
    if M2 < 1e-9:
        r12 = -1.0  # no end moments: treat as single curvature (conservative)
    else:
        same = (Mt * Mb) >= 0
        r12 = -(M1 / M2) if same else (M1 / M2)
    slr = klu_in / r
    limit = min(34 + 12 * r12, 40.0)
    out = {"klu/r": slr, "limit": limit, "slender": slr > limit and Pu > 0, "delta": 1.0,
           "Mc": M2, "Mc_nomin": M2, "Cm": math.nan, "Pc": math.nan, "M2min_used": False,
           "r12": r12}
    if not out["slender"]:
        return out
    EIeff = EI / (1 + beta_dns)                  # kip-in^2
    Pc = math.pi**2 * EIeff / klu_in**2
    M2min = Pu * (0.6 + 0.03 * h_dir) / 12      # kip-ft
    Cm = 0.6 - 0.4 * r12

    def dlt(cm_):
        if Pu >= 0.75 * Pc:
            return math.inf
        return max(cm_ / (1 - Pu / (0.75 * Pc)), 1.0)

    # moment without the minimum (used for the other axis, R6.6.4.5.4)
    Mc_nomin = dlt(Cm) * M2 if M2 > 0 else 0.0
    if M2 < M2min:
        M2, Cm = M2min, 1.0
        out["M2min_used"] = True
    delta = dlt(Cm)
    out.update(delta=delta, Mc=delta * M2, Mc_nomin=Mc_nomin, Cm=Cm, Pc=Pc)
    return out


def _pm_ratio(Pu, Mux, Muy, circ, dx, dy, phiPo, fc, Ag):
    """Demand/capacity for one load point (moments in kip-ft)."""
    if math.isinf(Mux) or math.isinf(Muy):
        return math.inf, "unstable: Pu >= 0.75 Pc"
    if circ:
        Mr = math.hypot(Mux, Muy)
        return cm.ray_ratio(Mr * 12, Pu, dx["phiMn"], dx["phiPn"])[0], "resultant M (circular)"
    if Muy <= 1e-6 * max(Mux, 1.0):
        return cm.ray_ratio(Mux * 12, Pu, dx["phiMn"], dx["phiPn"])[0], "uniaxial x"
    if Mux <= 1e-6 * max(Muy, 1.0):
        return cm.ray_ratio(Muy * 12, Pu, dy["phiMn"], dy["phiPn"])[0], "uniaxial y"
    if Pu >= 0.1 * fc * Ag / 1000:
        _, _, Pnx = cm.ray_ratio(Mux * 12, Pu, dx["phiMn"], dx["phiPn"])
        _, _, Pny = cm.ray_ratio(Muy * 12, Pu, dy["phiMn"], dy["phiPn"])
        if Pnx <= 0 or Pny <= 0:
            return math.inf, "Bresler: no capacity"
        Pni = 1 / (1 / Pnx + 1 / Pny - 1 / phiPo)
        return (Pu / Pni if Pni > 0 else math.inf), "Bresler 1/Pn"
    mxc = cm.M_at_P(Pu, dx["phiMn"], dx["phiPn"])
    myc = cm.M_at_P(Pu, dy["phiMn"], dy["phiPn"])
    ratio = (Mux * 12 / mxc if mxc > 0 else math.inf) + (Muy * 12 / myc if myc > 0 else math.inf)
    return ratio, "linear Mx/Mnx + My/Mny"


def run(shape="Rectangular tied", b=16.0, h=16.0, D=18.0, fc=4000.0, fy=60000.0,
        fyt=60000.0, bar=9, nx=3, ny=3, n_circ=8, cover=1.5, tie=3, s_tie=16.0, agg=0.75,
        lu=12.0, kx=1.0, ky=1.0, ei_method=EI_METHODS[0], case_loads=None):
    """case_loads: {case: (P, Mx_top, Mx_bot, My_top, My_bot)} in kip / kip-ft."""
    case_loads = {k: v for k, v in (case_loads or {}).items() if any(abs(x) > 0 for x in v)}
    bar, tie = int(bar), int(tie)
    spiral = "spiral" in shape
    trans = "spiral" if spiral else "tied"
    circ = shape.startswith("Circular")
    db, Ab = C.bar(bar)
    tie_db = C.bar(tie)[0]
    if circ:
        b = h = D
    pts = bar_layout(shape, b, h, D, cover, tie_db, db, nx, ny, n_circ)
    nbars = len(pts)
    outline_x = cm.CircleShape(D) if circ else cm.rect(b, h)
    outline_y = cm.CircleShape(D) if circ else cm.rect(h, b)
    sec_x = cm.Section(outline_x, [(h / 2 - y, Ab) for x, y in pts], fc, fy, trans)
    sec_y = cm.Section(outline_y, [(b / 2 - x, Ab) for x, y in pts], fc, fy, trans)
    Ag = sec_x.Ag
    Ast = nbars * Ab
    rho = Ast / Ag
    Ec = C.Ec(fc)
    dx = sec_x.interaction()
    dy = sec_y.interaction()
    phi_c = 0.75 if spiral else 0.65
    phiPo = phi_c * sec_x.Po

    r = Report(f"Concrete Column - {shape} " + (f"D={D:g} in" if circ else f"{b:g}x{h:g} in"),
               code="ACI 318-19, ASCE 7-22 (LRFD)")
    r.input("Shape", shape)
    r.input("Size", f"D = {D:g}" if circ else f"b = {b:g} (x), h = {h:g} (y)", "in")
    r.input("f'c / fy", f"{fc:g} / {fy:g}", "psi")
    r.input("Bars", f"{nbars} #{bar}" + ("" if circ else f" ({nx} per x-face, {ny} per y-face)"))
    r.input("Clear cover to tie/spiral", cover, "in")
    r.input("Ties / spiral", f"#{tie} @ {s_tie:g} in")
    r.input("Unsupported length lu", lu, "ft")
    r.input("k (x-axis / y-axis bending)", f"{kx:g} / {ky:g}")
    r.input("Frame", "Nonsway (braced)")
    for c, v in case_loads.items():
        r.input(f"Loads {c}: P, Mxt, Mxb, Myt, Myb", ", ".join(f"{x:g}" for x in v), "kip, kip-ft")

    r.result("Ag", Ag, "in^2")
    r.result("Ast", Ast, "in^2")
    r.result("rho_g", rho)
    r.result("Po", sec_x.Po, "kip", "ACI 22.4.2.2")
    r.result("phi Pn,max", dx["phiPn_max"], "kip", "ACI 22.4.2.1 / Table 21.2.2")
    M0x = cm.M_at_P(0.0, dx["phiMn"], dx["phiPn"])
    M0y = cm.M_at_P(0.0, dy["phiMn"], dy["phiPn"])
    r.result("phi Mnx at P = 0", M0x / 12, "kip-ft")
    r.result("phi Mny at P = 0", M0y / 12, "kip-ft")

    # -- reinforcement limits & detailing --------------------------------------
    r.check("rho_g >= 0.01", 0.01, rho, "", "ACI 10.6.1.1")
    r.check("rho_g <= 0.08", rho, 0.08, "", "ACI 10.6.1.1", note="0.04 practical at splices")
    r.check("Number of bars", 6 if spiral else 4, nbars, "", "ACI 10.7.3.1")
    clear_req = max(1.5, 1.5 * db, 4 / 3 * agg)
    if circ:
        Rb = D / 2 - cover - tie_db - db / 2
        clear = 2 * Rb * math.sin(math.pi / nbars) - db
    else:
        xb = b / 2 - cover - tie_db - db / 2
        yb = h / 2 - cover - tie_db - db / 2
        clear = min(2 * xb / (nx - 1), 2 * yb / (ny - 1)) - db
    r.check("Clear bar spacing", clear_req, clear, "in", "ACI 25.2.3")
    if spiral:
        Dc = (D if circ else min(b, h)) - 2 * cover
        Ach = math.pi * Dc**2 / 4
        Asp = C.bar(tie)[1]
        rho_s = 4 * Asp / (Dc * s_tie)
        rho_s_min = 0.45 * (Ag / Ach - 1) * fc / fyt
        r.check("Spiral ratio rho_s", rho_s_min, rho_s, "", "ACI 25.7.3.3 (Eq. 25.7.3.3)")
        r.check("Spiral bar >= 3/8 in", 0.375, tie_db, "in", "ACI 25.7.3.2")
        clear_sp = s_tie - tie_db
        r.check("Spiral clear pitch <= 3 in", clear_sp, 3.0, "in", "ACI 25.7.3.1")
        r.check("Spiral clear pitch >= 1 in", max(1.0, 4 / 3 * agg), clear_sp, "in", "ACI 25.7.3.1")
    else:
        min_tie = 4 if bar >= 11 else 3
        r.check("Tie size", min_tie, tie, "bar #", "ACI 25.7.2.2")
        least = D if circ else min(b, h)
        smax = min(16 * db, 48 * tie_db, least)
        r.check("Tie spacing", s_tie, smax, "in", "ACI 25.7.2.1")

    # -- load combinations and slenderness --------------------------------------
    rx = 0.25 * D if circ else 0.3 * h
    ry = 0.25 * D if circ else 0.3 * b
    Igx = outline_x.Ig
    Igy = outline_y.Ig
    Isx = sum(Ab * y**2 for x, y in pts)
    Isy = sum(Ab * x**2 for x, y in pts)
    Es = C.ES / 1000  # ksi
    Ec_k = Ec / 1000
    if ei_method.startswith("0.4"):
        EIx, EIy = 0.4 * Ec_k * Igx, 0.4 * Ec_k * Igy
    else:
        EIx, EIy = 0.2 * Ec_k * Igx + Es * Isx, 0.2 * Ec_k * Igy + Es * Isy
    r.result("r (x-axis / y-axis)", f"{rx:.2f} / {ry:.2f}", "in", "ACI 6.2.5.2")
    r.result("EI basis", ei_method, "", "ACI 6.6.4.4.4")

    combos = loads.active(loads.lrfd(), case_loads.keys()) if case_loads else []
    rows, worst = [], None
    for name, f in combos:
        v = np.zeros(5)
        for c, fac in f.items():
            if c in case_loads:
                v += fac * np.array(case_loads[c])
        Pu, Mxt, Mxb, Myt, Myb = (float(x) for x in v)
        PD = f.get("D", 0.0) * case_loads.get("D", (0,) * 5)[0]
        beta = min(max(PD / Pu, 0.0), 1.0) if Pu > 0 else 0.0
        mx = _magnify(Pu, Mxt, Mxb, kx * lu * 12, rx, h, EIx, beta, "x")
        my = _magnify(Pu, Myt, Myb, ky * lu * 12, ry, b, EIy, beta, "y")
        # M2,min need not act about both axes at once (ACI R6.6.4.5.4):
        # check (x with minimum, y without) and (x without, y with minimum).
        pairs = {(mx["Mc"], my["Mc_nomin"]), (mx["Mc_nomin"], my["Mc"])}
        if circ:
            # circular: resultant of the actual moments, but not less than the
            # magnified minimum moment in any one direction
            Mr = max(math.hypot(mx["Mc_nomin"], my["Mc_nomin"]), mx["Mc"], my["Mc"])
            pairs = {(Mr, 0.0)}
        ratio, how, Mux, Muy = -1.0, "", 0.0, 0.0
        for mxx, myy in pairs:
            rt, hw = _pm_ratio(Pu, mxx, myy, circ, dx, dy, phiPo, fc, Ag)
            if rt > ratio:
                ratio, how, Mux, Muy = rt, hw, mxx, myy
        row = {"Combination": name, "Pu (k)": Pu, "Mux (k-ft)": Mux, "Muy (k-ft)": Muy,
               "delta x": mx["delta"], "delta y": my["delta"], "Ratio": ratio, "Method": how}
        rows.append(row)
        if worst is None or ratio > worst[0]["Ratio"]:
            worst = (row, mx, my)

    if worst:
        row, mx, my = worst
        r.result("Governing combination", row["Combination"])
        for ax, m in (("x", mx), ("y", my)):
            txt = f"klu/r = {m['klu/r']:.1f}, limit = {m['limit']:.1f}"
            if m["slender"]:
                txt += (f", slender: Cm = {m['Cm']:.2f}, Pc = {m['Pc']:.0f} k, "
                        f"delta = {m['delta']:.3f}" + (", M2,min used" if m["M2min_used"] else ""))
            else:
                txt += ", slenderness neglected"
            r.result(f"Slenderness about {ax}", txt, "", "ACI 6.2.5.1 / 6.6.4.5")
            r.check(f"klu/r about {ax} <= 100", m["klu/r"], 100, "", "ACI 6.2.6",
                    note="above 100 a nonlinear second-order analysis is required")
        r.check("P-M interaction", row["Ratio"], 1.0, "ratio", "ACI 22.4 / 10.5.1",
                note=f"[{row['Combination']}] {row['Method']}")
        r.section("Load combinations (magnified moments)", rows)
    r.note("Nonsway frame assumed. Sway frames need 6.6.4.6 (not implemented).")
    r.note("Biaxial: Bresler reciprocal load method when Pu >= 0.1 f'c Ag, else the "
           "conservative linear sum Mux/phiMnx + Muy/phiMny.")
    r.note("Design curve is cut at phi Pn,max with a horizontal line (ACI 22.4.2.1).")
    if shape == "Rectangular spiral":
        r.note("Rectangular spiral: circular spiral with core diameter = least side - 2 cover.")
    return {"report": r, "diagram_x": dx, "diagram_y": dy, "rows": rows, "bars": pts,
            "sec_x": sec_x, "sec_y": sec_y, "circular": circ,
            "loads": dict(H=lu, case_loads=case_loads)}


def plot(out):
    import matplotlib.pyplot as plt

    from ..load_plots import draw_column_loads

    ld = out["loads"]
    cl = ld["case_loads"]
    pm_axes = [("x", out["diagram_x"], "Mux (k-ft)")]
    if not out["circular"]:
        pm_axes.append(("y", out["diagram_y"], "Muy (k-ft)"))
    n_tot = 1 + len(pm_axes)
    ratios = [0.9] + [1.2] * len(pm_axes)
    fig, ax = plt.subplots(1, n_tot, figsize=(4.5 * n_tot + 0.5, 5.0),
                           gridspec_kw={"width_ratios": ratios})
    draw_column_loads(
        ax[0],
        H=ld["H"],
        P_cases={c: v[0] for c, v in cl.items()},
        Mx_top={c: v[1] for c, v in cl.items()},
        Mx_bot={c: v[2] for c, v in cl.items()},
        My_top={c: v[3] for c, v in cl.items()},
        My_bot={c: v[4] for c, v in cl.items()},
        title="Applied service loads",
    )
    for a, (nm, d, key) in zip(ax[1:], pm_axes):
        a.plot(d["Mn"] / 12, d["Pn"], "--", color="gray", label="Nominal (Pn,max cap)")
        a.plot(d["phiMn"] / 12, d["phiPn"], color="tab:blue", lw=2, label="Design (phi)")
        pts = [(row[key], row["Pu (k)"]) for row in out["rows"]]
        if pts:
            a.scatter(*zip(*pts), color="tab:red", zorder=3, label="Factored loads")
        a.axhline(0, color="k", lw=0.6)
        a.axvline(0, color="k", lw=0.6)
        a.set_xlabel(f"M about {nm}-axis (kip-ft)")
        a.set_ylabel("P (kip)")
        a.grid(alpha=0.3)
        a.legend(fontsize=8)
        a.set_title(f"P–M interaction, bending about {nm}")
    fig.tight_layout()
    return fig
