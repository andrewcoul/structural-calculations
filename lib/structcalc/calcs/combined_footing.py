"""Combined rectangular footing with two columns along its length.

Soil (ASCE 7-22 2.4 ASD): linear pressure along the length with partial
contact (no tension), resultant location.
Concrete (ASCE 7-22 2.3 LRFD, ACI 318-19): shear and moment diagrams along
the length from net factored soil pressure, longitudinal bottom/top steel,
one-way shear (22.5.5.1, size effect), punching at each column (22.6,
8.4.4.2) with the perimeter clipped at the footing ends, transverse strip
steel under each column, bearing (22.8), minimum steel (7.6.1.1).

x is measured from the left end of the footing (ft). Column moments M are
+ clockwise (they move the soil resultant toward +x).
"""

from __future__ import annotations

import math

import numpy as np

from .. import concrete as aci
from .. import foundation as fd
from .. import loads
from ..report import Report

CASE_ORDER = ("D", "L", "Lr", "S", "R", "W", "E")


def _trapz(y, x):
    fn = getattr(np, "trapezoid", None) or np.trapz
    return float(fn(y, x))


def parse_pm(s):
    """'P, M' (kip, kip-ft) -> (P, M)."""
    vals = [float(v) for v in str(s or "").split(",") if v.strip()]
    if len(vals) > 2:
        raise ValueError(f"load '{s}': give P, M")
    return tuple(vals + [0.0] * (2 - len(vals)))


def _diagram(L, B, x_cols, c_cols, P_cols, M_cols, w_down, R_up, e_up, n=3000):
    """Return x, q_up(x) k/ft, V(x) kip, M(x) kip-ft (bottom tension +)."""
    p = fd.pressure_1d(R_up, e_up, L)
    if p is None:
        return None
    ql, qr, cont = p
    pts = [0.0, L] + [xc - cc / 2 for xc, cc in zip(x_cols, c_cols)] + \
          [xc + cc / 2 for xc, cc in zip(x_cols, c_cols)] + list(x_cols)
    if cont < L:
        pts.append(L - cont if qr > 0 else cont)
    x = np.unique(np.concatenate([np.linspace(0, L, n), np.clip(pts, 0, L)]))
    q = fd.pressure_1d_at(x, ql, qr, cont, L)
    # V: upward pressure (trapezoid exact for piecewise linear q with the kink in the grid)
    V = np.concatenate([[0.0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(x))])
    V -= w_down * x
    for xc, cc, P in zip(x_cols, c_cols, P_cols):
        V -= P * np.clip((x - (xc - cc / 2)) / cc, 0.0, 1.0)
    # M: integrate V (Simpson-like accuracy via fine grid) + concentrated moments
    M = np.concatenate([[0.0], np.cumsum(0.5 * (V[1:] + V[:-1]) * np.diff(x))])
    for xc, Mc in zip(x_cols, M_cols):
        M += np.where(x >= xc, Mc, 0.0)
    return x, q, V, M, (ql, qr, cont)


def run(L=20.0, B=6.0, h=30.0, x1=2.0, c1x=16.0, c1y=16.0, x2=16.0, c2x=18.0, c2y=18.0,
        soil_cover=1.0, gamma_c=150.0, gamma_s=110.0, sur_D=0.0, sur_L=0.0,
        qa=3.0, qa_increase_WE=1.0, fc=3000.0, fy=60000.0,
        bar_bot=7, s_bot=8.0, bar_top=6, s_top=10.0, bar_tr=6, s_tr=8.0,
        cover=3.0, top_cover=2.0, strip_ext=0.75, col1=None, col2=None, lam=1.0):
    """Units: L, B, x1, x2, soil_cover ft; h, column dims, spacings, covers in;
    surcharges psf; unit weights pcf; qa ksf; fc, fy psi.
    col1/col2: {case: (P kip, M kip-ft)} service loads.
    strip_ext: transverse strip extends strip_ext*d beyond each column face."""
    col1 = {c: tuple(v) for c, v in (col1 or {}).items() if any(abs(t) > 0 for t in v)}
    col2 = {c: tuple(v) for c, v in (col2 or {}).items() if any(abs(t) > 0 for t in v)}
    present = [c for c in CASE_ORDER if c in col1 or c in col2 or (c == "L" and sur_L)]
    A = L * B
    hf = h / 12
    W_self = gamma_c * hf * A / 1000 + gamma_s * soil_cover * A / 1000   # kip (D)
    x_cols = (x1, x2)
    c_cols = (c1x / 12, c2x / 12)
    for xc, cc, nm in ((x1, c1x / 12, "1"), (x2, c2x / 12, "2")):
        if xc - cc / 2 < -1e-9 or xc + cc / 2 > L + 1e-9:
            raise ValueError(f"column {nm} extends beyond the footing")

    r = Report("Combined Footing", code="ACI 318-19, ASCE 7-22 (ASD soil / LRFD concrete), IBC 2024")
    r.input("Footing L x B x h", f"{L:g} ft x {B:g} ft x {h:g} in")
    r.input("Column 1", f"x = {x1:g} ft, {c1x:g} x {c1y:g} in")
    r.input("Column 2", f"x = {x2:g} ft, {c2x:g} x {c2y:g} in")
    r.input("Soil cover / unit weights", f"{soil_cover:g} ft / {gamma_c:g}, {gamma_s:g} pcf")
    r.input("Surcharge D / L", f"{sur_D:g} / {sur_L:g}", "psf")
    r.input("Allowable soil pressure", qa, "ksf")
    r.input("f'c / fy", f"{fc:g} / {fy:g}", "psi")
    r.input("Longitudinal bottom / top", f"#{bar_bot} @ {s_bot:g} / #{bar_top} @ {s_top:g} in")
    r.input("Transverse (bottom, lowest layer)", f"#{bar_tr} @ {s_tr:g} in")
    for c in present:
        r.input(f"Load {c}: col 1 (P, M) / col 2 (P, M)",
                f"{col1.get(c, (0, 0))[0]:g}, {col1.get(c, (0, 0))[1]:g} / "
                f"{col2.get(c, (0, 0))[0]:g}, {col2.get(c, (0, 0))[1]:g}", "kip, kip-ft")

    def combo_loads(f):
        P1 = sum(fac * col1.get(c, (0, 0))[0] for c, fac in f.items())
        M1 = sum(fac * col1.get(c, (0, 0))[1] for c, fac in f.items())
        P2 = sum(fac * col2.get(c, (0, 0))[0] for c, fac in f.items())
        M2 = sum(fac * col2.get(c, (0, 0))[1] for c, fac in f.items())
        w = (f.get("D", 0) * (W_self / L + sur_D * B / 1000) + f.get("L", 0) * sur_L * B / 1000)
        R = P1 + P2 + w * L
        Mc = P1 * (x1 - L / 2) + P2 * (x2 - L / 2) + M1 + M2
        return P1, M1, P2, M2, w, R, Mc

    # ------------------------------------------------------------- ASD soil
    asd_rows, gov = [], None
    for name, f in loads.active(loads.asd(), present):
        P1, M1, P2, M2, w, R, Mc = combo_loads(f)
        wind_eq = any(f.get(k, 0) != 0 for k in ("W", "E"))
        qall = qa * (qa_increase_WE if wind_eq else 1.0)
        if R <= 0:
            asd_rows.append({"Combination": name, "R (k)": R, "e (ft)": float("nan"),
                             "q left (ksf)": float("nan"), "q right (ksf)": float("nan"),
                             "Contact (ft)": 0.0})
            gov = (math.inf, name, math.inf, qall) if gov is None or gov[0] < math.inf else gov
            continue
        e = Mc / R
        p = fd.pressure_1d(R, e, L)
        if p is None:
            asd_rows.append({"Combination": name, "R (k)": R, "e (ft)": e, "q left (ksf)": float("inf"),
                             "q right (ksf)": float("inf"), "Contact (ft)": 0.0})
            gov = (math.inf, name, math.inf, qall)
            continue
        ql, qr, cont = p
        qmax = max(ql, qr) / B
        asd_rows.append({"Combination": name, "R (k)": R, "e (ft)": e, "q left (ksf)": ql / B,
                         "q right (ksf)": qr / B, "Contact (ft)": cont})
        if gov is None or qmax / qall > gov[0]:
            gov = (qmax / qall, name, qmax, qall, (ql / B, qr / B, cont))
    if gov:
        r.result("Governing bearing combination", gov[1])
        r.check("Soil bearing", gov[2], gov[3], "ksf", "IBC 1806 / ASCE 7-22 2.4", note=f"[{gov[1]}]")
        r.result("Footing + soil weight", W_self, "kip")

    # ------------------------------------------------------------ LRFD concrete
    dbb, Abb = aci.bar(bar_bot)
    dbt, Abt = aci.bar(bar_top)
    dbr, Abr = aci.bar(bar_tr)
    d_tr = h - cover - dbr / 2
    d_bot = h - cover - dbr - dbb / 2
    d_top = h - top_cover - dbt / 2
    d_avg = (d_tr + d_bot) / 2
    b_in = B * 12
    As_bot = Abb * b_in / s_bot
    As_top = Abt * b_in / s_top
    r.result("d: longitudinal bottom / top / transverse", f"{d_bot:.2f} / {d_top:.2f} / {d_tr:.2f}", "in")

    G = {"M+": (0.0, ""), "M-": (0.0, ""), "V": (0.0, "", 0.0, 0.0), "P1": (0.0, ""), "P2": (0.0, ""),
         "v1": (0.0, "", None), "v2": (0.0, "", None)}
    diags = []
    lrfd_rows = []
    fl_b = aci.flexure_rect(As_bot, b_in, d_bot, fc, fy)
    fl_t = aci.flexure_rect(As_top, b_in, d_top, fc, fy)
    phiVc_b = aci.PHI_SHEAR * aci.Vc_oneway(b_in, d_bot, fc, As=As_bot, lam=lam)
    phiVc_t = aci.PHI_SHEAR * aci.Vc_oneway(b_in, d_top, fc, As=As_top, lam=lam)
    for name, f in loads.active(loads.lrfd(), present):
        P1, M1, P2, M2, w, R, Mc = combo_loads(f)
        if R <= 0:
            continue
        dg = _diagram(L, B, x_cols, c_cols, (P1, P2), (M1, M2), w, R, Mc / R)
        if dg is None:
            continue
        x, q, V, M, _ = dg
        diags.append((name, x, V, M))
        inside = np.zeros_like(x, dtype=bool)
        for xc, cc in zip(x_cols, c_cols):
            inside |= (x > xc - cc / 2 + 1e-9) & (x < xc + cc / 2 - 1e-9)
        Mout = M[~inside]
        Mpos, Mneg = float(max(Mout.max(), 0)), float(max(-Mout.min(), 0))
        # one-way shear outside the zones within d of the column faces
        crit = np.ones_like(x, dtype=bool)
        for xc, cc in zip(x_cols, c_cols):
            crit &= ~((x > xc - cc / 2 - d_bot / 12) & (x < xc + cc / 2 + d_bot / 12))
        Vmax_ratio, Vmax, Vcap = 0.0, 0.0, 1.0
        if crit.any():
            for xi, Vi, Mi in zip(x[crit], V[crit], M[crit]):
                cap = phiVc_b if Mi >= 0 else phiVc_t
                if abs(Vi) / cap > Vmax_ratio:
                    Vmax_ratio, Vmax, Vcap = abs(Vi) / cap, abs(float(Vi)), cap
        # punching at each column
        vus = []
        for i, (xc, cxi, cyi, Pi, Mi) in enumerate(((x1, c1x, c1y, P1, M1), (x2, c2x, c2y, P2, M2)), 1):
            ps = fd.punching_section(xc * 12, 0.0, cxi, cyi, d_avg, (0.0, L * 12), (-B * 6, B * 6))
            xa, xb = ps["x1"] / 12, ps["x2"] / 12
            m = (x >= xa) & (x <= xb)
            up = _trapz(q[m], x[m]) if m.sum() > 1 else 0.0
            frac = (ps["y2"] - ps["y1"]) / b_in
            Vu = Pi - (up - w * (xb - xa)) * frac
            Munb = abs(Mi) * 12 + abs(Vu) * abs(ps["xbar"] - xc * 12)
            gv = fd.gamma_v(ps["b1"], ps["b2"])
            vu = Vu * 1000 / (ps["bo"] * d_avg) + gv * Munb * 1000 * ps["c"] / ps["Jc"]
            vus.append(vu)
            if vu > G[f"v{i}"][0]:
                G[f"v{i}"] = (vu, name, ps)
            if Pi > G[f"P{i}"][0]:
                G[f"P{i}"] = (Pi, name)
        lrfd_rows.append({"Combination": name, "P1 (k)": P1, "P2 (k)": P2, "M+ (k-ft)": Mpos,
                          "M- (k-ft)": -Mneg, "V crit (k)": Vmax, "vu col 1 (psi)": vus[0],
                          "vu col 2 (psi)": vus[1], "Closure M(L)": float(M[-1])})
        if Mpos > G["M+"][0]:
            G["M+"] = (Mpos, name)
        if Mneg > G["M-"][0]:
            G["M-"] = (Mneg, name)
        if Vmax_ratio > G["V"][0]:
            G["V"] = (Vmax_ratio, name, Vmax, Vcap)

    if lrfd_rows:
        Mu, cn = G["M+"]
        As_req = aci.As_required(Mu * 12, b_in, d_bot, fc, fy)
        As_min = aci.As_min_slab(b_in, h, fy)
        r.result("Max positive moment (bottom tension)", Mu, "kip-ft", note=f"[{cn}]")
        r.result("Bottom As required / min / provided", f"{As_req:.2f} / {As_min:.2f} / {As_bot:.2f}", "in²")
        r.check("Flexure, bottom longitudinal", Mu, fl_b["phiMn"] / 12, "kip-ft", "ACI 22.2", note=f"[{cn}]")
        r.check("Minimum steel, bottom longitudinal", As_min, As_bot, "in²", "ACI 7.6.1.1 / 24.4.3.2")
        Mu, cn = G["M-"]
        if Mu > 0:
            As_req = aci.As_required(Mu * 12, b_in, d_top, fc, fy)
            r.result("Max negative moment (top tension)", -Mu, "kip-ft", note=f"[{cn}]")
            r.result("Top As required / provided", f"{As_req:.2f} / {As_top:.2f}", "in²")
            r.check("Flexure, top longitudinal", Mu, fl_t["phiMn"] / 12, "kip-ft", "ACI 22.2", note=f"[{cn}]")
        for nm_, s_ in (("bottom", s_bot), ("top", s_top), ("transverse", s_tr)):
            r.check(f"Bar spacing ({nm_})", s_, min(3 * h, 18.0), "in", "ACI 7.7.2.3")
        ratio, cn, Vu, cap = G["V"]
        r.check("One-way shear (longitudinal)", Vu, cap, "kip", "ACI 22.5.5.1 (Av < Av,min, size effect)",
                note=f"[{cn}] beyond d from column faces")
        for i, (cxi, cyi) in enumerate(((c1x, c1y), (c2x, c2y)), 1):
            vu, cn, ps = G[f"v{i}"]
            if ps is None:
                continue
            beta_c = max(cxi, cyi) / min(cxi, cyi)
            vc = aci.vc_twoway(d_avg, ps["bo"], beta_c, ps["alpha_s"], fc, lam)
            r.result(f"Column {i} punching: bo, sides, alpha_s",
                     f"{ps['bo']:.1f} in, {ps['n_sides']}, {ps['alpha_s']}", "", "ACI 22.6.4 / 22.6.5.3")
            r.check(f"Two-way shear, column {i}", vu, aci.PHI_SHEAR * vc, "psi", "ACI 22.6.5.2 / 8.4.4.2",
                    note=f"[{cn}]")
        # transverse strips
        for i, (xc, cxi, cyi) in enumerate(((x1, c1x, c1y), (x2, c2x, c2y)), 1):
            Pu, cn = G[f"P{i}"]
            if Pu <= 0:
                continue
            xa = max(0.0, xc - cxi / 24 - strip_ext * d_tr / 12)
            xb = min(L, xc + cxi / 24 + strip_ext * d_tr / 12)
            bw = (xb - xa) * 12
            proj = (B - cyi / 12) / 2
            Mu_t = Pu / B * proj**2 / 2
            As_t = Abr * bw / s_tr
            fl = aci.flexure_rect(As_t, bw, d_tr, fc, fy)
            As_req = aci.As_required(Mu_t * 12, bw, d_tr, fc, fy)
            As_min = aci.As_min_slab(bw, h, fy)
            r.result(f"Transverse strip {i}: width, Mu", f"{bw:.1f} in, {Mu_t:.1f} kip-ft", "",
                     note=f"[{cn}] strip = column + {strip_ext:g}d each side")
            r.result(f"Transverse strip {i}: As req / min / prov", f"{As_req:.2f} / {As_min:.2f} / {As_t:.2f}",
                     "in²")
            r.check(f"Transverse flexure, strip {i}", Mu_t, fl["phiMn"] / 12, "kip-ft", "ACI 22.2",
                    note=f"[{cn}]")
            r.check(f"Transverse minimum steel, strip {i}", As_min, As_t, "in²", "ACI 7.6.1.1")
            Vu_t = Pu / B * max(proj - d_tr / 12, 0)
            phiVc = aci.PHI_SHEAR * aci.Vc_oneway(bw, d_tr, fc, As=As_t, lam=lam)
            r.check(f"Transverse one-way shear, strip {i}", Vu_t, phiVc, "kip", "ACI 22.5.5.1",
                    note="strip carries the full column load")
            A1 = cxi * cyi
            avail_x = max(2 * min(xc, L - xc) * 12, cxi)  # concentric frustum within the footing
            k = fd.bearing_A2_factor(cxi, cyi, h, avail_x, B * 12)
            r.check(f"Bearing, column {i}", Pu, aci.PHI_BEARING * 0.85 * fc * A1 * k / 1000, "kip",
                    "ACI 22.8.3.2", note=f"sqrt(A2/A1)={k:.2f}")
        r.section("ASD combinations (soil)", asd_rows)
        r.section("LRFD combinations (concrete)", lrfd_rows)
    else:
        r.section("ASD combinations (soil)", asd_rows)
        r.note("No LRFD combination keeps the footing in compression: concrete checks not done.")
    r.note("Rigid footing: linear soil pressure along the length, no tension, uniform across the width.")
    r.note("Column loads are spread over the column width; column moments act at the column centerline.")
    r.note("Maximum moments exclude the zone inside each column (face moments govern).")
    r.note("Punching: unbalanced moment = column moment + Vu x (perimeter centroid offset); "
           "perimeter clipped at footing ends (3-sided for edge columns).")
    r.note("Transverse steel: each column load is carried by a strip under the column "
           "(column width + strip_ext x d each side), net pressure uniform across B.")
    r.note("Horizontal column shears are not included. Development lengths and hooks are not checked.")
    gov_soil = (gov[1], *gov[4]) if (gov and len(gov) > 4) else None
    return {"report": r, "diagrams": diags, "L": L, "B": B, "h": h, "x_cols": x_cols, "c_cols": c_cols,
            "col1": col1, "col2": col2, "gov_soil": gov_soil}


def plot(out):
    """Footing elevation with service loads & governing bearing pressure above LRFD V and M envelopes."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    from ..load_plots import _case_summary, _draw_moment_arc

    L, B = out["L"], out["B"]
    hf = out.get("h", 30.0) / 12.0
    fig, ax = plt.subplots(3, 1, figsize=(9.5, 7.8), sharex=True,
                           gridspec_kw={"height_ratios": [1.35, 1.0, 1.0]})

    # --- Top panel: Footing FBD + service loads + governing bearing pressure ---
    a0 = ax[0]
    a0.add_patch(Rectangle((0, 0), L, hf, fc="#cbd5e1", ec="#1e293b", lw=1.5, zorder=3))
    col_h = 0.65 * hf
    for idx, (xc, cc, cld) in enumerate(zip(out["x_cols"], out["c_cols"], (out["col1"], out["col2"])), 1):
        a0.add_patch(Rectangle((xc - cc / 2, hf), cc, col_h, fc="#94a3b8", ec="#1e293b", lw=1.3, zorder=3))
        y_top = hf + col_h
        P_txt = _case_summary({c: v[0] for c, v in cld.items()}, "k")
        M_txt = _case_summary({c: v[1] for c, v in cld.items()}, "k-ft")
        if P_txt:
            a0.annotate("", xy=(xc, y_top + 0.04), xytext=(xc, y_top + 0.85 * hf),
                        arrowprops=dict(arrowstyle="-|>,head_length=0.5,head_width=0.3", color="#dc2626", lw=2.2),
                        zorder=7)
            a0.text(xc, y_top + 0.92 * hf, f"Col {idx}: {P_txt}", ha="center", va="bottom",
                    fontsize=8.2, fontweight="bold", color="#dc2626",
                    bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="#dc2626", lw=0.7), zorder=8)
        if M_txt:
            _draw_moment_arc(a0, xc, y_top + 0.1 * hf, 0.032 * L, 0.34 * hf,
                             clockwise=(sum(v[1] for v in cld.values()) >= 0), color="#b45309")
            side_x = xc + (0.042 * L if xc < L / 2 else -0.042 * L)
            ha = "left" if xc < L / 2 else "right"
            a0.text(side_x, y_top + 0.25 * hf, f"M={M_txt}", ha=ha, va="center",
                    fontsize=7.8, fontweight="semibold", color="#b45309")

    # Governing service soil pressure underneath
    gs = out.get("gov_soil")
    y_bot = -0.9 * hf
    if gs is not None:
        gname, ql, qr, cont = gs
        xs = np.linspace(0, L, 200)
        q_arr = fd.pressure_1d_at(xs, ql, qr, cont, L)
        qmax = max(ql, qr, 1e-9)
        sc = 0.72 * hf / qmax
        a0.fill_between(xs, -q_arr * sc - 0.12 * hf, -0.12 * hf, color="#f59e0b", alpha=0.35, zorder=2)
        a0.plot(xs, -q_arr * sc - 0.12 * hf, color="#b45309", lw=1.4, zorder=3)
        for xa in np.linspace(0, L, 11):
            qv = float(fd.pressure_1d_at(np.array([xa]), ql, qr, cont, L)[0])
            if qv > 1e-6:
                a0.annotate("", xy=(xa, -0.12 * hf), xytext=(xa, -qv * sc - 0.12 * hf),
                            arrowprops=dict(arrowstyle="-|>,head_length=0.25,head_width=0.15", color="#b45309", lw=0.9))
        a0.text(0.2, -ql * sc - 0.18 * hf, f"{ql:.2f} ksf", ha="left", va="top", fontsize=7.8, color="#92400e")
        a0.text(L - 0.2, -qr * sc - 0.18 * hf, f"{qr:.2f} ksf", ha="right", va="top", fontsize=7.8, color="#92400e")
        a0.text(L / 2, -0.96 * hf, f"Service bearing ({gname})", ha="center", va="top",
                fontsize=7.8, color="#78350f", style="italic")
        y_bot = -1.25 * hf

    a0.set_xlim(-0.05 * L, 1.05 * L)
    a0.set_ylim(y_bot, hf + col_h + 1.35 * hf)
    a0.set_yticks([])
    a0.set_ylabel("Loads")
    a0.spines["top"].set_visible(False)
    a0.spines["right"].set_visible(False)
    a0.spines["left"].set_visible(False)
    a0.grid(axis="x", alpha=0.25)
    a0.set_title(f"Combined footing ({L:g} × {B:g} ft, h = {out.get('h', 30):g}\") — loads & LRFD envelopes",
                 fontsize=10, fontweight="semibold")

    for name, x, V, M in out["diagrams"]:
        ax[1].plot(x, V, lw=0.9, alpha=0.65)
        ax[2].plot(x, M, lw=0.9, alpha=0.65, label=name)
    for a in ax[1:]:
        a.axhline(0, color="k", lw=0.6)
        for xc, cc in zip(out["x_cols"], out["c_cols"]):
            a.axvspan(xc - cc / 2, xc + cc / 2, color="0.6", alpha=0.3)
        a.grid(alpha=0.3)
    ax[1].set_ylabel("Vu (kip)")
    ax[2].set_ylabel("Mu (kip-ft)\n+ = bottom tension")
    ax[2].set_xlabel("x from left end (ft)")
    if len(out["diagrams"]) <= 12:
        ax[2].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    return fig
