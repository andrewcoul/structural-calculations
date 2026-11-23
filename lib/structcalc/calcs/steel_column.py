"""Steel column: AISC 360-22 Ch. E, F, G, H (H1.1) with optional nonsway
B1 amplification (Appendix 8), ASCE 7-22 combinations.

Sign conventions
----------------
* Axial load P: + compression, - tension (kip).
* End moments M_top, M_bot (kip-ft) are INTERNAL bending moments at the
  column ends with one sign convention: equal values = single curvature
  (uniform moment), opposite signs = reverse (double) curvature.
* Eccentricity ex (in) produces strong-axis moment P*ex at the top;
  ey (in) produces weak-axis moment P*ey at the top.
* Uniform lateral loads wx / wy (kip/ft) along the height produce
  strong-axis / weak-axis bending (simply supported between ends).
"""

from __future__ import annotations

import math

import numpy as np

from .. import loads, steel
from ..report import Report

CASES = ("D", "L", "Lr", "S", "R", "W", "E")


def _moment_diagram(H, Mbot, Mtop, w, n=41):
    z = np.linspace(0.0, H, n)
    M = Mbot + (Mtop - Mbot) * z / H + w * z * (H - z) / 2.0
    V = (Mtop - Mbot) / H + w * (H / 2.0 - z)
    return z, M, V


def _cm(Mbot, Mtop, has_transverse):
    """AISC App. 8 Eq. A-8-4 (Cm = 0.6 - 0.4 M1/M2, M1/M2 + for reverse
    curvature). With the internal-moment convention used here this is
    Cm = 0.6 + 0.4 * (Msmall / Mlarge)."""
    if has_transverse:
        return 1.0
    a, b = (Mbot, Mtop) if abs(Mbot) <= abs(Mtop) else (Mtop, Mbot)
    if abs(b) < 1e-12:
        return 1.0
    return 0.6 + 0.4 * (a / b)


def analyze(section, Fy=50.0, method="LRFD", H=12.0, Lcx=12.0, Lcy=12.0, Lcz=12.0, Lb=None,
            case_loads=None, ex=0.0, ey=0.0, wx=None, wy=None, self_weight=True,
            amplify_B1=True, E=steel.E_STEEL):
    """case_loads: {case: dict(P=, Mx_top=, Mx_bot=, My_top=, My_bot=)}.
    wx, wy: {case: kip/ft}. Lengths in ft."""
    p = steel.shape(section)
    method = method.upper()
    case_loads = {c: dict(v) for c, v in (case_loads or {}).items()}
    wx = dict(wx or {})
    wy = dict(wy or {})
    if self_weight:
        d = case_loads.setdefault("D", {})
        d["P"] = d.get("P", 0.0) + p["weight"] / 1000.0 * H
    present = [c for c in CASES if any(abs(v) > 0 for v in case_loads.get(c, {}).values())
               or wx.get(c) or wy.get(c)]
    combos = loads.lrfd() if method == "LRFD" else loads.asd()
    combos = loads.active(combos, present)
    Lb = Lcy if not Lb else Lb

    comp = steel.compression(p, Fy, Lcx * 12, Lcy * 12, Lcz * 12, E=E)
    fc = steel.factor("c", method)
    fb = steel.factor("b", method)
    Pc = fc * comp.Pn
    Pt, _ = steel.tension_yield(p, Fy, method)
    fy_minor = steel.flexure(p, Fy, "y")
    Mcy = fb * fy_minor.Mn / 12.0
    Vcx, _, vnx = steel.shear(p, Fy, "x", method)
    Vcy, _, vny = steel.shear(p, Fy, "y", method)
    alpha = 1.0 if method == "LRFD" else 1.6
    Ix, Iy = p["Ix"], p["Iy"]
    Pe1x = math.pi**2 * E * Ix / (Lcx * 12) ** 2 if Lcx > 0 else math.inf
    Pe1y = math.pi**2 * E * Iy / (Lcy * 12) ** 2 if Lcy > 0 else math.inf

    rows = []
    gov = None
    for name, f in combos:
        def comb(key, f=f):
            return sum(fac * case_loads.get(c, {}).get(key, 0.0) for c, fac in f.items())
        P = comb("P")
        Mxt = comb("Mx_top") + P * ex / 12.0
        Mxb = comb("Mx_bot")
        Myt = comb("My_top") + P * ey / 12.0
        Myb = comb("My_bot")
        wxx = sum(fac * wx.get(c, 0.0) for c, fac in f.items())
        wyy = sum(fac * wy.get(c, 0.0) for c, fac in f.items())
        z, Mx, Vx = _moment_diagram(H, Mxb, Mxt, wxx)
        _, My, Vy = _moment_diagram(H, Myb, Myt, wyy)
        Mx1 = float(np.abs(Mx).max())
        My1 = float(np.abs(My).max())
        B1x = B1y = 1.0
        if amplify_B1 and P > 0 and (Mx1 > 0 or My1 > 0):
            Cmx = _cm(Mxb, Mxt, abs(wxx) > 0)
            Cmy = _cm(Myb, Myt, abs(wyy) > 0)
            B1x = Cmx / (1 - alpha * P / Pe1x) if alpha * P < Pe1x else math.inf
            B1y = Cmy / (1 - alpha * P / Pe1y) if alpha * P < Pe1y else math.inf
            B1x = max(B1x, 1.0) if Mx1 > 0 else 1.0
            B1y = max(B1y, 1.0) if My1 > 0 else 1.0
        Mrx = B1x * Mx1 if Mx1 > 0 else 0.0
        Mry = B1y * My1 if My1 > 0 else 0.0

        # Cb from the moment diagram (quarter points of the unbraced length)
        def at(frac, M=Mx):
            return float(np.interp(frac * H, z, M))
        Cb = steel.Cb_factor(at(0.25), at(0.5), at(0.75), Mx1) if Mx1 > 0 else 1.0
        fx = steel.flexure(p, Fy, "x", Lb * 12, Cb)
        Mcx = fb * fx.Mn / 12.0
        if P >= 0:
            ratio, eq = steel.interaction_H1(P, Pc, Mrx, Mcx, Mry, Mcy)
            Pcap = Pc
        else:
            ratio, eq = steel.interaction_H1(-P, Pt, Mrx, Mcx, Mry, Mcy)
            Pcap = Pt
            eq = eq + " (tension, H1.2)"
        Vux = float(np.abs(Vx).max())
        Vuy = float(np.abs(Vy).max())
        row = {"Combination": name, "P (k)": P, "Mrx (k-ft)": Mrx, "Mry (k-ft)": Mry,
               "B1x": B1x, "B1y": B1y, "Cb": Cb, "H1 ratio": ratio, "Eq.": eq,
               "Vx/Vc": Vux / Vcx, "Vy/Vc": Vuy / Vcy}
        rows.append(row)
        if gov is None or ratio > gov[0]["H1 ratio"]:
            gov = (row, dict(P=P, Pcap=Pcap, Mrx=Mrx, Mcx=Mcx, Mry=Mry, Mcy=Mcy, fx=fx,
                             Vux=Vux, Vuy=Vuy, z=z,
                             Mx=Mx * (B1x if math.isfinite(B1x) else 1.0),
                             My=My * (B1y if math.isfinite(B1y) else 1.0)))

    word = "phi" if method == "LRFD" else "1/Omega"
    r = Report(f"Steel Column - {p['name']}", code=f"AISC 360-22 ({method}), ASCE 7-22")
    r.input("Section", p["name"])
    r.input("Fy", Fy, "ksi")
    r.input("Column height H", H, "ft")
    r.input("Lcx / Lcy / Lcz", f"{Lcx:g} / {Lcy:g} / {Lcz:g}", "ft")
    r.input("Lb (lateral-torsional)", Lb, "ft")
    for c in CASES:
        v = case_loads.get(c)
        if v and any(abs(x) > 0 for x in v.values()):
            r.input(f"Loads {c}", ", ".join(f"{k}={x:g}" for k, x in v.items() if x), "kip, kip-ft")
        if wx.get(c) or wy.get(c):
            r.input(f"Lateral {c} wx/wy", f"{wx.get(c, 0) * 1000:g} / {wy.get(c, 0) * 1000:g}", "plf")
    if ex or ey:
        r.input("Eccentricity ex / ey", f"{ex:g} / {ey:g}", "in")
    r.input("Self weight included", bool(self_weight))
    r.input("Nonsway B1 amplification", bool(amplify_B1))

    r.result("Ag", p["area"], "in^2")
    r.result("KL/r governing", comp.KL_r, "", "AISC E2")
    r.result("Governing buckling axis", comp.governing_axis)
    r.result("Fe", comp.Fe, "ksi", "AISC E3-4 / E4")
    r.result("Fcr", comp.Fcr, "ksi", "AISC E3-2/E3-3")
    r.result("Ae", comp.Ae, "in^2", "AISC E7")
    r.result("Element slenderness", comp.detail)
    r.result(f"{word}*Pn (compression)", Pc, "kip", "AISC E3/E7")
    r.result(f"{word}*Pn (tension yield)", Pt, "kip", "AISC D2")
    r.result(f"{word}*Mny", Mcy, "kip-ft", "AISC F6/F7/F8")
    row, g = gov
    r.result("Governing combination", row["Combination"])
    r.result(f"{word}*Mnx (Cb={row['Cb']:.2f})", g["Mcx"], "kip-ft", "AISC F2-F8")
    r.result("Strong-axis limit state", g["fx"].limit_state)
    r.result("B1x / B1y", f"{row['B1x']:.3f} / {row['B1y']:.3f}", "", "AISC App. 8 A-8-3")
    r.check("Axial", abs(g["P"]), g["Pcap"], "kip", "AISC E3 / D2", note=f"[{row['Combination']}]")
    r.check("Combined axial + flexure", row["H1 ratio"], 1.0, "ratio", f"AISC {row['Eq.']}",
            note=f"[{row['Combination']}]")
    vx_g = max(rows, key=lambda d: d["Vx/Vc"])
    vy_g = max(rows, key=lambda d: d["Vy/Vc"])
    r.check("Shear (strong axis)", vx_g["Vx/Vc"] * Vcx, Vcx, "kip", "AISC G", note=vnx)
    r.check("Shear (weak axis)", vy_g["Vy/Vc"] * Vcy, Vcy, "kip", "AISC G", note=vny)
    r.section("Load combinations", rows)
    r.note("Effective length method: Lc = K*L entered directly. B1 uses Pe1 with EI (no 0.8 tau_b "
           "reduction); for the direct analysis method enter K = 1 and reduce stiffness separately.")
    r.note("B2 (sway) amplification is not computed. Enter second-order moments if the frame sways.")
    r.note("Lateral loads wx/wy act on a column pinned at both ends between the end moments.")
    return {"report": r, "section": p, "rows": rows, "governing": gov, "compression": comp,
            "loads": dict(H=H, case_loads=case_loads, ex=ex, ey=ey, wx=wx, wy=wy)}


def run(section="W10X33", auto_select=False, depth_filter="", **kw):
    if not auto_select:
        return analyze(section, **kw)
    names = steel.lightest("W", lambda p: True, nominal_depth=depth_filter or None)
    tried = []
    for n in names:
        out = analyze(n, **kw)
        tried.append((n, out["report"].max_ratio))
        if out["report"].ok:
            out["report"].note(f"Auto-select: lightest passing W shape is {n} ({len(tried)} checked).")
            out["tried"] = tried
            return out
    out = analyze(names[-1], **kw)
    out["report"].note("Auto-select: no W shape passed. Heaviest shape shown.")
    out["tried"] = tried
    return out


def plot(out):
    import matplotlib.pyplot as plt

    from ..load_plots import draw_column_loads

    row, g = out["governing"]
    ld = out["loads"]
    cl = ld["case_loads"]
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 4.4), sharey=True,
                           gridspec_kw={"width_ratios": [1.15, 1.0, 1.0]})
    draw_column_loads(
        ax[0],
        H=ld["H"],
        P_cases={c: v.get("P", 0.0) for c, v in cl.items()},
        Mx_top={c: v.get("Mx_top", 0.0) for c, v in cl.items()},
        Mx_bot={c: v.get("Mx_bot", 0.0) for c, v in cl.items()},
        My_top={c: v.get("My_top", 0.0) for c, v in cl.items()},
        My_bot={c: v.get("My_bot", 0.0) for c, v in cl.items()},
        w_lat_x=ld["wx"],
        w_lat_y=ld["wy"],
        ex=ld["ex"],
        ey=ld["ey"],
        title="Applied service loads",
    )
    ax[1].plot(g["Mx"], g["z"], color="tab:red")
    ax[1].fill_betweenx(g["z"], g["Mx"], color="tab:red", alpha=0.25)
    ax[1].set_title("Mx × B1 (kip-ft)")
    ax[2].plot(g["My"], g["z"], color="tab:blue")
    ax[2].fill_betweenx(g["z"], g["My"], color="tab:blue", alpha=0.25)
    ax[2].set_title("My × B1 (kip-ft)")
    for a in ax[1:]:
        a.axvline(0, color="k", lw=0.6)
        a.grid(alpha=0.3)
    fig.suptitle(f"{out['section']['name']} — {row['Combination']}")
    fig.tight_layout()
    return fig
