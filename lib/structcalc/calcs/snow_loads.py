"""Snow loads: ASCE 7-22 Chapter 7.

pg is the strength-level ground snow load for the Risk Category (ASCE 7
Hazard Tool). No snow importance factor in ASCE 7-22.

    pf = 0.7 Ce Ct pg                                  (7.3-1)
    ps = Cs pf                                         (7.4-1)
    pm = pg <= pm,max ; pm,max = 25/30/35/40 psf (RC I-IV)   (7.3.3, Table 7.3-4)
    gamma = 0.13 pg + 14 <= 30 pcf                     (7.7-1)
    hd = 1.5 sqrt(pg^0.74 lu^0.70 W2^1.7 / gamma)      (7.6-1)
    Leeward step drift: hd (lu = upper roof), slope 1:4, w = 4 hd
    Windward step drift: 0.75 hd (lu = lower roof), slope 1:8, w = 8 (0.75 hd)
    Truncated leeward drift: h = hc, w = 4 hd^2 / hc <= 8 hc
    Rain-on-snow: +8 psf where pg <= pm,max and slope < W/50 deg (7.10)
    Sliding snow: 0.4 pf W per ft of eave, over 15 ft (7.9)
Units: psf, ft, pcf, degrees.
"""

from __future__ import annotations

import math

import numpy as np

from ..report import Report

# Table 7.3-1 exposure factor Ce: terrain -> (fully, partially, sheltered)
CE = {
    "B": (0.9, 1.0, 1.2),
    "C": (0.9, 1.0, 1.1),
    "D": (0.8, 0.9, 1.0),
    "Above tree line (windswept mountains)": (0.7, 0.8, None),
    "Alaska, no trees within 2 mi": (0.7, 0.8, None),
}
EXPOSURES = ("Fully exposed", "Partially exposed", "Sheltered")

THERMAL = {
    "Heated, unventilated roof (Table 7.3-3)": None,
    "Kept just above freezing / cold ventilated roof": 1.2,
    "Unheated or open air": 1.2,
    "Freezer building": 1.3,
    "Continuously heated greenhouse": 0.85,
}

# Table 7.3-3: Ct for heated structures with unventilated roofs (R-value x pg)
_CT_R = [20, 30, 40, 50]
_CT_PG = [10, 20, 30, 40, 50, 60, 70]
_CT = [
    [1.20, 1.11, 1.05, 1.01, 1.00, 1.00, 1.00],
    [1.20, 1.17, 1.14, 1.13, 1.12, 1.11, 1.10],
    [1.20, 1.19, 1.17, 1.16, 1.16, 1.15, 1.15],
    [1.20, 1.20, 1.19, 1.19, 1.19, 1.18, 1.18],
]
PM_MAX = {"I": 25.0, "II": 30.0, "III": 35.0, "IV": 40.0}


def Ce(terrain, exposure):
    v = CE[terrain][EXPOSURES.index(exposure)]
    if v is None:
        raise ValueError("Sheltered exposure is not defined for this terrain (Table 7.3-1)")
    return v


def Ct_heated(R_value, pg):
    """Table 7.3-3 with bilinear interpolation; R > 50 -> 1.2."""
    if R_value > 50:
        return 1.2
    R = max(R_value, 20.0)
    p = min(max(pg, 10.0), 70.0)
    row = [float(np.interp(p, _CT_PG, r)) for r in _CT]
    return float(np.interp(R, _CT_R, row))


def Ct(condition, R_value, pg):
    v = THERMAL[condition]
    return Ct_heated(R_value, pg) if v is None else v


def Cs_slope(theta, Ct_val, slippery):
    """Fig. 7.4-1: curve (a) Ct <= 1.1, (b) 1.1 < Ct < 1.2, (c) Ct >= 1.2."""
    if Ct_val <= 1.1:
        t0, t1 = 30.0, 70.0          # (a): one line for all surfaces
        curve = "(a)"
    elif Ct_val < 1.2:
        t0, t1 = (10.0, 70.0) if slippery else (37.5, 70.0)
        curve = "(b)"
    else:
        t0, t1 = (15.0, 70.0) if slippery else (45.0, 70.0)
        curve = "(c)"
    if theta <= t0:
        return 1.0, curve
    if theta >= t1:
        return 0.0, curve
    return 1.0 - (theta - t0) / (t1 - t0), curve


def gamma(pg):
    return min(0.13 * pg + 14.0, 30.0)


def hd_eq(pg, lu, W2):
    """Eq. 7.6-1 drift height (ft)."""
    if lu <= 0:
        return 0.0
    return 1.5 * math.sqrt(pg**0.74 * lu**0.70 * W2**1.7 / gamma(pg))


def step_drift(pg, W2, pf_low, upper_len, lower_len, step_h):
    """Roof step drifts (7.7.1). Returns dict with leeward and windward drifts."""
    g = gamma(pg)
    hb = pf_low / g
    hc = step_h - hb
    out = {"gamma": g, "hb": hb, "hc": hc, "required": hc > 0 and hc / hb >= 0.2 if hb > 0 else False}
    for kind, lu, ratio, slope in (("leeward", upper_len, 1.0, 4.0), ("windward", lower_len, 0.75, 8.0)):
        h_full = ratio * hd_eq(pg, lu, W2)
        if hc <= 0:
            h, w, trunc = 0.0, 0.0, False
        elif h_full <= hc:
            h, w, trunc = h_full, slope * h_full, False
        else:
            h, w, trunc = hc, min(slope * h_full**2 / hc, 2 * slope * hc), True
        out[kind] = {"lu": lu, "h_full": h_full, "h": h, "w": w, "pd": g * h, "truncated": trunc,
                     "w_exceeds_roof": w > lower_len}
    return out


def drift_profile(pf_low, d, lower_len, n=200):
    x = np.linspace(0, lower_len, n)
    prof = {}
    for kind in ("leeward", "windward"):
        dd = d[kind]
        sur = np.where(x < dd["w"], dd["pd"] * (1 - x / dd["w"]) if dd["w"] > 0 else 0.0, 0.0)
        prof[kind] = pf_low + sur
    return x, prof


def run(pg=30.0, W2=0.55, risk_category="II", terrain="C", exposure="Partially exposed",
        thermal="Heated, unventilated roof (Table 7.3-3)", R_value=30.0, roof_type="Gable",
        theta=0.0, slippery=False, eave_to_ridge=20.0, prismatic_short=False,
        drift=True, upper_len=100.0, lower_len=40.0, step_h=10.0,
        sliding=False, upper_theta=26.57, upper_eave_to_ridge=20.0, upper_slippery=False):
    ce = Ce(terrain, exposure)
    ct = Ct(thermal, R_value, pg)
    pf = 0.7 * ce * ct * pg
    cs, curve = Cs_slope(theta, ct, slippery)
    ps = cs * pf
    pmax = PM_MAX[risk_category]
    pm_applies = theta < 15.0 and roof_type in ("Flat", "Monoslope", "Gable", "Hip")
    pm = min(pg, pmax) if pm_applies else 0.0
    ros = pg <= pmax and theta < eave_to_ridge / 50.0
    g = gamma(pg)

    r = Report("Snow Loads", code="ASCE 7-22 Ch. 7 (strength-level pg), IBC 2024")
    r.input("Ground snow load pg", pg, "psf", note=f"Risk Category {risk_category}, ASCE 7 Hazard Tool")
    r.input("Winter wind parameter W2", W2)
    r.input("Terrain / exposure", f"{terrain} / {exposure}")
    r.input("Thermal condition", thermal)
    if THERMAL[thermal] is None:
        r.input("Roof R-value", R_value, "h·ft²·°F/Btu")
    r.input("Roof type", roof_type)
    r.input("Roof slope", theta, "deg")
    r.input("Slippery unobstructed surface", bool(slippery))
    r.input("Eave to ridge distance W", eave_to_ridge, "ft")

    r.result("Ce", ce, "", "Table 7.3-1")
    r.result("Ct", ct, "", "Table 7.3-2 / 7.3-3")
    r.result("Flat roof snow load pf = 0.7 Ce Ct pg", pf, "psf", "Eq. 7.3-1")
    r.result(f"Slope factor Cs, Fig. 7.4-1{curve}", cs, "", "7.4")
    r.result("Sloped roof (balanced) snow load ps", ps, "psf", "Eq. 7.4-1")
    if pm_applies:
        r.result("Minimum roof snow load pm (separate case)", pm, "psf", "7.3.3, Table 7.3-4")
    if ros:
        r.result("Rain-on-snow surcharge", 8.0, "psf", "7.10")
        r.result("Balanced load with rain-on-snow", ps + 8.0, "psf", "7.10")
    design_uniform = max(ps + (8.0 if ros else 0.0), pm)
    r.result("Governing uniform roof snow load", design_uniform, "psf")
    r.result("Snow density gamma", g, "pcf", "Eq. 7.7-1")
    r.result("Snow in seismic weight (pf > 45 psf)", 0.15 * pf if pf > 45 else 0.0, "psf", "12.7.2")

    rows = []
    # Unbalanced loads for hip and gable roofs (7.6.1)
    unb = None
    if roof_type in ("Gable", "Hip"):
        if 2.38 <= theta <= 30.2:
            if prismatic_short and eave_to_ridge <= 20.0:
                unb = {"windward": 0.0, "leeward": pg, "surcharge": 0.0, "length": 0.0}
                r.result("Unbalanced: windward / leeward (W <= 20 ft)", f"0 / {pg:.1f}", "psf", "7.6.1")
            else:
                S = 1.0 / math.tan(math.radians(theta))
                hdu = hd_eq(pg, eave_to_ridge, W2)
                sur = hdu * g / math.sqrt(S)
                ln = 8.0 * math.sqrt(S) * hdu / 3.0
                unb = {"windward": 0.3 * ps, "leeward": ps, "surcharge": sur, "length": ln, "hd": hdu}
                r.result("Unbalanced windward load 0.3 ps", 0.3 * ps, "psf", "7.6.1")
                r.result("Unbalanced leeward hd (lu = W)", hdu, "ft", "Eq. 7.6-1")
                r.result("Leeward surcharge hd gamma / sqrt(S)", sur, "psf", "7.6.1")
                r.result("Surcharge length from ridge 8 sqrt(S) hd / 3", ln, "ft", "7.6.1")
        else:
            r.note("Unbalanced gable/hip load not required: slope outside 2.38 to 30.2 deg (7.6.1).")

    d = None
    if drift:
        pf_low = pf  # lower roof assumed flat with the same Ce, Ct
        d = step_drift(pg, W2, pf_low, upper_len, lower_len, step_h)
        r.input("Upper roof length (leeward fetch)", upper_len, "ft")
        r.input("Lower roof length (windward fetch)", lower_len, "ft")
        r.input("Step height (top of lower roof to top of upper roof)", step_h, "ft")
        r.result("Balanced snow depth hb", d["hb"], "ft", "7.7.1")
        r.result("Clear height hc", d["hc"], "ft", "7.7.1")
        if not d["required"]:
            r.note("hc/hb < 0.2: roof step drift load not required (7.7.1).")
        else:
            for kind in ("leeward", "windward"):
                dd = d[kind]
                rows.append({"Drift": kind.capitalize(), "lu (ft)": dd["lu"],
                             "h untruncated (ft)": dd["h_full"], "h design (ft)": dd["h"],
                             "w (ft)": dd["w"], "pd surcharge (psf)": dd["pd"],
                             "max load pf+pd (psf)": pf_low + dd["pd"],
                             "truncated": "yes" if dd["truncated"] else "no"})
                if dd["w_exceeds_roof"]:
                    r.note(f"{kind.capitalize()} drift width exceeds the lower roof length: the drift is "
                           "truncated at the far end of the lower roof (load not reduced to zero).")
            r.section("Roof step drifts (7.7.1), on the balanced lower roof load pf", rows)

    if sliding:
        cs_up, _ = Cs_slope(upper_theta, ct, upper_slippery)
        trigger = (upper_slippery and upper_theta > math.degrees(math.atan(0.25 / 12))) or \
                  (not upper_slippery and upper_theta > math.degrees(math.atan(2 / 12)))
        if trigger:
            Ws = 0.4 * pf * upper_eave_to_ridge
            dist = min(15.0, lower_len)
            r.result("Sliding snow load 0.4 pf W (per ft of eave)", Ws, "plf", "7.9")
            r.result("Distributed over", dist, "ft", "7.9")
            r.result("Sliding surcharge on lower roof", Ws / dist, "psf", "7.9")
        else:
            r.note("Sliding snow not required: upper roof slope below the 7.9 threshold.")

    if pg <= 10:
        r.note("Low-snow site: check the 7.2 exception (pg <= 10 psf with lu <= 100 ft, or pg <= 5 psf "
               "with lu <= 300 ft, or no drift potential and pg below the factored roof live load).")
    r.note("ASCE 7-22 combinations: LRFD 1.0S (principal) / 0.3S (companion), ASD 0.7S.")
    r.note("pm is a separate uniform load case; it is not combined with drift, sliding, unbalanced, "
           "or partial loads, and rain-on-snow is not added to it.")
    r.note("Windward drift truncation uses w = 8 h^2/hc <= 16 hc (same snow area as the 1:8 triangle); "
           "verify against 7.7.1. Leeward drift: hd not limited by lower-roof length (conservative).")
    r.note("Parapet/rooftop-unit drifts (7.8), partial loads (7.5), and ponding (7.11) are not computed.")
    geom = dict(roof_type=roof_type, theta=theta, eave_to_ridge=eave_to_ridge,
                upper_len=upper_len, lower_len=lower_len, step_h=step_h,
                design_uniform=design_uniform, ros=ros)
    return {"report": r, "pf": pf, "ps": ps, "pm": pm, "Ct": ct, "Cs": cs, "drift": d,
            "unbalanced": unb, "gamma": g, "geom": geom}


def plot(out, lower_len=None):
    """Roof snow geometry sketch (left) + lower-roof drift pressure profile (right)."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Rectangle

    g = out.get("geom", {})
    lower_len = lower_len if lower_len is not None else g.get("lower_len", 40.0)
    d = out["drift"]
    unb = out["unbalanced"]

    fig, (axs, axp) = plt.subplots(1, 2, figsize=(11.8, 4.1), gridspec_kw={"width_ratios": [1.2, 1.0]})

    # --- Left: Roof geometry & snow depth/profile sketch ---
    if d and d["required"]:
        step_h = g.get("step_h", 10.0)
        hb, hc = d["hb"], d["hc"]
        u_show = min(g.get("upper_len", 100.0) * 0.25, max(lower_len * 0.45, 12.0))
        # Upper and lower roof structures
        axs.add_patch(Rectangle((-u_show, 0), u_show, step_h, fc="#cbd5e1", ec="#1e293b", lw=1.5, zorder=2))
        axs.plot([0, lower_len], [0, 0], color="#1e293b", lw=2.5, zorder=3)
        # Balanced snow layer hb on both roofs
        axs.add_patch(Rectangle((-u_show, step_h), u_show, hb, fc="#bfdbfe", ec="#3b82f6", lw=1.0, zorder=3))
        axs.add_patch(Rectangle((0, 0), lower_len, hb, fc="#bfdbfe", ec="#3b82f6", lw=1.0,
                                label=f"Balanced hb={hb:.2f} ft (pf={out['pf']:.1f} psf)", zorder=3))
        # Leeward drift triangle
        dl = d["leeward"]
        wl = min(dl["w"], lower_len)
        hl_end = dl["h"] * max(1.0 - wl / dl["w"], 0.0) if dl["w"] > 0 else 0.0
        axs.add_patch(Polygon([(0, hb), (0, hb + dl["h"]), (wl, hb + hl_end), (wl, hb)],
                              closed=True, fc="#60a5fa", alpha=0.55, ec="#1d4ed8", lw=1.5,
                              label=f"Leeward hd={dl['h']:.2f} ft, w={dl['w']:.1f} ft", zorder=4))
        # Windward drift dashed outline
        dw = d["windward"]
        ww = min(dw["w"], lower_len)
        hw_end = dw["h"] * max(1.0 - ww / dw["w"], 0.0) if dw["w"] > 0 else 0.0
        axs.plot([0, ww], [hb + dw["h"], hb + hw_end], color="#ea580c", ls="--", lw=1.8,
                 label=f"Windward hd={dw['h']:.2f} ft, w={dw['w']:.1f} ft", zorder=5)
        # Step height & hc annotations
        axs.annotate("", xy=(-0.08 * u_show, step_h), xytext=(-0.08 * u_show, 0),
                     arrowprops=dict(arrowstyle="<->", color="#334155", lw=1.0))
        axs.text(-0.14 * u_show, step_h / 2, f"step\n{step_h:g} ft", ha="right", va="center", fontsize=8)
        axs.text(lower_len * 0.52, hb + dl["h"] * 0.65,
                 f"hc = {hc:.2f} ft\nmax = {out['pf'] + max(dl['pd'], dw['pd']):.1f} psf",
                 ha="left", va="center", fontsize=8.2, fontweight="semibold", color="#1e3a8a",
                 bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#93c5fd"))
        axs.set_xlim(-u_show - 2.0, lower_len + 2.0)
        axs.set_ylim(-0.6, max(step_h + hb + 1.5, hb + dl["h"] + 1.5))
        axs.set_xlabel("Distance from roof step (ft)")
        axs.set_ylabel("Height above lower roof (ft)")
        axs.set_title("Roof step drift geometry", fontsize=9.6, fontweight="semibold")
        axs.legend(loc="upper right", fontsize=7.4)
        axs.grid(alpha=0.25)
    elif unb is not None:
        W = g.get("eave_to_ridge", 20.0)
        rise = W * math.tan(math.radians(g.get("theta", 15.0)))
        axs.plot([-W, 0, W], [0, rise, 0], color="#1e293b", lw=2.4, zorder=3)
        # Windward & leeward unbalanced load bands above roof
        sc = max(rise, 6.0) * 0.38 / max(unb["leeward"] + unb["surcharge"], out["ps"], 1.0)
        axs.fill([-W, 0, 0, -W], [0, rise, rise + unb["windward"] * sc, unb["windward"] * sc],
                 fc="#93c5fd", alpha=0.5, ec="#2563eb", label=f"Windward: {unb['windward']:.1f} psf")
        axs.fill([0, W, W, 0], [rise, 0, unb["leeward"] * sc, rise + unb["leeward"] * sc],
                 fc="#60a5fa", alpha=0.5, ec="#1d4ed8", label=f"Leeward: {unb['leeward']:.1f} psf")
        if unb["surcharge"] > 0 and unb["length"] > 0:
            x_s = min(unb["length"], W)
            y_s = rise * (1 - x_s / W)
            axs.fill([0, x_s, x_s, 0],
                     [rise + unb["leeward"] * sc, y_s + unb["leeward"] * sc,
                      y_s + (unb["leeward"] + unb["surcharge"]) * sc,
                      rise + (unb["leeward"] + unb["surcharge"]) * sc],
                     fc="#f97316", alpha=0.55, ec="#c2410c",
                     label=f"+Surcharge: {unb['surcharge']:.1f} psf ({unb['length']:.1f} ft)")
        axs.set_xlabel("Distance from ridge (ft)")
        axs.set_ylabel("Roof elevation (ft)")
        axs.set_title(f"Unbalanced {g.get('roof_type', 'gable').lower()} roof snow ({g.get('theta', 0):g}°)",
                      fontsize=9.6, fontweight="semibold")
        axs.legend(fontsize=7.6)
        axs.grid(alpha=0.25)
    else:
        W = g.get("eave_to_ridge", 20.0)
        axs.plot([0, W], [0, 0], color="#1e293b", lw=2.5)
        axs.add_patch(Rectangle((0, 0), W, out["pf"] / out["gamma"], fc="#bfdbfe", ec="#2563eb", lw=1.2))
        axs.text(W / 2, out["pf"] / out["gamma"] + 0.15, f"pf = {out['pf']:.1f} psf", ha="center", va="bottom",
                 fontweight="bold", color="#1d4ed8")
        axs.set_ylim(-0.3, out["pf"] / out["gamma"] * 2.0 + 0.5)
        axs.set_xlabel("Roof width (ft)")
        axs.set_ylabel("Snow depth (ft)")
        axs.set_title("Balanced roof snow load", fontsize=9.6, fontweight="semibold")
        axs.grid(alpha=0.25)

    # --- Right: Lower roof snow load profile (psf) ---
    if d and d["required"]:
        x, prof = drift_profile(out["pf"], d, lower_len)
        axp.plot(x, prof["leeward"], label="Leeward drift + pf", color="tab:blue", lw=1.8)
        axp.plot(x, prof["windward"], label="Windward drift + pf", color="tab:orange", ls="--", lw=1.8)
    axp.axhline(out["pf"], color="k", lw=0.8, label="Balanced pf")
    axp.set_xlabel("Distance from roof step (ft)")
    axp.set_ylabel("Snow load (psf)")
    axp.set_xlim(0, lower_len)
    axp.set_ylim(bottom=0)
    axp.grid(alpha=0.3)
    axp.legend(fontsize=8)
    axp.set_title("Lower roof snow load at step (psf)", fontsize=9.6, fontweight="semibold")
    fig.tight_layout()
    return fig
