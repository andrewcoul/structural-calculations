"""Seismic loads: ASCE 7-22 Equivalent Lateral Force procedure (12.8),
ELF Method 1 (two-period spectrum: SDS, SD1, TL).

Equations:
    Cs = SDS / (R/Ie)                                   (12.8-2)
    Cs <= SD1 / (T (R/Ie))           for T <= TL         (12.8-3)
    Cs <= SD1 TL / (T^2 (R/Ie))      for T >  TL         (12.8-4)
    Cs >= max(0.044 SDS Ie, 0.01)                        (12.8-5)
    Cs >= 0.5 S1 / (R/Ie)            where S1 >= 0.6 g   (12.8-6)
    V = Cs W                                             (12.8-1)
    Ta = Ct hn^x (Table 12.8-2), T <= Cu Ta (Table 12.8-1)
    Fx = Cvx V, Cvx = wx hx^k / sum(wi hi^k)              (12.8-11, 12.8-12)
    Ev = 0.2 SDS D                                       (12.4-4a)
Optional 12.8.1.3: SDS for Cs and Ev may be taken as max(1.0, 0.7 SDS)
when all conditions of 12.8.1.3 are met (user confirms).
"""

from __future__ import annotations

import numpy as np

from ..report import Report

CT_X = {
    "Steel moment-resisting frame": (0.028, 0.8),
    "Concrete moment-resisting frame": (0.016, 0.9),
    "Steel EBF / BRBF": (0.03, 0.75),
    "All other structural systems": (0.02, 0.75),
}
IE = {"I": 1.0, "II": 1.0, "III": 1.25, "IV": 1.5}


def Cu(SD1):
    """Table 12.8-1 (linear interpolation permitted)."""
    return float(np.interp(SD1, [0.1, 0.15, 0.2, 0.3, 0.4], [1.7, 1.6, 1.5, 1.4, 1.4]))


def _sdc_val(v, limits, rc):
    a, b, c = limits
    if v < a:
        return "A"
    if v < b:
        return "C" if rc == "IV" else "B"
    if v < c:
        return "D" if rc == "IV" else "C"
    return "D"


def sdc(SDS, SD1, S1, rc):
    """Tables 11.6-1 and 11.6-2, plus 11.6 (S1 >= 0.75)."""
    if S1 >= 0.75:
        return "F" if rc == "IV" else "E", "S1 >= 0.75g"
    a = _sdc_val(SDS, (0.167, 0.33, 0.50), rc)
    b = _sdc_val(SD1, (0.067, 0.133, 0.20), rc)
    return max(a, b), f"Table 11.6-1: {a}, Table 11.6-2: {b}"


def k_exp(T):
    return float(np.interp(T, [0.5, 2.5], [1.0, 2.0]))


def _floats(text):
    return [float(x) for x in str(text).replace(";", ",").split(",") if x.strip()]


def run(SDS=1.0, SD1=0.6, S1=0.6, TL=8.0, risk_category="II", R=6.5, Omega0=3.0, Cd=4.0,
        rho=1.0, system="All other structural systems", T_computed=0.0,
        story_heights="10; 10", story_weights="400; 300", use_12813=False, site_class="D"):
    """story_heights: story heights (ft) from the base up; story_weights: seismic weight
    at each level (kip), same order (level 1 = first floor above base)."""
    hs = _floats(story_heights)
    ws = _floats(story_weights)
    if len(hs) != len(ws) or not hs:
        raise ValueError("story_heights and story_weights must have the same number of values")
    Ie = IE[risk_category]
    hx = np.cumsum(hs)
    hn = float(hx[-1])
    W = float(sum(ws))
    Ct, x = CT_X[system]
    Ta = Ct * hn**x
    cu = Cu(SD1)
    T = min(T_computed, cu * Ta) if T_computed > 0 else Ta
    sds_cs = max(1.0, 0.7 * SDS) if use_12813 else SDS
    sds_cs = min(sds_cs, SDS)
    RI = R / Ie
    Cs_2 = sds_cs / RI
    Cs_max = SD1 / (T * RI) if T <= TL else SD1 * TL / (T**2 * RI)
    Cs_min = max(0.044 * sds_cs * Ie, 0.01)
    Cs_min6 = 0.5 * S1 / RI if S1 >= 0.6 else 0.0
    Cs = max(min(Cs_2, Cs_max), Cs_min, Cs_min6)
    if Cs == Cs_2:
        gov = "12.8-2"
    elif Cs == Cs_max:
        gov = "12.8-3" if T <= TL else "12.8-4"
    elif Cs == Cs_min6:
        gov = "12.8-6"
    else:
        gov = "12.8-5"
    V = Cs * W
    k = k_exp(T)
    whk = np.array(ws) * hx**k
    Cvx = whk / whk.sum()
    F = Cvx * V
    Vx = np.cumsum(F[::-1])[::-1]
    Mot = [float(sum(F[j] * (hx[j] - (hx[i - 1] if i else 0.0)) for j in range(i, len(F))))
           for i in range(len(F))]
    cat, why = sdc(SDS, SD1, S1, risk_category)

    r = Report("Seismic Loads - Equivalent Lateral Force", code="ASCE 7-22 Ch. 11, 12 (ELF Method 1)")
    r.input("SDS", SDS, "g")
    r.input("SD1", SD1, "g")
    r.input("S1", S1, "g")
    r.input("TL", TL, "s")
    r.input("Site class", site_class)
    r.input("Risk category", risk_category)
    r.input("Seismic force-resisting system (period)", system)
    r.input("R / Omega0 / Cd", f"{R:g} / {Omega0:g} / {Cd:g}", "", note="Table 12.2-1")
    r.input("rho", rho, "", note="12.3.4")
    r.input("Story heights", story_heights, "ft")
    r.input("Level weights", story_weights, "kip")
    if T_computed:
        r.input("Computed period T", T_computed, "s")
    r.input("Use 12.8.1.3 SDS cap", bool(use_12813))

    r.result("Ie", Ie, "", "Table 1.5-2")
    r.result("Seismic design category", cat, "", why)
    r.result("hn", hn, "ft")
    r.result("Ct, x", f"{Ct:g}, {x:g}", "", "Table 12.8-2")
    r.result("Ta = Ct hn^x", Ta, "s", "Eq. 12.8-8")
    r.result("Cu", cu, "", "Table 12.8-1")
    r.result("Cu Ta", cu * Ta, "s", "12.8.2")
    r.result("T used", T, "s", "12.8.2")
    r.result("Ts = SD1/SDS", SD1 / SDS if SDS else 0.0, "s")
    if use_12813:
        r.result("SDS used for Cs and Ev", sds_cs, "g", "12.8.1.3")
    r.result("Cs = SDS/(R/Ie)", Cs_2, "", "Eq. 12.8-2")
    r.result("Cs maximum", Cs_max, "", "Eq. 12.8-3/12.8-4")
    r.result("Cs minimum", max(Cs_min, Cs_min6), "", "Eq. 12.8-5/12.8-6")
    r.result("Cs (governs: Eq. " + gov + ")", Cs, "", "12.8.1.1")
    r.result("Effective seismic weight W", W, "kip", "12.7.2")
    r.result("Base shear V = Cs W", V, "kip", "Eq. 12.8-1")
    r.result("Base shear (ASD, 0.7V)", 0.7 * V, "kip", "2.4.5")
    r.result("Eh = rho QE at base", rho * V, "kip", "Eq. 12.4-3")
    r.result("Ev = 0.2 SDS D", f"{0.2 * sds_cs:.3f} D", "", "Eq. 12.4-4a")
    r.result("Overstrength base shear Omega0 V", Omega0 * V, "kip", "12.4.3")
    r.result("Distribution exponent k", k, "", "12.8.3")
    rows = []
    for i in range(len(hs) - 1, -1, -1):
        rows.append({"Level": i + 1, "hx (ft)": float(hx[i]), "wx (k)": ws[i],
                     "wx hx^k": float(whk[i]), "Cvx": float(Cvx[i]), "Fx (k)": float(F[i]),
                     "Story shear Vx (k)": float(Vx[i]), "Overturning below level (k-ft)": Mot[i]})
    r.section("Vertical distribution of seismic forces (strength level, E = rho QE)", rows)
    r.note("Forces are strength level QE. Multiply by rho for Eh; ASD combinations use 0.7E.")
    r.note("Effective seismic weight per 12.7.2 includes 15% of the flat roof snow load where "
           "pf > 45 psf (ASCE 7-22). Include it in the level weights.")
    r.note("ELF Method 2 (multi-period spectrum, 12.8.1.1) is not computed. Verify ELF is permitted "
           "(Table 12.6-1) and check the system height limits of Table 12.2-1.")
    r.note("Story drift (12.8.6, Cd) and P-delta (12.8.7) require a structural model.")
    return {"report": r, "V": V, "Cs": Cs, "T": T, "k": k, "F": F, "hx": hx, "Vx": Vx, "Mot": Mot, "SDC": cat}


def plot(out):
    import matplotlib.pyplot as plt

    hx, F, Vx = out["hx"], out["F"], out["Vx"]
    Mot = out.get("Mot") or [float(sum(F[j] * (hx[j] - (hx[i - 1] if i else 0.0))
                                       for j in range(i, len(F)))) for i in range(len(F))]
    fig, ax = plt.subplots(1, 3, figsize=(10.8, 4.2), sharey=True)
    ax[0].barh(hx, F, height=max(hx) * 0.055, color="tab:orange")
    for h_i, f_i in zip(hx, F):
        ax[0].text(f_i + 0.03 * max(F), h_i, f"{f_i:.1f}", va="center", ha="left", fontsize=8)
    ax[0].set_xlim(0, max(F) * 1.24)
    ax[0].set_xlabel("Story force Fx (kip)")
    ax[0].set_ylabel("Height (ft)")

    lv = np.concatenate([[0.0], hx])
    xs, ys = [], []
    for i, v in enumerate(Vx):
        xs += [v, v]
        ys += [lv[i], lv[i + 1]]
    ax[1].fill_betweenx(ys, 0, xs, color="tab:blue", alpha=0.22)
    ax[1].plot(xs, ys, color="tab:blue", lw=1.8)
    ax[1].set_xlim(left=0)
    ax[1].set_xlabel("Story shear Vx (kip)")

    m_vals = list(Mot) + [0.0]
    ax[2].fill_betweenx(lv, 0, m_vals, color="tab:red", alpha=0.22)
    ax[2].plot(m_vals, lv, color="tab:red", lw=1.8)
    ax[2].set_xlim(left=0)
    ax[2].set_xlabel("Overturning Mot (kip-ft)")

    for a in ax:
        a.grid(alpha=0.3)
    fig.suptitle(f"ELF (SDC {out['SDC']}): V = {out['V']:.1f} kip, Cs = {out['Cs']:.4f}, T = {out['T']:.3f} s")
    fig.tight_layout()
    return fig
