"""Concrete beam: ACI 318-19 flexure (strain compatibility, doubly reinforced,
rectangular or T), shear with stirrups, crack control, bar fit, and
deflection with Ie (Table 24.2.3.5) plus long-term multiplier (24.2.4.1).
Analysis by structcalc.beam with ASCE 7-22 LRFD combinations.

Units at the interface: dims in in, spans in ft, fc/fy in psi,
loads in kip/ft and kip (+ down).
"""

from __future__ import annotations

import math

import numpy as np

from .. import concrete as C
from .. import concrete_member as cm
from .. import loads
from ..beam import Beam
from ..beaminput import apply_loads, make_supports, parse_loads, regions
from ..report import Report

SECTION_TYPES = ["Rectangular", "T-beam"]
LT_LIMITS = {"L/240 (nonstructural elements not likely damaged)": 240,
             "L/480 (nonstructural elements likely damaged)": 480}


def _layer_spacing(b, cover, ds, n, db):
    """(center-to-center spacing, clear spacing) of n bars in one layer."""
    if n <= 1:
        return math.inf, math.inf
    cc = (b - 2 * (cover + ds) - db) / (n - 1)
    return cc, cc - db


def run(section_type="Rectangular", b=12.0, h=24.0, bf=48.0, hf=4.0,
        fc=4000.0, fy=60000.0, fyt=60000.0, lam=1.0, wc=150.0,
        n_bot=3, bar_bot=8, n_top=2, bar_top=5, cover=1.5, stirrup=3, legs=2,
        s_stirrup=10.0, agg=0.75,
        span_type="Simple span", span=20.0, overhang_left=0.0, overhang_right=0.0,
        custom_supports="", w=None, extra_loads="", self_weight=True,
        deflection_method="Effective Ie (ACI 24.2.3.5)", sustained_live=0.25, xi=2.0,
        live_limit=360, lt_limit=240):
    w = dict(w or {})
    bar_bot, bar_top, stirrup, lt_limit = int(bar_bot), int(bar_top), int(stirrup), int(lt_limit)
    fyt_eff = min(fyt, 60000.0)  # ACI 20.2.2.4 limit for shear reinforcement
    ds = C.bar(stirrup)[0] if (stirrup and s_stirrup > 0) else 0.0
    db_b, Ab = C.bar(bar_bot)
    db_t, At = C.bar(bar_top) if n_top else (0.0, 0.0)
    As_bot, As_top = n_bot * Ab, n_top * At
    d_bot = h - cover - ds - db_b / 2           # depth of bottom bars from top
    d_top = cover + ds + db_t / 2               # depth of top bars from top
    tbeam = section_type == "T-beam"
    shape = cm.tee(b, h, bf, hf) if tbeam else cm.rect(b, h)
    bars = [(d_bot, As_bot)] + ([(d_top, As_top)] if As_top else [])
    sec_pos = cm.Section(shape, bars, fc, fy)
    sec_neg = sec_pos.flipped()                 # compression at bottom
    Ag, Ig, ybar = shape.area, shape.Ig, shape.ybar
    Ec = C.Ec(fc, wc)                           # psi

    # -- analysis ---------------------------------------------------------
    length, supports = make_supports(span_type, span, overhang_left, overhang_right, custom_supports)
    beam = Beam(length, Ec / 1000.0, Ig, supports)
    apply_loads(beam, w, parse_loads(extra_loads))
    if self_weight:
        beam.self_weight = wc * Ag / 144 / 1000  # kip/ft
    per_case = beam.solve()
    results = beam.combos(loads.lrfd())
    Mpos = max(max(r.Mmax for *_, r in results), 0.0)
    Mneg = max(max(-r.Mmin for *_, r in results), 0.0)
    Vu = max(r.Vabs for *_, r in results)
    gp = max(results, key=lambda t: t[2].Mmax)[0]
    gn = max(results, key=lambda t: -t[2].Mmin)[0]
    gv = max(results, key=lambda t: t[2].Vabs)[0]

    r = Report(f"Concrete Beam - {section_type} {b:g}x{h:g} in",
               code="ACI 318-19, ASCE 7-22 (LRFD), IBC 2024")
    r.input("Section", section_type + (f", bf={bf:g} in, hf={hf:g} in" if tbeam else ""))
    r.input("Web width bw", b, "in")
    r.input("Depth h", h, "in")
    r.input("f'c", fc, "psi")
    r.input("fy / fyt", f"{fy:g} / {fyt:g}", "psi")
    r.input("Bottom bars", f"{n_bot} #{bar_bot}")
    r.input("Top bars", f"{n_top} #{bar_top}" if n_top else "none")
    r.input("Clear cover to stirrup", cover, "in")
    r.input("Stirrups", f"#{stirrup} x {legs} legs @ {s_stirrup:g} in" if s_stirrup > 0 else "none")
    r.input("Span type", span_type)
    r.input("Total length", length, "ft")
    r.input("Supports", ", ".join(f"{x:g} ft {t}" for x, t in supports))
    for c, v in w.items():
        if v:
            r.input(f"Uniform load {c}", v * 1000, "plf")
    if extra_loads.strip():
        r.input("Extra loads", extra_loads.strip())
    r.input("Self weight included", bool(self_weight), note=f"wc = {wc:g} pcf")

    r.result("d (bottom bars)", d_bot, "in")
    r.result("d' (top bars)", d_top if As_top else "-", "in")
    r.result("Ag", Ag, "in^2")
    r.result("Ig (gross, steel ignored)", Ig, "in^4")
    r.result("Ec", Ec / 1000, "ksi", "ACI 19.2.2.1")
    r.result("beta1", C.beta1(fc), "", "ACI Table 22.2.2.4.3")

    # -- flexure ------------------------------------------------------------
    flex_pos = sec_pos.flexure()
    r.result("Max +Mu", Mpos, "kip-ft", f"[{gp}]")
    r.result("+M: c / a", f"{flex_pos['c']:.2f} / {flex_pos['a']:.2f}", "in")
    r.result("+M: eps_t / phi", f"{flex_pos['eps_t']:.5f} / {flex_pos['phi']:.3f}", "",
             "ACI Table 21.2.2")
    r.check("Flexure +M", Mpos, flex_pos["phi"] * flex_pos["Mn"] / 12, "kip-ft",
            "ACI 22.2 / 9.5.1.1", note="strain compatibility")
    r.check("Net tensile strain +M (eps_t >= 0.004)", 0.004, flex_pos["eps_t"], "",
            "ACI 9.3.3.1")
    As_req = C.As_required(Mpos * 12, b, d_bot, fc, fy) if Mpos > 0 else 0.0  # web width: conservative
    Asmin = C.As_min_beam(b, d_bot, fc, fy)
    if Mpos > 0:
        need = min(Asmin, 4 / 3 * As_req)
        r.check("As,min bottom", need, As_bot, "in^2", "ACI 9.6.1.2 / 9.6.1.3",
                note=f"As,min={Asmin:.2f}, 4/3 As,req={4 / 3 * As_req:.2f}")
    flex_neg = None
    if Mneg > 1e-6:
        r.result("Max -Mu", Mneg, "kip-ft", f"[{gn}]")
        if As_top > 0:
            flex_neg = sec_neg.flexure()
            r.result("-M: eps_t / phi", f"{flex_neg['eps_t']:.5f} / {flex_neg['phi']:.3f}")
            r.check("Flexure -M", Mneg, flex_neg["phi"] * flex_neg["Mn"] / 12, "kip-ft",
                    "ACI 22.2 / 9.5.1.1")
            r.check("Net tensile strain -M (eps_t >= 0.004)", 0.004, flex_neg["eps_t"], "",
                    "ACI 9.3.3.1")
            dneg = h - d_top
            As_req_n = C.As_required(Mneg * 12, b, dneg, fc, fy)
            Asmin_n = C.As_min_beam(b, dneg, fc, fy)
            r.check("As,min top", min(Asmin_n, 4 / 3 * As_req_n), As_top, "in^2",
                    "ACI 9.6.1.2 / 9.6.1.3")
        else:
            r.check("Flexure -M", Mneg, 0.0, "kip-ft", "ACI 22.2",
                    note="negative moment with no top bars")

    # -- shear ----------------------------------------------------------------
    d = d_bot
    sqfc = min(math.sqrt(fc), 100.0)
    Av = legs * C.bar(stirrup)[1] if s_stirrup > 0 else 0.0
    if s_stirrup > 0:
        Avmin = max(0.75 * sqfc * b * s_stirrup / fyt_eff, 50 * b * s_stirrup / fyt_eff)
        has_min = Av >= Avmin
        Vs = Av * fyt_eff * d / s_stirrup / 1000
    else:
        Avmin, has_min, Vs = 0.0, False, 0.0
    Vc = C.Vc_oneway(b, d, fc, As=As_bot, lam=lam, Av_ge_min=has_min)
    Vs_max = 8 * sqfc * b * d / 1000
    phiVn = C.PHI_SHEAR * (Vc + min(Vs, Vs_max))
    r.result("Max Vu", Vu, "kip", f"[{gv}] at support (conservative, not at d)")
    r.result("Vc", Vc, "kip", "ACI Table 22.5.5.1 " + ("(a)" if has_min else "(c), size effect"))
    r.result("Vs = Av fyt d / s", Vs, "kip", "ACI 22.5.8.5.3")
    r.check("Shear", Vu, phiVn, "kip", "ACI 22.5.1.1, phi=0.75")
    r.check("Shear section limit Vs <= 8 sqrt(f'c) bw d", Vs, Vs_max, "kip", "ACI 22.5.1.2")
    if Vu > C.PHI_SHEAR * lam * sqfc * b * d / 1000:
        r.check("Av,min required (Vu > phi lam sqrt(f'c) bw d)", Avmin if s_stirrup > 0 else 1.0,
                Av if s_stirrup > 0 else 0.0, "in^2", "ACI 9.6.3.1 / 9.6.3.4")
    if s_stirrup > 0:
        smax = min(d / 2, 24.0) if Vs <= 4 * sqfc * b * d / 1000 else min(d / 4, 12.0)
        r.check("Stirrup spacing", s_stirrup, smax, "in", "ACI 9.7.6.2.2")

    # -- detailing --------------------------------------------------------------
    clear_req = max(1.0, db_b, 4 / 3 * agg)
    cc_b, clr_b = _layer_spacing(b, cover, ds, n_bot, db_b)
    if n_bot > 1:
        r.check("Bar fit bottom (clear spacing)", clear_req, clr_b, "in", "ACI 25.2.1",
                note="single layer")
    fs = 2 / 3 * fy
    c_c = cover + ds
    s_crack = min(15 * (40000 / fs) - 2.5 * c_c, 12 * (40000 / fs))
    if n_bot > 1:
        r.check("Crack control spacing bottom", cc_b, s_crack, "in", "ACI 24.3.2",
                note=f"fs = 2/3 fy, cc = {c_c:.2f} in")
    if n_top > 1:
        cc_t, clr_t = _layer_spacing(b, cover, ds, n_top, db_t)
        r.check("Bar fit top (clear spacing)", max(1.0, db_t, 4 / 3 * agg), clr_t, "in", "ACI 25.2.1")
        if Mneg > 1e-6:
            r.check("Crack control spacing top", cc_t, s_crack, "in", "ACI 24.3.2")

    # -- deflection ------------------------------------------------------------
    use_ie = deflection_method.startswith("Effective")
    fr = C.fr(fc, lam)
    Mcr_pos = fr * Ig / (h - ybar) / 1000 / 12       # kip-ft
    Mcr_neg = fr * Ig / ybar / 1000 / 12
    kd_p, Icr_p = sec_pos.cracked(Ec)
    neg_for_icr = sec_neg if As_top > 0 else None
    kd_n, Icr_n = neg_for_icr.cracked(Ec) if neg_for_icr else (0.0, 0.0)
    r.result("Mcr (+M / -M)", f"{Mcr_pos:.1f} / {Mcr_neg:.1f}", "kip-ft", "ACI 24.2.3.5, fr=7.5 lam sqrt(f'c)")
    r.result("Icr (+M)", Icr_p, "in^4", f"kd = {kd_p:.2f} in")
    if neg_for_icr:
        r.result("Icr (-M)", Icr_n, "in^4", f"kd = {kd_n:.2f} in")

    def svc(factors):
        res = None
        for cse, f in factors.items():
            if cse in per_case:
                res = per_case[cse] * f if res is None else res + per_case[cse] * f
        return res

    live = {"L": 1.0, "Lr": 1.0, "S": 1.0}
    rD = svc({"D": 1.0})
    rDL = svc({"D": 1.0, **live})
    rows = []
    for x1, x2, Lref, rname in regions(length, supports):
        def region_vals(res):
            if res is None:
                return 0.0, 0.0
            m = (res.x >= x1 - 1e-9) & (res.x <= x2 + 1e-9)
            Mr = res.M[m]
            # ACI 24.2.3.7: midspan (positive) Ie for spans, support Ie for
            # cantilevers and overhangs
            if rname.startswith("Span") and Mr.max() > 0:
                Ma = Mr.max()
            else:
                Ma = Mr.min() if abs(Mr.min()) > 1e-9 else Mr.max()
            return float(Ma), float(np.abs(res.defl[m]).max())

        MaD, dD = region_vals(rD)
        MaDL, dDL = region_vals(rDL)

        def ie(Ma):
            if not use_ie:
                return Ig
            if Ma >= 0:
                return cm.Ie_aci(Ma, Mcr_pos, Ig, Icr_p)
            if neg_for_icr is None:
                return cm.Ie_aci(Ma, Mcr_neg, Ig, Icr_p)  # no top steel: approx
            return cm.Ie_aci(Ma, Mcr_neg, Ig, Icr_n)

        IeD, IeDL = ie(MaD), ie(MaDL)
        DD = dD * Ig / IeD
        DDL = dDL * Ig / IeDL
        DL = max(DDL - DD, 0.0)
        # compression steel ratio at the region's critical section (24.2.4.1.2)
        if MaDL >= 0:
            rho_p = As_top / (b * d_bot)
        else:
            rho_p = As_bot / (b * (h - d_top))
        lam_d = xi / (1 + 50 * rho_p)
        Dsus = DD + sustained_live * DL
        D_lt = lam_d * Dsus + DL
        rows.append({"Region": rname, "Ma D (k-ft)": MaD, "Ma D+L (k-ft)": MaDL,
                     "Ie D (in4)": IeD, "Ie D+L (in4)": IeDL, "Delta D (in)": DD,
                     "Delta L (in)": DL, "lambda_delta": lam_d, "LT + L (in)": D_lt})
        r.check(f"Deflection immediate live - {rname}", DL, Lref * 12 / live_limit, "in",
                f"ACI Table 24.2.2 (L/{live_limit})",
                note=f"L/{Lref * 12 / DL:.0f}" if DL > 1e-9 else "")
        r.check(f"Deflection long-term + live - {rname}", D_lt, Lref * 12 / lt_limit, "in",
                f"ACI Table 24.2.2 (L/{lt_limit})",
                note=f"lambda_delta={lam_d:.2f}, sustained L {sustained_live:.0%}")
    r.section("Deflection by region (service)", rows)

    reac_rows = []
    for x in sorted(per_case[next(iter(per_case))].reactions):
        row = {"Support x (ft)": x}
        for c, res in per_case.items():
            row[f"R {c} (k)"] = res.reactions[x][0]
        reac_rows.append(row)
    r.section("Service reactions by load case (kip, + up)", reac_rows)
    r.section("Load combinations", [{"Combination": n, "M+ (k-ft)": round(res.Mmax, 4),
                                     "M- (k-ft)": round(res.Mmin, 4), "V (k)": res.Vabs}
                                    for n, _, res in results])
    r.note("Analysis uses gross Ig for all spans (relative stiffness); deflections are "
           "scaled by Ig/Ie for each region and load level (ACI 24.2.3.7 midspan/support Ie).")
    r.note("Live deflection = Delta(D+L) - Delta(D); L includes L, Lr and S.")
    r.note("Shear demand is the maximum shear (at the support face or load), not at d "
           "from the support: conservative. Torsion is not checked.")
    r.note("Each bar group is one layer. Development and splice lengths are not checked.")
    if tbeam:
        r.note("T-beam: flange width bf must satisfy ACI 6.3.2.1. For negative moment the "
               "top bars are assumed within the web width for spacing checks.")
    return {"report": r, "per_case": per_case, "combos": results, "sec_pos": sec_pos,
            "sec_neg": sec_neg, "length": length, "supports": supports,
            "governing": {"pos": gp, "neg": gn, "shear": gv}}
