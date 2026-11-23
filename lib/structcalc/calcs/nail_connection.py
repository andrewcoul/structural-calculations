"""Nailed / wood-screwed connection: NDS 2024 Chapter 12, single shear.

Lateral (yield limit equations), withdrawal, combined lateral + withdrawal
(NDS 12.4.1), toe-nail, end grain and diaphragm factors, and the number of
fasteners required.
"""

from __future__ import annotations

import math

from .. import wood_conn as wc
from ..report import Report

FASTENERS = ["Common nail", "Box nail", "Sinker nail", "Wood screw"]
INSTALL = ["Side grain", "End grain", "Toe-nail"]


def fastener_props(fastener, size, screw_length=2.5):
    """Return (L, D, head, Fyb, kind)."""
    if fastener == "Wood screw":
        if size not in wc.WOOD_SCREWS:
            raise ValueError(f"wood screw size must be one of {list(wc.WOOD_SCREWS)}, got {size}")
        D = wc.WOOD_SCREWS[size]
        return screw_length, D, None, wc.screw_fyb(D), "screw"
    table = wc.NAILS[fastener.split()[0].lower()]
    if size not in table:
        raise ValueError(f"{fastener} size {size} not in NDS Table L4: {list(table)}")
    L, D, H = table[size]
    return L, D, H, wc.nail_fyb(D), "nail"


def run(fastener="Common nail", size="16d", screw_length=2.5, D_override=0.0, Fyb_override=0.0,
        side_member="Douglas Fir-Larch", main_member="Douglas Fir-Larch",
        G_side=None, G_main=None, t_side=1.5, t_main=3.5, install="Side grain",
        diaphragm=False, CD=1.0, CM=1.0, Ct=1.0, method="ASD", lam=1.0,
        n_fasteners=4, lateral_load=0.0, withdrawal_load=0.0):
    """Lengths in inches, loads in lb on the whole group (ASD level for ASD,
    factored for LRFD). Returns dict with report and key values."""
    L, D, H, Fyb, kind = fastener_props(fastener, size, screw_length)
    if D_override:
        D = D_override
    if Fyb_override:
        Fyb = Fyb_override
    Gs = G_side if G_side else wc.SPECIES_G[side_member]
    Gm = G_main if G_main else wc.SPECIES_G[main_member]

    toe = install == "Toe-nail"
    end_grain = install == "End grain"
    if toe:
        # NDS 12.1.6.4 / Figure 12B: started L/3 from member end at 30 deg
        ls = L / 3.0
        p = L * math.cos(math.radians(30)) - L / 3.0
    else:
        ls = t_side
        p = min(L - t_side, t_main)
    Fes = wc.fe_wood(Gs, D)
    Fem = wc.fe_wood(Gm, D)
    lm = max(p, 1e-6)
    modes = wc.yield_modes(D, lm, ls, Fem, Fes, Fyb, "single", 0.0)
    Z, mode = wc.governing(modes)

    Ceg = wc.CEG_LATERAL if end_grain else 1.0
    Ctn = wc.CTN_LATERAL if toe else 1.0
    Cdi = wc.CDI if (diaphragm and kind == "nail") else 1.0
    lrfd = method.upper() == "LRFD"
    time = wc.KF_CONN * wc.PHI_Z * lam if lrfd else CD
    Zp = Z * time * CM * Ct * Ceg * Cdi * Ctn

    # withdrawal (per inch of penetration into the main member)
    if kind == "nail":
        W = wc.withdrawal_nail(Gm, D)
        p_w = p
    else:
        W = wc.withdrawal_wood_screw(Gm, D)
        p_w = min(p, 2.0 / 3.0 * L)  # threaded length about 2/3 L (NDS Table L3)
    if end_grain:
        W = 0.0  # NDS 12.2.2 / 12.2.3: no withdrawal from end grain for nails/screws
    Ctn_w = wc.CTN_WITHDRAWAL if toe else 1.0
    Wp = W * time * CM * Ct * Ctn_w * p_w  # lb per fastener

    # combined loading
    load_res = math.hypot(lateral_load, withdrawal_load)
    alpha = math.degrees(math.atan2(withdrawal_load, lateral_load)) if load_res else 0.0
    if withdrawal_load > 0 and lateral_load > 0:
        Za = wc.combined_lateral_withdrawal(Zp, Wp, alpha, "nail" if kind == "nail" else "screw")
    elif withdrawal_load > 0:
        Za = Wp
    else:
        Za = Zp
    cap = n_fasteners * Za
    n_req = math.ceil(load_res / Za) if Za > 0 and load_res > 0 else 0

    title = f"{'Nailed' if kind == 'nail' else 'Wood Screw'} Connection - {size} {fastener.lower()}"
    r = Report(title, code=f"NDS 2024 Ch. 12 ({method.upper()})")
    r.input("Fastener", f"{size} {fastener.lower()}: L={L:g} in, D={D:g} in"
            + (f", head {H:g} in" if H else ""))
    r.input("Fyb", Fyb, "psi", note="NDS Table 12N/12L footnote (verify for product)")
    r.input("Side member", f"{side_member if not G_side else 'G=' + str(G_side)}, t={t_side:g} in"
            if not toe else f"{side_member} (toe-nailed)")
    r.input("Main member", f"{main_member if not G_main else 'G=' + str(G_main)}, t={t_main:g} in")
    r.input("Installation", install + (", diaphragm" if diaphragm else ""))
    if lrfd:
        r.input("lambda / CM / Ct", f"{lam:g} / {CM:g} / {Ct:g}")
    else:
        r.input("CD / CM / Ct", f"{CD:g} / {CM:g} / {Ct:g}")
    r.input("Number of fasteners", n_fasteners)
    r.input("Lateral load on group", lateral_load, "lb")
    r.input("Withdrawal load on group", withdrawal_load, "lb")

    r.result("Gs / Gm", f"{Gs:g} / {Gm:g}", "", "NDS Table 12.3.3A")
    r.result("Fe side = Fe main (D < 1/4 in)" if Gs == Gm else "Fes / Fem",
             f"{Fes:.0f}" if Gs == Gm else f"{Fes:.0f} / {Fem:.0f}", "psi", "NDS 12.3.3: 16,600 G^1.84")
    r.result("Side member bearing length ls", ls, "in", "NDS Fig. 12B" if toe else "")
    r.result("Penetration into main member p", p, "in")
    r.result("p / D", p / D, "", "NDS 12.1.6: >= 6 required, 10 for full value")
    r.result("Reference Z (governing mode)", Z, "lb", f"NDS Table 12.3.1A, mode {mode}")
    r.result("Ceg / Ctn / Cdi", f"{Ceg:g} / {Ctn:g} / {Cdi:g}", "", "NDS 12.5.2-12.5.4")
    r.result("Z' per fastener", Zp, "lb", f"KF*phi*lambda={time:.3f}" if lrfd else f"CD={CD:g}")
    r.result("Reference W per inch", W, "lb/in",
             "NDS Eq. 12.2-3" if kind == "nail" else "NDS Eq. 12.2-2")
    r.result("W' x p per fastener", Wp, "lb", f"p = {p_w:.3f} in")
    if withdrawal_load > 0 and lateral_load > 0:
        r.result("Load angle alpha to wood surface", alpha, "deg")
        r.result("Z'alpha per fastener", Za, "lb",
                 "NDS Eq. 12.4-2" if kind == "nail" else "NDS Eq. 12.4-1")
    r.result("Fasteners required", n_req)
    r.section("Yield limit equations (per fastener, lb)",
              [{"Mode": m, "Z (lb)": v, "Governs": "<-" if m == mode else ""}
               for m, v in modes.items()])

    r.check("Connection capacity", load_res, cap, "lb", "NDS 12.3 / 12.2 / 12.4",
            note=f"{n_fasteners} x {Za:.0f} lb")
    r.check("Penetration p >= 6D", 6 * D, p, "in", "NDS 12.1.6 (lateral)")
    if withdrawal_load > 0 and end_grain:
        r.check("Withdrawal from end grain not permitted", 1, 0, "", "NDS 12.2.2/12.2.3")
    if toe:
        r.note("Toe-nail geometry: nail started L/3 from the member end at 30 deg to the face "
               "(NDS 12.1.6, Fig. 12B); ls = L/3, p = L cos30 - L/3.")
    if p / D < 10:
        r.note("p < 10D: the yield equations use lm = p, which reduces the lateral value. Many "
               "NDS tabulated values assume p >= 10D.")
    if kind == "screw":
        r.note("Wood screw: shank diameter used for D (unthreaded shank at the shear plane). "
               "Enter the root diameter Dr in D override if threads cross the shear plane. "
               "Threaded length taken as 2/3 L for withdrawal.")
    r.note("Spacing, end and edge distances for nails and screws: use spacing that prevents "
           "splitting (NDS 12.1.6.2); no C_Delta applies for D < 1/4 in.")
    r.note("Cg = 1.0 for D < 1/4 in (NDS 12.3.6). Fyb values per NDS footnotes; verify the "
           "product.")
    return {"report": r, "Z": Z, "mode": mode, "modes": modes, "Zp": Zp, "Wp": Wp,
            "Za": Za, "n_required": n_req, "p": p, "ls": ls, "D": D, "L": L}
