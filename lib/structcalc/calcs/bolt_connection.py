"""Bolted wood connection: NDS 2024 Chapter 12 (yield limit equations),
group action Cg, geometry factor C_Delta, ASD or LRFD.
"""

from __future__ import annotations

import math

from .. import wood_conn as wc
from ..report import Report

SHEAR = ["Single shear", "Double shear"]
SIDE = ["Wood", "Steel plate"]


def run(D=0.625, shear="Single shear", side_material="Wood",
        species_main="Douglas Fir-Larch", species_side="Douglas Fir-Larch",
        G_main=None, G_side=None, steel_grade="ASTM A36 (Fu=58 ksi)",
        t_main=3.5, w_main=7.25, t_side=1.5, w_side=7.25,
        E_main=1.6e6, E_side=1.6e6,
        theta_main=0.0, theta_side=0.0,
        n_rows=1, n_per_row=2, spacing=2.5, row_spacing=0.0,
        end_main=7.0, edge_main=1.5, stress_main="tension", loaded_edge_main=True,
        end_side=7.0, edge_side=1.5, stress_side="tension", loaded_edge_side=True,
        softwood=True, Fyb=wc.FYB_BOLT,
        CD=1.0, CM=1.0, Ct=1.0, method="ASD", lam=1.0, load=0.0):
    """All lengths in inches, stresses psi, angles degrees, load in lb.

    load: ASD-level load (ASD) or factored load (LRFD) on the whole connection.
    For double shear, t_side and w_side are for ONE side member.
    """
    steel_side = side_material == "Steel plate"
    double = shear == "Double shear"
    Gm = G_main if G_main else wc.SPECIES_G[species_main]
    Gs = None if steel_side else (G_side if G_side else wc.SPECIES_G[species_side])

    Fem = wc.fe_wood(Gm, D, theta_main)
    if steel_side:
        Fes = wc.STEEL_FE[steel_grade]
        theta_max = theta_main
        E_side = 29.0e6
    else:
        Fes = wc.fe_wood(Gs, D, theta_side)
        theta_max = max(theta_main, theta_side)
    lm, ls = t_main, t_side
    modes = wc.yield_modes(D, lm, ls, Fem, Fes, Fyb, "double" if double else "single", theta_max)
    Z, mode = wc.governing(modes)

    # group action (each row has n_per_row bolts)
    Am = t_main * w_main
    As = t_side * w_side * (2 if double else 1)
    Cg = wc.group_action(n_per_row, spacing, D, Am, As, E_main, E_side, wood_to_steel=steel_side)

    # geometry factor and minimums for each wood member (NDS 12.5.1)
    def dirn(theta):
        return "parallel" if theta < 45 else "perpendicular"

    l_over_D = (lm / D) if steel_side else min(lm, ls) / D
    members = [("Main", theta_main, end_main, edge_main, stress_main, loaded_edge_main)]
    if not steel_side:
        members.append(("Side", theta_side, end_side, edge_side, stress_side, loaded_edge_side))
    geo = []
    for nm, th, e_end, e_edge, stress, loaded in members:
        cdel, det = wc.geometry_factor(D, dirn(th), e_end, spacing, stress, softwood, n_per_row)
        geo.append((nm, th, e_end, e_edge, loaded, cdel, det))
    CDelta = min(g[5] for g in geo)

    n_total = n_rows * n_per_row
    if method.upper() == "LRFD":
        f_time = wc.KF_CONN * wc.PHI_Z * lam
        Zp = Z * CM * Ct * Cg * CDelta * f_time
        fnote = f"KF={wc.KF_CONN}, phi={wc.PHI_Z}, lambda={lam}"
    else:
        Zp = Z * CD * CM * Ct * Cg * CDelta
        fnote = f"CD={CD}"
    Zp_total = n_total * Zp

    r = Report(f"Bolted Wood Connection - {D:g} in bolts, {shear.lower()}",
               code=f"NDS 2024 Ch. 12 ({method.upper()})")
    r.input("Bolt diameter D", D, "in")
    r.input("Bolt Fyb", Fyb, "psi")
    r.input("Connection", f"{shear}, {side_material.lower()} side member(s)")
    r.input("Main member", f"{species_main if not G_main else 'G=' + str(G_main)}, "
                           f"{t_main:g} x {w_main:g} in")
    if steel_side:
        r.input("Side plate(s)", f"{steel_grade}, {t_side:g} x {w_side:g} in")
    else:
        r.input("Side member(s)", f"{species_side if not G_side else 'G=' + str(G_side)}, "
                                  f"{t_side:g} x {w_side:g} in")
    r.input("Load angle to grain: main / side", f"{theta_main:g} / {theta_side if not steel_side else '-'}", "deg")
    r.input("Bolts", f"{n_rows} row(s) x {n_per_row} per row = {n_total}")
    r.input("Spacing in row / between rows", f"{spacing:g} / {row_spacing:g}", "in")
    r.input("Main member end / edge distance", f"{end_main:g} / {edge_main:g}", "in")
    if not steel_side:
        r.input("Side member end / edge distance", f"{end_side:g} / {edge_side:g}", "in")
    if method.upper() == "LRFD":
        r.input("lambda / CM / Ct", f"{lam:g} / {CM:g} / {Ct:g}")
    else:
        r.input("CD / CM / Ct", f"{CD:g} / {CM:g} / {Ct:g}")
    r.input("Applied load", load, "lb")

    r.result("G main", Gm, "", "NDS Table 12.3.3A")
    if Gs is not None:
        r.result("G side", Gs, "", "NDS Table 12.3.3A")
    r.result("Fem (at angle)", Fem, "psi", "NDS 12.3.3, Eq. 12.3-11")
    r.result("Fes" + (" (steel)" if steel_side else " (at angle)"), Fes, "psi",
             "NDS Table 12.3.3B" if steel_side else "NDS 12.3.3")
    r.result("Re = Fem/Fes", Fem / Fes)
    r.result("Reference Z (governing mode)", Z, "lb", f"NDS Table 12.3.1A, mode {mode}")
    r.result("Group action factor Cg", Cg, "", "NDS Eq. 12.3-1")
    r.result("Geometry factor C_Delta", CDelta, "", "NDS 12.5.1")
    r.result("Z' per bolt", Zp, "lb", fnote)
    r.result("Z' connection", Zp_total, "lb")
    r.section("Yield limit equations (per bolt, lb)",
              [{"Mode": m, "Z (lb)": v, "Governs": "<-" if m == mode else ""}
               for m, v in modes.items()])

    r.check("Lateral capacity", load, Zp_total, "lb", "NDS 12.3, 11.3")
    # geometry minimums (demand = minimum required, capacity = provided)
    if D >= 0.25:
        rs = row_spacing if n_rows > 1 else 0.0
        for nm, th, e_end, e_edge, loaded, cdel, det in geo:
            d = dirn(th)
            r.result(f"{nm} member: C_Delta / end distance for C_Delta=1.0",
                     f"{cdel:.3f} / {det['end_full']:.3g} in", "", "NDS Table 12.5.1A")
            end_note = f"{d} to grain" + (", loaded edge" if d == "perpendicular" and loaded else "")
            r.check(f"{nm} end distance >= minimum", det["end_min"], e_end, "in",
                    "NDS Table 12.5.1A", note=end_note)
            e_min = wc.edge_distance_min(D, d, l_over_D, rs, loaded)
            r.check(f"{nm} edge distance >= minimum", e_min, e_edge, "in", "NDS Table 12.5.1C")
            if n_per_row > 1:
                r.check(f"{nm} spacing in row >= minimum", det["spacing_min"], spacing, "in",
                        "NDS Table 12.5.1B")
            if n_rows > 1:
                r.check(f"{nm} spacing between rows >= minimum",
                        wc.row_spacing_min(D, d, l_over_D), row_spacing, "in", "NDS Table 12.5.1D")
                if d == "parallel" and row_spacing * (n_rows - 1) > 5.0:
                    r.note(f"{nm} member: outer rows more than 5 in. apart with load parallel to "
                           "grain - use separate splice plates or check shrinkage (NDS 12.5.1.3).")
    else:
        r.note("D < 1/4 in: C_Delta and Cg do not apply.")
    r.note("Loads at 45 deg or more to grain use the perpendicular-to-grain spacing and edge "
           "rules; less than 45 deg uses the parallel rules (simplification of NDS 12.5.1).")
    r.note("Not checked: net section and shear tearout of members (NDS Appendix E, 11.1.2), "
           "steel side plate strength (AISC 360), bolt bearing on steel, and fire design.")
    r.note("Cg uses gross member areas without division among rows (NDS 12.3.6).")
    return {"report": r, "Z": Z, "mode": mode, "modes": modes, "Zp": Zp, "Zp_total": Zp_total,
            "Cg": Cg, "CDelta": CDelta, "Fem": Fem, "Fes": Fes}
