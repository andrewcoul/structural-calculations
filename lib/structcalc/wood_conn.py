"""NDS 2024 Chapter 12 dowel-type fasteners (bolts, lag screws, nails, wood screws).

Units: in, lb, psi (connection values in NDS are in lb). Calc modules
convert to kip only for display where useful.

Contents
--------
* Specific gravity table (NDS Table 12.3.3A, common species).
* Dowel bearing strength Fe (NDS 12.3.3, Table 12.3.3B for steel).
* Yield limit equations (NDS Table 12.3.1A) with Rd (Table 12.3.1B).
* Adjustment factors: Cg (Eq. 12.3-1), C_Delta (12.5.1), Ceg, Ctn, Cdi.
* Withdrawal (NDS 12.2): nails, wood screws, lag screws.
* Combined lateral and withdrawal loads (NDS 12.4.1).
* Fastener dimension tables (NDS Appendix L).

A licensed engineer must verify the tables against the current NDS.
"""

from __future__ import annotations

import math

# ---------------------------------------------------------------------------
# Specific gravity, NDS Table 12.3.3A (assigned G for dowel bearing)
# ---------------------------------------------------------------------------
SPECIES_G = {
    "Douglas Fir-Larch": 0.50,
    "Douglas Fir-Larch (North)": 0.49,
    "Douglas Fir-South": 0.46,
    "Hem-Fir": 0.43,
    "Hem-Fir (North)": 0.46,
    "Southern Pine": 0.55,
    "Spruce-Pine-Fir": 0.42,
    "Spruce-Pine-Fir (South)": 0.36,
    "Western Woods": 0.36,
    "Redwood (close grain)": 0.44,
    "Eastern Softwoods": 0.36,
    "Mixed Maple": 0.55,
    "Red Oak": 0.67,
    "Mixed Oak": 0.68,
    "Yellow Poplar": 0.43,
    # Wood structural panels (NDS 12.3.3.3 / Table 12.3.3A footnote)
    "Plywood - Structural I / Marine": 0.50,
    "Plywood - other grades": 0.42,
    "OSB - all grades": 0.50,
}

# Dowel bearing strength of steel members, NDS Table 12.3.3B (psi), = 1.5*Fu
# for hot-rolled plate; light-gauge sheet per the same table. VERIFY.
STEEL_FE = {
    "ASTM A36 (Fu=58 ksi)": 87000.0,
    "ASTM A572 Gr 50 (Fu=65 ksi)": 97500.0,
    "ASTM A992 (Fu=65 ksi)": 97500.0,
    "ASTM A653 Gr 33 sheet (Fu=45 ksi)": 61850.0,
    "ASTM A653 Gr 40 sheet (Fu=55 ksi)": 75600.0,
}

# ---------------------------------------------------------------------------
# Fasteners
# ---------------------------------------------------------------------------
BOLT_DIAMETERS = [0.5, 0.625, 0.75, 0.875, 1.0]
FYB_BOLT = 45000.0  # psi, NDS Table 12A footnote

# NDS Table L4 (nails): pennyweight -> (L in, D in, head dia in)
NAILS = {
    "common": {
        "6d": (2.0, 0.113, 0.266), "8d": (2.5, 0.131, 0.281), "10d": (3.0, 0.148, 0.312),
        "12d": (3.25, 0.148, 0.312), "16d": (3.5, 0.162, 0.344), "20d": (4.0, 0.192, 0.406),
        "30d": (4.5, 0.207, 0.438), "40d": (5.0, 0.225, 0.469), "50d": (5.5, 0.244, 0.500),
        "60d": (6.0, 0.263, 0.531),
    },
    "box": {
        "6d": (2.0, 0.099, 0.266), "8d": (2.5, 0.113, 0.297), "10d": (3.0, 0.128, 0.312),
        "12d": (3.25, 0.128, 0.312), "16d": (3.5, 0.135, 0.344), "20d": (4.0, 0.148, 0.375),
        "30d": (4.5, 0.148, 0.375), "40d": (5.0, 0.162, 0.406),
    },
    "sinker": {
        "6d": (1.875, 0.092, 0.234), "8d": (2.375, 0.113, 0.266), "10d": (2.875, 0.120, 0.281),
        "12d": (3.125, 0.135, 0.312), "16d": (3.25, 0.148, 0.344), "20d": (3.75, 0.177, 0.375),
        "30d": (4.25, 0.192, 0.406), "40d": (4.75, 0.207, 0.438), "60d": (5.75, 0.244, 0.500),
    },
}

# Wood screws, NDS Table L3: gauge -> D (shank) in
WOOD_SCREWS = {"#6": 0.138, "#8": 0.164, "#10": 0.190, "#12": 0.216, "#14": 0.242}


def nail_fyb(D):
    """Nail bending yield strength (psi) by diameter, NDS Table 12N footnote /
    Appendix I. VERIFY for the specific fastener."""
    if D <= 0.142:
        return 100000.0
    if D <= 0.177:
        return 90000.0
    if D <= 0.236:
        return 80000.0
    if D <= 0.273:
        return 70000.0
    if D <= 0.344:
        return 60000.0
    return 45000.0


def screw_fyb(D):
    """Wood screw Fyb (psi), NDS Table 12L footnote. VERIFY."""
    if D <= 0.142:
        return 100000.0
    if D <= 0.177:
        return 90000.0
    if D <= 0.236:
        return 80000.0
    return 70000.0


# ---------------------------------------------------------------------------
# Dowel bearing strength
# ---------------------------------------------------------------------------
def fe_wood(G, D, theta_deg=0.0):
    """Dowel bearing strength (psi) of wood, NDS 12.3.3 + Hankinson (Eq. 12.3-11).

    D < 1/4 in: Fe = 16,600 G^1.84 (no grain-angle effect).
    D >= 1/4 in: Fe_par = 11,200 G; Fe_perp = 6,100 G^1.45 / sqrt(D).
    """
    if D < 0.25:
        return 16600.0 * G**1.84
    fpar = 11200.0 * G
    fperp = 6100.0 * G**1.45 / math.sqrt(D)
    return hankinson(fpar, fperp, theta_deg)


def hankinson(Fpar, Fperp, theta_deg):
    t = math.radians(theta_deg)
    return Fpar * Fperp / (Fpar * math.sin(t) ** 2 + Fperp * math.cos(t) ** 2)


# ---------------------------------------------------------------------------
# Yield limit equations
# ---------------------------------------------------------------------------
def rd_factors(D, theta_max_deg=0.0):
    """Reduction term Rd by mode, NDS Table 12.3.1B."""
    if D <= 0.17:
        kd = 2.2
        return {m: kd for m in ("Im", "Is", "II", "IIIm", "IIIs", "IV")}
    if D < 0.25:
        kd = 10 * D + 0.5
        return {m: kd for m in ("Im", "Is", "II", "IIIm", "IIIs", "IV")}
    kt = 1 + 0.25 * min(max(theta_max_deg, 0.0), 90.0) / 90.0
    return {"Im": 4 * kt, "Is": 4 * kt, "II": 3.6 * kt, "IIIm": 3.2 * kt, "IIIs": 3.2 * kt,
            "IV": 3.2 * kt}


def yield_modes(D, lm, ls, Fem, Fes, Fyb, shear="single", theta_max_deg=0.0):
    """Reference lateral design value Z for each yield mode (lb), NDS Table 12.3.1A.

    D: diameter (in); lm, ls: dowel bearing lengths in main and side member (in);
    Fem, Fes: dowel bearing strengths (psi); Fyb: bending yield strength (psi).
    shear: "single" or "double" (double = symmetric side members of thickness ls).
    Returns dict mode -> Z (lb).
    """
    Rd = rd_factors(D, theta_max_deg)
    Re = Fem / Fes
    Rt = lm / ls
    k3 = -1 + math.sqrt(2 * (1 + Re) / Re + 2 * Fyb * (2 + Re) * D**2 / (3 * Fem * ls**2))
    mode_iv_core = D**2 * math.sqrt(2 * Fem * Fyb / (3 * (1 + Re)))
    if shear == "double":
        return {
            "Im": D * lm * Fem / Rd["Im"],
            "Is": 2 * D * ls * Fes / Rd["Is"],
            "IIIs": 2 * k3 * D * ls * Fem / ((2 + Re) * Rd["IIIs"]),
            "IV": 2 * mode_iv_core / Rd["IV"],
        }
    k1 = (math.sqrt(Re + 2 * Re**2 * (1 + Rt + Rt**2) + Rt**2 * Re**3) - Re * (1 + Rt)) / (1 + Re)
    k2 = -1 + math.sqrt(2 * (1 + Re) + 2 * Fyb * (1 + 2 * Re) * D**2 / (3 * Fem * lm**2))
    return {
        "Im": D * lm * Fem / Rd["Im"],
        "Is": D * ls * Fes / Rd["Is"],
        "II": k1 * D * ls * Fes / Rd["II"],
        "IIIm": k2 * D * lm * Fem / ((1 + 2 * Re) * Rd["IIIm"]),
        "IIIs": k3 * D * ls * Fem / ((2 + Re) * Rd["IIIs"]),
        "IV": mode_iv_core / Rd["IV"],
    }


def governing(modes):
    m = min(modes, key=modes.get)
    return modes[m], m


# ---------------------------------------------------------------------------
# Adjustment factors
# ---------------------------------------------------------------------------
def group_action(n, s, D, Am, As, Em, Es, wood_to_steel=False):
    """Group action factor Cg, NDS Eq. 12.3-1 (D >= 1/4 in).

    n: fasteners in a row; s: spacing in the row (in); Am, As: gross areas
    of main and side members (in^2; As = sum of side members for double shear);
    Em, Es: moduli (psi). gamma = 180,000 D^1.5 (wood-wood) or 270,000 D^1.5
    (wood-steel) lb/in.
    """
    if n <= 1 or D < 0.25:
        return 1.0
    gamma = (270000.0 if wood_to_steel else 180000.0) * D**1.5
    EmAm, EsAs = Em * Am, Es * As
    REA = min(EsAs / EmAm, EmAm / EsAs)
    u = 1 + gamma * s / 2 * (1 / EmAm + 1 / EsAs)
    m = u - math.sqrt(u**2 - 1)
    num = m * (1 - m ** (2 * n))
    den = n * ((1 + REA * m**n) * (1 + m) - 1 + m ** (2 * n))
    return num / den * (1 + REA) / (1 - m)


def geometry_factor(D, load_dir, end_dist, spacing, member="tension", softwood=True, n_in_row=1):
    """Geometry factor C_Delta, NDS 12.5.1 (D >= 1/4 in).

    load_dir: "parallel" or "perpendicular" to grain (of the member checked).
    member: "tension" (load toward the end) or "compression" (load away from end).
    Returns (C_Delta, details dict with minimums and full-value distances).
    """
    if D < 0.25:
        return 1.0, {}
    det = {}
    # end distance, Table 12.5.1A
    if load_dir == "perpendicular":
        e_min, e_full = 2 * D, 4 * D
    elif member == "compression":
        e_min, e_full = 2 * D, 4 * D
    else:
        e_min, e_full = (3.5 * D, 7 * D) if softwood else (2.5 * D, 5 * D)
    det.update(end_min=e_min, end_full=e_full)
    c_end = 1.0 if end_dist >= e_full else (end_dist / e_full if end_dist >= e_min else 0.0)
    # spacing in a row, Table 12.5.1B
    s_min, s_full = 3 * D, 4 * D
    det.update(spacing_min=s_min, spacing_full=s_full)
    if n_in_row > 1:
        c_sp = 1.0 if spacing >= s_full else (spacing / s_full if spacing >= s_min else 0.0)
    else:
        c_sp = 1.0
    det.update(c_end=c_end, c_spacing=c_sp)
    return min(c_end, c_sp), det


def edge_distance_min(D, load_dir, l_over_D, row_spacing=0.0, loaded_edge=True):
    """Minimum edge distance (in), NDS Table 12.5.1C."""
    if load_dir == "parallel":
        if l_over_D <= 6:
            return 1.5 * D
        return max(1.5 * D, row_spacing / 2)
    return 4 * D if loaded_edge else 1.5 * D


def row_spacing_min(D, load_dir, l_over_D):
    """Minimum spacing between rows (in), NDS Table 12.5.1D."""
    if load_dir == "parallel":
        return 1.5 * D
    if l_over_D <= 2:
        return 2.5 * D
    if l_over_D < 6:
        return (5 * l_over_D * D + 10 * D) / 8
    return 5 * D


CEG_LATERAL = 0.67   # NDS 12.5.2 end grain factor, lateral
CTN_LATERAL = 0.83   # NDS 12.5.4 toe-nail factor, lateral
CTN_WITHDRAWAL = 0.67
CDI = 1.1            # NDS 12.5.3 diaphragm factor (nails in diaphragms)

# LRFD conversion (NDS Appendix N)
KF_CONN = 3.32
PHI_Z = 0.65


# ---------------------------------------------------------------------------
# Withdrawal, NDS 12.2 (lb per inch of thread / shank penetration)
# ---------------------------------------------------------------------------
def withdrawal_nail(G, D):
    """Smooth-shank nails and spikes, NDS Eq. 12.2-3."""
    return 1380.0 * G**2.5 * D


def withdrawal_wood_screw(G, D):
    """Wood screws, NDS Eq. 12.2-2 (D = shank diameter)."""
    return 2850.0 * G**2 * D


def withdrawal_lag(G, D):
    """Lag screws, NDS Eq. 12.2-1."""
    return 1800.0 * G**1.5 * D**0.75


def combined_lateral_withdrawal(Zp, Wp, alpha_deg, fastener="nail"):
    """Adjusted capacity at angle alpha between wood surface and load (NDS 12.4.1).

    Zp: adjusted lateral Z'; Wp: adjusted withdrawal W' * p (lb).
    Nails/spikes: Eq. 12.4-2 (cos, sin); others: Eq. 12.4-1 (cos^2, sin^2).
    """
    a = math.radians(alpha_deg)
    if fastener == "nail":
        den = Wp * math.cos(a) + Zp * math.sin(a)
    else:
        den = Wp * math.cos(a) ** 2 + Zp * math.sin(a) ** 2
    return Wp * Zp / den if den > 0 else 0.0
