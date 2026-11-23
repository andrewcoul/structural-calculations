"""NDS 2024 wood member design helpers (sawn lumber, glulam, SCL).

UNITS: stresses and E in psi, dimensions in inches, forces in lb unless a
function says otherwise.

Reference design values
-----------------------
* Table 4A (visually graded dimension lumber, all species except SP):
  DF-L, Hem-Fir, SPF. Values are base values; the size factor CF applies.
* Table 4B (Southern Pine dimension lumber): size-specific values per
  SPIB 2013 revision (still current in NDS 2024 Supplement); CF applies
  only to 4x lumber 8 in. and wider (Fb x1.1) and widths > 12 in. (x0.9).
* Table 4D (timbers 5x5 and larger): DF-L Beams & Stringers and Posts &
  Timbers. Other species: enter values with the override fields.
* Table 5A (glulam, DF): 24F-V4, 24F-V8 bending combinations.
* LVL / PSL: GENERIC typical values only - verify with the manufacturer's
  evaluation report (ICC-ES / APA) before use.

Every value can be overridden by the user (``override`` dicts).
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

# --------------------------------------------------------------------------
# Reference design values
# Keys: Fb, Ft, Fv, Fc_perp, Fc, E, Emin (psi), G (specific gravity)
# --------------------------------------------------------------------------


def _v(Fb, Ft, Fv, Fcp, Fc, E, Emin, G):
    return {"Fb": Fb, "Ft": Ft, "Fv": Fv, "Fc_perp": Fcp, "Fc": Fc, "E": E, "Emin": Emin, "G": G}


# NDS Supplement 2024 Table 4A (2"-4" thick, 2" & wider). Base values; use CF.
TABLE_4A = {
    "DF-L": {
        "Select Structural": _v(1500, 1000, 180, 625, 1700, 1_900_000, 690_000, 0.50),
        "No.1": _v(1000, 675, 180, 625, 1500, 1_700_000, 620_000, 0.50),
        "No.2": _v(900, 575, 180, 625, 1350, 1_600_000, 580_000, 0.50),
        "No.3": _v(525, 325, 180, 625, 775, 1_400_000, 510_000, 0.50),
        "Stud": _v(700, 450, 180, 625, 850, 1_400_000, 510_000, 0.50),
    },
    "Hem-Fir": {
        "Select Structural": _v(1400, 925, 150, 405, 1500, 1_600_000, 580_000, 0.43),
        "No.1": _v(975, 625, 150, 405, 1350, 1_500_000, 550_000, 0.43),
        "No.2": _v(850, 525, 150, 405, 1300, 1_300_000, 470_000, 0.43),
        "No.3": _v(500, 300, 150, 405, 725, 1_200_000, 440_000, 0.43),
        "Stud": _v(675, 400, 150, 405, 800, 1_200_000, 440_000, 0.43),
    },
    "SPF": {
        "Select Structural": _v(1250, 700, 135, 425, 1400, 1_500_000, 550_000, 0.42),
        "No.1": _v(875, 450, 135, 425, 1150, 1_400_000, 510_000, 0.42),   # No.1/No.2
        "No.2": _v(875, 450, 135, 425, 1150, 1_400_000, 510_000, 0.42),   # No.1/No.2
        "No.3": _v(500, 250, 135, 425, 650, 1_200_000, 440_000, 0.42),
        "Stud": _v(675, 350, 135, 425, 725, 1_200_000, 440_000, 0.42),
    },
}

# NDS Supplement Table 4B, Southern Pine (SPIB 2013 values), by nominal width.
# Width groups: 4 (2"-4" wide), 6 (5"-6"), 8, 10, 12.  "No.3"/"Stud" = "No.3 and Stud".
_SP = {
    4: {"Select Structural": (2350, 1650, 1900, 1_800_000), "No.1": (1500, 1000, 1650, 1_600_000),
        "No.2": (1100, 675, 1450, 1_400_000), "No.3": (650, 400, 850, 1_300_000)},
    6: {"Select Structural": (2100, 1450, 1800, 1_800_000), "No.1": (1350, 875, 1550, 1_600_000),
        "No.2": (1000, 600, 1400, 1_400_000), "No.3": (575, 350, 800, 1_300_000)},
    8: {"Select Structural": (1950, 1350, 1700, 1_800_000), "No.1": (1250, 800, 1500, 1_600_000),
        "No.2": (925, 550, 1350, 1_400_000), "No.3": (525, 325, 775, 1_300_000)},
    10: {"Select Structural": (1700, 1150, 1650, 1_800_000), "No.1": (1050, 700, 1450, 1_600_000),
         "No.2": (800, 475, 1300, 1_400_000), "No.3": (475, 275, 750, 1_300_000)},
    12: {"Select Structural": (1600, 1100, 1650, 1_800_000), "No.1": (1000, 650, 1400, 1_600_000),
         "No.2": (750, 450, 1250, 1_400_000), "No.3": (450, 250, 725, 1_300_000)},
}
_EMIN = {1_800_000: 660_000, 1_600_000: 580_000, 1_400_000: 510_000, 1_300_000: 470_000}

# NDS Supplement Table 4D, DF-L timbers (5"x5" and larger).
TABLE_4D = {
    ("DF-L", "Beams & Stringers"): {
        "Select Structural": _v(1600, 950, 170, 625, 1100, 1_600_000, 580_000, 0.50),
        "No.1": _v(1350, 675, 170, 625, 925, 1_600_000, 580_000, 0.50),
        "No.2": _v(875, 425, 170, 625, 600, 1_300_000, 470_000, 0.50),
    },
    ("DF-L", "Posts & Timbers"): {
        "Select Structural": _v(1500, 1000, 170, 625, 1150, 1_600_000, 580_000, 0.50),
        "No.1": _v(1200, 825, 170, 625, 1000, 1_600_000, 580_000, 0.50),
        "No.2": _v(750, 475, 170, 625, 700, 1_300_000, 470_000, 0.50),
    },
}

# NDS Supplement Table 5A, Douglas Fir glulam (bending about x-x unless noted).
GLULAM = {
    "24F-V4": {"Fb": 2400, "Fb_neg": 1850, "Ft": 1100, "Fv": 265, "Fc_perp": 650,
               "Fc": 1650, "E": 1_800_000, "Emin": 950_000, "Ey": 1_600_000,
               "Ey_min": 850_000, "Fby": 1450, "G": 0.50},
    "24F-V8": {"Fb": 2400, "Fb_neg": 2400, "Ft": 1100, "Fv": 265, "Fc_perp": 650,
               "Fc": 1650, "E": 1_800_000, "Emin": 950_000, "Ey": 1_600_000,
               "Ey_min": 850_000, "Fby": 1450, "G": 0.50},
}

# GENERIC SCL values: VERIFY WITH MANUFACTURER (ICC-ES/APA report).
SCL = {
    "LVL 1.9E (generic)": {"Fb": 2600, "Ft": 1555, "Fv": 285, "Fc_perp": 750, "Fc": 2510,
                           "E": 1_900_000, "Emin": 966_000, "G": 0.50, "depth_exp": 0.136},
    "LVL 2.0E (generic)": {"Fb": 2900, "Ft": 2000, "Fv": 285, "Fc_perp": 750, "Fc": 3200,
                           "E": 2_000_000, "Emin": 1_016_000, "G": 0.50, "depth_exp": 0.143},
    "PSL 2.2E (generic)": {"Fb": 2900, "Ft": 2025, "Fv": 290, "Fc_perp": 625, "Fc": 2900,
                           "E": 2_200_000, "Emin": 1_120_000, "G": 0.50, "depth_exp": 0.111},
}

SPECIES = ["DF-L", "Hem-Fir", "SPF", "Southern Pine"]
GRADES = ["Select Structural", "No.1", "No.2", "No.3", "Stud"]
MATERIALS = ["Sawn lumber", "Glulam", "LVL 1.9E (generic)", "LVL 2.0E (generic)",
             "PSL 2.2E (generic)"]

SPECIFIC_GRAVITY = {"DF-L": 0.50, "Hem-Fir": 0.43, "SPF": 0.42, "Southern Pine": 0.55}

# --------------------------------------------------------------------------
# Sections
# --------------------------------------------------------------------------
_DRESSED = {2: 1.5, 3: 2.5, 4: 3.5, 5: 4.5, 6: 5.5, 8: 7.25, 10: 9.25, 12: 11.25,
            14: 13.25, 16: 15.25}


def dressed(nominal: float, timber: bool = False) -> float:
    """Dressed (actual) size of a nominal dimension, in.

    Dimension lumber (2"-4" thick): 2->1.5, 4->3.5, 6->5.5, 8->7.25, 10->9.25...
    Timbers (5"x5" and larger): nominal - 0.5 in.
    """
    if timber:
        return nominal - 0.5
    n = int(round(nominal))
    if n in _DRESSED and abs(n - nominal) < 1e-9:
        return _DRESSED[n]
    return nominal - (0.5 if nominal <= 6 else 0.75)


@dataclass
class Section:
    name: str
    b: float            # total width (in), all plies
    d: float            # depth (in)
    plies: int = 1
    nominal_b: float = 0.0
    nominal_d: float = 0.0
    kind: str = "sawn"  # sawn | timber | glulam | scl

    @property
    def A(self):
        return self.b * self.d

    @property
    def Sx(self):
        return self.b * self.d**2 / 6

    @property
    def Ix(self):
        return self.b * self.d**3 / 12

    @property
    def Sy(self):
        return self.d * self.b**2 / 6

    @property
    def Iy(self):
        return self.d * self.b**3 / 12

    @property
    def ply_b(self):
        return self.b / self.plies


def parse_size(size: str, material: str = "Sawn lumber") -> Section:
    """Parse '2x10', '(3)2x12', '3-2x12', '6x12', '5.125x18', '1.75x11.875'.

    Sawn lumber uses nominal sizes (dressed automatically). Glulam and SCL
    use actual dimensions; '(2)1.75x11.875' gives a 2-ply LVL beam.
    """
    s = size.strip().lower().replace(" ", "").replace("×", "x")
    plies = 1
    m = re.match(r"^\((\d+)\)(.*)$", s) or re.match(r"^(\d+)-(.*x.*)$", s)
    if m:
        plies, s = int(m.group(1)), m.group(2)
    parts = s.split("x")
    if len(parts) != 2:
        raise ValueError(f"size '{size}' must look like 2x10, (2)2x12 or 5.125x18")
    nb, nd = (float(eval_frac(p)) for p in parts)
    if material == "Sawn lumber":
        timber = nb >= 5 and nd >= 5
        b1, d = dressed(nb, timber), dressed(nd, timber)
        kind = "timber" if timber else "sawn"
    else:
        b1, d = nb, nd
        kind = "glulam" if material == "Glulam" else "scl"
    return Section(name=size.strip(), b=b1 * plies, d=d, plies=plies, nominal_b=nb,
                   nominal_d=nd, kind=kind)


def eval_frac(t: str) -> float:
    t = t.strip()
    if "-" in t and "/" in t:  # 1-3/4
        w, f = t.split("-")
        a, b = f.split("/")
        return float(w) + float(a) / float(b)
    if "/" in t:
        a, b = t.split("/")
        return float(a) / float(b)
    return float(t)


# --------------------------------------------------------------------------
# Reference values lookup
# --------------------------------------------------------------------------
def timber_class(sec: Section) -> str:
    """Beams & Stringers when width (depth) > thickness + 2 in. (nominal)."""
    return "Beams & Stringers" if sec.nominal_d > sec.nominal_b + 2 else "Posts & Timbers"


def reference_values(material, species, grade, sec: Section, override=None):
    """Return (values dict, table note). ``override``: dict of non-zero values."""
    note = ""
    if material == "Sawn lumber":
        if sec.kind == "timber":
            cls = timber_class(sec)
            tbl = TABLE_4D.get((species, cls))
            if tbl is None or grade not in tbl:
                if not override:
                    raise ValueError(f"NDS Table 4D values for {species} {cls} {grade} are not in "
                                     "this library: enter all values with the override fields")
                vals = {}
            else:
                vals = dict(tbl[grade])
            note = f"NDS Supplement Table 4D ({species} {cls})"
        elif species == "Southern Pine":
            w = sec.nominal_d
            grp = 4 if w <= 4 else 6 if w <= 6 else 8 if w <= 8 else 10 if w <= 10 else 12
            g = "No.3" if grade in ("No.3", "Stud") else grade
            Fb, Ft, Fc, E = _SP[grp][g]
            vals = _v(Fb, Ft, 175, 565, Fc, E, _EMIN[E], 0.55)
            note = f"NDS Supplement Table 4B (Southern Pine, {grp}\" width group)"
        else:
            vals = dict(TABLE_4A[species][grade])
            note = f"NDS Supplement Table 4A ({species})"
    elif material == "Glulam":
        g = grade if grade in GLULAM else "24F-V4"
        vals = dict(GLULAM[g])
        note = f"NDS Supplement Table 5A (DF glulam {g})"
        if grade not in GLULAM:
            note += f" - grade '{grade}' is not a glulam combination, 24F-V4 used"
    else:
        vals = dict(SCL[material])
        note = "GENERIC SCL values - verify with manufacturer evaluation report"
    for k, v in (override or {}).items():
        if v:
            vals[k] = float(v)
    for k in ("Fb", "Ft", "Fv", "Fc_perp", "Fc", "E", "Emin"):
        if k not in vals:
            raise ValueError(f"reference value {k} missing: enter it in the override fields")
    return vals, note


# --------------------------------------------------------------------------
# Adjustment factors
# --------------------------------------------------------------------------
CD_BY_CASE = {"D": 0.9, "L": 1.0, "S": 1.15, "Lr": 1.25, "R": 1.25, "W": 1.6, "E": 1.6}


def CD_for_combo(factors: dict) -> float:
    """Load duration factor of an ASD combination: the shortest-duration
    load present governs (NDS 2.3.2, Appendix B)."""
    present = [k for k, v in factors.items() if v]
    return max(CD_BY_CASE.get(k, 1.0) for k in present) if present else 0.9


def lambda_for_combo(factors: dict) -> float:
    """NDS Table N3 time effect factor for an ASCE 7 strength combination."""
    if factors.get("W", 0) or factors.get("E", 0):
        return 1.0
    keys = {k for k, v in factors.items() if v}
    if keys == {"D"}:
        return 0.6
    return 0.8  # occupancy live / roof live / snow principal (0.7 for storage live)


def size_factor(sec: Section, species: str, grade: str, material="Sawn lumber"):
    """Return dict CF for Fb, Ft, Fc."""
    one = {"Fb": 1.0, "Ft": 1.0, "Fc": 1.0}
    if material != "Sawn lumber":
        return one
    if sec.kind == "timber":
        d = sec.d
        return {"Fb": (12 / d) ** (1 / 9) if d > 12 else 1.0, "Ft": 1.0, "Fc": 1.0}
    w, t = sec.nominal_d, sec.nominal_b
    if species == "Southern Pine":
        cf = 1.0
        if t >= 4 and w >= 8:
            cf = 1.1
        if w > 12:
            return {"Fb": 0.9 * cf, "Ft": 0.9, "Fc": 0.9}
        return {"Fb": cf, "Ft": 1.0, "Fc": 1.0}
    if grade == "Stud":
        if w <= 4:
            return {"Fb": 1.1, "Ft": 1.1, "Fc": 1.05}
        if w <= 6:
            return {"Fb": 1.0, "Ft": 1.0, "Fc": 1.0}
        raise ValueError("Stud grade not tabulated wider than 6 in. (use No.3)")
    # Select Structural, No.1, No.2, No.3 (Table 4A adjustment factors)
    table = [  # width, Fb(2-3 thick), Fb(4 thick), Ft, Fc
        (4, 1.5, 1.5, 1.5, 1.15), (5, 1.4, 1.4, 1.4, 1.1), (6, 1.3, 1.3, 1.3, 1.1),
        (8, 1.2, 1.3, 1.2, 1.05), (10, 1.1, 1.2, 1.1, 1.0), (12, 1.0, 1.1, 1.0, 1.0),
        (99, 0.9, 1.0, 0.9, 0.9)]
    for wmax, fb23, fb4, ft, fc in table:
        if w <= wmax:
            return {"Fb": fb4 if t >= 4 else fb23, "Ft": ft, "Fc": fc}
    return one


def wet_service(material, sec: Section, vals, CF):
    """CM factors for wet service (MC > 19 %)."""
    if material == "Glulam":
        return {"Fb": 0.8, "Ft": 0.8, "Fv": 0.875, "Fc_perp": 0.53, "Fc": 0.73, "E": 0.833}
    if material != "Sawn lumber":
        return {"Fb": 1.0, "Ft": 1.0, "Fv": 1.0, "Fc_perp": 1.0, "Fc": 1.0, "E": 1.0}
    if sec.kind == "timber":
        return {"Fb": 1.0, "Ft": 1.0, "Fv": 1.0, "Fc_perp": 0.67, "Fc": 0.91, "E": 1.0}
    return {"Fb": 1.0 if vals["Fb"] * CF["Fb"] <= 1150 else 0.85, "Ft": 1.0, "Fv": 0.97,
            "Fc_perp": 0.67, "Fc": 1.0 if vals["Fc"] * CF["Fc"] <= 750 else 0.8, "E": 0.9}


def temperature(temp="T <= 100F", wet=False):
    """NDS Table 2.3.3 Ct."""
    if temp == "100F < T <= 125F":
        o = 0.7 if wet else 0.8
        return {"Ft": 0.9, "E": 0.9, "Fb": o, "Fv": o, "Fc": o, "Fc_perp": o}
    if temp == "125F < T <= 150F":
        o = 0.5 if wet else 0.7
        return {"Ft": 0.9, "E": 0.9, "Fb": o, "Fv": o, "Fc": o, "Fc_perp": o}
    return {k: 1.0 for k in ("Ft", "E", "Fb", "Fv", "Fc", "Fc_perp")}


def volume_factor(L_ft, sec: Section, southern_pine=False):
    """Glulam CV (NDS 5.3.6)."""
    x = 20 if southern_pine else 10
    return min((21 / L_ft) ** (1 / x) * (12 / sec.d) ** (1 / x) * (5.125 / sec.b) ** (1 / x), 1.0)


def bearing_area_factor(lb_in, at_end=False):
    """NDS 3.10.4 Cb (bearing < 6 in. long, not nearer than 3 in. to the end)."""
    if at_end or lb_in >= 6:
        return 1.0
    return (lb_in + 0.375) / lb_in


def le_beam(lu_in, d_in, load="general"):
    """NDS Table 3.3.3 effective length (in).

    load: 'general' (single span or cantilever, any load - conservative),
          'uniform' (single span uniform), 'cantilever_uniform', 'center_point'.
    """
    r = lu_in / d_in if d_in else 0.0
    if load == "uniform":
        return 2.06 * lu_in if r < 7 else 1.63 * lu_in + 3 * d_in
    if load == "center_point":
        return 1.80 * lu_in if r < 7 else 1.37 * lu_in + 3 * d_in
    if load == "cantilever_uniform":
        return 1.33 * lu_in if r < 7 else 0.90 * lu_in + 3 * d_in
    if r < 7:
        return 2.06 * lu_in
    if r <= 14.3:
        return 1.63 * lu_in + 3 * d_in
    return 1.84 * lu_in


def beam_stability(Fb_star, Emin_p, le_in, d, b):
    """NDS 3.3.3: return (CL, RB, FbE)."""
    if le_in <= 0:
        return 1.0, 0.0, math.inf
    RB = math.sqrt(le_in * d / b**2)
    FbE = 1.20 * Emin_p / RB**2
    a = (1 + FbE / Fb_star) / 1.9
    CL = a - math.sqrt(a**2 - (FbE / Fb_star) / 0.95)
    return min(CL, 1.0), RB, FbE


def column_stability(Fc_star, Emin_p, le_over_d, c=0.8):
    """NDS 3.7.1: return (CP, FcE)."""
    if le_over_d <= 0:
        return 1.0, math.inf
    FcE = 0.822 * Emin_p / le_over_d**2
    a = (1 + FcE / Fc_star) / (2 * c)
    CP = a - math.sqrt(a**2 - (FcE / Fc_star) / c)
    return CP, FcE


def column_c(material):
    return 0.8 if material == "Sawn lumber" else 0.9


# LRFD format conversion and resistance factors (NDS Table 4.3.1 / N1, N2)
KF = {"Fb": 2.54, "Ft": 2.70, "Fv": 2.88, "Fc": 2.40, "Fc_perp": 1.67, "Emin": 1.76}
PHI = {"Fb": 0.85, "Ft": 0.80, "Fv": 0.75, "Fc": 0.90, "Fc_perp": 0.90, "Emin": 0.85}


def time_factor(prop, method, factors):
    """ASD: CD (none for Fc_perp, E, Emin). LRFD: KF*phi*lambda (lambda not on
    Fc_perp/Emin)."""
    if method.upper() == "ASD":
        return CD_for_combo(factors) if prop in ("Fb", "Ft", "Fv", "Fc") else 1.0
    lam = lambda_for_combo(factors) if prop in ("Fb", "Ft", "Fv", "Fc") else 1.0
    return KF.get(prop, 1.0) * PHI.get(prop, 1.0) * lam


@dataclass
class Member:
    """Material + section + service conditions; computes adjusted values."""
    material: str
    species: str
    grade: str
    sec: Section
    wet: bool = False
    temp: str = "T <= 100F"
    incised: bool = False
    repetitive: bool = False
    flat_use: bool = False
    override: dict = field(default_factory=dict)

    def __post_init__(self):
        self.ref, self.table_note = reference_values(self.material, self.species, self.grade,
                                                     self.sec, self.override)
        self.CF = size_factor(self.sec, self.species, self.grade, self.material)
        if self.material.startswith(("LVL", "PSL")):
            ex = self.ref.get("depth_exp", 0.136)
            self.CF = {"Fb": min((12 / self.sec.d) ** ex, 1.15) if self.sec.d else 1.0,
                       "Ft": 1.0, "Fc": 1.0}
        self.CM = wet_service(self.material, self.sec, self.ref, self.CF) if self.wet else \
            {k: 1.0 for k in ("Fb", "Ft", "Fv", "Fc_perp", "Fc", "E")}
        self.Ct = temperature(self.temp, self.wet)
        inc = self.incised and self.material == "Sawn lumber"
        self.Ci = {"Fb": 0.8 if inc else 1.0, "Ft": 0.8 if inc else 1.0,
                   "Fv": 0.8 if inc else 1.0, "Fc": 0.8 if inc else 1.0,
                   "Fc_perp": 1.0, "E": 0.95 if inc else 1.0}
        self.Cr = 1.15 if (self.repetitive and self.material == "Sawn lumber"
                           and self.sec.kind == "sawn") else 1.0
        self.Cfu = 1.0  # flat use not used for edgewise beams

    @property
    def G(self):
        return self.ref.get("G", 0.5)

    def E_adj(self):
        return self.ref["E"] * self.CM["E"] * self.Ct["E"] * self.Ci["E"]

    def Emin_adj(self, method="ASD"):
        f = KF["Emin"] * PHI["Emin"] if method.upper() == "LRFD" else 1.0
        return self.ref["Emin"] * self.CM["E"] * self.Ct["E"] * self.Ci["E"] * f

    def Fb_star(self, method, factors):
        """Fb with all factors except CL and CV (for CL computation)."""
        r = self.ref
        return (r["Fb"] * time_factor("Fb", method, factors) * self.CM["Fb"] * self.Ct["Fb"]
                * self.CF["Fb"] * self.Ci["Fb"] * self.Cr)

    def Fv_adj(self, method, factors):
        return (self.ref["Fv"] * time_factor("Fv", method, factors) * self.CM["Fv"]
                * self.Ct["Fv"] * self.Ci["Fv"])

    def Ft_adj(self, method, factors):
        return (self.ref["Ft"] * time_factor("Ft", method, factors) * self.CM["Ft"]
                * self.Ct["Ft"] * self.CF["Ft"] * self.Ci["Ft"])

    def Fc_star(self, method, factors):
        return (self.ref["Fc"] * time_factor("Fc", method, factors) * self.CM["Fc"]
                * self.Ct["Fc"] * self.CF["Fc"] * self.Ci["Fc"])

    def Fc_perp_adj(self, method="ASD", Cb=1.0):
        f = KF["Fc_perp"] * PHI["Fc_perp"] if method.upper() == "LRFD" else 1.0
        return self.ref["Fc_perp"] * self.CM["Fc_perp"] * self.Ct["Fc_perp"] * Cb * f

    def factor_table(self):
        return [
            {"Factor": "CF (size) Fb/Ft/Fc", "Value": f"{self.CF['Fb']:.3f} / {self.CF['Ft']:.3f} / {self.CF['Fc']:.3f}"},
            {"Factor": "CM (wet service) Fb", "Value": f"{self.CM['Fb']:.3f}"},
            {"Factor": "Ct (temperature) Fb", "Value": f"{self.Ct['Fb']:.3f}"},
            {"Factor": "Ci (incising) Fb", "Value": f"{self.Ci['Fb']:.3f}"},
            {"Factor": "Cr (repetitive)", "Value": f"{self.Cr:.3f}"},
        ]


def reference_inputs(report, m: Member):
    """Add reference values to a report."""
    r = m.ref
    report.result("Reference values", m.table_note)
    report.result("Fb / Ft / Fv", f"{r['Fb']:.0f} / {r['Ft']:.0f} / {r['Fv']:.0f}", "psi")
    report.result("Fc_perp / Fc", f"{r['Fc_perp']:.0f} / {r['Fc']:.0f}", "psi")
    report.result("E / Emin", f"{r['E']:,.0f} / {r['Emin']:,.0f}", "psi")
