"""Hand-check tests: wood.py, 03 wood beam, 04 wood column, 10 wood shear wall."""

import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))

import pytest  # noqa: E402

from structcalc import wood  # noqa: E402
from structcalc.calcs import wood_beam, wood_column, wood_shear_wall  # noqa: E402

approx = pytest.approx


# --- wood.py ----------------------------------------------------------------
def test_dressed_sizes():
    assert (wood.parse_size("2x10").b, wood.parse_size("2x10").d) == (1.5, 9.25)
    s = wood.parse_size("(3)2x12")
    assert (s.b, s.d, s.plies) == (4.5, 11.25, 3)
    s = wood.parse_size("6x12")
    assert (s.b, s.d, s.kind) == (5.5, 11.5, "timber")
    assert wood.timber_class(s) == "Beams & Stringers"
    assert wood.timber_class(wood.parse_size("8x8")) == "Posts & Timbers"


def test_table_values():
    v, _ = wood.reference_values("Sawn lumber", "DF-L", "No.2", wood.parse_size("2x10"))
    assert (v["Fb"], v["Fv"], v["Fc"], v["E"], v["Emin"]) == (900, 180, 1350, 1_600_000, 580_000)
    v, _ = wood.reference_values("Sawn lumber", "Southern Pine", "No.2", wood.parse_size("2x10"))
    assert (v["Fb"], v["Ft"], v["Fc"], v["E"]) == (800, 475, 1300, 1_400_000)  # SPIB 2013 / Table 4B
    v, _ = wood.reference_values("Sawn lumber", "Southern Pine", "No.2", wood.parse_size("2x4"))
    assert (v["Fb"], v["Emin"]) == (1100, 510_000)


def test_size_factor():
    assert wood.size_factor(wood.parse_size("2x10"), "DF-L", "No.2")["Fb"] == 1.1
    assert wood.size_factor(wood.parse_size("4x10"), "DF-L", "No.2")["Fb"] == 1.2
    assert wood.size_factor(wood.parse_size("2x4"), "DF-L", "Stud")["Fc"] == 1.05
    assert wood.size_factor(wood.parse_size("6x16"), "DF-L", "No.1")["Fb"] == approx((12 / 15.5) ** (1 / 9))


def test_CD_combo():
    assert wood.CD_for_combo({"D": 1.0}) == 0.9
    assert wood.CD_for_combo({"D": 1.0, "S": 0.7}) == 1.15
    assert wood.CD_for_combo({"D": 0.6, "W": -0.6}) == 1.6


def test_CP_hand():
    # 6x6 DF-L No.1 P&T, le = 12 ft: Fc* = 1000, Emin = 580000, le/d = 26.18
    CP, FcE = wood.column_stability(1000, 580_000, 144 / 5.5, 0.8)
    assert FcE == approx(0.822 * 580_000 / (144 / 5.5) ** 2)
    a = (1 + FcE / 1000) / 1.6
    assert CP == approx(a - math.sqrt(a * a - FcE / 1000 / 0.8))
    assert CP == approx(0.5561, abs=1e-4)


# --- 03 wood beam -------------------------------------------------------------
def test_wood_beam_2x10():
    o = wood_beam.run("2x10", w={"D": 0.1, "L": 0.2}, span=12, self_weight=False, repetitive=True)
    ck = {c.name: c for c in o["report"].checks}
    assert ck["Bending fb"].demand == approx(5.4 * 12000 / (1.5 * 9.25**2 / 6), rel=1e-5)  # 3029 psi
    assert ck["Bending fb"].capacity == approx(900 * 1.0 * 1.1 * 1.15, rel=1e-9)          # CD CF Cr
    assert ck["Shear fv = 1.5V/A"].demand == approx(1.5 * 1800 / 13.875, rel=1e-5)
    dL = 5 * (0.2 / 12) * 144**4 / (384 * 1600 * (1.5 * 9.25**3 / 12))
    assert ck["Deflection Live (L) - Span 1"].demand == approx(dL, rel=1e-5)


def test_wood_beam_glulam_CL():
    o = wood_beam.run("5.125x16.5", material="Glulam", grade="24F-V4", w={"D": 0.3, "S": 0.5},
                      span=24, bracing="At supports only", self_weight=False)
    seg = o["governing"]["flex"][1]
    # hand: lu/d = 17.45 > 14.3 -> le = 1.84 lu; Fb* = 2400*1.15
    le = 1.84 * 288
    RB2 = le * 16.5 / 5.125**2
    FbE = 1.2 * 950_000 / RB2
    x = FbE / 2760
    a = (1 + x) / 1.9
    CL = a - math.sqrt(a * a - x / 0.95)
    assert seg["CL"] == approx(CL, rel=1e-6)
    assert seg["Fbp"] == approx(2760 * CL, rel=1e-6)  # CL < CV governs


def test_wood_beam_autoselect():
    o = wood_beam.run(auto_select=True, w={"D": 0.1, "L": 0.2}, span=12)
    assert o["report"].ok


def test_wood_beam_lrfd():
    o = wood_beam.run("2x12", method="LRFD", w={"D": 0.1, "L": 0.2}, span=12, self_weight=False)
    flex = next(c for c in o["report"].checks if c.name == "Bending fb")
    # 1.2D+1.6L governs, lambda=0.8, KF*phi = 2.54*0.85, CF = 1.0
    assert flex.capacity == approx(900 * 2.54 * 0.85 * 0.8 * 1.0, rel=1e-6)
    assert flex.demand == approx((1.2 * 0.1 + 1.6 * 0.2) * 144 / 8 * 12000 / (1.5 * 11.25**2 / 6), rel=1e-5)


# --- 04 wood column ----------------------------------------------------------
def test_column_6x6():
    o = wood_column.run(size="6x6", grade="No.1", stud_wall=False, weak_bracing="Same as strong axis",
                        height=12, P={"D": 5, "L": 10})
    g = o["governing"]
    assert g["CP"] == approx(0.5561, abs=1e-4)
    assert g["Fcp"] == approx(1000 * 0.5561, rel=1e-3)


def test_stud_wind_interaction():
    o = wood_column.run(size="2x6", grade="Stud", height=10, P={"D": 0.4, "L": 0.8, "S": 0.6},
                        wind_psf=30, trib_width=16)
    g = o["governing"]
    # D+0.6W, CD=1.6: Fc* = 850*1.6, CP=0.5288; fb1 = 0.6*40*10^2/8*12/Sx
    fc = 400 / 8.25
    FcE = 0.822 * 510_000 / (120 / 5.5) ** 2
    Fcp = 850 * 1.6 * 0.52883
    fb = 0.6 * 40 * 100 / 8 * 12 / (1.5 * 5.5**2 / 6)
    expect = (fc / Fcp) ** 2 + fb / (700 * 1.6 * 1.15 * (1 - fc / FcE))
    assert g["ratio"] == approx(expect, rel=1e-3)


# --- 10 wood shear wall ------------------------------------------------------
def test_shear_wall_capacity():
    o = wood_shear_wall.run(Eh=4.0, W=5.0, segments="8, 4", height=9)
    ck = {c.name: c for c in o["report"].checks}
    cap = 980 * (8 + 4 * (1.25 - 0.125 * 9 / 4)) / 2.8 / 1000
    assert ck["Line shear capacity (seismic)"].capacity == approx(cap, rel=1e-9)
    assert ck["Line shear capacity (wind)"].capacity == approx(cap * 2.8 / 2.0, rel=1e-9)


def test_shear_wall_published_example():
    """StructSuite SDPWS 2021 guide (Nov 2026): (2)2x4 DF No.2 end post,
    le = 103.5 in -> Pallow = 5,438 lb; drift terms 0.021/0.137/0.168 in."""
    o = wood_shear_wall.run(sheathing="Sheathing 7/16 in. 8d", spacing=3, height=9, segments="6",
                            Eh=426.1 * 6 / 1000, W=0, rho=1.0, chord="4x4", chord_grade="No.2",
                            hd_defl=0.112, Cd=1.0, wall_D=0)
    res = {i.name: i.value for i in o["report"].results}
    d1, d2, d3 = (float(t) for t in res["delta bending / shear / anchorage"].split("/"))
    assert (d1, d2, d3) == (approx(0.021, abs=1e-3), approx(0.137, abs=1e-3), approx(0.168, abs=1e-3))
    o2 = wood_shear_wall.run(Eh=4.0, W=0, segments="8", height=9, chord="(2)2x4", wall_D=0)
    chord = next(c for c in o2["report"].checks if c.name.startswith("Compression chord"))
    assert chord.capacity == approx(5.438, abs=0.002)


def test_table_43a_values():
    assert wood_shear_wall.table_43a("Sheathing 7/16 in. 8d", 6) == (670, 15)
    assert wood_shear_wall.table_43a("Sheathing 7/16 in. 8d", 2) == (1640, 42)
    assert wood_shear_wall.table_43a("Sheathing 15/32 in. 8d", 4)[0] == 1065
    assert wood_shear_wall.table_43a("Sheathing 15/32 in. 10d", 3)[0] == 1680
