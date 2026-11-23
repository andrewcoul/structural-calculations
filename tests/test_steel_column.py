"""Hand checks for the steel column module (AISC 360-22)."""

import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))

import pytest  # noqa: E402

from structcalc import steel  # noqa: E402
from structcalc.calcs import steel_column  # noqa: E402

approx = pytest.approx


def test_axial_only_matches_E3():
    # W14X90, Lc = 14 ft, Fy = 50: AISC Manual Table 4-1a gives phiPn = 1030 kips.
    o = steel_column.run("W14X90", H=14, Lcx=14, Lcy=14, Lcz=14,
                         case_loads={"D": {"P": 100}, "L": {"P": 200}}, self_weight=False)
    ax = next(c for c in o["report"].checks if c.name == "Axial")
    assert ax.capacity == approx(1030, abs=8)
    assert ax.demand == approx(1.2 * 100 + 1.6 * 200)


def test_hand_calc_Fcr():
    p = steel.shape("W14X90")
    KLr = 14 * 12 / p["ry"]
    Fe = math.pi**2 * 29000 / KLr**2
    Fcr = 0.658 ** (50 / Fe) * 50
    c = steel.compression(p, 50, 14 * 12, 14 * 12)
    assert c.Fcr == approx(Fcr, rel=1e-9)


def test_H1_1a_with_eccentricity_and_B1():
    # W10X33, H = Lc = 14 ft, P_D=60, P_L=80, ex=6 in at top only (bottom pinned):
    o = steel_column.run("W10X33", H=14, Lcx=14, Lcy=14, Lcz=14, Lb=14,
                         case_loads={"D": {"P": 60}, "L": {"P": 80}}, ex=6.0,
                         self_weight=False, amplify_B1=True)
    row, g = o["governing"]
    Pu = 1.2 * 60 + 1.6 * 80
    assert row["P (k)"] == approx(Pu)
    Mu1 = Pu * 6 / 12  # moment at top, zero at bottom -> Cm = 0.6
    p = steel.shape("W10X33")
    Pe1 = math.pi**2 * 29000 * p["Ix"] / (14 * 12) ** 2
    B1 = max(0.6 / (1 - Pu / Pe1), 1.0)
    assert row["B1x"] == approx(B1, rel=1e-9)
    assert g["Mrx"] == approx(B1 * Mu1, rel=1e-9)
    # Cb for a linear diagram from M to 0: 12.5/(2.5+3*0.25... ) -> 1.67
    assert row["Cb"] == approx(12.5 / (2.5 + 3 * 0.75 + 4 * 0.5 + 3 * 0.25), rel=1e-6)
    comp = steel.compression(p, 50, 168, 168, 168)
    Pc = 0.9 * comp.Pn
    Mcx = 0.9 * steel.flexure(p, 50, "x", 168, row["Cb"]).Mn / 12
    expect = Pu / Pc + 8 / 9 * (B1 * Mu1 / Mcx) if Pu / Pc >= 0.2 else Pu / (2 * Pc) + B1 * Mu1 / Mcx
    assert row["H1 ratio"] == approx(expect, rel=1e-9)


def test_lateral_load_moment():
    # 0.2 klf wind on 12 ft pinned column: M = wL^2/8 = 3.6 kip-ft, 1.0W LRFD
    o = steel_column.run("HSS6X6X3/8", H=12, Lcx=12, Lcy=12, Lcz=12,
                         case_loads={"D": {"P": 10}}, wx={"W": 0.2}, self_weight=False,
                         amplify_B1=False)
    mx = max(r["Mrx (k-ft)"] for r in o["rows"])
    assert mx == approx(3.6, rel=1e-3)


def test_tension_case():
    o = steel_column.run("W8X31", H=10, Lcx=10, Lcy=10, Lcz=10,
                         case_loads={"D": {"P": 5}, "W": {"P": -60}}, self_weight=False)
    p = steel.shape("W8X31")
    ax = next(c for c in o["report"].checks if c.name == "Axial")
    # governing is tension: 0.9D - 1.0W -> -55.5 kip (or compression 1.2D+W)
    rows = {r["Combination"]: r for r in o["rows"]}
    assert any(r["P (k)"] < 0 for r in rows.values())
    assert 0.9 * 50 * p["area"] == approx(steel.tension_yield(p, 50)[0])
    assert ax.ok


def test_autoselect():
    o = steel_column.run(auto_select=True, depth_filter="10", H=14, Lcx=14, Lcy=14, Lcz=14,
                         case_loads={"D": {"P": 100}, "L": {"P": 150}})
    assert o["report"].ok and o["section"]["name"].startswith("W10X")
