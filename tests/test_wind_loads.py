"""Hand checks for ASCE 7-22 wind loads."""

import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))

import pytest  # noqa: E402

from structcalc.calcs import wind_loads as w  # noqa: E402

approx = pytest.approx


def test_Kz_asce722_constants():
    # Table 26.10-1 (ASCE 7-22): 0-15 ft C = 0.85, D = 1.03; formula 2.41 (z/zg)^(2/alpha)
    assert w.Kz(15, "C") == approx(0.85, abs=0.005)
    assert w.Kz(10, "D") == approx(1.03, abs=0.01)
    assert w.Kz(30, "C") == approx(2.41 * (30 / 2460) ** (2 / 9.8), rel=1e-12)
    assert w.Kz(60, "B") == approx(2.41 * (60 / 3280) ** (2 / 7.5), rel=1e-12)


def test_qz_has_no_Kd():
    q = w.qz(30, 115, "C")
    assert q == approx(0.00256 * w.Kz(30, "C") * 115**2, rel=1e-12)
    assert w.Ke(5000) == approx(math.exp(-0.0000362 * 5000))


def test_leeward_wall_cp():
    assert w.cp_leeward_wall(0.5) == -0.5
    assert w.cp_leeward_wall(2.0) == -0.3
    assert w.cp_leeward_wall(3.0) == approx(-0.25)
    assert w.cp_leeward_wall(6.0) == -0.2


def test_roof_cp_table_points():
    assert w.cp_roof_windward(20, 0.25) == approx((-0.3, 0.2))
    assert w.cp_roof_windward(30, 1.0) == approx((-0.3, 0.2))
    assert w.cp_roof_leeward(10, 0.5) == approx(-0.5)
    assert w.cp_roof_windward(70, 0.5) == approx((0.7, 0.7))
    z = w.flat_roof_zones(20, 100)  # h/L = 0.2
    assert [round(c[2], 3) for c in z] == [-0.9, -0.9, -0.5, -0.3]


def test_base_shear_flat_box_hand_calc():
    # V=100 mph, Exp C, flat roof, eave 15 ft, B=100, L=50 (L/B=0.5 -> leeward -0.5)
    o = w.run(V=100, exposure="C", roof_type="Flat", h_eave=15, B_ridge=100, L_across=50,
              use_Ke=False)
    q = 0.00256 * w.Kz(15, "C") * 100**2
    p_ww = q * 0.85 * 0.85 * 0.8
    p_lw = q * 0.85 * 0.85 * 0.5
    F = (p_ww + p_lw) * 15 * 100 / 1000
    d1 = o["directions"][0]
    assert d1["Base shear (k)"] == approx(F, rel=1e-6)
    assert d1["Overturning (k-ft)"] == approx(F * 7.5, rel=1e-6)
    assert d1["Minimum 27.1.5 (k)"] == approx(16 * 1500 / 1000)


def test_cc_wall_zone5():
    # A = 10 ft^2, theta <= 10 -> 0.9 * (-1.4) = -1.26 ; p = qh Kd (GCp - GCpi)
    o = w.run(V=115, roof_type="Flat", h_eave=30, cc_area=10, use_Ke=False)
    qh = 0.00256 * w.Kz(30, "C") * 115**2
    z5 = o["cc"][1]
    assert z5["GCp -"] == approx(-1.26)
    assert z5["p - (psf)"] == approx(qh * 0.85 * (-1.26 - 0.18), rel=1e-9)
