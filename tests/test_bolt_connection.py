"""Hand checks for NDS bolt connections (wood_conn + calcs.bolt_connection)."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))

import pytest  # noqa: E402

from structcalc import wood_conn as wc  # noqa: E402
from structcalc.calcs import bolt_connection  # noqa: E402

approx = pytest.approx


def test_dowel_bearing():
    assert wc.fe_wood(0.5, 0.625, 0) == approx(5600)
    assert wc.fe_wood(0.5, 0.625, 90) == approx(2824, abs=1)      # Caltrans D-13
    assert wc.fe_wood(0.5, 0.625, 50) == approx(3551, abs=1)      # Hankinson
    assert wc.fe_wood(0.5, 0.162) == approx(16600 * 0.5**1.84)


def test_rd():
    assert wc.rd_factors(0.162)["IV"] == 2.2
    assert wc.rd_factors(0.20)["Im"] == approx(2.5)
    rd = wc.rd_factors(0.625, 50)
    assert rd["Im"] == approx(4.556, abs=1e-3) and rd["II"] == approx(4.1, abs=1e-3)
    assert rd["IV"] == approx(3.644, abs=1e-3)


def test_caltrans_D13_double_shear():
    """Caltrans Falsework Manual App. D Example 13 (2020): 6 x 5/8 in bolts,
    12x12 DF-L post at 50 deg, 2x8 DF-L braces. Z = 1389 lb (IIIs), Cg = 0.99,
    Z' = 13,201 lb (with Cg rounded to 0.99)."""
    o = bolt_connection.run(D=0.625, shear="Double shear", t_main=12, w_main=12, t_side=1.5,
                            w_side=7.25, E_main=1.3e6, E_side=1.6e6, theta_main=50,
                            theta_side=0, n_rows=2, n_per_row=3, spacing=2.5, row_spacing=4.75,
                            end_main=6, edge_main=3.6, end_side=5, edge_side=1.25, CD=1.6,
                            load=13000)
    assert o["mode"] == "IIIs"
    assert o["Z"] == approx(1389, abs=1)
    assert o["modes"]["Im"] == approx(5846, abs=2)
    assert o["modes"]["Is"] == approx(2305, abs=1)
    assert o["modes"]["IV"] == approx(1731, abs=1)
    assert o["Cg"] == approx(0.99, abs=0.005)
    assert o["Zp_total"] == approx(6 * 1389 * 1.6 * o["Cg"], rel=1e-3)
    assert o["report"].ok


def test_single_shear_table_value():
    """1/2 in bolt, 1-1/2 in DF-L side and main, parallel: yield equations
    give 483 lb (mode II); NDS Table 12A tabulates rounded values."""
    Fe = wc.fe_wood(0.5, 0.5, 0)
    m = wc.yield_modes(0.5, 1.5, 1.5, Fe, Fe, 45000)
    z, mode = wc.governing(m)
    assert mode == "II"
    assert z == approx(483, abs=1)


def test_geometry_factor():
    c, d = wc.geometry_factor(0.625, "parallel", 3.0, 2.5, "tension", True, 3)
    assert d["end_full"] == approx(4.375) and d["end_min"] == approx(2.1875)
    assert c == approx(3.0 / 4.375)
    c, _ = wc.geometry_factor(0.625, "parallel", 2.0, 2.5, "tension", True, 3)
    assert c == 0.0  # below minimum
    assert wc.edge_distance_min(0.625, "perpendicular", 2.4, loaded_edge=True) == approx(2.5)
    assert wc.row_spacing_min(0.625, "perpendicular", 8) == approx(3.125)


def test_lrfd_conversion():
    o_asd = bolt_connection.run(n_per_row=1, CD=1.0, load=100)
    o_lrfd = bolt_connection.run(n_per_row=1, method="LRFD", lam=0.8, load=100)
    assert o_lrfd["Zp"] / o_asd["Zp"] == approx(3.32 * 0.65 * 0.8, rel=1e-9)


def test_steel_side_plate():
    o = bolt_connection.run(D=0.75, side_material="Steel plate", t_side=0.25, t_main=3.5,
                            n_per_row=1, load=1000)
    assert o["Fes"] == 87000
    assert o["report"].ok
