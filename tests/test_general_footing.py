"""Hand checks: soil pressure solver (foundation.py) and general footing (07)."""

import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))

import pytest  # noqa: E402

from structcalc import foundation as fd  # noqa: E402
from structcalc.calcs import general_footing as gf  # noqa: E402

approx = pytest.approx


# --- soil pressure ---------------------------------------------------------
def test_kern_biaxial_closed_form():
    p = fd.bearing_pressure(100, 0.5, 0.3, 8, 6)
    assert p.full_contact
    assert p.qmax == approx(100 / 48 * (1 + 6 * 0.5 / 8 + 6 * 0.3 / 6))
    assert p.qmin == approx(100 / 48 * (1 - 6 * 0.5 / 8 - 6 * 0.3 / 6))


def test_one_way_partial_contact():
    # q = 2P / (3 B (L/2 - e)), contact length 3 (L/2 - e)
    p = fd.bearing_pressure(100, 2.0, 0.0, 8, 6)
    assert p.qmax == approx(2 * 100 / (3 * 6 * (4 - 2.0)), rel=1e-9)
    assert p.contact_ratio == approx(3 * (4 - 2.0) / 8, rel=1e-9)


def test_corner_contact_tetrahedron():
    # triangular corner contact with legs a=6, b=8: resultant at a/4, b/4 from
    # the corner, qmax = 6P/(a b)
    P = 100
    p = fd.bearing_pressure(P, 5 - 6 / 4, 5 - 8 / 4, 10, 10)
    assert p.qmax == approx(6 * P / 48, rel=1e-9)
    assert p.contact_ratio == approx(0.5 * 6 * 8 / 100, rel=1e-9)


def test_biaxial_partial_equilibrium():
    p = fd.bearing_pressure(50, 1.5, 1.0, 8, 6)
    Q, Qx, Qy = p.integrals()
    assert not p.full_contact
    assert (Q, Qx, Qy) == (approx(50), approx(75), approx(50))


def test_resultant_outside_footing():
    assert fd.bearing_pressure(10, 4.1, 0, 8, 6) is None
    assert fd.bearing_pressure(-5, 0, 0, 8, 6) is None


def test_pressure_1d():
    assert fd.pressure_1d(120, 1.0, 12) == approx((5.0, 15.0, 12))
    ql, qr, c = fd.pressure_1d(120, -3.0, 12)
    assert ql == approx(2 * 120 / (3 * 3)) and qr == 0 and c == approx(9)


# --- general footing --------------------------------------------------------
def _ck(r, name):
    return next(c for c in r.checks if c.name.startswith(name))


def test_concentric_footing_hand_calc():
    # 8x8 ft, 24 in, 18x18 col, Pu = 1.2*80 + 1.6*60 = 192 k -> qu = 3.0 ksf net
    o = gf.run(Lx=8, Ly=8, h=24, col_x=18, col_y=18, bar_x=5, s_x=12, bar_y=5, s_y=12,
               col_loads={"D": (80, 0, 0, 0, 0), "L": (60, 0, 0, 0, 0)})
    r = o["report"]
    assert _ck(r, "Flexure, bars along x").demand == approx(3.0 * 8 * 3.25**2 / 2, rel=1e-6)
    dx = 24 - 3 - 0.625 / 2
    assert _ck(r, "One-way shear, x").demand == approx(3.0 * 8 * (3.25 - dx / 12), rel=1e-6)
    davg = (dx + (24 - 3 - 0.625 - 0.3125)) / 2
    b = (18 + davg) / 12
    Vu2 = 192 - 3.0 * b * b
    assert _ck(r, "Two-way").demand == approx(Vu2 * 1000 / (4 * (18 + davg) * davg), rel=1e-6)
    # service bearing: (80+60 + footing 19.2 + soil 7.04) / 64
    assert _ck(r, "Soil bearing").demand == approx((140 + 19.2 + 7.04) / 64, rel=1e-9)
    assert _ck(r, "Minimum steel (x)").demand == approx(0.0018 * 96 * 24)


def test_footing_partial_contact_service():
    o = gf.run(Lx=8, Ly=6, h=18, soil_cover=0, col_loads={"D": (50, 0, 0, 0, 0)})
    W = o["W"]
    Ptot = 50 + W
    o = gf.run(Lx=8, Ly=6, h=18, soil_cover=0, col_loads={"D": (50, 0, 2.0 * Ptot, 0, 0)})
    p = o["pressures"]["1: D"]
    assert p.ex == approx(2.0)
    assert p.qmax == approx(2 * Ptot / (3 * 6 * (4 - 2.0)), rel=1e-9)


def test_shear_adds_base_moment():
    o = gf.run(Lx=8, Ly=8, h=24, soil_cover=0, ped_height=12, col_loads={"D": (100, 0, 0, 10, 0)})
    p = o["pressures"]["1: D"]
    assert p.ex == approx(10 * 3.0 / (100 + o["W"]))  # V * (pedestal 1 ft + footing 2 ft)


def test_uplift_and_overturning_flags():
    o = gf.run(Lx=6, Ly=6, h=18, col_loads={"D": (5, 0, 0, 0, 0), "W": (-40, 0, 0, 0, 0)})
    up = _ck(o["report"], "Uplift")
    assert not up.ok
