"""Hand checks for earth pressure coefficients and the retaining wall (09)."""

import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))

import pytest  # noqa: E402

from structcalc import foundation as fd  # noqa: E402
from structcalc.calcs import retaining_wall as rw  # noqa: E402

approx = pytest.approx


def test_coefficients():
    assert fd.ka_rankine(30) == approx(1 / 3)
    assert fd.ka_rankine(30, 15) == approx(0.373, abs=5e-4)     # Das, Table (phi 30, slope 15)
    assert fd.ka_coulomb(30, 0, 0) == approx(1 / 3)
    assert fd.ka_coulomb(30, 20, 0) == approx(0.297, abs=5e-4)  # Das, Coulomb table
    assert fd.kp_rankine(30) == approx(3.0)
    assert fd.kae_mononobe_okabe(30, 0, 0, 0.0) == approx(1 / 3)
    th = math.atan(0.2)
    num = math.cos(math.radians(30) - th) ** 2
    den = math.cos(th) ** 2 * (1 + math.sqrt(0.5 * math.sin(math.radians(30) - th) / math.cos(th))) ** 2
    assert fd.kae_mononobe_okabe(30, 0, 0, 0.2) == approx(num / den)


def _ck(r, name):
    return next(c for c in r.checks if c.name.startswith(name))


HAND = dict(H_stem=10, H_retained=10, t_top=12, t_bot=12, toe=2, heel=4, h_ftg=12,
            gamma=120, phi=30, mu=0.5, passive_efp=0, include_toe_soil=False)


def test_stability_hand_calc():
    # Pa = 0.5(1/3)(0.12)(11^2) = 2.42 k at 11/3; stem 1.5 k @ 2.5; footing 1.05 @ 3.5;
    # soil 4.8 @ 5.0 -> Mr = 31.425, Mo = 8.873
    r = rw.run(**HAND)["report"]
    Mo = 2.42 * 11 / 3
    assert _ck(r, "Overturning FS (Static)").capacity == approx(31.425 / Mo, rel=1e-4)
    assert _ck(r, "Sliding FS (Static)").capacity == approx(0.5 * 7.35 / 2.42, rel=1e-4)
    xbar = (31.425 - Mo) / 7.35
    e = 3.5 - xbar
    assert _ck(r, "Soil bearing (Static)").demand == approx(7.35 / 7 * (1 + 6 * e / 7), rel=1e-4)


def test_stem_and_heel_hand_calc():
    r = rw.run(**HAND)["report"]
    assert _ck(r, "Stem flexure at base").demand == approx(1.6 * (1 / 3) * 0.12 * 1000 / 6, rel=1e-4)
    # 0.9D + 1.6H: partial contact (e > B/6); heel moment by hand = 7.618 k-ft/ft
    N = 0.9 * 7.35
    xbar = (0.9 * 31.425 - 1.6 * 2.42 * 11 / 3) / N
    qt = 2 * N / (3 * xbar)
    c = 3 * xbar
    Mq = qt / c * (c - 3) ** 3 / 6           # integral of the triangular pressure over the heel
    Mh = 0.9 * (0.15 + 1.2) * 4**2 / 2 - Mq
    assert _ck(r, "Heel flexure").demand == approx(Mh, rel=2e-3)


def test_seed_whitman_increment():
    o = rw.run(**HAND, seismic="Seed-Whitman", kh=0.2, wall_inertia=False)
    f = [x for x in o["forces"] if x[0] == "Seismic earth increment"][0]
    assert f[3] == approx(3 / 8 * 0.2 * 120 * 11**2 / 1000)
    assert f[5] == approx(0.6 * 11)


def test_water_and_coulomb_run():
    o = rw.run(**HAND, h_water=4, method="Coulomb", delta=15, beta=10, surcharge=100, strip_q=400)
    names = [x[0] for x in o["forces"]]
    assert "Water uplift" in names and "Active earth (vertical comp.)" in names
    w = [x for x in o["forces"] if x[0] == "Water pressure"][0]
    assert w[3] == approx(0.5 * 62.4 * 16 / 1000, rel=1e-3)
