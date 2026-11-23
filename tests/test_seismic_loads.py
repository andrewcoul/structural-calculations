"""Hand checks for ASCE 7-22 ELF (12.8)."""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))

import pytest  # noqa: E402

from structcalc.calcs import seismic_loads as s  # noqa: E402

approx = pytest.approx


def test_short_period_cs_12_8_2():
    # Two-story wood house example: SDS=1.58, SD1=1.10, R=6.5, Ie=1 -> Cs = 0.2431
    o = s.run(SDS=1.58, SD1=1.10, S1=0.72, R=6.5, story_heights="10;10",
              story_weights="50;31.1")
    assert o["Cs"] == approx(1.58 / 6.5, rel=1e-9)
    assert o["V"] == approx(1.58 / 6.5 * 81.1, rel=1e-9)


def test_period_and_cs_max():
    # 8-story steel MRF, hn = 105 ft, SDS=1.0, SD1=0.6, R=8, RC II
    o = s.run(SDS=1.0, SD1=0.6, S1=0.6, R=8, system="Steel moment-resisting frame",
              story_heights="14;" + ";".join(["13"] * 7), story_weights=";".join(["1000"] * 8))
    Ta = 0.028 * 105**0.8
    assert o["T"] == approx(Ta, rel=1e-9)
    Cs_max = 0.6 / (Ta * 8)
    assert o["Cs"] == approx(max(Cs_max, 0.044 * 1.0, 0.5 * 0.6 / 8), rel=1e-9)
    assert o["k"] == approx(1 + (Ta - 0.5) / 2, rel=1e-9)


def test_cu_cap_and_vertical_distribution():
    o = s.run(SDS=0.5, SD1=0.2, S1=0.2, R=3, T_computed=2.0, story_heights="12;12;12",
              story_weights="100;100;80")
    Ta = 0.02 * 36**0.75
    assert o["T"] == approx(1.5 * Ta)  # Cu = 1.5 at SD1 = 0.2
    F = o["F"]
    assert F.sum() == approx(o["V"])
    k = o["k"]
    w = [100 * 12**k, 100 * 24**k, 80 * 36**k]
    assert F[2] / o["V"] == approx(w[2] / sum(w))


def test_sdc_and_minimums():
    assert s.sdc(0.4, 0.1, 0.1, "II")[0] == "C"
    assert s.sdc(0.2, 0.25, 0.2, "II")[0] == "D"
    assert s.sdc(0.2, 0.1, 0.8, "IV")[0] == "F"
    o = s.run(SDS=0.1, SD1=0.05, S1=0.04, R=8, story_heights="200", story_weights="1000",
              system="Steel moment-resisting frame")
    assert o["Cs"] == approx(0.01)


def test_12813_cap():
    o = s.run(SDS=1.6, SD1=1.0, S1=0.7, R=6.5, use_12813=True, story_heights="10",
              story_weights="100")
    assert o["Cs"] == approx(max(1.0, 0.7 * 1.6) / 6.5)
