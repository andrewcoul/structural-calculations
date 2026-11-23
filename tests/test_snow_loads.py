"""Hand checks for ASCE 7-22 snow loads (published worked examples)."""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))

import pytest  # noqa: E402

from structcalc.calcs import snow_loads as s  # noqa: E402

approx = pytest.approx


def test_ct_table_733():
    assert s.Ct_heated(30, 30) == approx(1.14)
    assert s.Ct_heated(38, 60) == approx(1.142, abs=0.001)   # Baltimore example: 1.14
    assert s.Ct_heated(15, 80) == approx(1.00)
    assert s.Ct_heated(60, 30) == approx(1.2)


def test_kansas_city_roof_step():
    # pg = 30, W2 = 0.55, Ce = 1.0, R-30 -> Ct = 1.14; upper 100 ft, lower 40 ft, step 10 ft
    o = s.run(pg=30, W2=0.55, terrain="C", exposure="Partially exposed", R_value=30,
              roof_type="Flat", theta=0, upper_len=100, lower_len=40, step_h=10)
    assert o["pf"] == approx(23.94, abs=0.01)
    assert o["pm"] == approx(30.0)
    assert o["gamma"] == approx(17.9)
    d = o["drift"]
    assert d["hc"] == approx(8.66, abs=0.01)
    assert d["leeward"]["h"] == approx(3.763, abs=0.002)
    assert d["leeward"]["w"] == approx(15.05, abs=0.02)
    assert d["leeward"]["pd"] == approx(67.4, abs=0.1)
    assert d["windward"]["h"] == approx(2.048, abs=0.002)
    assert d["windward"]["w"] == approx(16.4, abs=0.02)
    assert d["windward"]["pd"] == approx(36.7, abs=0.1)


def test_baltimore_drift():
    o = s.run(pg=60, W2=0.45, R_value=38, roof_type="Flat", theta=0, upper_len=80, lower_len=80,
              step_h=6)
    assert o["pf"] == approx(47.9, abs=0.1)
    assert o["drift"]["leeward"]["h"] == approx(3.44, abs=0.01)
    assert o["drift"]["windward"]["h"] == approx(2.58, abs=0.01)


def test_slope_factor_curves():
    # R-38 metal roof at 6:12, pg = 60 -> Ct = 1.14, Fig. 7.4-1(b) slippery: Cs = 0.724
    cs, curve = s.Cs_slope(26.565, 1.14, True)
    assert curve == "(b)" and cs == approx(0.724, abs=0.001)
    assert s.Cs_slope(26.565, 1.10, True) == (1.0, "(a)")
    assert s.Cs_slope(50, 1.2, False)[0] == approx(1 - 5 / 25)


def test_truncated_leeward_drift():
    d = s.step_drift(50, 0.6, 30, 200, 60, 4)
    lw = d["leeward"]
    assert lw["truncated"]
    assert lw["h"] == approx(d["hc"])
    assert lw["w"] == approx(min(4 * lw["h_full"] ** 2 / d["hc"], 8 * d["hc"]))


def test_unbalanced_gable_and_rain_on_snow():
    o = s.run(pg=20, W2=0.5, R_value=20, roof_type="Gable", theta=18.43, eave_to_ridge=30,
              drift=False)
    u = o["unbalanced"]
    assert u["windward"] == approx(0.3 * o["ps"])
    assert u["surcharge"] == approx(s.hd_eq(20, 30, 0.5) * s.gamma(20) / 3**0.5, rel=1e-3)
    # rain-on-snow: pg <= 30 and slope 18.4 > W/50 = 0.6 deg -> not applied
    assert "Rain-on-snow surcharge" not in [i.name for i in o["report"].results]
    o2 = s.run(pg=20, W2=0.5, roof_type="Flat", theta=0.5, eave_to_ridge=50, drift=False)
    assert "Rain-on-snow surcharge" in [i.name for i in o2["report"].results]
