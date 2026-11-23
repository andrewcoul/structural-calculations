"""Hand checks for NDS nail / wood screw connections."""
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))

import pytest  # noqa: E402

from structcalc import wood_conn as wc  # noqa: E402
from structcalc.calcs import nail_connection  # noqa: E402

approx = pytest.approx


def test_16d_common_table_12N():
    """NDS Table 12N: 16d common (D=0.162), 1-1/2 in side, G=0.50 -> Z = 141 lb."""
    o = nail_connection.run(size="16d", t_side=1.5, t_main=3.5)
    assert o["mode"] == "IV"
    assert round(o["Z"]) == 141


def test_withdrawal_16d():
    """NDS Table 12.2C: 16d common, G=0.50 -> 40 lb/in."""
    assert wc.withdrawal_nail(0.5, 0.162) == approx(39.5, abs=0.1)
    assert round(wc.withdrawal_nail(0.5, 0.162)) == 40


def test_combined_nail():
    Zp, Wp = 140.0, 80.0
    a = 30.0
    exp = Wp * Zp / (Wp * math.cos(math.radians(a)) + Zp * math.sin(math.radians(a)))
    assert wc.combined_lateral_withdrawal(Zp, Wp, a, "nail") == approx(exp)
    assert wc.combined_lateral_withdrawal(Zp, Wp, 0, "nail") == approx(Zp)
    assert wc.combined_lateral_withdrawal(Zp, Wp, 90, "screw") == approx(Wp)


def test_toe_nail():
    o = nail_connection.run(size="16d", install="Toe-nail")
    assert o["ls"] == approx(3.5 / 3)
    assert o["p"] == approx(3.5 * math.cos(math.radians(30)) - 3.5 / 3)
    assert o["Zp"] == approx(o["Z"] * 0.83)


def test_end_grain_withdrawal_ng():
    o = nail_connection.run(size="16d", install="End grain", withdrawal_load=50)
    assert not o["report"].ok
    o = nail_connection.run(size="16d", install="End grain", lateral_load=100)
    assert o["Zp"] == approx(o["Z"] * 0.67)


def test_short_penetration_ng():
    o = nail_connection.run(size="6d", t_side=1.5, t_main=3.5, lateral_load=10, n_fasteners=1)
    # p = 0.5 in < 6D = 0.678 in
    assert not o["report"].ok


def test_n_required():
    o = nail_connection.run(size="16d", lateral_load=1000, n_fasteners=4)
    assert o["n_required"] == math.ceil(1000 / o["Zp"])


def test_wood_screw():
    o = nail_connection.run(fastener="Wood screw", size="#10", screw_length=3.0,
                            withdrawal_load=200, n_fasteners=2)
    assert o["Wp"] == approx(2850 * 0.25 * 0.190 * 1.5)
    with pytest.raises(ValueError):
        nail_connection.run(fastener="Wood screw", size="16d")
