"""Hand checks for the combined footing (08)."""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))

import pytest  # noqa: E402

from structcalc.calcs import combined_footing as cf  # noqa: E402

approx = pytest.approx


def _run(**kw):
    base = dict(L=20, B=6, h=30, x1=2, c1x=16, c1y=16, x2=16, c2x=18, c2y=18,
                col1={"D": (100, 0), "L": (60, 0)}, col2={"D": (150, 0), "L": (90, 0)})
    base.update(kw)
    return cf.run(**base)


def test_negative_moment_hand_calc():
    # 1.2D+1.6L: P1=216, P2=324; net pressure linear 23.76 -> 30.24 k/ft (self weight cancels)
    # V=0 at 0.162x^2 + 23.76x - 216 = 0; M = 11.88x^2 + 0.054x^3 - 216(x-2)
    o = _run()
    import math
    x = (-23.76 + math.sqrt(23.76**2 + 4 * 0.162 * 216)) / (2 * 0.162)
    M = 11.88 * x**2 + 0.054 * x**3 - 216 * (x - 2)
    ck = next(c for c in o["report"].checks if c.name == "Flexure, top longitudinal")
    assert ck.demand == approx(-M, rel=1e-4)


def test_closure_and_bearing():
    o = _run()
    for name, x, V, M in o["diagrams"]:
        assert abs(M[-1]) < 1e-3 and abs(V[-1]) < 1e-6
    # D+L: R = 400 + 58.2; e = (160*(2-10) + 240*(16-10)) / R
    R = 400 + 58.2
    e = (160 * -8 + 240 * 6) / R
    ck = next(c for c in o["report"].checks if c.name == "Soil bearing")
    assert ck.demand == approx(R / 20 * (1 + 6 * e / 20) / 6, rel=1e-9)


def test_column_moment_jump_and_edge_punching():
    o = _run(x1=0.75, col2={"D": (150, 40), "L": (90, 0)})
    name, x, V, M = o["diagrams"][0]
    i = (abs(x - 16)).argmin()
    assert M[i + 1] - M[i - 1] == approx(1.4 * 40, abs=0.5)  # 1.4D combo, jump 56 k-ft
    res = [i for i in o["report"].results if i.name.startswith("Column 1 punching")][0]
    assert ", 3, 30" in res.value  # 3-sided perimeter at the footing end


def test_transverse_strip_moment():
    o = _run()
    Pu = 1.2 * 150 + 1.6 * 90
    ck = next(c for c in o["report"].checks if c.name == "Transverse flexure, strip 2")
    assert ck.demand == approx(Pu / 6 * ((6 - 1.5) / 2) ** 2 / 2, rel=1e-9)
