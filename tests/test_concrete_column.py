"""Hand-check tests: ACI 318-19 concrete column (module 06).

Reference column: 16 x 16 in, 8 #9 (3 bars per face), f'c = 4 ksi, fy = 60 ksi,
1.5 in clear cover to #3 ties -> d' = 1.5 + 0.375 + 1.128/2 = 2.439 in.
"""

import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))

import pytest  # noqa: E402

from structcalc import concrete as C  # noqa: E402
from structcalc import concrete_member as cm  # noqa: E402
from structcalc.calcs import concrete_column as cc  # noqa: E402

approx = pytest.approx
DP = 1.5 + 0.375 + 1.128 / 2
FC, FY, ES = 4.0, 60.0, 29000.0


def ref_section():
    return cm.Section(cm.rect(16, 16), [(DP, 3.0), (8.0, 2.0), (16 - DP, 3.0)], 4000, 60000)


def test_pure_axial():
    Po = 0.85 * FC * (256 - 8) + FY * 8        # 1323.2 kip
    s = ref_section()
    assert s.Po == approx(Po)
    assert s.interaction()["phiPn_max"] == approx(0.65 * 0.80 * Po)


def test_balanced_point():
    d = 16 - DP
    c = 0.003 / (0.003 + FY / ES) * d
    a = 0.85 * c
    Cc = 0.85 * FC * 16 * a
    eps_top = 0.003 * (c - DP) / c
    fs_top = min(ES * eps_top, FY)
    eps_mid = 0.003 * (c - 8) / c
    fs_mid = ES * eps_mid                       # ~0, elastic
    Ftop = 3 * (fs_top - 0.85 * FC)             # inside block: displaced concrete
    Fmid = 2 * fs_mid
    Fbot = -3 * FY
    P = Cc + Ftop + Fmid + Fbot
    M = Cc * (8 - a / 2) + Ftop * (8 - DP) - Fbot * (16 - DP - 8)
    st = ref_section().state(c)
    assert st["Pn"] == approx(P, rel=1e-6)      # ~361.5 kip
    assert st["Mn"] == approx(M, rel=1e-6)      # ~304.0 kip-ft
    assert st["phi"] == 0.65                    # eps_t = eps_ty -> compression controlled
    assert P == approx(361.5, abs=0.1) and M / 12 == approx(304.0, abs=0.1)


def test_pure_bending():
    # top bars elastic, middle and bottom yield in tension:
    # 46.24 c + 3(87(c - d')/c - 3.4) - 300 = 0
    k = 0.85 * FC * 16 * 0.85
    A, B, Cq = k, 3 * (87 - 3.4) - 5 * FY, -3 * 87 * DP
    c = (-B + math.sqrt(B * B - 4 * A * Cq)) / (2 * A)
    a = 0.85 * c
    fs_top = 87 * (c - DP) / c
    M = k * c * (8 - a / 2) + 3 * (fs_top - 3.4) * (8 - DP) + 3 * FY * (16 - DP - 8)
    f = ref_section().flexure()
    assert f["c"] == approx(c, rel=1e-6)
    assert f["Mn"] == approx(M, rel=1e-6)       # ~232.7 kip-ft
    assert f["phi"] == 0.9


def test_ray_ratio_on_curve():
    d = ref_section().interaction()
    k = len(d["phiMn"]) // 2
    r, _, _ = cm.ray_ratio(d["phiMn"][k], d["phiPn"][k], d["phiMn"], d["phiPn"])
    assert r == approx(1.0, abs=1e-9)
    r2, _, _ = cm.ray_ratio(0.5 * d["phiMn"][k], 0.5 * d["phiPn"][k], d["phiMn"], d["phiPn"])
    assert r2 == approx(0.5, abs=1e-9)


def test_module_layout_and_limits():
    o = cc.run(case_loads={"D": (200, 30, -15, 0, 0)})
    assert len(o["bars"]) == 8
    r = o["report"]
    rho = next(c for c in r.checks if c.name == "rho_g >= 0.01")
    assert rho.capacity == approx(8 / 256)


def test_circular_slenderness_hand():
    """D = 20 in, 8 #9, lu = 20 ft, k = 1, Pu = 1.2*300 + 1.6*150 = 600 kip,
    M2 = 28 k-ft < M2,min = 600(0.6 + 0.03*20)/12 = 60 k-ft -> Cm = 1.0."""
    o = cc.run(shape="Circular spiral", D=20, n_circ=8, s_tie=2.0, lu=20,
               case_loads={"D": (300, 10, 10, 0, 0), "L": (150, 10, 10, 0, 0)})
    row = next(rw for rw in o["rows"] if rw["Combination"].startswith("2:"))
    Ec = C.Ec(4000) / 1000
    Ig = math.pi * 20**4 / 64
    Rb = 10 - 1.5 - 0.375 - 1.128 / 2
    Ise = 8 * 1.0 * Rb**2 / 2
    beta = 1.2 * 300 / 600
    EI = (0.2 * Ec * Ig + 29000 * Ise) / (1 + beta)
    Pc = math.pi**2 * EI / 240**2
    delta = 1.0 / (1 - 600 / (0.75 * Pc))
    assert row["delta x"] == approx(delta, rel=1e-9)
    assert row["Mux (k-ft)"] == approx(delta * 60, rel=1e-9)


def test_short_column_not_magnified():
    o = cc.run(case_loads={"D": (100, 20, 20, 0, 0)}, lu=8)   # klu/r = 96/4.8 = 20 < 22
    row = o["rows"][0]
    assert row["delta x"] == 1.0
    assert row["Mux (k-ft)"] == approx(1.4 * 20)


def test_overload_fails():
    o = cc.run(case_loads={"D": (600, 50, 50, 0, 0), "L": (300, 0, 0, 0, 0)})
    assert not o["report"].ok
