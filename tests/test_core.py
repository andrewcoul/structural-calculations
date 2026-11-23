"""Hand-check tests for the shared engine and the steel beam module.

Run: .venv/bin/python -m pytest -q tests
"""

import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))

import pytest  # noqa: E402

from structcalc import concrete, loads, steel  # noqa: E402
from structcalc.beam import Beam  # noqa: E402

approx = pytest.approx


# --- beam solver vs closed form -------------------------------------------
def test_simple_span_udl():
    b = Beam(24, 29000, 510, [(0, "pin"), (24, "pin")]).add_udl(1.0, "D")
    r = b.solve()["D"]
    assert r.Mmax == approx(72.0, rel=1e-6)
    assert r.dmax == approx(5 * 24**4 * 1728 / (384 * 29000 * 510), rel=1e-6)


def test_cantilever_point():
    r = Beam(10, 29000, 100, [(0, "fixed")]).add_point(10, 5, "D").solve()["D"]
    assert r.Mmin == approx(-50, rel=1e-6)
    assert r.dmax == approx(5 * 1000 * 1728 / (3 * 29000 * 100), rel=1e-6)


def test_two_span_continuous():
    r = Beam(40, 29000, 500, [(0, "pin"), (20, "pin"), (40, "pin")]).add_udl(1, "D").solve()["D"]
    assert r.Mmin == approx(-50, rel=1e-6)
    assert r.reactions[20.0][0] == approx(25, rel=1e-6)


def test_fixed_fixed():
    r = Beam(12, 29000, 100, [(0, "fixed"), (12, "fixed")]).add_udl(2, "D").solve()["D"]
    assert r.Mmin == approx(-24, rel=1e-6)
    assert r.Mmax == approx(12, rel=1e-6)


def test_cases_add():
    b = Beam(20, 29000, 200, [(0, "pin"), (20, "pin")])
    b.add_udl(1, "D").add_point(5, 3, "L")
    res = b.combos([("x", {"D": 1.2, "L": 1.6})])
    assert res[0][2].reactions[0.0][0] == approx(1.2 * 10 + 1.6 * 3 * 15 / 20, rel=1e-6)


# --- load combinations ----------------------------------------------------
def test_asce722_snow_factors():
    names = {n: f for n, f in loads.lrfd()}
    assert any(f.get("S") == 1.0 and f.get("D") == 1.2 for f in names.values())
    assert any(f.get("S") == 0.3 and f.get("L") == 1.6 for f in names.values())
    assert any(f.get("S") == 0.7 and len(f) == 2 for _, f in loads.asd())


# --- AISC 360-22 vs Manual tables -----------------------------------------
def test_W18X50_Lp_Lr_phiMp():
    f = steel.flexure_W_major(steel.shape("W18X50"), 50, 0)
    assert f.Lp / 12 == approx(5.83, abs=0.01)
    assert f.Lr / 12 == approx(16.9, abs=0.05)
    assert 0.9 * f.Mn / 12 == approx(379, abs=0.5)


def test_W14X90_noncompact_flange():
    f = steel.flexure_W_major(steel.shape("W14X90"), 50, 0)
    assert 0.9 * f.Mn / 12 == approx(574, abs=1)  # AISC Table 3-2


def test_W_shear():
    assert steel.shear(steel.shape("W18X50"), 50)[0] == approx(192, abs=0.5)  # Table 3-2
    assert steel.shear(steel.shape("W12X26"), 50)[0] == approx(84.2, abs=0.2)


def test_W14X90_compression():
    p = steel.shape("W14X90")
    assert 0.9 * steel.compression(p, 50, 1e-3, 1e-3).Pn == approx(1190, abs=3)  # Table 4-1a, Lc=0


def test_HSS_flexure_compact():
    f = steel.flexure(steel.shape("HSS6X6X1/4"), 50)
    assert 0.9 * f.Mn / 12 == approx(42.0, abs=0.1)


def test_shape_names():
    assert steel.shape("hss6.625x0.280")["OD"] == approx(6.625, abs=0.01)
    with pytest.raises(KeyError):
        steel.shape("W99X1")


# --- ACI 318-19 core --------------------------------------------------------
def test_beta1():
    assert concrete.beta1(4000) == 0.85
    assert concrete.beta1(5000) == approx(0.80)
    assert concrete.beta1(9000) == 0.65


def test_flexure_rect():
    r = concrete.flexure_rect(2.37, 12, 21.5, 4000, 60000)
    assert r["a"] == approx(3.485, abs=0.002)
    assert r["Mn"] == approx(2.37 * 60 * (21.5 - 3.485 / 2), rel=1e-3)
    assert r["phi"] == 0.9
    As = concrete.As_required(r["phiMn"], 12, 21.5, 4000, 60000)
    assert As == approx(2.37, rel=1e-6)


def test_vc():
    assert concrete.Vc_oneway(12, 20, 4000, Av_ge_min=True) == approx(2 * math.sqrt(4000) * 240 / 1000)


# --- steel beam module ----------------------------------------------------
def test_steel_beam_module():
    from structcalc.calcs import steel_beam

    o = steel_beam.run("W18X35", w={"D": 0.4, "L": 0.8}, extra_loads="P, 12, 5, L", span=24,
                       bracing="Continuous (both flanges)")
    r = o["report"]
    flex = next(c for c in r.checks if c.name == "Flexure")
    Mu = 1.2 * (0.4 + 0.035) * 24**2 / 8 + 1.6 * (0.8 * 24**2 / 8 + 5 * 24 / 4)
    assert flex.demand == approx(Mu, rel=1e-4)
    assert flex.capacity == approx(0.9 * 50 * 66.5 / 12, rel=1e-6)


def test_steel_beam_autoselect():
    from structcalc.calcs import steel_beam

    o = steel_beam.run(auto_select=True, depth_filter="16", w={"D": 0.4, "L": 0.8}, span=24,
                       bracing="Top flange continuous")
    assert o["report"].ok
    assert o["section"]["name"].startswith("W16X")


def test_load_schematic_plots():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from structcalc.beam import plot as beam_plot
    from structcalc.calcs import steel_beam, steel_column, wood_column, concrete_column

    sb = steel_beam.run("W18X35", span_type="Overhang one end (right)", span=18, overhang_right=6,
                        w={"D": 0.4, "L": 0.8}, extra_loads="P, 12, 5, L; t, 0, 10, 0.2, 0.6, S; M, 6, 15, W")
    _, _, res = sb["governing"]["flex"]
    fig = beam_plot(res, title="Beam loads test", beam=sb["beam"])
    assert len(fig.axes) == 4
    plt.close(fig)

    sc = steel_column.run(section="W10X49", H=14, Lcx=14, Lcy=14, Lcz=14, Lb=14,
                          case_loads={"D": {"P": 60, "Mx_top": 20, "Mx_bot": -10}, "L": {"P": 80}},
                          ex=6.0, wx={"W": 0.2})
    fig = steel_column.plot(sc)
    assert len(fig.axes) == 3
    plt.close(fig)

    wc = wood_column.run(size="2x6", P={"D": 0.4, "L": 0.8}, wind_psf=25.0, trib_width=16.0)
    fig = wood_column.plot(wc)
    assert len(fig.axes) == 2
    plt.close(fig)

    cc = concrete_column.run(case_loads={"D": (200, 30, -15, 10, 10), "L": (120, 20, -10, 5, 5)})
    fig = concrete_column.plot(cc)
    assert len(fig.axes) == 4
    plt.close(fig)

    from structcalc.calcs import (
        bolt_connection, combined_footing, concrete_beam, general_footing,
        nail_connection, retaining_wall, seismic_loads, snow_loads,
        wind_loads, wood_shear_wall,
    )
    for mod, kw in [
        (concrete_beam, dict(section_type="T-beam", w={"D": 1.0, "L": 1.2})),
        (concrete_column, dict(shape="Circular spiral", case_loads={"D": (200, 20, -10, 0, 0)})),
        (general_footing, dict(col_loads={"D": (80, 0, 0, 0, 0), "W": (0, 0, 40, 5, 0)})),
        (combined_footing, dict(col1={"D": (80, 0)}, col2={"D": (120, 0), "W": (0, 60)})),
        (retaining_wall, dict(surcharge=100, strip_q=200, P_D=2.0)),
        (wood_shear_wall, dict()),
        (bolt_connection, dict(n_rows=2, n_per_row=3, theta_main=45.0, load=12000)),
        (nail_connection, dict(install="Toe-nail", lateral_load=200, withdrawal_load=50)),
        (wind_loads, dict()),
        (seismic_loads, dict()),
        (snow_loads, dict(roof_type="Gable", theta=20.0, drift=False)),
    ]:
        fig = mod.plot(mod.run(**kw))
        assert len(fig.axes) >= 1
        plt.close(fig)
