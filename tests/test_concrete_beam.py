"""Hand-check tests: ACI 318-19 concrete beam (module 05) and section helpers."""

import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))

import pytest  # noqa: E402

from structcalc import concrete as C  # noqa: E402
from structcalc import concrete_member as cm  # noqa: E402
from structcalc.calcs import concrete_beam  # noqa: E402

approx = pytest.approx


def test_doubly_reinforced_hand():
    """b=12, h=24, As=4.0 in2 at d=21.5, As'=1.2 in2 at d'=2.5, f'c=4 ksi, fy=60 ksi.
    Compression steel does not yield; solve c from the quadratic (ksi units):
    0.85 fc b beta1 c + As'(87(c-d')/c - 0.85fc) = As fy  (Es*0.003 = 87 ksi)."""
    fc, fy, b, d, dp, As, Asp = 4.0, 60.0, 12.0, 21.5, 2.5, 4.0, 1.2
    k1 = 0.85 * fc * b * 0.85                       # 34.68
    # k1 c^2 + (Asp*(87 - 3.4) - As*fy) c - 87*Asp*dp = 0
    A, B, Cq = k1, Asp * (87 - 0.85 * fc) - As * fy, -87 * Asp * dp
    c = (-B + math.sqrt(B * B - 4 * A * Cq)) / (2 * A)
    fsp = 87 * (c - dp) / c
    assert fsp < fy
    a = 0.85 * c
    Mn = 0.85 * fc * b * a * (d - a / 2) + Asp * (fsp - 0.85 * fc) * (d - dp)  # kip-in
    eps_t = 0.003 * (d - c) / c
    sec = cm.Section(cm.rect(b, 24), [(d, As), (dp, Asp)], 4000, 60000)
    f = sec.flexure()
    assert f["c"] == approx(c, rel=1e-6)
    assert f["Mn"] == approx(Mn, rel=1e-6)
    assert f["eps_t"] == approx(eps_t, rel=1e-6)
    assert f["phi"] == 0.9  # eps_t = 0.0087 > 0.005


def test_singly_matches_core_helper():
    sec = cm.Section(cm.rect(12, 24), [(21.5, 2.37)], 4000, 60000)
    f = sec.flexure()
    ref = C.flexure_rect(2.37, 12, 21.5, 4000, 60000)
    assert f["Mn"] == approx(ref["Mn"], rel=1e-6)


def test_tbeam_flange_block():
    sec = cm.Section(cm.tee(12, 24, 48, 4), [(21.5, 2.37)], 4000, 60000)
    a = 2.37 * 60 / (0.85 * 4 * 48)
    assert sec.flexure()["Mn"] == approx(2.37 * 60 * (21.5 - a / 2), rel=1e-6)


def test_icr_rectangular_closed_form():
    Ec = 57000 * math.sqrt(4000)
    n = C.ES / Ec
    b, d, As = 12, 21.5, 2.37
    rho = As / (b * d)
    k = math.sqrt((rho * n) ** 2 + 2 * rho * n) - rho * n
    Icr = b * (k * d) ** 3 / 3 + n * As * (d - k * d) ** 2
    kd, I = cm.Section(cm.rect(b, 24), [(d, As)], 4000, 60000).cracked(Ec)
    assert kd == approx(k * d, rel=1e-6)
    assert I == approx(Icr, rel=1e-6)


def test_Ie_aci_318_19():
    Ig, Icr, Mcr = 13824, 5000, 45.0
    assert cm.Ie_aci(0.6 * Mcr, Mcr, Ig, Icr) == Ig            # Ma <= 2/3 Mcr
    Ma = 120.0
    expect = Icr / (1 - ((2 / 3) * Mcr / Ma) ** 2 * (1 - Icr / Ig))
    assert cm.Ie_aci(Ma, Mcr, Ig, Icr) == approx(expect)


def _check(r, name):
    return next(c for c in r.checks if c.name == name)


def test_module_simple_span_flexure_shear():
    o = concrete_beam.run(w={"D": 1.0, "L": 1.2}, span=20)
    r = o["report"]
    wsw = 0.150 * 12 * 24 / 144                  # 0.30 kip/ft
    Mu = (1.2 * (1.0 + wsw) + 1.6 * 1.2) * 20**2 / 8
    assert _check(r, "Flexure +M").demand == approx(Mu, rel=1e-6)
    d = 24 - 1.5 - 0.375 - 0.5
    Vc = 2 * math.sqrt(4000) * 12 * d / 1000     # Av >= Av,min -> Eq. (a)
    Vs = 0.22 * 60 * d / 10
    assert _check(r, "Shear").capacity == approx(0.75 * (Vc + Vs), rel=1e-6)
    assert _check(r, "Shear").demand == approx(Mu * 8 / 20**2 * 20 / 2, rel=1e-6)
    assert r.ok


def test_module_no_stirrups_size_effect():
    o = concrete_beam.run(w={"D": 0.3, "L": 0.3}, span=16, s_stirrup=0)
    r = o["report"]
    d = 24 - 1.5 - 0.5  # no stirrup
    rho = 2.37 / (12 * d)
    ls = math.sqrt(2 / (1 + d / 10))
    Vc = 8 * ls * rho ** (1 / 3) * math.sqrt(4000) * 12 * d / 1000
    assert _check(r, "Shear").capacity == approx(0.75 * Vc, rel=1e-6)


def test_module_gross_deflection():
    o = concrete_beam.run(w={"D": 0.0, "L": 1.0}, span=24, self_weight=False,
                          deflection_method="Gross Ig")
    r = o["report"]
    Ec = C.Ec(4000) / 1000
    Ig = 12 * 24**3 / 12
    dL = 5 * (1.0 / 12) * (24 * 12) ** 4 / (384 * Ec * Ig)
    assert _check(r, "Deflection immediate live - Span 1").demand == approx(dL, rel=1e-5)


def test_module_negative_moment_needs_top_bars():
    o = concrete_beam.run(span_type="Cantilever", span=8, w={"D": 1.0, "L": 0.5}, n_top=0)
    r = o["report"]
    assert not _check(r, "Flexure -M").ok
    o = concrete_beam.run(span_type="Cantilever", span=8, w={"D": 1.0, "L": 0.5},
                          n_top=4, bar_top=7)
    neg = _check(o["report"], "Flexure -M")
    Mu = (1.2 * (1.0 + 0.3) + 1.6 * 0.5) * 8**2 / 2
    assert neg.demand == approx(Mu, rel=1e-6)
    assert neg.ok


def test_crack_control_formula():
    o = concrete_beam.run(w={"D": 1.0, "L": 1.2})
    c = _check(o["report"], "Crack control spacing bottom")
    fs = 40000
    cc = 1.5 + 0.375
    assert c.capacity == approx(min(15 * 40000 / fs - 2.5 * cc, 12 * 40000 / fs))
