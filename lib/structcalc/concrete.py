"""ACI 318-19 core helpers shared by concrete members and foundations.

UNITS (read carefully):
    fc, fy, fyt  -> psi            (e.g. fc=4000, fy=60000)
    dimensions   -> in
    forces       -> kip            moments -> kip-in
    Functions that take sqrt(fc) use fc in psi and convert to kip.
"""

from __future__ import annotations

import math

ES = 29_000_000.0  # psi, reinforcement modulus (ACI 20.2.2.2)
EPS_CU = 0.003     # ACI 22.2.2.1

# Bar: (diameter in, area in^2) - ASTM A615 / ACI Appendix
BARS = {
    3: (0.375, 0.11), 4: (0.500, 0.20), 5: (0.625, 0.31), 6: (0.750, 0.44),
    7: (0.875, 0.60), 8: (1.000, 0.79), 9: (1.128, 1.00), 10: (1.270, 1.27),
    11: (1.410, 1.56), 14: (1.693, 2.25), 18: (2.257, 4.00),
}


def bar(n):
    """Return (diameter in, area in^2) of bar size #n."""
    return BARS[int(n)]


def Ec(fc, wc=150.0):
    """ACI 19.2.2.1: modulus of elasticity, psi."""
    return wc**1.5 * 33 * math.sqrt(fc)


def fr(fc, lam=1.0):
    """ACI 19.2.3.1: modulus of rupture, psi."""
    return 7.5 * lam * math.sqrt(fc)


def beta1(fc):
    """ACI Table 22.2.2.4.3."""
    if fc <= 4000:
        return 0.85
    return max(0.85 - 0.05 * (fc - 4000) / 1000, 0.65)


def eps_ty(fy):
    return fy / ES


def phi_flexure_axial(eps_t, fy=60000, transverse="tied"):
    """ACI Table 21.2.2 strength reduction factor from net tensile strain."""
    ety = eps_ty(fy)
    phc = 0.75 if transverse == "spiral" else 0.65
    if eps_t <= ety:
        return phc
    if eps_t >= ety + 0.003:
        return 0.90
    return phc + (0.90 - phc) * (eps_t - ety) / 0.003


PHI_SHEAR = 0.75   # ACI 21.2.1
PHI_BEARING = 0.65


def flexure_rect(As, b, d, fc, fy, dt=None):
    """Singly reinforced rectangular section, ACI 22.2.

    Returns dict(Mn kip-in, a, c, eps_t, phi, phiMn kip-in, tension_controlled).
    """
    dt = d if dt is None else dt
    a = As * fy / (0.85 * fc * b)
    c = a / beta1(fc)
    eps_t = EPS_CU * (dt - c) / c if c > 0 else math.inf
    Mn = As * fy * (d - a / 2) / 1000.0  # kip-in
    phi = phi_flexure_axial(eps_t, fy)
    return {"Mn": Mn, "a": a, "c": c, "eps_t": eps_t, "phi": phi, "phiMn": phi * Mn,
            "tension_controlled": eps_t >= eps_ty(fy) + 0.003}


def As_required(Mu, b, d, fc, fy, phi=0.9):
    """Required tension steel (in^2) for factored moment Mu (kip-in).
    Returns math.inf when the section is too small (no real root)."""
    Mu_lb = abs(Mu) * 1000.0
    k = 0.85 * fc * b
    disc = d**2 - 2 * Mu_lb / (phi * k)
    if disc < 0:
        return math.inf
    a = d - math.sqrt(disc)
    return k * a / fy


def As_min_beam(b, d, fc, fy):
    """ACI 9.6.1.2."""
    return max(3 * math.sqrt(fc) / fy, 200 / fy) * b * d


def As_min_slab(b, h, fy):
    """ACI 7.6.1.1 / 24.4.3.2 shrinkage and temperature minimum (per width b)."""
    if fy < 60000:
        rho = 0.0020
    else:
        rho = max(0.0018 * 60000 / fy, 0.0014)
    return rho * b * h


def lambda_s(d):
    """ACI 22.5.5.1.3 size effect factor (d in inches)."""
    return min(math.sqrt(2 / (1 + d / 10)), 1.0)


def Vc_oneway(b, d, fc, As=None, lam=1.0, Nu=0.0, Ag=None, Av_ge_min=False):
    """ACI 318-19 Table 22.5.5.1 one-way shear strength Vc (kip).

    Av_ge_min=True  -> Vc = 2*lam*sqrt(fc)*b*d (+ axial term) (Eq. a)
    Av_ge_min=False -> Vc = 8*lam_s*lam*rho_w^(1/3)*sqrt(fc)*b*d (Eq. c);
                       requires As (longitudinal tension steel).
    Nu in kip (+ compression). Upper limit 5*lam*sqrt(fc)*b*d.
    """
    sq = min(math.sqrt(fc), 100.0)  # ACI 22.5.3.1
    axial = 0.0
    if Nu and Ag:
        axial = Nu * 1000 / (6 * Ag)
    if Av_ge_min:
        v = 2 * lam * sq + axial
    else:
        if As is None:
            raise ValueError("As required for Vc when Av < Av,min")
        rho = As / (b * d)
        v = 8 * lambda_s(d) * lam * rho ** (1 / 3) * sq + axial
    v = max(v, 0.0)
    v = min(v, 5 * lam * sq)
    return v * b * d / 1000.0


def vc_twoway(d, bo, beta, alpha_s, fc, lam=1.0):
    """ACI 318-19 Table 22.6.5.2 two-way shear stress vc (psi), no shear reinf.

    beta: long/short side of the column; alpha_s: 40 interior, 30 edge, 20 corner.
    """
    sq = min(math.sqrt(fc), 100.0)
    ls = lambda_s(d)
    return min(4 * ls * lam * sq,
               (2 + 4 / beta) * ls * lam * sq,
               (2 + alpha_s * d / bo) * ls * lam * sq)


def gamma_v(b1, b2):
    """ACI 8.4.4.2.2: fraction of unbalanced moment transferred by shear.
    b1 = critical section dimension parallel to moment span, b2 perpendicular."""
    gf = 1 / (1 + 2 / 3 * math.sqrt(b1 / b2))
    return 1 - gf


def ld_simplified(db, fc, fy, lam=1.0, psi_t=1.0, psi_e=1.0, psi_g=None, clear_ok=True):
    """ACI 318-19 Table 25.4.2.3 tension development length (in), simplified.

    clear_ok: clear spacing >= 2db and cover >= db (or min stirrups) case.
    psi_g: grade factor (1.0 Gr 60, 1.15 Gr 80); default from fy.
    """
    if psi_g is None:
        psi_g = 1.0 if fy <= 60000 else (1.15 if fy <= 80000 else 1.3)
    sq = min(math.sqrt(fc), 100.0)
    small = db <= 0.75
    if clear_ok:
        k = 25 if small else 20
    else:
        k = 50 / 3 if small else 40 / 3
    ld = fy * psi_t * psi_e * psi_g / (k * lam * sq) * db
    return max(ld, 12.0)
