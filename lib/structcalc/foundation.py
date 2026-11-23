"""Foundation helpers: soil pressure under rigid footings, earth pressure
coefficients, and exact polygon integration used by the footing and
retaining wall modules.

Units: ft, kip, ksf (pressures under footings), pcf/psf where named.

Soil pressure model: rigid footing on linear (Winkler-type, no tension)
soil.  q(x, y) = max(0, a + b*x + c*y).  The coefficients are found so
that the pressure resultant equals the applied load P at (ex, ey).
Inside the kern the solution is the closed form P/A*(1 +/- 6ex/Lx +/- 6ey/Ly);
outside it, the contact polygon is found by fixed-point iteration with
exact polygon integrals (Green's theorem), so the result is exact for
biaxial partial contact (no grid error).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# ---------------------------------------------------------------------------
# Polygon geometry (exact integrals)
# ---------------------------------------------------------------------------


def rect(x1, y1, x2, y2):
    """Counter-clockwise rectangle polygon."""
    return [(x1, y1), (x2, y1), (x2, y2), (x1, y2)]


def clip(poly, a, b, c):
    """Clip polygon to the half plane a + b*x + c*y >= 0 (Sutherland-Hodgman)."""
    if not poly:
        return []
    out = []
    n = len(poly)
    for i in range(n):
        p = poly[i]
        q = poly[(i + 1) % n]
        fp = a + b * p[0] + c * p[1]
        fq = a + b * q[0] + c * q[1]
        if fp >= 0:
            out.append(p)
        if (fp >= 0) != (fq >= 0):
            t = fp / (fp - fq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def moments(poly):
    """Exact area integrals of a simple polygon (CCW).

    Returns dict A=∫1, Sx=∫x, Sy=∫y, Ixx=∫x², Iyy=∫y², Ixy=∫xy.
    """
    A = Sx = Sy = Ixx = Iyy = Ixy = 0.0
    n = len(poly)
    if n < 3:
        return dict(A=0.0, Sx=0.0, Sy=0.0, Ixx=0.0, Iyy=0.0, Ixy=0.0)
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        cr = x0 * y1 - x1 * y0
        A += cr
        Sx += (x0 + x1) * cr
        Sy += (y0 + y1) * cr
        Ixx += (x0 * x0 + x0 * x1 + x1 * x1) * cr
        Iyy += (y0 * y0 + y0 * y1 + y1 * y1) * cr
        Ixy += (x0 * y1 + 2 * x0 * y0 + 2 * x1 * y1 + x1 * y0) * cr
    return dict(A=A / 2, Sx=Sx / 6, Sy=Sy / 6, Ixx=Ixx / 12, Iyy=Iyy / 12, Ixy=Ixy / 24)


def area_centroid(poly):
    m = moments(poly)
    if abs(m["A"]) < 1e-12:
        return 0.0, 0.0, 0.0
    return m["A"], m["Sx"] / m["A"], m["Sy"] / m["A"]


# ---------------------------------------------------------------------------
# Soil pressure under a rectangular footing (biaxial, no tension)
# ---------------------------------------------------------------------------


@dataclass
class Pressure:
    """Linear soil pressure q = a + b x + c y (ksf) over the contact polygon.

    Origin at footing center; x along Lx, y along Ly (ft)."""

    a: float
    b: float
    c: float
    Lx: float
    Ly: float
    P: float
    ex: float
    ey: float
    contact: list
    iterations: int = 0

    def q(self, x, y):
        return np.maximum(self.a + self.b * np.asarray(x) + self.c * np.asarray(y), 0.0)

    @property
    def corners(self):
        hx, hy = self.Lx / 2, self.Ly / 2
        return [(-hx, -hy), (hx, -hy), (hx, hy), (-hx, hy)]

    @property
    def qmax(self):
        return float(max(self.a + self.b * x + self.c * y for x, y in self.corners))

    @property
    def qmin(self):
        return float(max(min(self.a + self.b * x + self.c * y for x, y in self.corners), 0.0))

    @property
    def contact_ratio(self):
        return moments(self.contact)["A"] / (self.Lx * self.Ly)

    @property
    def full_contact(self):
        return self.contact_ratio > 1 - 1e-9

    def integrals(self, halfplanes=()):
        """∫q, ∫q·x, ∫q·y over contact ∩ all half planes (a+bx+cy >= 0)."""
        poly = self.contact
        for hp in halfplanes:
            poly = clip(poly, *hp)
        m = moments(poly)
        a, b, c = self.a, self.b, self.c
        Q = a * m["A"] + b * m["Sx"] + c * m["Sy"]
        Qx = a * m["Sx"] + b * m["Ixx"] + c * m["Ixy"]
        Qy = a * m["Sy"] + b * m["Ixy"] + c * m["Iyy"]
        return Q, Qx, Qy


def in_kern(ex, ey, Lx, Ly):
    return abs(ex) / Lx + abs(ey) / Ly <= 1 / 6 + 1e-12


def bearing_pressure(P, ex, ey, Lx, Ly, tol=1e-10, max_iter=500):
    """Soil pressure for vertical load P (kip, > 0) at eccentricity (ex, ey) ft.

    Returns Pressure, or None when P <= 0 (net uplift) or when the resultant
    is outside the footing (overturning: no equilibrium possible).
    """
    if P <= 0:
        return None
    if abs(ex) >= Lx / 2 - 1e-9 or abs(ey) >= Ly / 2 - 1e-9:
        return None
    full = rect(-Lx / 2, -Ly / 2, Lx / 2, Ly / 2)
    A = Lx * Ly
    Iy = Ly * Lx**3 / 12  # ∫x²
    Ix = Lx * Ly**3 / 12  # ∫y²
    a, b, c = P / A, P * ex / Iy, P * ey / Ix
    if in_kern(ex, ey, Lx, Ly):
        return Pressure(a, b, c, Lx, Ly, P, ex, ey, full, 0)
    t = np.array([P, P * ex, P * ey])
    v = np.array([a, b, c])
    it = 0
    for it in range(1, max_iter + 1):
        poly = clip(full, *v)
        m = moments(poly)
        if m["A"] <= 1e-12:
            v = v * 0.5 + np.array([P / A, 0, 0]) * 0.5
            continue
        M = np.array([[m["A"], m["Sx"], m["Sy"]],
                      [m["Sx"], m["Ixx"], m["Ixy"]],
                      [m["Sy"], m["Ixy"], m["Iyy"]]])
        try:
            vn = np.linalg.solve(M, t)
        except np.linalg.LinAlgError:
            vn = np.linalg.lstsq(M, t, rcond=None)[0]
        if np.max(np.abs(vn - v)) <= tol * max(1.0, np.max(np.abs(vn))):
            v = vn
            break
        # mild damping keeps the iteration stable for corner contact
        v = 0.5 * v + 0.5 * vn if it > 50 else vn
    poly = clip(full, *v)
    return Pressure(float(v[0]), float(v[1]), float(v[2]), Lx, Ly, P, ex, ey, poly, it)


def pressure_1d(R, e, L):
    """Uniaxial pressure (force per unit length, e.g. kip/ft) under a rigid
    strip of length L with resultant R at eccentricity e from the center
    (+ toward the right end).

    Returns (q_left, q_right, contact_length). Partial contact: triangular
    pressure, q_max = 2R / (3*(L/2 - |e|)) over 3*(L/2 - |e|).
    """
    if R <= 0 or abs(e) >= L / 2:
        return None
    if abs(e) <= L / 6 + 1e-12:
        return R / L * (1 - 6 * e / L), R / L * (1 + 6 * e / L), L
    xc = L / 2 - abs(e)
    qm = 2 * R / (3 * xc)
    if e > 0:
        return 0.0, qm, 3 * xc
    return qm, 0.0, 3 * xc


def pressure_1d_at(x, q_left, q_right, contact, L):
    """Evaluate the 1D pressure at positions x (0 at the left end)."""
    x = np.asarray(x, dtype=float)
    if contact >= L - 1e-9:
        return q_left + (q_right - q_left) * x / L
    if q_right > 0:  # contact at right end
        x0 = L - contact
        return np.where(x >= x0, q_right * (x - x0) / contact, 0.0)
    return np.where(x <= contact, q_left * (contact - x) / contact, 0.0)


# ---------------------------------------------------------------------------
# Earth pressure coefficients
# ---------------------------------------------------------------------------


def ka_rankine(phi_deg, beta_deg=0.0):
    """Rankine active coefficient for a backfill slope beta (resultant
    parallel to the slope)."""
    p, b = math.radians(phi_deg), math.radians(beta_deg)
    if beta_deg >= phi_deg:
        raise ValueError("backfill slope must be less than the friction angle")
    cb, cp = math.cos(b), math.cos(p)
    r = math.sqrt(cb**2 - cp**2)
    return cb * (cb - r) / (cb + r)


def kp_rankine(phi_deg):
    return math.tan(math.radians(45 + phi_deg / 2)) ** 2


def ka_coulomb(phi_deg, delta_deg=0.0, beta_deg=0.0, omega_deg=0.0):
    """Coulomb active coefficient. omega = back-face inclination from the
    vertical (0 = vertical virtual back); resultant at delta from the normal."""
    p, d, b, w = (math.radians(v) for v in (phi_deg, delta_deg, beta_deg, omega_deg))
    if beta_deg >= phi_deg:
        raise ValueError("backfill slope must be less than the friction angle")
    num = math.cos(p - w) ** 2
    den = math.cos(w) ** 2 * math.cos(d + w) * (
        1 + math.sqrt(math.sin(p + d) * math.sin(p - b) / (math.cos(d + w) * math.cos(b - w)))) ** 2
    return num / den


def kae_mononobe_okabe(phi_deg, delta_deg=0.0, beta_deg=0.0, kh=0.0, kv=0.0, omega_deg=0.0):
    """Mononobe-Okabe seismic active coefficient KAE (AASHTO A11.3.1 form)."""
    theta = math.atan(kh / (1 - kv))
    p, d, b, w = (math.radians(v) for v in (phi_deg, delta_deg, beta_deg, omega_deg))
    if p - theta - b < 0:
        raise ValueError("phi - theta - beta < 0: M-O has no solution (reduce kh or slope)")
    num = math.cos(p - theta - w) ** 2
    den = math.cos(theta) * math.cos(w) ** 2 * math.cos(d + w + theta) * (
        1 + math.sqrt(math.sin(p + d) * math.sin(p - theta - b)
                      / (math.cos(d + w + theta) * math.cos(b - w)))) ** 2
    return num / den


def strip_lateral(q, a, b, z):
    """Lateral pressure (same unit as q) at depth z on a rigid wall from a
    strip load q of width b at distance a from the wall (elastic solution
    with the factor 2 for a rigid wall, e.g. Das, Principles of Foundation
    Engineering, Eq. for strip surcharge): sigma = 2q/pi (beta - sin(beta)cos(2 alpha))."""
    z = np.maximum(np.asarray(z, dtype=float), 1e-9)
    t1 = np.arctan(a / z)
    t2 = np.arctan((a + b) / z)
    beta = t2 - t1
    alpha = t1 + beta / 2
    return 2 * q / math.pi * (beta - np.sin(beta) * np.cos(2 * alpha))


# ---------------------------------------------------------------------------
# Small helpers for the concrete checks
# ---------------------------------------------------------------------------


def bearing_A2_factor(c1, c2, h, avail1, avail2):
    """sqrt(A2/A1) for ACI 22.8.3.2 with a 1:2 frustum.

    c1, c2: loaded area dims (in); h: support depth (in); avail1/2: max
    available dims of the supporting surface (in), concentric load assumed.
    """
    k = min((c1 + 4 * h) / c1, (c2 + 4 * h) / c2, avail1 / c1, avail2 / c2)
    return min(max(k, 1.0), 2.0)


def bars_area(bar_no, spacing_in, width_in):
    """Area of bars at given spacing across a width (in²)."""
    from .concrete import bar

    return bar(bar_no)[1] * width_in / spacing_in


def punching_section(xc, yc, cx, cy, d, x_lim, y_lim):
    """Critical perimeter at d/2 from a rectangular column (ACI 22.6.4.1),
    clipped at footing edges.

    xc, yc: column center; cx, cy: column dims; d: effective depth; all in
    the same unit (ft). x_lim, y_lim: (min, max) footing extents.
    Returns dict with rectangle bounds, bo, list of present sides, alpha_s,
    b1 (along x), b2 (along y), centroid xbar, Jc about y-axis (moment in
    x-z plane), c_max.
    """
    x1, x2 = xc - (cx + d) / 2, xc + (cx + d) / 2
    y1, y2 = yc - (cy + d) / 2, yc + (cy + d) / 2
    sides = []
    if x1 > x_lim[0] + 1e-9:
        sides.append(("x1", x1, y1, x1, y2))
    if x2 < x_lim[1] - 1e-9:
        sides.append(("x2", x2, y1, x2, y2))
    if y1 > y_lim[0] + 1e-9:
        sides.append(("y1", max(x1, x_lim[0]), y1, min(x2, x_lim[1]), y1))
    if y2 < y_lim[1] - 1e-9:
        sides.append(("y2", max(x1, x_lim[0]), y2, min(x2, x_lim[1]), y2))
    x1c, x2c = max(x1, x_lim[0]), min(x2, x_lim[1])
    y1c, y2c = max(y1, y_lim[0]), min(y2, y_lim[1])
    bo = sum(math.hypot(s[3] - s[1], s[4] - s[2]) for s in sides)
    n_sides = len(sides)
    alpha_s = {4: 40, 3: 30}.get(n_sides, 20)
    # centroid and Jc (moment about the y-axis, i.e. varying stress along x)
    Lsum = 0.0
    xs = 0.0
    for _, ax, ay, bx, by in sides:
        ln = math.hypot(bx - ax, by - ay)
        Lsum += ln
        xs += ln * (ax + bx) / 2
    xbar = xs / Lsum if Lsum else xc
    Jc = 0.0
    for _, ax, ay, bx, by in sides:
        ln = math.hypot(bx - ax, by - ay)
        xm = (ax + bx) / 2
        if abs(by - ay) < 1e-12:  # side parallel to x (spans in moment direction)
            Jc += d * ln**3 / 12 + ln * d**3 / 12 + ln * d * (xm - xbar) ** 2
        else:  # side parallel to y
            Jc += ln * d * (xm - xbar) ** 2
    cmax = max(abs(x1c - xbar), abs(x2c - xbar)) if sides else 0.0
    return dict(x1=x1c, x2=x2c, y1=y1c, y2=y2c, bo=bo, n_sides=n_sides, alpha_s=alpha_s,
                b1=x2c - x1c, b2=y2c - y1c, xbar=xbar, Jc=Jc, c=cmax)


def gamma_v(b1, b2):
    """ACI 318-19 Eq. 8.4.2.2.2: gamma_v = 1 - 1/(1 + (2/3)sqrt(b1/b2))."""
    return 1 - 1 / (1 + 2 / 3 * math.sqrt(b1 / b2))
