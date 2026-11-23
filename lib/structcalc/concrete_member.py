"""ACI 318-19 member helpers: strain-compatibility sections, P-M interaction,
cracked section properties, and effective moment of inertia.

UNITS follow structcalc.concrete: fc, fy in psi; dimensions in in;
forces in kip; moments in kip-in.

Depth coordinate ``y`` is measured DOWN from the compression (top) fiber.
Bars are (y, As). Compression is positive for P, strain and bar stress.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from . import concrete as C


# ---------------------------------------------------------------------------
# Concrete outlines
# ---------------------------------------------------------------------------
@dataclass
class RectShape:
    """Stack of rectangles: pieces = [(y_top, y_bot, width)], y down from top."""

    pieces: list

    @property
    def depth(self):
        return max(p[1] for p in self.pieces)

    def block(self, a):
        """Area and first moment about the top fiber of the region 0 <= y <= a."""
        A = Q = 0.0
        for y1, y2, w in self.pieces:
            lo, hi = y1, min(y2, a)
            if hi > lo:
                A += w * (hi - lo)
                Q += w * (hi**2 - lo**2) / 2
        return A, Q

    def width_at(self, y):
        return sum(w for y1, y2, w in self.pieces if y1 <= y <= y2)

    @property
    def area(self):
        return self.block(self.depth)[0]

    @property
    def ybar(self):
        A, Q = self.block(self.depth)
        return Q / A

    @property
    def Ig(self):
        yb = self.ybar
        I = 0.0
        for y1, y2, w in self.pieces:
            h = y2 - y1
            I += w * h**3 / 12 + w * h * ((y1 + y2) / 2 - yb) ** 2
        return I

    def flipped(self):
        D = self.depth
        return RectShape(sorted((D - y2, D - y1, w) for y1, y2, w in self.pieces))

    def compression_integrals(self, kd):
        """(int b*(kd-y) dy, int b*(kd-y)^2 dy) over 0..kd (cracked section)."""
        s1 = s2 = 0.0
        for y1, y2, w in self.pieces:
            lo, hi = y1, min(y2, kd)
            if hi > lo:
                s1 += w * ((kd - lo) ** 2 - (kd - hi) ** 2) / 2
                s2 += w * ((kd - lo) ** 3 - (kd - hi) ** 3) / 3
        return s1, s2


@dataclass
class CircleShape:
    D: float

    @property
    def depth(self):
        return self.D

    def block(self, a):
        R = self.D / 2
        a = min(max(a, 0.0), self.D)
        if a <= 0:
            return 0.0, 0.0
        if a >= self.D:
            A = math.pi * R**2
            return A, A * R
        th = 2 * math.acos((R - a) / R)           # central angle
        A = R**2 / 2 * (th - math.sin(th))
        ybar_c = 4 * R * math.sin(th / 2) ** 3 / (3 * (th - math.sin(th)))  # from center
        return A, A * (R - ybar_c)

    @property
    def area(self):
        return math.pi * self.D**2 / 4

    @property
    def ybar(self):
        return self.D / 2

    @property
    def Ig(self):
        return math.pi * self.D**4 / 64

    def flipped(self):
        return self


def rect(b, h):
    return RectShape([(0.0, h, b)])


def tee(bw, h, bf, hf):
    if bf <= bw or hf <= 0:
        return rect(bw, h)
    return RectShape([(0.0, hf, bf), (hf, h, bw)])


# ---------------------------------------------------------------------------
# Section with reinforcement
# ---------------------------------------------------------------------------
@dataclass
class Section:
    shape: object
    bars: list          # [(y from top, As)]
    fc: float           # psi
    fy: float           # psi
    transverse: str = "tied"   # tied | spiral (phi for compression-controlled)
    Es: float = C.ES
    yref: float = None  # moment reference (in from top); default = gross centroid

    def __post_init__(self):
        if self.yref is None:
            self.yref = self.shape.ybar
        self.bars = [(float(y), float(a)) for y, a in self.bars]

    @property
    def h(self):
        return self.shape.depth

    @property
    def Ast(self):
        return sum(a for _, a in self.bars)

    @property
    def Ag(self):
        return self.shape.area

    @property
    def Po(self):
        """ACI 22.4.2.2, kip."""
        return (0.85 * self.fc * (self.Ag - self.Ast) + self.fy * self.Ast) / 1000

    @property
    def Pn_max(self):
        """ACI 22.4.2.1 (0.80 Po tied, 0.85 Po spiral), kip."""
        return (0.85 if self.transverse == "spiral" else 0.80) * self.Po

    @property
    def dt(self):
        return max(y for y, _ in self.bars)

    def flipped(self):
        H = self.h
        return Section(self.shape.flipped(), [(H - y, a) for y, a in self.bars], self.fc,
                       self.fy, self.transverse, self.Es, H - self.yref)

    def state(self, c):
        """Forces at neutral axis depth c (in). Returns dict with Pn (kip),
        Mn (kip-in about yref), eps_t, phi, a, bar stresses (psi)."""
        b1 = C.beta1(self.fc)
        a = min(b1 * c, self.h)
        Ac, Qc = self.shape.block(a)
        Cc = 0.85 * self.fc * Ac                    # lb
        yc = Qc / Ac if Ac > 0 else 0.0
        P = Cc
        M = Cc * (self.yref - yc)
        fs_list = []
        for y, As in self.bars:
            eps = C.EPS_CU * (c - y) / c
            fs = max(min(self.Es * eps, self.fy), -self.fy)
            F = As * fs
            if fs > 0 and y <= a:
                F -= As * 0.85 * self.fc        # displaced concrete
            P += F
            M += F * (self.yref - y)
            fs_list.append(fs)
        eps_t = C.EPS_CU * (self.dt - c) / c
        phi = C.phi_flexure_axial(eps_t, self.fy, self.transverse)
        return {"c": c, "a": a, "Pn": P / 1000, "Mn": M / 1000, "eps_t": eps_t, "phi": phi,
                "fs": fs_list}

    # -- flexure (P = 0) ------------------------------------------------------
    def flexure(self, P=0.0):
        """Solve c for axial load P (kip, compression +). Returns state dict."""
        lo, hi = 1e-4, 50 * self.h
        if self.state(hi)["Pn"] < P:
            raise ValueError("axial load exceeds section capacity")
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if self.state(mid)["Pn"] > P:
                hi = mid
            else:
                lo = mid
            if hi - lo < 1e-9:
                break
        return self.state(0.5 * (lo + hi))

    # -- interaction ------------------------------------------------------------
    def interaction(self, n=160):
        """Nominal and design P-M curve, compression face at y = 0.

        Returns dict of numpy arrays ordered from pure tension to pure
        compression: Pn, Mn, phi, phiPn, phiMn (design curve capped at
        phi*Pn,max per ACI 22.4.2.1).
        """
        h = self.h
        cs = np.concatenate([np.geomspace(0.01 * h, 1.0 * h, int(n * 0.6)),
                             np.geomspace(1.0 * h, 30 * h, int(n * 0.4))[1:]])
        pts = [self.state(c) for c in cs]
        Pt = -self.fy * self.Ast / 1000
        Pn = np.array([Pt] + [p["Pn"] for p in pts] + [self.Po])
        Mn = np.array([0.0] + [p["Mn"] for p in pts] + [0.0])
        phi = np.array([0.9] + [p["phi"] for p in pts] + [0.75 if self.transverse == "spiral" else 0.65])
        phiPn = phi * Pn
        phiMn = phi * Mn
        cap = phi[-1] * self.Pn_max
        over = phiPn > cap
        if over.any():
            # cut the curve at phi*Pn,max: interpolate the crossing, then a
            # horizontal line to the axis (ACI 22.4.2.1)
            i = int(np.argmax(over))
            p1, p2, m1, m2 = phiPn[i - 1], phiPn[i], phiMn[i - 1], phiMn[i]
            mc = m1 + (m2 - m1) * (cap - p1) / (p2 - p1)
            dP = np.concatenate([phiPn[:i], [cap, cap]])
            dM = np.concatenate([phiMn[:i], [mc, 0.0]])
        else:
            dP, dM = phiPn, phiMn
        nP = np.where(Pn > self.Pn_max, self.Pn_max, Pn)
        return {"Pn": nP, "Mn": Mn, "phi": phi, "phiPn": dP, "phiMn": dM,
                "phiPn_max": cap, "Po": self.Po}

    # -- elastic / cracked properties (beams) ---------------------------------
    def cracked(self, Ec):
        """Cracked transformed section, compression face at top.

        Returns (kd in, Icr in^4). Ec in psi. Bars below the NA are
        transformed with n, bars above with (n - 1).
        """
        n = self.Es / Ec
        if not isinstance(self.shape, RectShape):
            raise NotImplementedError("cracked section only for rectangular outlines")

        def f(kd):
            s1, _ = self.shape.compression_integrals(kd)
            m = s1
            for y, As in self.bars:
                m += ((n - 1) if y < kd else n) * As * (kd - y)
            return m

        lo, hi = 1e-6, self.h
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        kd = 0.5 * (lo + hi)
        _, s2 = self.shape.compression_integrals(kd)
        Icr = s2 + sum(((n - 1) if y < kd else n) * As * (kd - y) ** 2 for y, As in self.bars)
        return kd, Icr


def Ie_aci(Ma, Mcr, Ig, Icr):
    """ACI 318-19 Table 24.2.3.5 effective moment of inertia.

    Ma, Mcr in the same units. Ie = Ig when Ma <= (2/3)Mcr, else
    Ie = Icr / (1 - ((2/3)Mcr/Ma)^2 (1 - Icr/Ig)).
    """
    Ma = abs(Ma)
    if Ma <= 2 / 3 * Mcr:
        return Ig
    return min(Icr / (1 - ((2 / 3) * Mcr / Ma) ** 2 * (1 - Icr / Ig)), Ig)


def ray_ratio(Mu, Pu, curveM, curveP):
    """Demand/capacity along the ray from the origin through (|Mu|, Pu).

    curveM/curveP: design curve ordered from tension to compression with
    M >= 0. Returns (ratio, M_cap, P_cap).
    """
    Mu = abs(Mu)
    r = math.hypot(Mu, Pu)
    if r < 1e-12:
        return 0.0, 0.0, 0.0
    th = math.atan2(Pu, Mu)
    ang = np.arctan2(curveP, curveM)
    # curve angles run from -pi/2 (pure tension) to +pi/2 (pure compression)
    for i in range(len(ang) - 1):
        a1, a2 = ang[i], ang[i + 1]
        if min(a1, a2) - 1e-12 <= th <= max(a1, a2) + 1e-12:
            # intersect ray with segment
            x1, y1, x2, y2 = curveM[i], curveP[i], curveM[i + 1], curveP[i + 1]
            dx, dy = math.cos(th), math.sin(th)
            den = dx * (y2 - y1) - dy * (x2 - x1)
            if abs(den) < 1e-14:
                continue
            t = (x1 * (y2 - y1) - y1 * (x2 - x1)) / den
            if t <= 0:
                continue
            return r / t, t * dx, t * dy
    return math.inf, 0.0, 0.0


def M_at_P(P, curveM, curveP):
    """Largest design moment on the curve at axial load P (kip)."""
    best = 0.0
    for i in range(len(curveP) - 1):
        p1, p2 = curveP[i], curveP[i + 1]
        if min(p1, p2) - 1e-9 <= P <= max(p1, p2) + 1e-9:
            if abs(p2 - p1) < 1e-12:
                m = max(curveM[i], curveM[i + 1])
            else:
                m = curveM[i] + (curveM[i + 1] - curveM[i]) * (P - p1) / (p2 - p1)
            best = max(best, m)
    return best
