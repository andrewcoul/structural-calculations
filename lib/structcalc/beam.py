"""Beam analysis by the direct stiffness method (Euler-Bernoulli elements).

Handles simple spans, cantilevers, overhangs, fixed ends, and continuous
multi-span beams with any mix of loads.

Sign convention (user side)
---------------------------
* x is measured from the left end, in ft.
* Loads are positive DOWNWARD (gravity). Uniform loads in kip/ft, point
  loads in kip, applied moments in kip-ft (positive = clockwise).
* Results: shear V in kip (beam convention), moment M in kip-ft
  (positive = sagging, tension at bottom), deflection in inches
  (positive = DOWNWARD), reactions in kip (positive = UPWARD).

Example
-------
    from structcalc.beam import Beam
    b = Beam(length=24, E=29000, I=510, supports=[(0, "pin"), (24, "pin")])
    b.add_udl(1.2, case="D")
    b.add_point(10, 6.0, case="L")
    res = b.solve()                 # dict of case -> BeamResult
    env = b.envelope(combos)        # combined
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import loads as _loads


@dataclass
class BeamResult:
    x: np.ndarray            # ft
    V: np.ndarray            # kip
    M: np.ndarray            # kip-ft
    defl: np.ndarray         # in, + down
    reactions: dict          # x(ft) -> (Ry kip up, Mz kip-ft) at supports

    def __add__(self, other):
        return BeamResult(self.x, self.V + other.V, self.M + other.M, self.defl + other.defl,
                          {k: (self.reactions[k][0] + other.reactions[k][0],
                               self.reactions[k][1] + other.reactions[k][1])
                           for k in self.reactions})

    def __mul__(self, f):
        return BeamResult(self.x, self.V * f, self.M * f, self.defl * f,
                          {k: (v[0] * f, v[1] * f) for k, v in self.reactions.items()})

    __rmul__ = __mul__

    @property
    def Mmax(self):
        return float(self.M.max())

    @property
    def Mmin(self):
        return float(self.M.min())

    @property
    def Mabs(self):
        return float(np.abs(self.M).max())

    @property
    def Vabs(self):
        return float(np.abs(self.V).max())

    @property
    def dmax(self):
        """Largest downward deflection (in)."""
        return float(self.defl.max())

    @property
    def dabs(self):
        return float(np.abs(self.defl).max())


@dataclass
class Beam:
    length: float                     # ft
    E: float = 29000.0                # ksi
    I: float = 100.0                  # in^4
    supports: list = field(default_factory=list)   # [(x_ft, "pin"|"fixed"), ...]
    n_elem: int = 240
    udls: list = field(default_factory=list)       # (x1, x2, w1, w2, case)
    points: list = field(default_factory=list)     # (x, P, case)
    moments: list = field(default_factory=list)    # (x, M, case)
    self_weight: float = 0.0                       # kip/ft, added to case D

    # -- load input -------------------------------------------------------
    def add_udl(self, w, case="D", x1=0.0, x2=None, w2=None):
        """Uniform (or trapezoidal when w2 given) load in kip/ft, + down."""
        x2 = self.length if x2 is None else x2
        w2 = w if w2 is None else w2
        if x2 <= x1:
            raise ValueError("x2 must be greater than x1")
        self.udls.append((float(x1), float(x2), float(w), float(w2), case))
        return self

    def add_point(self, x, P, case="D"):
        self.points.append((float(x), float(P), case))
        return self

    def add_moment(self, x, M, case="D"):
        self.moments.append((float(x), float(M), case))
        return self

    @property
    def cases(self):
        cs = {c for *_, c in self.udls} | {c for *_, c in self.points} | {c for *_, c in self.moments}
        if self.self_weight:
            cs.add("D")
        return sorted(cs, key=lambda c: _loads.CASES.index(c) if c in _loads.CASES else 99)

    # -- mesh -------------------------------------------------------------
    def _mesh(self):
        L = self.length
        pts = set(np.linspace(0.0, L, self.n_elem + 1).round(9))
        for x, _ in self.supports:
            pts.add(round(float(x), 9))
        for x1, x2, *_ in self.udls:
            pts.update([round(x1, 9), round(x2, 9)])
        for x, *_ in self.points + self.moments:
            pts.add(round(x, 9))
        xs = np.array(sorted(p for p in pts if -1e-9 <= p <= L + 1e-9))
        # merge near-duplicates
        keep = [xs[0]]
        for p in xs[1:]:
            if p - keep[-1] > 1e-6:
                keep.append(p)
        return np.array(keep)

    # -- solve --------------------------------------------------------------
    def _solve_case(self, case, xs):
        n = len(xs)
        ndof = 2 * n
        EI = self.E * self.I / 144.0  # kip-ft^2 (E ksi * I in^4 -> kip*in^2 / 144)
        K = np.zeros((ndof, ndof))
        F = np.zeros(ndof)         # equivalent nodal loads, y up
        Fe_all = np.zeros((n - 1, 4))
        h = np.diff(xs)
        for e in range(n - 1):
            le = h[e]
            k = EI / le**3 * np.array([
                [12, 6 * le, -12, 6 * le],
                [6 * le, 4 * le**2, -6 * le, 2 * le**2],
                [-12, -6 * le, 12, -6 * le],
                [6 * le, 2 * le**2, -6 * le, 4 * le**2]])
            idx = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
            K[np.ix_(idx, idx)] += k

        def wload(x):
            w = 0.0
            for x1, x2, w1, w2, c in self.udls:
                if c != case:
                    continue
                if x1 - 1e-9 <= x <= x2 + 1e-9:
                    w += w1 + (w2 - w1) * (x - x1) / (x2 - x1)
            if case == "D" and self.self_weight:
                w += self.self_weight
            return w

        for e in range(n - 1):
            le = h[e]
            # sample at element ends (just inside, so partial loads are exact)
            qa = -wload(xs[e] + 1e-7)
            qb = -wload(xs[e + 1] - 1e-7)
            fe = qa * np.array([7 * le / 20, le**2 / 20, 3 * le / 20, -le**2 / 30]) + \
                qb * np.array([3 * le / 20, le**2 / 30, 7 * le / 20, -le**2 / 20])
            Fe_all[e] = fe
            F[2 * e:2 * e + 4] += fe

        def node(x):
            return int(np.argmin(np.abs(xs - x)))

        for x, P, c in self.points:
            if c == case:
                F[2 * node(x)] += -P
        for x, M, c in self.moments:
            if c == case:
                F[2 * node(x) + 1] += -M  # clockwise positive (user) -> ccw positive (internal)

        fixed = []
        for x, typ in self.supports:
            i = node(x)
            fixed.append(2 * i)
            if typ == "fixed":
                fixed.append(2 * i + 1)
        fixed = sorted(set(fixed))
        free = [d for d in range(ndof) if d not in fixed]
        d = np.zeros(ndof)
        d[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
        R = K @ d - F

        V = np.zeros(n)
        M = np.zeros(n)
        # internal forces from element end forces
        Vl = np.zeros(n - 1)
        Ml = np.zeros(n - 1)
        Vr = np.zeros(n - 1)
        Mr = np.zeros(n - 1)
        for e in range(n - 1):
            le = h[e]
            k = EI / le**3 * np.array([
                [12, 6 * le, -12, 6 * le],
                [6 * le, 4 * le**2, -6 * le, 2 * le**2],
                [-12, -6 * le, 12, -6 * le],
                [6 * le, 2 * le**2, -6 * le, 4 * le**2]])
            f = k @ d[2 * e:2 * e + 4] - Fe_all[e]
            Vl[e], Ml[e], Vr[e], Mr[e] = f[0], -f[1], -f[2], f[3]
        # node values: use left-element right-end for M (continuous) and
        # report shear just right of each node (and last node from left)
        M[:-1] = Ml
        M[-1] = Mr[-1]
        V[:-1] = Vl
        V[-1] = Vr[-1]
        # Build a step-friendly shear array: duplicate x at nodes with jumps
        # Nodes where shear can jump (supports, point loads, load ends) get a
        # duplicate x entry; the set is the same for every load case, so
        # results of different cases have the same length and can be added.
        special = {int(np.argmin(np.abs(xs - x))) for x, _ in self.supports}
        special |= {int(np.argmin(np.abs(xs - x))) for x, *_ in self.points + self.moments}
        xx, VV, MM, DD = [], [], [], []
        for i in range(n):
            if 0 < i < n - 1 and i in special:
                xx.append(xs[i]); VV.append(Vr[i - 1]); MM.append(Mr[i - 1]); DD.append(-d[2 * i] * 12)
            xx.append(xs[i]); VV.append(V[i]); MM.append(M[i]); DD.append(-d[2 * i] * 12)
        reactions = {}
        for x, typ in self.supports:
            i = node(x)
            reactions[round(float(x), 6)] = (float(R[2 * i]),
                                             float(-R[2 * i + 1]) if typ == "fixed" else 0.0)
        return BeamResult(np.array(xx), np.array(VV), np.array(MM), np.array(DD), reactions)

    def solve(self):
        if not self.supports:
            raise ValueError("define at least one support")
        xs = self._mesh()
        out = {c: self._solve_case(c, xs) for c in self.cases}
        if not out:
            out["D"] = self._solve_case("D", xs)
        return out

    def combos(self, combos):
        """Return list of (name, factors, BeamResult) for each combination."""
        per = self.solve()
        res = []
        for name, f in _loads.active(combos, per.keys()):
            r = None
            for c, fac in f.items():
                if c in per:
                    r = per[c] * fac if r is None else r + per[c] * fac
            if r is not None:
                res.append((name, f, r))
        return res


def plot(result: BeamResult, title="", show_defl=True):
    import matplotlib.pyplot as plt

    n = 3 if show_defl else 2
    fig, ax = plt.subplots(n, 1, figsize=(9, 2.2 * n), sharex=True)
    ax[0].fill_between(result.x, result.V, color="tab:blue", alpha=0.3)
    ax[0].plot(result.x, result.V, color="tab:blue")
    ax[0].set_ylabel("V (kip)")
    ax[1].fill_between(result.x, result.M, color="tab:red", alpha=0.3)
    ax[1].plot(result.x, result.M, color="tab:red")
    ax[1].set_ylabel("M (kip-ft)")
    if show_defl:
        ax[2].plot(result.x, -result.defl, color="tab:green")
        ax[2].set_ylabel("Defl (in)")
    for a in ax:
        a.axhline(0, color="k", lw=0.6)
        a.grid(alpha=0.3)
    ax[-1].set_xlabel("x (ft)")
    if title:
        ax[0].set_title(title)
    fig.tight_layout()
    return fig
