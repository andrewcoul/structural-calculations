"""ASCE 7-22 load combinations (Sections 2.3 and 2.4).

Load case keys:
    D  dead          L  live          Lr roof live     S  snow
    R  rain          W  wind (strength-level, ASCE 7-22)
    E  horizontal seismic effect Eh (= rho*QE)

The vertical seismic effect Ev = 0.2*SDS*D (ASCE 7-22 12.4.2.2) is
included when ``sds`` > 0.

Each combination is returned as (name, {case: factor}). Wind and seismic
combinations are generated with both signs of W/E where ``reversible``
is True, because wind uplift and seismic overturning can act in either
direction.
"""

from __future__ import annotations

CASES = ("D", "L", "Lr", "S", "R", "W", "E")


def _combo(name, **f):
    return name, {k: v for k, v in f.items() if v}


def lrfd(sds: float = 0.0, reversible: bool = True, roof_live_as_L: bool = False):
    """ASCE 7-22 2.3.1 / 2.3.6 strength design combinations."""
    c = []
    c.append(_combo("1: 1.4D", D=1.4))
    for nm, kw in (("0.5Lr", {"Lr": 0.5}), ("0.3S", {"S": 0.3}), ("0.5R", {"R": 0.5})):
        c.append(_combo(f"2: 1.2D+1.6L+{nm}", D=1.2, L=1.6, **kw))
    for nm, kw in (("1.6Lr", {"Lr": 1.6}), ("1.0S", {"S": 1.0}), ("1.6R", {"R": 1.6})):
        c.append(_combo(f"3: 1.2D+{nm}+L", D=1.2, L=1.0, **kw))
    signs = (1, -1) if reversible else (1,)
    for s in signs:
        sg = "+" if s > 0 else "-"
        for nm, kw in (("1.6Lr", {"Lr": 1.6}), ("1.0S", {"S": 1.0}), ("1.6R", {"R": 1.6})):
            c.append(_combo(f"3: 1.2D+{nm}{sg}0.5W", D=1.2, W=0.5 * s, **kw))
        for nm, kw in (("0.5Lr", {"Lr": 0.5}), ("0.3S", {"S": 0.3}), ("0.5R", {"R": 0.5})):
            c.append(_combo(f"4: 1.2D{sg}1.0W+L+{nm}", D=1.2, W=1.0 * s, L=1.0, **kw))
        c.append(_combo(f"5: 0.9D{sg}1.0W", D=0.9, W=1.0 * s))
        ev = 0.2 * sds
        c.append(_combo(f"6: (1.2+0.2SDS)D{sg}1.0E+L+0.15S", D=1.2 + ev, E=1.0 * s, L=1.0, S=0.15))
        c.append(_combo(f"7: (0.9-0.2SDS)D{sg}1.0E", D=0.9 - ev, E=1.0 * s))
    return _dedupe(c)


def asd(sds: float = 0.0, reversible: bool = True):
    """ASCE 7-22 2.4.1 / 2.4.5 allowable stress design combinations."""
    c = []
    c.append(_combo("1: D", D=1.0))
    c.append(_combo("2: D+L", D=1.0, L=1.0))
    for nm, kw in (("Lr", {"Lr": 1.0}), ("0.7S", {"S": 0.7}), ("R", {"R": 1.0})):
        c.append(_combo(f"3: D+{nm}", D=1.0, **kw))
    for nm, kw in (("0.75Lr", {"Lr": 0.75}), ("0.525S", {"S": 0.525}), ("0.75R", {"R": 0.75})):
        c.append(_combo(f"4: D+0.75L+{nm}", D=1.0, L=0.75, **kw))
    signs = (1, -1) if reversible else (1,)
    for s in signs:
        sg = "+" if s > 0 else "-"
        c.append(_combo(f"5: D{sg}0.6W", D=1.0, W=0.6 * s))
        for nm, kw in (("0.75Lr", {"Lr": 0.75}), ("0.525S", {"S": 0.525}), ("0.75R", {"R": 0.75})):
            c.append(_combo(f"6: D+0.75L{sg}0.45W+{nm}", D=1.0, L=0.75, W=0.45 * s, **kw))
        c.append(_combo(f"7: 0.6D{sg}0.6W", D=0.6, W=0.6 * s))
        ev = 0.2 * sds
        c.append(_combo(f"8: (1.0+0.14SDS)D{sg}0.7E", D=1.0 + 0.7 * ev, E=0.7 * s))
        c.append(_combo(f"9: (1.0+0.105SDS)D{sg}0.525E+0.75L+0.525S",
                        D=1.0 + 0.525 * ev, E=0.525 * s, L=0.75, S=0.525))
        c.append(_combo(f"10: (0.6-0.14SDS)D{sg}0.7E", D=0.6 - 0.7 * ev, E=0.7 * s))
    return _dedupe(c)


def service():
    """Serviceability combinations used for deflection checks."""
    return [("D", {"D": 1.0}), ("L", {"L": 1.0}), ("D+L", {"D": 1.0, "L": 1.0}),
            ("Lr", {"Lr": 1.0}), ("S", {"S": 1.0}), ("D+L+Lr+S", {"D": 1, "L": 1, "Lr": 1, "S": 1}),
            ("W", {"W": 1.0}), ("0.6W", {"W": 0.6})]


def _dedupe(combos):
    seen, out = set(), []
    for name, f in combos:
        key = tuple(sorted(f.items()))
        if key in seen:
            continue
        seen.add(key)
        out.append((name, f))
    return out


def active(combos, present):
    """Keep only combinations whose variable loads are all present.

    ``present`` is an iterable of case keys with non-zero load. A combo is
    kept if it contains at least one present case; combos that differ only
    by an absent case collapse onto the same effect, so they are removed.
    """
    present = set(present) | {"D"}
    out, seen = [], set()
    for name, f in combos:
        g = {k: v for k, v in f.items() if k in present}
        if not g:
            continue
        key = tuple(sorted(g.items()))
        if key in seen:
            continue
        seen.add(key)
        num = name.split(":", 1)[0]
        out.append((f"{num}: {label(g)}", g))
    return out


def label(factors):
    """Readable label, e.g. {'D': 1.2, 'L': 1.6} -> '1.2D+1.6L'."""
    parts = []
    for k in CASES:
        if k not in factors:
            continue
        v = factors[k]
        mag = abs(v)
        txt = (k if abs(mag - 1) < 1e-9 else f"{mag:.3g}{k}")
        parts.append(("-" if v < 0 else "+") + txt)
    s = "".join(parts)
    return s[1:] if s.startswith("+") else s


def combine(effects: dict, factors: dict):
    """Linear combination: sum(factor * effects[case])."""
    total = None
    for k, f in factors.items():
        if k not in effects:
            continue
        v = effects[k] * f
        total = v if total is None else total + v
    return total if total is not None else 0.0
