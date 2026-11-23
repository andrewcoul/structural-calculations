"""Shared input handling for beam modules (steel, wood, concrete).

Span types (Colab dropdown values):
    "Simple span", "Cantilever", "Overhang one end (right)",
    "Overhang both ends", "Fixed-fixed", "Fixed-pinned (propped)",
    "Custom supports"

Custom supports text: "0:pin, 20:pin, 40:pin" (x in ft; pin or fixed).

Extra loads text, one load per line or separated by ";" (x in ft, loads + down):
    P, x, kip, case                      point load
    w, x1, x2, kip/ft, case              partial uniform load
    t, x1, x2, w1, w2, case              trapezoidal load (kip/ft)
    M, x, kip-ft, case                   applied moment (clockwise +)
Case: D, L, Lr, S, R, W, E.  Lines starting with # are ignored.
"""

from __future__ import annotations

from .beam import Beam
from .loads import CASES

SPAN_TYPES = [
    "Simple span",
    "Cantilever",
    "Overhang one end (right)",
    "Overhang both ends",
    "Fixed-fixed",
    "Fixed-pinned (propped)",
    "Custom supports",
]


def make_supports(span_type, span, overhang_left=0.0, overhang_right=0.0, custom=""):
    """Return (total_length_ft, supports list)."""
    t = span_type
    if t == "Simple span":
        return span, [(0.0, "pin"), (span, "pin")]
    if t == "Cantilever":
        return span, [(0.0, "fixed")]
    if t == "Overhang one end (right)":
        L = span + overhang_right
        return L, [(0.0, "pin"), (span, "pin")]
    if t == "Overhang both ends":
        L = overhang_left + span + overhang_right
        return L, [(overhang_left, "pin"), (overhang_left + span, "pin")]
    if t == "Fixed-fixed":
        return span, [(0.0, "fixed"), (span, "fixed")]
    if t == "Fixed-pinned (propped)":
        return span, [(0.0, "fixed"), (span, "pin")]
    if t == "Custom supports":
        sup = []
        for tok in custom.replace(";", ",").split(","):
            tok = tok.strip()
            if not tok:
                continue
            x, _, typ = tok.partition(":")
            typ = (typ or "pin").strip().lower()
            if typ not in ("pin", "fixed"):
                raise ValueError(f"support type must be pin or fixed, got '{typ}'")
            sup.append((float(x), typ))
        if not sup:
            raise ValueError("custom supports text is empty")
        sup.sort()
        L = max(span, sup[-1][0])
        return L, sup
    raise ValueError(f"unknown span type '{span_type}'")


def parse_loads(text):
    """Parse the extra-loads text box. Returns list of tuples for apply_loads."""
    out = []
    import re

    for n, raw in enumerate(re.split(r"[;\n]", text or ""), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(",")]
        kind = parts[0].upper()
        try:
            case = parts[-1]
            if case not in CASES:
                raise ValueError(f"case '{case}' not one of {CASES}")
            nums = [float(p) for p in parts[1:-1]]
            if kind == "P" and len(nums) == 2:
                out.append(("P", nums[0], nums[1], case))
            elif kind == "W" and len(nums) == 3:
                out.append(("w", nums[0], nums[1], nums[2], nums[2], case))
            elif kind == "T" and len(nums) == 4:
                out.append(("w", nums[0], nums[1], nums[2], nums[3], case))
            elif kind == "M" and len(nums) == 2:
                out.append(("M", nums[0], nums[1], case))
            else:
                raise ValueError("wrong number of values")
        except ValueError as e:
            raise ValueError(f"extra loads line {n} '{raw.strip()}': {e}") from None
    return out


def apply_loads(beam: Beam, uniform: dict, extra: list):
    """uniform: {case: kip/ft over full length}; extra: parse_loads output."""
    for case, w in uniform.items():
        if w:
            beam.add_udl(w, case)
    for ld in extra:
        if ld[0] == "P":
            beam.add_point(ld[1], ld[2], ld[3])
        elif ld[0] == "w":
            beam.add_udl(ld[3], ld[5], x1=ld[1], x2=ld[2], w2=ld[4])
        elif ld[0] == "M":
            beam.add_moment(ld[1], ld[2], ld[3])
    return beam


def regions(length, supports):
    """Split the beam into spans and overhangs for deflection limits.

    Returns list of (x1, x2, reference_length_ft, label). Cantilever /
    overhang regions use twice their length (IBC Table 1604.3 note).
    """
    xs = sorted(x for x, _ in supports)
    out = []
    if len(xs) == 1:
        x0 = xs[0]
        if x0 > 1e-9:
            out.append((0.0, x0, 2 * x0, "Cantilever (left)"))
        if length - x0 > 1e-9:
            out.append((x0, length, 2 * (length - x0), "Cantilever"))
        return out
    if xs[0] > 1e-9:
        out.append((0.0, xs[0], 2 * xs[0], "Left overhang"))
    for i in range(len(xs) - 1):
        out.append((xs[i], xs[i + 1], xs[i + 1] - xs[i], f"Span {i + 1}"))
    if length - xs[-1] > 1e-9:
        out.append((xs[-1], length, 2 * (length - xs[-1]), "Right overhang"))
    return out


def deflection_checks(report, beam_cases, length, supports, limits, Ief_note=""):
    """Add deflection checks to ``report``.

    beam_cases: dict case -> BeamResult (service, unfactored).
    limits: list of (label, {case: factor}, ratio_denominator), e.g.
            [("Live", {"L": 1}, 360), ("Total", {"D": 1, "L": 1}, 240)]
    """
    import numpy as np

    regs = regions(length, supports)
    for label, combo, denom in limits:
        if not any(c in beam_cases for c in combo):
            continue
        r = None
        for c, f in combo.items():
            if c in beam_cases:
                r = beam_cases[c] * f if r is None else r + beam_cases[c] * f
        for x1, x2, Lref, rname in regs:
            m = (r.x >= x1 - 1e-9) & (r.x <= x2 + 1e-9)
            d = float(np.abs(r.defl[m]).max()) if m.any() else 0.0
            allow = Lref * 12 / denom
            ratio_txt = f"L/{Lref * 12 / d:.0f}" if d > 1e-9 else "L/inf"
            report.check(f"Deflection {label} - {rname}", d, allow, "in",
                         ref=f"IBC 1604.3 (L/{denom})", note=f"{ratio_txt}{Ief_note}")
