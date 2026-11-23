"""AISC 360-22 member design (LRFD and ASD) for W, rectangular HSS and round HSS.

Shape properties: AISC Shapes Database v16.0 (W, HSS rect/square, HSS round).
Units: kip, in, ksi. Moments returned in kip-in unless named *_kft.
"""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass

E_STEEL = 29000.0  # ksi
G_STEEL = 11200.0  # ksi

_DATA = os.path.join(os.path.dirname(__file__), "data")
_FILES = {"W": "W.csv", "HSS_rect": "HSS_rect.csv", "HSS_round": "HSS_round.csv"}
_CACHE: dict = {}

# phi (LRFD) and Omega (ASD)
PHI = {"b": 0.90, "c": 0.90, "v": 0.90, "v_rolled_I": 1.00, "t": 0.90}
OMEGA = {"b": 1.67, "c": 1.67, "v": 1.67, "v_rolled_I": 1.50, "t": 1.67}


def _load(kind):
    if kind not in _CACHE:
        rows = {}
        with open(os.path.join(_DATA, _FILES[kind]), newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                name = r.pop("shape")
                d = {}
                for k, v in r.items():
                    v = (v or "").strip()
                    try:
                        d[k] = float(v)
                    except ValueError:
                        d[k] = None
                d["kind"] = kind
                d["name"] = name
                rows[name] = d
        _CACHE[kind] = rows
    return _CACHE[kind]


def normalize_name(name: str) -> str:
    n = name.strip().upper().replace(" ", "").replace("×", "X")
    n = n.replace("/", "_")
    if n.startswith("HSS"):
        n = n.replace(".", "_")
    return n


def shape(name: str) -> dict:
    """Look up a shape by AISC name, e.g. 'W18X35', 'HSS6X6X1/4', 'HSS6.625X0.280'."""
    n = normalize_name(name)
    for kind in _FILES:
        rows = _load(kind)
        if n in rows:
            return dict(rows[n])
    # round HSS names in the database use 3 decimals: HSS6_625X0_280
    raise KeyError(f"shape '{name}' not found in AISC v16.0 W / HSS tables")


def shapes(kind="W"):
    """List shape names of one kind ('W', 'HSS_rect', 'HSS_round')."""
    return list(_load(kind).keys())


def lightest(kind, predicate, nominal_depth=None):
    """Return names that pass ``predicate(props)`` sorted by weight."""
    out = []
    for n, p in _load(kind).items():
        if nominal_depth and kind == "W" and not n.startswith(f"W{nominal_depth}X"):
            continue
        try:
            if predicate(p):
                out.append((p["weight"], n))
        except (TypeError, ZeroDivisionError):
            continue
    return [n for _, n in sorted(out)]


def factor(kind: str, method: str):
    """Return the resistance multiplier: phi (LRFD) or 1/Omega (ASD)."""
    if method.upper() == "LRFD":
        return PHI[kind]
    return 1.0 / OMEGA[kind]


# ---------------------------------------------------------------------------
# Flexure
# ---------------------------------------------------------------------------
@dataclass
class FlexureResult:
    Mn: float          # nominal, kip-in
    limit_state: str
    Mp: float
    Lp: float = math.nan   # in
    Lr: float = math.nan   # in
    detail: str = ""


def _w_props(p):
    d, bf, tf, tw = p["d"], p["bf"], p["tf"], p["tw"]
    kdes = p["k"]
    h = d - 2 * kdes
    return d, bf, tf, tw, h


def flexure_W_major(p, Fy, Lb_in, Cb=1.0, E=E_STEEL):
    """AISC 360-22 F2 / F3 strong-axis flexure of doubly symmetric W shapes."""
    d, bf, tf, tw, h = _w_props(p)
    Zx, Sx, ry, J, rts, ho = p["Zx"], p["Sx"], p["ry"], p["J"], p["rts"], p["ho"]
    Mp = Fy * Zx
    c = 1.0
    Lp = 1.76 * ry * math.sqrt(E / Fy)
    jc = J * c / (Sx * ho)
    Lr = 1.95 * rts * E / (0.7 * Fy) * math.sqrt(jc + math.sqrt(jc**2 + 6.76 * (0.7 * Fy / E) ** 2))
    # yielding
    Mn_y = Mp
    # LTB
    if Lb_in <= Lp:
        Mn_ltb = Mp
    elif Lb_in <= Lr:
        Mn_ltb = min(Cb * (Mp - (Mp - 0.7 * Fy * Sx) * (Lb_in - Lp) / (Lr - Lp)), Mp)
    else:
        Fcr = Cb * math.pi**2 * E / (Lb_in / rts) ** 2 * math.sqrt(1 + 0.078 * jc * (Lb_in / rts) ** 2)
        Mn_ltb = min(Fcr * Sx, Mp)
    # flange local buckling (F3)
    lam = bf / (2 * tf)
    lpf = 0.38 * math.sqrt(E / Fy)
    lrf = 1.0 * math.sqrt(E / Fy)
    if lam <= lpf:
        Mn_flb = Mp
        flange = "compact"
    elif lam <= lrf:
        Mn_flb = Mp - (Mp - 0.7 * Fy * Sx) * (lam - lpf) / (lrf - lpf)
        flange = "noncompact"
    else:
        kc = min(max(4 / math.sqrt(h / tw), 0.35), 0.76)
        Mn_flb = 0.9 * E * kc * Sx / lam**2
        flange = "slender"
    web_ok = h / tw <= 3.76 * math.sqrt(E / Fy)
    cands = [(Mn_y, "Yielding (F2.1)"), (Mn_ltb, "Lateral-torsional buckling (F2.2)"),
             (Mn_flb, f"Flange local buckling (F3.2, {flange} flange)")]
    Mn, ls = min(cands, key=lambda t: t[0])
    detail = f"bf/2tf={lam:.2f} ({flange}); h/tw={h / tw:.1f} ({'compact' if web_ok else 'NONCOMPACT - F4 required'})"
    return FlexureResult(Mn, ls, Mp, Lp, Lr, detail)


def flexure_W_minor(p, Fy, E=E_STEEL):
    """AISC 360-22 F6 weak-axis flexure of W shapes."""
    bf, tf = p["bf"], p["tf"]
    Zy, Sy = p["Zy"], p["Sy"]
    Mp = min(Fy * Zy, 1.6 * Fy * Sy)
    lam = bf / (2 * tf)
    lpf = 0.38 * math.sqrt(E / Fy)
    lrf = 1.0 * math.sqrt(E / Fy)
    if lam <= lpf:
        return FlexureResult(Mp, "Yielding (F6.1)", Mp)
    if lam <= lrf:
        Mn = Mp - (Mp - 0.7 * Fy * Sy) * (lam - lpf) / (lrf - lpf)
        return FlexureResult(min(Mn, Mp), "Flange local buckling (F6.2)", Mp)
    Fcr = 0.70 * E / lam**2
    return FlexureResult(Fcr * Sy, "Flange local buckling (F6.2, slender)", Mp)


def _box_Se(B, H, t, bflat, be, S_table):
    """Effective elastic section modulus of a box with the compression
    flange reduced from flat width ``bflat`` to ``be`` (AISC F7.2c).

    A sharp-corner plate model gives the ratio Se/S; the ratio is applied
    to the tabulated S so that corner radii are respected.
    """
    def s_top(loss):
        Ag = 2 * B * t + 2 * (H - 2 * t) * t
        A_rm = loss * t
        y_rm = H - t / 2
        Ig = (B * H**3 - (B - 2 * t) * (H - 2 * t) ** 3) / 12
        yc = H / 2
        A = Ag - A_rm
        ybar = (Ag * yc - A_rm * y_rm) / A
        I = Ig + Ag * (yc - ybar) ** 2 - (A_rm * t**2 / 12 + A_rm * (y_rm - ybar) ** 2)
        # extreme fiber with the larger stress (compression side governs)
        return min(I / (H - ybar), I / ybar)

    return S_table * s_top(max(bflat - be, 0.0)) / s_top(0.0)


def flexure_HSS_rect(p, Fy, axis="x", Lb_in=0.0, Cb=1.0, E=E_STEEL):
    """AISC 360-22 F7 flexure of rectangular/square HSS."""
    t = p["tdes"]
    if axis == "x":
        Z, S, B, H, bflat, hflat, ry = p["Zx"], p["Sx"], p["B"], p["Ht"], p["b"], p["h"], p["ry"]
    else:
        Z, S, B, H, bflat, hflat, ry = p["Zy"], p["Sy"], p["Ht"], p["B"], p["h"], p["b"], p["rx"]
    Mp = Fy * Z
    lam_f = bflat / t
    lam_w = hflat / t
    lpf, lrf = 1.12 * math.sqrt(E / Fy), 1.40 * math.sqrt(E / Fy)
    lpw, lrw = 2.42 * math.sqrt(E / Fy), 5.70 * math.sqrt(E / Fy)
    cands = [(Mp, "Yielding (F7.1)")]
    if lam_f > lpf:
        if lam_f <= lrf:
            Mn = min(Mp - (Mp - Fy * S) * (3.57 * lam_f * math.sqrt(Fy / E) - 4.0), Mp)
            cands.append((Mn, "Flange local buckling (F7.2b, noncompact)"))
        else:
            be = min(1.92 * t * math.sqrt(E / Fy) * (1 - 0.38 / lam_f * math.sqrt(E / Fy)), bflat)
            Se = _box_Se(B, H, t, bflat, be, S)
            cands.append((Fy * Se, "Flange local buckling (F7.2c, slender)"))
    if lam_w > lpw:
        if lam_w <= lrw:
            Mn = min(Mp - (Mp - Fy * S) * (0.305 * lam_w * math.sqrt(Fy / E) - 0.738), Mp)
            cands.append((Mn, "Web local buckling (F7.3b)"))
        else:
            cands.append((Fy * S * 0.9, "Web local buckling (F7.3c, slender web - approx.)"))
    # LTB (F7.4) - only relevant for deep, narrow sections about the major axis
    J, Ag = p["J"], p["area"]
    if Lb_in > 0 and axis == "x" and H > B:
        Lp = 0.13 * E * ry * math.sqrt(J * Ag) / Mp
        Lr = 2 * E * ry * math.sqrt(J * Ag) / (0.7 * Fy * S)
        if Lb_in > Lp:
            if Lb_in <= Lr:
                Mn = min(Cb * (Mp - (Mp - 0.7 * Fy * S) * (Lb_in - Lp) / (Lr - Lp)), Mp)
            else:
                Mn = min(2 * E * Cb * math.sqrt(J * Ag) / (Lb_in / ry), Mp)
            cands.append((Mn, "Lateral-torsional buckling (F7.4)"))
    Mn, ls = min(cands, key=lambda c: c[0])
    return FlexureResult(Mn, ls, Mp, detail=f"b/t={lam_f:.1f}, h/t={lam_w:.1f}")


def flexure_HSS_round(p, Fy, E=E_STEEL):
    """AISC 360-22 F8 flexure of round HSS."""
    D, t, Z, S = p["OD"], p["tdes"], p["Zx"], p["Sx"]
    Dt = D / t
    Mp = Fy * Z
    if Dt <= 0.07 * E / Fy:
        return FlexureResult(Mp, "Yielding (F8.1)", Mp, detail=f"D/t={Dt:.1f} compact")
    if Dt <= 0.31 * E / Fy:
        Mn = (0.021 * E / Dt + Fy) * S
        return FlexureResult(min(Mn, Mp), "Local buckling (F8.2a)", Mp, detail=f"D/t={Dt:.1f} noncompact")
    return FlexureResult(0.33 * E / Dt * S, "Local buckling (F8.2b)", Mp, detail=f"D/t={Dt:.1f} slender")


def flexure(p, Fy, axis="x", Lb_in=0.0, Cb=1.0):
    k = p["kind"]
    if k == "W":
        return flexure_W_major(p, Fy, Lb_in, Cb) if axis == "x" else flexure_W_minor(p, Fy)
    if k == "HSS_rect":
        return flexure_HSS_rect(p, Fy, axis, Lb_in, Cb)
    return flexure_HSS_round(p, Fy)


def Cb_factor(MA, MB, MC, Mmax):
    """AISC 360-22 Eq. F1-1 (Rm = 1.0)."""
    MA, MB, MC, Mmax = abs(MA), abs(MB), abs(MC), abs(Mmax)
    if Mmax == 0:
        return 1.0
    return min(12.5 * Mmax / (2.5 * Mmax + 3 * MA + 4 * MB + 3 * MC), 3.0)


# ---------------------------------------------------------------------------
# Shear
# ---------------------------------------------------------------------------
def shear(p, Fy, axis="x", method="LRFD", E=E_STEEL):
    """Return (design/allowable shear strength kip, nominal Vn, note)."""
    k = p["kind"]
    if k == "W":
        d, bf, tf, tw, h = _w_props(p)
        if axis == "x":
            Aw = d * tw
            htw = h / tw
            kv = 5.34
            if htw <= 2.24 * math.sqrt(E / Fy):
                Cv1 = 1.0
                Vn = 0.6 * Fy * Aw * Cv1
                f = factor("v_rolled_I", method)
                return f * Vn, Vn, f"G2.1(a): h/tw={htw:.1f}, Cv1=1.0"
            Cv1 = 1.0 if htw <= 1.10 * math.sqrt(kv * E / Fy) else 1.10 * math.sqrt(kv * E / Fy) / htw
            Vn = 0.6 * Fy * Aw * Cv1
            return factor("v", method) * Vn, Vn, f"G2.1(b): h/tw={htw:.1f}, Cv1={Cv1:.3f}"
        # weak axis G6: two flanges, kv = 1.2
        lam = bf / (2 * tf)
        kv = 1.2
        Cv2 = _cv2(lam, kv, Fy, E)
        Vn = 0.6 * Fy * (2 * bf * tf) * Cv2
        return factor("v", method) * Vn, Vn, f"G6: Cv2={Cv2:.3f}"
    if k == "HSS_rect":
        t = p["tdes"]
        h = p["h"] if axis == "x" else p["b"]
        Aw = 2 * h * t
        Cv2 = _cv2(h / t, 5.0, Fy, E)
        Vn = 0.6 * Fy * Aw * Cv2
        return factor("v", method) * Vn, Vn, f"G4: h/t={h / t:.1f}, Cv2={Cv2:.3f}"
    # round HSS, G5 (conservative: Lv term ignored)
    D, t, Ag = p["OD"], p["tdes"], p["area"]
    Fcr = min(0.78 * E / (D / t) ** 1.5, 0.6 * Fy)
    Vn = Fcr * Ag / 2
    return factor("v", method) * Vn, Vn, f"G5: Fcr={Fcr:.1f} ksi"


def _cv2(htw, kv, Fy, E):
    a = 1.10 * math.sqrt(kv * E / Fy)
    b = 1.37 * math.sqrt(kv * E / Fy)
    if htw <= a:
        return 1.0
    if htw <= b:
        return a / htw
    return 1.51 * kv * E / (htw**2 * Fy)


# ---------------------------------------------------------------------------
# Compression
# ---------------------------------------------------------------------------
@dataclass
class CompressionResult:
    Pn: float
    Fcr: float
    Fe: float
    KL_r: float
    Ae: float
    governing_axis: str
    detail: str = ""


def _fcr_flexural(Fy, Fe):
    if Fy / Fe <= 2.25:
        return 0.658 ** (Fy / Fe) * Fy
    return 0.877 * Fe


def _eff_width(b, t, lam_r, c1, c2, Fy, Fcr):
    lam = b / t
    if lam <= lam_r * math.sqrt(Fy / Fcr):
        return b
    Fel = (c2 * lam_r / lam) ** 2 * Fy
    r = math.sqrt(Fel / Fcr)
    return min(b * (1 - c1 * r) * r, b)


def compression(p, Fy, Lcx_in, Lcy_in, Lcz_in=None, E=E_STEEL, G=G_STEEL):
    """AISC 360-22 Chapter E (E3 flexural buckling, E4 torsional for W, E7 slender)."""
    k = p["kind"]
    Ag = p["area"]
    rx, ry = p["rx"], p["ry"]
    slx = Lcx_in / rx
    sly = Lcy_in / ry
    sl = max(slx, sly)
    axis = "x" if slx >= sly else "y"
    Fe = math.pi**2 * E / sl**2 if sl > 0 else math.inf
    detail = []
    if k == "W" and Lcz_in:
        # E4 torsional buckling of doubly symmetric members
        Cw, J, Ix, Iy = p["Cw"], p["J"], p["Ix"], p["Iy"]
        Fe_t = (math.pi**2 * E * Cw / Lcz_in**2 + G * J) / (Ix + Iy)
        if Fe_t < Fe:
            Fe = Fe_t
            axis = "torsional"
            detail.append(f"E4 torsional Fe={Fe_t:.1f} ksi governs")
    Fcr = _fcr_flexural(Fy, Fe) if Fe != math.inf else Fy
    # slender elements (E7)
    Ae = Ag
    if k == "W":
        d, bf, tf, tw, h = _w_props(p)
        bfl = bf / 2
        lr_f = 0.56 * math.sqrt(E / Fy)
        lr_w = 1.49 * math.sqrt(E / Fy)
        be_f = _eff_width(bfl, tf, lr_f, 0.22, 1.49, Fy, Fcr)
        he = _eff_width(h, tw, lr_w, 0.18, 1.31, Fy, Fcr)
        Ae = Ag - 4 * (bfl - be_f) * tf - (h - he) * tw
        detail.append(f"bf/2tf={bfl / tf:.2f} (λr={lr_f:.1f}), h/tw={h / tw:.1f} (λr={lr_w:.1f})")
    elif k == "HSS_rect":
        t = p["tdes"]
        lr = 1.40 * math.sqrt(E / Fy)
        be = _eff_width(p["b"], t, lr, 0.20, 1.38, Fy, Fcr)
        he = _eff_width(p["h"], t, lr, 0.20, 1.38, Fy, Fcr)
        Ae = Ag - 2 * (p["b"] - be) * t - 2 * (p["h"] - he) * t
        detail.append(f"b/t={p['b'] / t:.1f}, h/t={p['h'] / t:.1f} (λr={lr:.1f})")
    else:
        Dt = p["OD"] / p["tdes"]
        if Dt > 0.11 * E / Fy:
            Ae = min((0.038 * E / (Fy * Dt) + 2 / 3) * Ag, Ag)
        detail.append(f"D/t={Dt:.1f} (λr={0.11 * E / Fy:.1f})")
    if Ae < Ag - 1e-6:
        detail.append(f"slender elements: Ae={Ae:.2f} in² (E7)")
    Pn = Fcr * Ae
    return CompressionResult(Pn, Fcr, Fe, sl, Ae, axis, "; ".join(detail))


# ---------------------------------------------------------------------------
# Combined forces
# ---------------------------------------------------------------------------
def interaction_H1(Pr, Pc, Mrx, Mcx, Mry=0.0, Mcy=1.0):
    """AISC 360-22 H1.1. Returns (ratio, equation)."""
    Pr, Mrx, Mry = abs(Pr), abs(Mrx), abs(Mry)
    if Pc <= 0:
        return math.inf, "-"
    if Pr / Pc >= 0.2:
        return Pr / Pc + 8 / 9 * (Mrx / Mcx + Mry / Mcy), "H1-1a"
    return Pr / (2 * Pc) + (Mrx / Mcx + Mry / Mcy), "H1-1b"


def tension_yield(p, Fy, method="LRFD"):
    Pn = Fy * p["area"]
    return factor("t", method) * Pn, Pn
