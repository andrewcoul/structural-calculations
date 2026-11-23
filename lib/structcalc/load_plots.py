"""Structural free-body / load schematic graphics for beams and columns."""

from __future__ import annotations

import numpy as np

CASE_COLORS = {
    "D": "#475569",   # slate
    "L": "#2563eb",   # blue
    "Lr": "#0d9488",  # teal
    "S": "#0284c7",   # sky
    "R": "#0891b2",   # cyan
    "W": "#7c3aed",   # violet
    "E": "#ea580c",   # orange
}


def _fmt_num(v: float) -> str:
    av = abs(v)
    if av >= 100:
        return f"{v:.0f}"
    if av >= 10:
        return f"{v:.1f}".rstrip("0").rstrip(".")
    return f"{v:.2f}".rstrip("0").rstrip(".")


def _case_summary(vals: dict, unit: str = "", scale: float = 1.0) -> str:
    parts = []
    for c in ("D", "L", "Lr", "S", "R", "W", "E"):
        v = vals.get(c, 0.0) * scale
        if abs(v) > 1e-9:
            parts.append(f"{_fmt_num(v)} {c}")
    if not parts:
        return ""
    s = ", ".join(parts)
    return f"{s} {unit}".strip()


def _draw_pin(ax, x: float, y0: float, h: float, w: float):
    from matplotlib.patches import Polygon

    tri = Polygon([(x, y0), (x - w / 2, y0 - h), (x + w / 2, y0 - h)],
                  closed=True, facecolor="#cbd5e1", edgecolor="#1e293b", lw=1.3, zorder=5)
    ax.add_patch(tri)
    ax.plot([x - w * 0.72, x + w * 0.72], [y0 - h, y0 - h], color="#1e293b", lw=1.5, zorder=5)


def _draw_fixed(ax, x: float, y0: float, h: float, w: float, side: str = "left"):
    from matplotlib.patches import Rectangle

    dx = -w if side == "left" else (0.0 if side == "right" else -w / 2)
    rect = Rectangle((x + dx, y0 - h), w, 2 * h,
                     facecolor="#cbd5e1", edgecolor="#1e293b", lw=1.3, hatch="////", zorder=5)
    ax.add_patch(rect)
    ax.plot([x, x], [y0 - h, y0 + h], color="#1e293b", lw=2.0, zorder=6)


def _draw_moment_arc(ax, cx: float, cy: float, rx: float, ry: float, clockwise: bool, color: str):
    """Draw a clean curved moment arrow centered at (cx, cy) with radii (rx, ry)."""
    if clockwise:
        th = np.linspace(np.deg2rad(145), np.deg2rad(-35), 40)
    else:
        th = np.linspace(np.deg2rad(35), np.deg2rad(215), 40)
    xs = cx + rx * np.cos(th)
    ys = cy + ry * np.sin(th)
    ax.plot(xs[:-2], ys[:-2], color=color, lw=2.0, zorder=7)
    ax.annotate(
        "",
        xy=(xs[-1], ys[-1]),
        xytext=(xs[-5], ys[-5]),
        arrowprops=dict(arrowstyle="-|>,head_length=0.45,head_width=0.28", color=color, lw=1.8),
        zorder=8,
    )


def draw_beam_loads(ax, beam, title: str = "", brace_pts=None, bracing_note: str = ""):
    """Draw a free-body diagram of a Beam (supports, UDLs/trapezoids, point loads, moments)."""
    L = float(beam.length)
    ax.set_xlim(-0.06 * L, 1.06 * L)

    # Collect all distributed loads including self weight (in kip/ft)
    udls = list(beam.udls)
    if beam.self_weight:
        # Merge into an existing full-length D load if present, else add separately
        merged = False
        for i, (x1, x2, w1, w2, c) in enumerate(udls):
            if c == "D" and abs(x1) < 1e-9 and abs(x2 - L) < 1e-9 and abs(w1 - w2) < 1e-9:
                udls[i] = (x1, x2, w1 + beam.self_weight, w2 + beam.self_weight, "D+SW")
                merged = True
                break
        if not merged:
            udls.insert(0, (0.0, L, beam.self_weight, beam.self_weight, "D (SW)"))

    # Separate downward vs upward UDLs for clean stacking
    down_udls = [u for u in udls if (u[2] + u[3]) >= 0]
    up_udls = [u for u in udls if (u[2] + u[3]) < 0]

    # Normalized vertical layout coordinates (beam at y = 0)
    band_h = 0.72
    gap = 0.18
    y_cursor = 0.14
    max_udl = max((max(abs(u[2]), abs(u[3])) for u in udls), default=1.0) or 1.0

    def _udl_color(case: str) -> str:
        base = case.split("+")[0].split()[0]
        return CASE_COLORS.get(base, "#334155")

    # Draw downward distributed loads stacked neatly above the beam
    for x1, x2, w1, w2, case in down_udls:
        col = _udl_color(case)
        h1 = max(abs(w1) / max_udl, 0.28) * band_h
        h2 = max(abs(w2) / max_udl, 0.28) * band_h
        y_bot = y_cursor
        xs_poly = [x1, x2, x2, x1]
        ys_poly = [y_bot, y_bot, y_bot + h2, y_bot + h1]
        ax.fill(xs_poly, ys_poly, color=col, alpha=0.14, zorder=2)
        ax.plot([x1, x2], [y_bot + h1, y_bot + h2], color=col, lw=1.6, zorder=3)
        # Evenly spaced downward arrows across [x1, x2]
        span_frac = max((x2 - x1) / L, 0.05)
        n_arr = max(3, min(11, int(round(9 * span_frac)) + 1))
        for xa in np.linspace(x1, x2, n_arr):
            t = (xa - x1) / max(x2 - x1, 1e-9)
            ht = h1 + (h2 - h1) * t
            ax.annotate(
                "",
                xy=(xa, y_bot),
                xytext=(xa, y_bot + ht),
                arrowprops=dict(arrowstyle="-|>,head_length=0.32,head_width=0.18", color=col, lw=1.1),
                zorder=4,
            )
        # Magnitude label centered above the band
        if abs(w1 - w2) < 1e-6:
            mag = f"{w1 * 1000:.0f} plf" if abs(w1) < 1.5 else f"{_fmt_num(w1)} k/ft"
        else:
            mag = f"{w1 * 1000:.0f}–{w2 * 1000:.0f} plf" if max(abs(w1), abs(w2)) < 1.5 else f"{_fmt_num(w1)}–{_fmt_num(w2)} k/ft"
        xm = 0.5 * (x1 + x2)
        # Avoid placing the UDL label directly under a point-load shaft
        if beam.points and any(abs(xp - xm) < 0.14 * L for xp, *_ in beam.points):
            cands = [x1 + f * (x2 - x1) for f in (0.25, 0.75, 0.18, 0.82)]
            xm = max(cands, key=lambda xc: min(abs(xc - xp) for xp, *_ in beam.points))
        t_lbl = (xm - x1) / max(x2 - x1, 1e-9)
        ym = y_bot + (h1 + (h2 - h1) * t_lbl) + 0.06
        ax.text(xm, ym, f"{case}: {mag}", ha="center", va="bottom", fontsize=8.5,
                fontweight="semibold", color=col,
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none", alpha=0.82), zorder=6)
        y_cursor = y_bot + max(h1, h2) + 0.34

    # Draw upward distributed loads (if any) below beam
    y_below_cursor = -0.52
    for x1, x2, w1, w2, case in up_udls:
        col = _udl_color(case)
        h1 = max(abs(w1) / max_udl, 0.28) * band_h
        h2 = max(abs(w2) / max_udl, 0.28) * band_h
        y_top = y_below_cursor
        ax.fill([x1, x2, x2, x1], [y_top, y_top, y_top - h2, y_top - h1], color=col, alpha=0.14, zorder=2)
        ax.plot([x1, x2], [y_top - h1, y_top - h2], color=col, lw=1.6, zorder=3)
        n_arr = max(3, min(9, int(round(9 * max((x2 - x1) / L, 0.05))) + 1))
        for xa in np.linspace(x1, x2, n_arr):
            t = (xa - x1) / max(x2 - x1, 1e-9)
            ht = h1 + (h2 - h1) * t
            ax.annotate("", xy=(xa, y_top), xytext=(xa, y_top - ht),
                        arrowprops=dict(arrowstyle="-|>,head_length=0.32,head_width=0.18", color=col, lw=1.1),
                        zorder=4)
        mag = f"{abs(w1) * 1000:.0f} plf" if abs(w1) < 1.5 else f"{_fmt_num(abs(w1))} k/ft"
        ax.text(0.5 * (x1 + x2), y_top - max(h1, h2) - 0.08, f"{case}: {mag} ↑", ha="center", va="top",
                fontsize=8.5, fontweight="semibold", color=col, zorder=6)
        y_below_cursor = y_top - max(h1, h2) - 0.38

    # Point loads: bold crimson/case arrows pointing cleanly to the beam
    y_top_max = max(y_cursor, 0.95)
    if beam.points:
        p_top = y_cursor + 0.72
        y_top_max = p_top + 0.34
        for idx, (xp, P, case) in enumerate(beam.points):
            col = "#dc2626" if case in ("L", "D") else _udl_color(case)
            if P >= 0:
                ax.annotate(
                    "",
                    xy=(xp, 0.08),
                    xytext=(xp, p_top),
                    arrowprops=dict(arrowstyle="-|>,head_length=0.55,head_width=0.32", color=col, lw=2.3),
                    zorder=8,
                )
                ax.text(xp, p_top + 0.06, f"P = {_fmt_num(P)} k ({case})", ha="center", va="bottom",
                        fontsize=8.5, fontweight="bold", color=col,
                        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=col, lw=0.8, alpha=0.95), zorder=9)
            else:
                p_bot = y_below_cursor - 0.55
                y_below_cursor = min(y_below_cursor, p_bot - 0.3)
                ax.annotate(
                    "",
                    xy=(xp, -0.08),
                    xytext=(xp, p_bot),
                    arrowprops=dict(arrowstyle="-|>,head_length=0.55,head_width=0.32", color=col, lw=2.3),
                    zorder=8,
                )
                ax.text(xp, p_bot - 0.06, f"P = {_fmt_num(abs(P))} k ↑ ({case})", ha="center", va="top",
                        fontsize=8.5, fontweight="bold", color=col,
                        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=col, lw=0.8, alpha=0.95), zorder=9)

    # Applied moments
    if beam.moments:
        for xm, M, case in beam.moments:
            col = "#b45309"
            _draw_moment_arc(ax, xm, 0.0, 0.025 * L, 0.38, clockwise=(M >= 0), color=col)
            ax.text(xm, 0.44, f"M = {_fmt_num(abs(M))} k-ft ({case})", ha="center", va="bottom",
                    fontsize=8.5, fontweight="bold", color=col,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=col, lw=0.7, alpha=0.95), zorder=9)
            y_top_max = max(y_top_max, 0.95)

    # Beam bar
    ax.plot([0.0, L], [0.0, 0.0], color="#0f172a", lw=4.5, solid_capstyle="butt", zorder=5)

    # Optional lateral brace markers
    if brace_pts:
        for xb in brace_pts:
            ax.plot([xb], [0.0], marker="x", ms=6.5, mew=1.8, color="#0d9488", zorder=7)
    if bracing_note:
        ax.text(L, y_top_max + 0.08, bracing_note, ha="right", va="top", fontsize=7.8,
                color="#0f766e", style="italic",
                bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="#99f6e4", lw=0.7, alpha=0.9), zorder=8)

    # Supports
    sup_h = 0.34
    sup_w = 0.032 * L
    for xs, typ in beam.supports:
        if typ == "fixed":
            side = "left" if xs <= 1e-6 else ("right" if xs >= L - 1e-6 else "center")
            _draw_fixed(ax, xs, 0.0, sup_h * 0.95, sup_w * 0.7, side=side)
        else:
            _draw_pin(ax, xs, -0.04, sup_h, sup_w)

    # Span dimension ticks below supports
    dim_y = min(y_below_cursor, -0.58)
    pts = sorted({0.0, L} | {x for x, _ in beam.supports})
    for a, b in zip(pts[:-1], pts[1:]):
        if b - a < 1e-6:
            continue
        ax.annotate("", xy=(b, dim_y), xytext=(a, dim_y),
                    arrowprops=dict(arrowstyle="<->", color="#475569", lw=1.0), zorder=4)
        ax.text(0.5 * (a + b), dim_y - 0.06, f"{_fmt_num(b - a)} ft", ha="center", va="top",
                fontsize=8.2, color="#334155")

    ax.set_ylim(dim_y - 0.38, y_top_max + 0.14)
    ax.set_yticks([])
    ax.set_ylabel("Loads", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", alpha=0.25)
    if title:
        ax.set_title(title, fontsize=10.5, fontweight="semibold")


def draw_column_loads(
    ax,
    H: float,
    P_cases: dict,
    Mx_top: dict | None = None,
    Mx_bot: dict | None = None,
    My_top: dict | None = None,
    My_bot: dict | None = None,
    w_lat_x: dict | None = None,
    w_lat_y: dict | None = None,
    ex: float = 0.0,
    ey: float = 0.0,
    bracing_note: str = "",
    brace_z: list | None = None,
    title: str = "Applied service loads",
):
    """Draw a vertical column free-body / load schematic on ``ax`` (y-axis = height in ft)."""
    Mx_top, Mx_bot = Mx_top or {}, Mx_bot or {}
    My_top, My_bot = My_top or {}, My_bot or {}
    w_lat_x, w_lat_y = w_lat_x or {}, w_lat_y or {}

    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-0.22 * H, 1.32 * H)

    # Column shaft
    ax.plot([0.0, 0.0], [0.0, H], color="#0f172a", lw=5.5, solid_capstyle="butt", zorder=5)

    # End pin / support triangles (base & top lateral restraint)
    from matplotlib.patches import Polygon

    bh = 0.06 * H
    bw = 0.24
    ax.add_patch(Polygon([(0, 0), (-bw / 2, -bh), (bw / 2, -bh)],
                         closed=True, fc="#cbd5e1", ec="#1e293b", lw=1.2, zorder=4))
    ax.plot([-bw * 0.85, bw * 0.85], [-bh, -bh], color="#1e293b", lw=1.5, zorder=4)
    ax.plot([0.0], [H], marker="o", ms=5, mfc="white", mec="#1e293b", mew=1.4, zorder=6)

    # Height dimension on the right
    dim_x = 1.18
    ax.annotate("", xy=(dim_x, H), xytext=(dim_x, 0.0),
                arrowprops=dict(arrowstyle="<->", color="#475569", lw=1.0))
    ax.text(dim_x + 0.08, H / 2, f"H = {_fmt_num(H)} ft", rotation=90,
            ha="left", va="center", fontsize=8.5, color="#334155")

    # Optional intermediate weak-axis brace markers
    if brace_z:
        for zb in brace_z:
            if 0 < zb < H:
                ax.plot([0.0], [zb], marker="x", ms=7, mew=1.8, color="#0d9488", zorder=7)
    if bracing_note:
        ax.text(0.0, -0.17 * H, bracing_note, ha="center", va="top", fontsize=7.8,
                color="#475569", style="italic")

    # Axial load at top (offset horizontally if eccentric)
    P_txt = _case_summary(P_cases, "k")
    if P_txt:
        P_net = sum(P_cases.values())
        has_ecc = abs(ex) > 1e-9 or abs(ey) > 1e-9
        ax_x = 0.34 * (1 if (ex if abs(ex) >= abs(ey) else ey) >= 0 else -1) if has_ecc else 0.0
        col = "#dc2626" if P_net >= 0 else "#0284c7"
        if has_ecc:
            # Rigid bracket at top showing eccentricity
            ax.plot([0.0, ax_x], [H, H], color="#334155", lw=2.2, zorder=6)
            ecc_parts = []
            if abs(ex) > 1e-9:
                ecc_parts.append(f"ex={_fmt_num(ex)} in")
            if abs(ey) > 1e-9:
                ecc_parts.append(f"ey={_fmt_num(ey)} in")
            ax.text(ax_x / 2, H - 0.03 * H, ", ".join(ecc_parts),
                    ha="center", va="top", fontsize=7.8, color="#334155",
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85), zorder=7)
        if P_net >= 0:
            ax.annotate("", xy=(ax_x, H + 0.01 * H), xytext=(ax_x, H + 0.20 * H),
                        arrowprops=dict(arrowstyle="-|>,head_length=0.55,head_width=0.34", color=col, lw=2.5),
                        zorder=8)
        else:
            ax.annotate("", xy=(ax_x, H + 0.21 * H), xytext=(ax_x, H + 0.02 * H),
                        arrowprops=dict(arrowstyle="-|>,head_length=0.55,head_width=0.34", color=col, lw=2.5),
                        zorder=8)
        ax.text(0.0, H + 0.22 * H, f"P = {P_txt}", ha="center", va="bottom",
                fontsize=8.5, fontweight="bold", color=col,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=col, lw=0.8, alpha=0.95), zorder=9)

    # End moments (top and bottom)
    def _m_label(mx_d, my_d):
        lines = []
        tx = _case_summary(mx_d, "k-ft")
        ty = _case_summary(my_d, "k-ft")
        if tx:
            lines.append(f"Mx: {tx}")
        if ty:
            lines.append(f"My: {ty}")
        return "\n".join(lines)

    top_mlab = _m_label(Mx_top, My_top)
    if top_mlab:
        m_sign = sum(Mx_top.values()) or sum(My_top.values())
        _draw_moment_arc(ax, 0.0, H, 0.32, 0.065 * H, clockwise=(m_sign >= 0), color="#b45309")
        ax.text(-0.14, H - 0.09 * H, top_mlab, ha="right", va="top", fontsize=7.6,
                fontweight="semibold", color="#b45309",
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="#b45309", lw=0.7, alpha=0.92), zorder=9)

    bot_mlab = _m_label(Mx_bot, My_bot)
    if bot_mlab:
        m_sign = sum(Mx_bot.values()) or sum(My_bot.values())
        _draw_moment_arc(ax, 0.0, 0.0, 0.32, 0.065 * H, clockwise=(m_sign >= 0), color="#b45309")
        ax.text(-0.14, 0.09 * H, bot_mlab, ha="right", va="bottom", fontsize=7.6,
                fontweight="semibold", color="#b45309",
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="#b45309", lw=0.7, alpha=0.92), zorder=9)

    # Uniform lateral load(s) along height (arrows pushing from the left)
    lat_lines = []
    tx = _case_summary(w_lat_x, "plf", scale=1000.0)
    ty = _case_summary(w_lat_y, "plf", scale=1000.0)
    if tx:
        lat_lines.append(f"wx: {tx}")
    if ty:
        lat_lines.append(f"wy: {ty}")
    if lat_lines:
        col = "#7c3aed"
        x_tail = -0.82
        ax.fill([x_tail, -0.06, -0.06, x_tail], [0, 0, H, H], color=col, alpha=0.12, zorder=2)
        ax.plot([x_tail, x_tail], [0, H], color=col, lw=1.6, zorder=3)
        for za in np.linspace(0.06 * H, 0.94 * H, 7):
            ax.annotate("", xy=(-0.05, za), xytext=(x_tail, za),
                        arrowprops=dict(arrowstyle="-|>,head_length=0.35,head_width=0.20", color=col, lw=1.2),
                        zorder=4)
        ax.text(x_tail - 0.06, H / 2, "\n".join(lat_lines), rotation=90,
                ha="right", va="center", fontsize=8.0, fontweight="semibold", color=col,
                bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none", alpha=0.85), zorder=6)

    ax.set_xticks([])
    ax.set_ylabel("Height (ft)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.grid(axis="y", alpha=0.25)
    ax.set_title(title, fontsize=10, fontweight="semibold")
