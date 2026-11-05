"""Sketch of the two pipe problems in pipe_formulas.tex (plots/pipe_formulas.png).

(a) Pipe pushed sideways: only the front half touches the soil, and the surface
    in front falls away from the pipe at the angle theta.  Forces H and V.
(b) Pipe pushed down: the surfaces on both sides fall away at theta.  Force Q.

The slip-line nets are those of OneSidedDiskSolver (alpha = theta, contact from
the invert), with the interface shear of a horizontal push in (a) and of
vertical penetration in (b).
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from slip_line.boundary import DiskFoundationBoundary
from slip_line.lateral import OneSidedDiskSolver

PLOTS = Path(__file__).parent / "plots"
PLOTS.mkdir(exist_ok=True)

D, P_OVER_D, ROUGHNESS, THETA = 1.0, 0.2, 0.3, np.radians(15.0)
GAP = 0.06  # how far the pipe in (a) has moved away from the soil behind it
SOIL, PIPE, CONTACT, LOAD = "#f1e6cc", "#d9d9d9", "darkred", "darkgreen"
ALPHA_COLOR, BETA_COLOR = "royalblue", "tomato"
FS = 19


def solve_net(omega_sign: float) -> OneSidedDiskSolver:
    """Net with the given sign of omega (+1: horizontal push, -1: vertical)."""
    disk = DiskFoundationBoundary(cu=1.0, D=D, p=P_OVER_D, r=ROUGHNESS)
    omega = omega_sign * 0.5 * np.arcsin(ROUGHNESS)
    s = OneSidedDiskSolver(
        disk, phi_start=0.0, omega=omega, alpha=THETA, n_fan=6, n_disk=12
    )
    s.solve()
    return s


def draw_net(ax, s, mirror: bool) -> None:
    signs = (1.0, -1.0) if mirror else (1.0,)
    for lines, color, label in (
        (s.all_alpha_lines(), ALPHA_COLOR, r"$\alpha$-linjer"),
        (s.all_beta_lines(), BETA_COLOR, r"$\beta$-linjer"),
    ):
        for k, line in enumerate(lines):
            for sx in signs:
                ax.plot(
                    [sx * n.x for n in line],
                    [-n.y for n in line],
                    color=color,
                    lw=0.7,
                    zorder=3,
                    label=label if (k == 0 and sx > 0) else None,
                )


def net_extent(s) -> tuple[float, float]:
    nodes = [n for line in s.all_alpha_lines() + s.all_beta_lines() for n in line]
    return max(n.x for n in nodes), max(n.y for n in nodes)


def region_centres(s) -> dict[str, tuple[float, float]]:
    def centre(nodes):
        return np.mean([n.x for n in nodes]), -np.mean([n.y for n in nodes])

    flat = {
        "I": [n for row in s.region_I_mesh for n in row],
        "II": [n for row in s.region_II_nodes for n in row],
        "III": [n for row in s.region_III_nodes for n in row],
    }
    return {k: centre(v) for k, v in flat.items()}


def draw_pipe(ax, s) -> None:
    yc = -s.y_c
    ax.add_patch(
        mpatches.Circle(
            (0.0, yc), s.R, facecolor=PIPE, edgecolor="black", lw=1.8, zorder=5
        )
    )
    ax.plot(0.0, yc, "k.", ms=5, zorder=7)
    a = np.radians(140.0)
    ax.plot([0, s.R * np.cos(a)], [yc, yc + s.R * np.sin(a)], "k-", lw=1.0, zorder=7)
    ax.text(
        0.5 * s.R * np.cos(a) - 0.02,
        yc + 0.5 * s.R * np.sin(a) + 0.03,
        "$R$",
        fontsize=FS,
        zorder=8,
    )


def draw_contact(ax, s, t0: float, t1: float) -> None:
    t = np.linspace(t0, t1, 100)
    ax.plot(
        s.R * np.sin(t), -(s.y_c + s.R * np.cos(t)), color=CONTACT, lw=3.0, zorder=6
    )


def draw_phi(ax, s) -> None:
    yc, R, tm = -s.y_c, s.R, s.theta_max
    for t in (0.0, tm):
        ax.plot([0, R * np.sin(t)], [yc, yc - R * np.cos(t)], "k--", lw=0.8, zorder=7)
    rr = 0.38 * R
    ax.add_patch(
        mpatches.Arc(
            (0, yc),
            2 * rr,
            2 * rr,
            theta1=-90,
            theta2=-90 + np.degrees(tm),
            lw=1.0,
            zorder=7,
        )
    )
    a = -0.5 * np.pi + 0.5 * tm
    ax.text(
        0.55 * R * np.cos(a),
        yc + 0.55 * R * np.sin(a),
        r"$\varphi$",
        fontsize=FS,
        ha="center",
        va="center",
        zorder=8,
    )


def draw_theta(ax, x0: float, side: float) -> None:
    """theta between the horizontal and the slope at the contact point x0."""
    L = 0.45
    ax.plot([x0, x0 + side * L], [0, 0], "k--", lw=0.8, zorder=7)
    rr = 0.36
    t1, t2 = (
        (-np.degrees(THETA), 0.0) if side > 0 else (180.0, 180.0 + np.degrees(THETA))
    )
    ax.add_patch(
        mpatches.Arc((x0, 0), 2 * rr, 2 * rr, theta1=t1, theta2=t2, lw=1.0, zorder=7)
    )
    a = np.radians(0.5 * (t1 + t2))
    ax.text(
        x0 + 0.46 * np.cos(a),
        0.46 * np.sin(a),
        r"$\theta$",
        fontsize=FS,
        ha="center",
        va="center",
        zorder=8,
        bbox={"boxstyle": "round,pad=0.1", "fc": "white", "ec": "none", "alpha": 0.8},
    )


def draw_P(ax, s, x_left: float) -> None:
    """P from the contact level (dashed chord) down to the invert, inside the pipe."""
    P, x_d = s.p * s.D, -0.09
    ax.plot([x_left, s.B_s], [0, 0], "k--", lw=0.8, zorder=7)
    ax.annotate(
        "",
        xy=(x_d, -P),
        xytext=(x_d, 0),
        arrowprops={"arrowstyle": "<->", "lw": 1.0},
        zorder=8,
    )
    ax.text(x_d - 0.03, -0.5 * P, "$P$", fontsize=FS, ha="right", va="center", zorder=8)


def draw_load(ax, xy, xytext, label, text_offset) -> None:
    ax.annotate(
        "",
        xy=xy,
        xytext=xytext,
        arrowprops={"arrowstyle": "-|>", "color": LOAD, "lw": 2.2},
        zorder=9,
    )
    ax.text(
        xytext[0] + text_offset[0],
        xytext[1] + text_offset[1],
        label,
        color=LOAD,
        fontsize=FS + 1,
        ha="center",
        va="center",
    )


def label_regions(ax, s, mirror: bool) -> None:
    for name, (x, y) in region_centres(s).items():
        for sx in (1.0, -1.0) if mirror else (1.0,):
            ax.text(
                sx * x,
                y,
                name,
                fontsize=FS - 6,
                ha="center",
                va="center",
                zorder=8,
                bbox={
                    "boxstyle": "round,pad=0.1",
                    "fc": "white",
                    "ec": "0.6",
                    "alpha": 0.9,
                },
            )


def soil_patch(ax, xs, ys) -> None:
    ax.add_patch(
        mpatches.Polygon(
            np.column_stack([xs, ys]), closed=True, fc=SOIL, ec="none", zorder=1
        )
    )


def slope_y(x: float, x0: float) -> float:
    return -abs(x - x0) * np.tan(THETA)


def panel_push(ax) -> None:
    s = solve_net(+1.0)
    R, B, tm = s.R, s.B_s, s.theta_max
    yc = -s.y_c
    x_r, depth = net_extent(s)
    x_r += 0.05
    x_l = -(B + GAP + 0.75)
    y_bot = -(depth + 0.08)

    # Soil: flat surface behind, the trench left by the pipe, the front contact, the slope.
    t_back = np.linspace(-tm, 0.0, 60)
    t_front = np.linspace(0.0, tm, 60)
    xs = [x_l, -GAP - B]
    ys = [0.0, 0.0]
    xs += list(-GAP + R * np.sin(t_back))
    ys += list(-(s.y_c + R * np.cos(t_back)))
    xs += list(R * np.sin(t_front))
    ys += list(-(s.y_c + R * np.cos(t_front)))
    xs += [x_r, x_r, x_l]
    ys += [slope_y(x_r, B), y_bot, y_bot]
    soil_patch(ax, xs, ys)
    ax.plot([x_l, -GAP - B], [0, 0], color="saddlebrown", lw=1.2, zorder=2)
    ax.plot(
        -GAP + R * np.sin(t_back),
        -(s.y_c + R * np.cos(t_back)),
        color="saddlebrown",
        lw=1.2,
        zorder=2,
    )
    ax.plot([B, x_r], [0, slope_y(x_r, B)], color="saddlebrown", lw=1.2, zorder=2)

    draw_net(ax, s, mirror=False)
    draw_pipe(ax, s)
    draw_contact(ax, s, 0.0, tm)
    draw_phi(ax, s)
    draw_theta(ax, B, +1.0)
    draw_P(ax, s, -B)
    label_regions(ax, s, mirror=False)
    draw_load(ax, (-R, yc), (-R - 0.35, yc), "$H$", (0.0, 0.07))
    draw_load(ax, (0.0, yc + R), (0.0, yc + R + 0.25), "$V$", (0.07, 0.0))
    t_gap = -0.6 * tm
    gap_pt = (-0.5 * GAP + R * np.sin(t_gap), -(s.y_c + R * np.cos(t_gap)))
    ax.annotate(
        "ingen kontakt",
        xy=gap_pt,
        xytext=(-B - GAP - 0.42, 0.16),
        fontsize=FS - 3,
        color="0.25",
        ha="center",
        arrowprops={"arrowstyle": "-", "color": "0.35", "lw": 0.8},
        zorder=8,
    )
    ax.legend(loc="upper right", fontsize=FS - 3, framealpha=0.9)
    ax.set_title("(a) Sideveis skyv", fontsize=FS + 1)
    ax.set_xlim(x_l - 0.02, x_r + 0.02)
    ax.set_ylim(y_bot - 0.02, yc + R + 0.33)


def panel_vertical(ax) -> None:
    s = solve_net(-1.0)
    R, B, tm = s.R, s.B_s, s.theta_max
    yc = -s.y_c
    x_r, depth = net_extent(s)
    x_r += 0.05
    y_bot = -(depth + 0.08)

    t = np.linspace(-tm, tm, 120)
    xs = [-x_r, -B] + list(R * np.sin(t)) + [B, x_r, x_r, -x_r]
    ys = (
        [slope_y(-x_r, -B), 0.0]
        + list(-(s.y_c + R * np.cos(t)))
        + [0.0, slope_y(x_r, B), y_bot, y_bot]
    )
    soil_patch(ax, xs, ys)
    for sx in (1.0, -1.0):
        ax.plot(
            [sx * B, sx * x_r],
            [0, slope_y(x_r, B)],
            color="saddlebrown",
            lw=1.2,
            zorder=2,
        )

    draw_net(ax, s, mirror=True)
    draw_pipe(ax, s)
    draw_contact(ax, s, -tm, tm)
    draw_phi(ax, s)
    draw_theta(ax, B, +1.0)
    draw_P(ax, s, -B)
    label_regions(ax, s, mirror=True)
    draw_load(ax, (0.0, yc + R), (0.0, yc + R + 0.25), "$Q$", (0.07, 0.0))
    ax.set_title("(b) Vertikal nedtrykking", fontsize=FS + 1)
    ax.set_xlim(-x_r - 0.02, x_r + 0.02)
    ax.set_ylim(y_bot - 0.02, yc + R + 0.33)


fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5.0))
panel_push(ax_a)
panel_vertical(ax_b)
for ax in (ax_a, ax_b):
    ax.set_aspect("equal")
    ax.axis("off")
fig.tight_layout()
out = PLOTS / "pipe_formulas.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print(f"saved {out}")
