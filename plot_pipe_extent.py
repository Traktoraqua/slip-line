"""How far Region III reaches from the pipe (plots/pipe_extent.png).

(a) The net of a pipe pushed sideways (OneSidedDiskSolver, as in
    pipe_formulas.png) with the closed-form outer alpha-line drawn on top,
    the straight beta-line l from the contact edge E and the length L of
    Region III along the slope.
(b) L/R against P/R from the formulas in pipe_formulas.tex.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

from slip_line.boundary import DiskFoundationBoundary
from slip_line.lateral import OneSidedDiskSolver

PLOTS = Path(__file__).parent / "plots"
PLOTS.mkdir(exist_ok=True)

D, P_OVER_D, ROUGHNESS, THETA = 1.0, 0.2, 0.3, np.radians(15.0)
SOIL, PIPE, CONTACT = "#f1e6cc", "#d9d9d9", "darkred"
ALPHA_COLOR, BETA_COLOR = "royalblue", "tomato"
FS = 19


def extent(phi: float, omega: float, theta: float) -> float:
    """L/R: Region III along the surface, from the contact edge.

    omega > 0: lateral push (a), omega < 0: vertical penetration (b).
    """
    t0 = 0.5 * np.pi - omega - theta
    k = np.cos(omega) - np.sin(omega)
    if phi <= t0:
        return phi * k
    return (
        t0 * k
        + np.cos(omega)
        + np.sin(omega)
        - np.cos(phi + theta)
        - np.sin(phi + theta)
    )


def direction(a):
    """Unit vector (x, y_down) at the angle a from straight down toward +x."""
    return np.array([np.sin(a), np.cos(a)])


def outer_alpha_line(R: float, phi: float, omega: float, theta: float) -> np.ndarray:
    """Closed-form boundary of the plastic zone, from the invert to the surface.

    Region I: the straight beta-line from the arc point t makes the angle
    kappa = pi/4 - omega with the radius and has the length R*t*sin(kappa).
    Then the fan (a circle about the edge E) and Region III.
    """
    kappa = 0.25 * np.pi - omega
    yc = -R * np.cos(phi)
    t0 = 0.5 * np.pi - omega - theta
    t_end = min(phi, t0)
    pts = []
    for t in np.linspace(0.0, t_end, 80):
        p = np.array([R * np.sin(t), yc + R * np.cos(t)])
        pts.append(p + R * t * np.sin(kappa) * direction(t - kappa))
    E = np.array([R * np.sin(phi), 0.0])
    if phi <= t0:
        ell = R * phi * np.sin(kappa)
        for a in np.linspace(phi - kappa, 0.25 * np.pi - theta, 40):
            pts.append(E + ell * direction(a))
    L = R * extent(phi, omega, theta)
    pts.append(E + L * np.array([np.cos(theta), np.sin(theta)]))
    return np.array(pts)


def to_plot(xy) -> tuple:
    """(x, y_down) -> plot coordinates (y up)."""
    xy = np.asarray(xy)
    return xy[..., 0], -xy[..., 1]


def panel_geometry(ax) -> None:
    disk = DiskFoundationBoundary(cu=1.0, D=D, p=P_OVER_D, r=ROUGHNESS)
    omega = 0.5 * np.arcsin(ROUGHNESS)
    s = OneSidedDiskSolver(disk, 0.0, omega, alpha=THETA, n_fan=6, n_disk=12)
    s.solve()
    R, phi, B = s.R, s.theta_max, s.B_s
    kappa = 0.25 * np.pi - omega
    ell = R * phi * np.sin(kappa)
    L = R * extent(phi, omega, THETA)
    print(f"(a) L/R: code {s.L / R:.5f}, formula {L / R:.5f}")
    yc = R * np.cos(phi)  # pipe centre, plot coordinates
    E = np.array([B, 0.0])
    F = E + L * np.array([np.cos(THETA), -np.sin(THETA)])
    G = E + ell * np.array(
        [np.sin(0.25 * np.pi - THETA), -np.cos(0.25 * np.pi - THETA)]
    )

    # Soil and slope.
    x_l, x_r, y_bot = -0.2, F[0] + 0.12, -0.36
    t = np.linspace(-phi, phi, 120)
    xs = [x_l, -B] + list(R * np.sin(t)) + [B, x_r, x_r, x_l]
    ys = (
        [0.0, 0.0]
        + list(yc - R * np.cos(t))
        + [0.0, -(x_r - B) * np.tan(THETA), y_bot, y_bot]
    )
    ax.add_patch(
        mpatches.Polygon(np.column_stack([xs, ys]), fc=SOIL, ec="none", zorder=1)
    )
    ax.plot(
        [B, x_r], [0, -(x_r - B) * np.tan(THETA)], color="saddlebrown", lw=1.2, zorder=2
    )
    ax.add_patch(
        mpatches.Polygon(
            [E, G, F], closed=True, fc="#ffe9a8", ec="none", alpha=0.9, zorder=2
        )
    )

    # Net from the code, thin.
    for lines, color in (
        (s.all_alpha_lines(), ALPHA_COLOR),
        (s.all_beta_lines(), BETA_COLOR),
    ):
        for line in lines:
            ax.plot(
                [n.x for n in line],
                [-n.y for n in line],
                color=color,
                lw=0.6,
                alpha=0.7,
                zorder=3,
            )

    # Closed-form boundary and the beta-line from E, thick.
    x, y = to_plot(outer_alpha_line(R, phi, omega, THETA))
    ax.plot(x, y, color=ALPHA_COLOR, lw=2.6, zorder=4)
    ax.plot([E[0], G[0]], [E[1], G[1]], color=BETA_COLOR, lw=2.6, zorder=4)

    # Pipe and contact arc.
    ax.add_patch(mpatches.Circle((0, yc), R, fc=PIPE, ec="black", lw=1.8, zorder=5))
    tc = np.linspace(0.0, phi, 80)
    ax.plot(R * np.sin(tc), yc - R * np.cos(tc), color=CONTACT, lw=3.0, zorder=6)
    ax.plot(0, yc, "k.", ms=5, zorder=7)

    # kappa at E: between the radius (extended into the soil) and the beta-line.
    n_out = np.array([np.sin(phi), -np.cos(phi)])
    ax.plot([0, B + 0.13 * n_out[0]], [yc, 0.13 * n_out[1]], "k--", lw=0.8, zorder=7)
    a1 = np.degrees(np.arctan2(n_out[1], n_out[0]))
    a2 = np.degrees(np.arctan2(G[1] - E[1], G[0] - E[0]))
    rr = 0.09
    ax.add_patch(
        mpatches.Arc(E, 2 * rr, 2 * rr, theta1=a2, theta2=a1, lw=1.2, zorder=8)
    )
    am = np.radians(0.5 * (a1 + a2))
    ax.text(
        E[0] + 1.45 * rr * np.cos(am),
        E[1] + 1.45 * rr * np.sin(am),
        r"$\kappa$",
        fontsize=FS,
        ha="center",
        va="center",
        zorder=8,
        bbox={"boxstyle": "round,pad=0.05", "fc": "white", "ec": "none", "alpha": 0.8},
    )

    # ell along the beta-line (on the Region I side), L along the slope.
    mid = 0.5 * (E + G) + 0.035 * np.array([-np.cos(phi - kappa), -np.sin(phi - kappa)])
    ax.text(
        mid[0],
        mid[1],
        r"$\ell$",
        fontsize=FS,
        ha="center",
        va="center",
        zorder=8,
        bbox={"boxstyle": "round,pad=0.05", "fc": "white", "ec": "none", "alpha": 0.8},
    )
    up = 0.08 * np.array([np.sin(THETA), np.cos(THETA)])
    for p in (E, F):
        ax.plot(
            [p[0], p[0] + 1.15 * up[0]],
            [p[1], p[1] + 1.15 * up[1]],
            "k-",
            lw=0.8,
            zorder=8,
        )
    ax.annotate(
        "",
        xy=F + up,
        xytext=E + up,
        arrowprops={"arrowstyle": "<->", "lw": 1.2},
        zorder=8,
    )
    m = 0.5 * (E + F) + 1.55 * up
    ax.text(m[0], m[1], "$L$", fontsize=FS + 1, ha="center", va="center", zorder=8)
    for p, name, off in ((E, "E", (0.015, 0.015)), (F, "F", (0.02, -0.05))):
        ax.plot(*p, "ko", ms=5, zorder=9)
        ax.text(p[0] + off[0], p[1] + off[1], name, fontsize=FS - 2, zorder=9)

    # Region names (II is the narrow fan at E, named in pipe_formulas.png).
    for name, nodes in (("I", [n for row in s.region_I_mesh for n in row]),):
        ax.text(
            np.mean([n.x for n in nodes]),
            -np.mean([n.y for n in nodes]),
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
    c3 = (E + G + F) / 3
    ax.text(
        c3[0],
        c3[1],
        "III",
        fontsize=FS - 6,
        ha="center",
        va="center",
        zorder=8,
        bbox={"boxstyle": "round,pad=0.1", "fc": "white", "ec": "0.6", "alpha": 0.9},
    )

    ax.legend(
        handles=[
            Line2D([], [], color=ALPHA_COLOR, lw=2.6, label=r"ytre $\alpha$-linje"),
            Line2D([], [], color=BETA_COLOR, lw=2.6, label=r"$\beta$-linje fra E"),
        ],
        loc="upper right",
        fontsize=FS - 4,
        framealpha=0.9,
    )
    ax.set_title("(a) Sone III ved sideveis skyv", fontsize=FS + 1)
    ax.set_xlim(x_l, x_r)
    ax.set_ylim(y_bot, yc + 0.06)
    ax.set_aspect("equal")
    ax.axis("off")


def panel_curves(ax) -> None:
    pr = np.linspace(0.0, 1.0, 401)
    phis = np.arccos(1.0 - pr)
    for r, color in ((0.0, "black"), (0.5, "tab:blue"), (1.0, "tab:red")):
        omega = 0.5 * np.arcsin(r)
        ok = phis <= 0.5 * np.pi + omega - THETA + 1e-12
        for sign, ls in ((1.0, "-"), (-1.0, "--")):
            L = np.array([extent(p, sign * omega, THETA) for p in phis[ok]])
            ax.plot(pr[ok], L, color=color, ls=ls, lw=2.0)
        t0 = 0.5 * np.pi - omega - THETA
        if t0 < phis[ok][-1]:  # the parallel zone starts inside the range
            ax.plot(1.0 - np.cos(t0), extent(t0, omega, THETA), "o", color=color, ms=7)
    comma = FuncFormatter(lambda v, _: f"{v:g}".replace(".", ","))
    ax.xaxis.set_major_formatter(comma)
    ax.yaxis.set_major_formatter(comma)
    ax.set_xlabel("$P/R$", fontsize=FS)
    ax.set_ylabel("$L/R$", fontsize=FS)
    ax.tick_params(labelsize=FS - 4)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, None)
    ax.grid(alpha=0.3)
    ax.legend(
        handles=[
            Line2D([], [], color="black", lw=2.0, label="$r = 0$"),
            Line2D([], [], color="tab:blue", lw=2.0, label="$r = 0{,}5$"),
            Line2D([], [], color="tab:red", lw=2.0, label="$r = 1$"),
            Line2D([], [], color="0.4", lw=2.0, ls="-", label="(a) sideveis"),
            Line2D([], [], color="0.4", lw=2.0, ls="--", label="(b) vertikal"),
            Line2D(
                [], [], color="0.4", marker="o", ls="none", ms=7, label="parallellsone"
            ),
        ],
        loc="upper left",
        fontsize=FS - 5,
        framealpha=0.9,
    )
    ax.set_title(r"(b) Lengden $L$ ved $\theta$ = 15°", fontsize=FS + 1)


fig, (ax_a, ax_b) = plt.subplots(
    1, 2, figsize=(14, 5.6), gridspec_kw={"width_ratios": [1.35, 1]}
)
panel_geometry(ax_a)
panel_curves(ax_b)
fig.tight_layout()
out = PLOTS / "pipe_extent.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print(f"saved {out}")
