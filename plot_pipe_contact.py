"""The general lateral push in pipe_formulas.tex: contact from lambda to phi
(plots/pipe_contact.png).

(a) lambda < 0: the contact goes past the invert onto the back of the pipe.
(b) lambda > 0: the contact ends before the invert.

The nets are those of OneSidedDiskSolver(phi_start=lambda) with the interface
shear of a horizontal push.  The thick line is the outer alpha-line from the
closed forms: it starts at t = lambda and ends at F,
L = R(phi - lambda)sqrt(1 - r) from the contact edge E.
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

D, P_OVER_D, ROUGHNESS, THETA = 1.0, 0.25, 0.3, np.radians(15.0)
LAMBDAS = (np.radians(-25.0), np.radians(20.0))
GAP = 0.06  # radial gap between the pipe and the soil at t = -phi
SOIL, PIPE, CONTACT, LOAD = "#f1e6cc", "#d9d9d9", "darkred", "darkgreen"
ALPHA_COLOR, BETA_COLOR, OUTER = "royalblue", "tomato", "navy"
FS = 19


def on_pipe(s, t, gap=0.0):
    """Screen (x, y up) of the pipe point t, moved gap out along the radius."""
    return (s.R + gap) * np.sin(t), -(s.y_c + (s.R + gap) * np.cos(t))


def outer_alpha_line(s, lam: float) -> np.ndarray:
    """Closed form, screen coordinates: Region I from t = lam, the fan about E,
    then Region III to F (no parallel zone)."""
    R, phi, kappa = s.R, s.theta_max, 0.25 * np.pi - s.omega_signed

    def direction(b):
        return np.array([np.sin(b), np.cos(b)])  # (x, y down)

    pts = []
    for t in np.linspace(lam, phi, 80):
        p = np.array([R * np.sin(t), s.y_c + R * np.cos(t)])
        pts.append(p + R * (t - lam) * np.sin(kappa) * direction(t - kappa))
    E, ell = np.array([s.B_s, 0.0]), R * (phi - lam) * np.sin(kappa)
    for b in np.linspace(phi - kappa, 0.25 * np.pi - THETA, 40):
        pts.append(E + ell * direction(b))
    pts.append(E + np.sqrt(2.0) * ell * np.array([np.cos(THETA), np.sin(THETA)]))
    pts = np.array(pts)
    return np.column_stack((pts[:, 0], -pts[:, 1]))


def angle_arc(
    ax, s, t1: float, radius: float, label: str, label_r: float, at: float = 0.5
) -> None:
    """Arc at the pipe centre from the invert (t = 0) to t1, with a label at
    the fraction at of the arc."""
    yc = -s.y_c
    lo, hi = sorted((-90.0, -90.0 + np.degrees(t1)))
    ax.add_patch(
        mpatches.Arc(
            (0, yc), 2 * radius, 2 * radius, theta1=lo, theta2=hi, lw=1.0, zorder=7
        )
    )
    b = -0.5 * np.pi + at * t1
    ax.text(
        label_r * np.cos(b),
        yc + label_r * np.sin(b),
        label,
        fontsize=FS - 1,
        ha="center",
        va="center",
        zorder=8,
    )


def solve(lam: float) -> OneSidedDiskSolver:
    disk = DiskFoundationBoundary(cu=1.0, D=D, p=P_OVER_D, r=ROUGHNESS)
    s = OneSidedDiskSolver(
        disk, lam, 0.5 * np.arcsin(ROUGHNESS), alpha=THETA, n_fan=6, n_disk=12
    )
    s.solve()
    assert s.complete() and not s.use_parallel_zone
    return s


def net_extent(s) -> tuple[float, float]:
    nodes = [n for line in s.all_alpha_lines() + s.all_beta_lines() for n in line]
    return max(n.x for n in nodes), max(n.y for n in nodes)


def panel(ax, s, title: str, first: bool, x_l: float, x_r: float, y_bot: float):
    """One panel; the first one also gets the legend and the region labels."""
    lam, R, B, phi, yc = s.phi_start, s.R, s.B_s, s.theta_max, -s.y_c

    # Soil: flat behind, a gap that closes at t = lam, the contact, the slope.
    t_gap = np.linspace(-phi - 0.3, lam, 300)
    gx, gy = on_pipe(s, t_gap, GAP * (lam - t_gap) / (lam + phi))
    below = gy <= 0.0  # the part of the gap line under the flat surface
    t_gap, gx, gy = t_gap[below], gx[below], gy[below]
    cx, cy = on_pipe(s, np.linspace(lam, phi, 60))
    y_r = -(x_r - B) * np.tan(THETA)
    xs = np.r_[x_l, gx, cx, x_r, x_r, x_l]
    ys = np.r_[0.0, gy, cy, y_r, y_bot, y_bot]
    ax.add_patch(
        mpatches.Polygon(np.column_stack([xs, ys]), fc=SOIL, ec="none", zorder=1)
    )
    edge = {"color": "saddlebrown", "lw": 1.2, "zorder": 2}
    ax.plot(np.r_[x_l, gx], np.r_[0.0, gy], **edge)
    ax.plot([B, x_r], [0.0, y_r], **edge)

    # Net, the closed-form outer alpha-line, the pipe and the contact.
    for lines, color, label in (
        (s.all_alpha_lines(), ALPHA_COLOR, r"$\alpha$-linjer"),
        (s.all_beta_lines(), BETA_COLOR, r"$\beta$-linjer"),
    ):
        for k, line in enumerate(lines):
            ax.plot(
                [n.x for n in line],
                [-n.y for n in line],
                color=color,
                lw=0.7,
                zorder=3,
                label=label if k == 0 else None,
            )
    o = outer_alpha_line(s, lam)
    ax.plot(
        o[:, 0], o[:, 1], color=OUTER, lw=2.6, zorder=4, label=r"ytre $\alpha$-linje"
    )
    ax.add_patch(mpatches.Circle((0, yc), R, fc=PIPE, ec="black", lw=1.8, zorder=5))
    ax.plot(*on_pipe(s, np.linspace(lam, phi, 100)), color=CONTACT, lw=3.0, zorder=6)

    # Angles from the invert: phi to E, lambda to the start of the contact.
    ax.plot(0.0, yc, "k.", ms=5, zorder=7)
    for t in (0.0, lam, phi):
        ax.plot(*np.array([[0.0, yc], on_pipe(s, t)]).T, "k--", lw=0.8, zorder=7)
    if lam < 0.0:
        angle_arc(ax, s, phi, 0.40 * R, r"$\varphi$", 0.56 * R)
        angle_arc(ax, s, lam, 0.40 * R, r"$\lambda$", 0.56 * R)
    else:
        angle_arc(ax, s, phi, 0.58 * R, r"$\varphi$", 0.72 * R, at=0.68)
        angle_arc(ax, s, lam, 0.36 * R, r"$\lambda$", 0.50 * R)
    ang = np.radians(140.0)
    ax.plot([0, R * np.cos(ang)], [yc, yc + R * np.sin(ang)], "k-", lw=1.0, zorder=7)
    ax.text(
        0.5 * R * np.cos(ang) - 0.02,
        yc + 0.5 * R * np.sin(ang) + 0.03,
        "$R$",
        fontsize=FS,
        zorder=8,
    )

    # theta at E, the points E and F, and the region labels.
    ax.plot([B, B + 0.42], [0, 0], "k--", lw=0.8, zorder=7)
    ax.add_patch(
        mpatches.Arc(
            (B, 0), 0.68, 0.68, theta1=-np.degrees(THETA), theta2=0.0, lw=1.0, zorder=7
        )
    )
    b = -0.5 * THETA
    ax.text(
        B + 0.44 * np.cos(b),
        0.44 * np.sin(b),
        r"$\theta$",
        fontsize=FS - 1,
        ha="center",
        va="center",
        zorder=8,
    )
    F = o[-1]
    ax.plot(*np.array([[B, 0.0], F]).T, "o", color="black", ms=6, zorder=9)
    ax.text(B - 0.03, 0.05, "E", fontsize=FS, ha="center", va="bottom", zorder=9)
    ax.text(F[0] + 0.02, F[1] - 0.04, "F", fontsize=FS, ha="left", va="top", zorder=9)
    regions = (("I", s.region_I_mesh), ("II", s.region_II_nodes))
    regions += (("III", s.region_III_nodes),)
    for name, rows in regions if first else ():  # too cramped in panel (b)
        flat = [n for row in rows for n in row]
        ax.text(
            np.mean([n.x for n in flat]),
            -np.mean([n.y for n in flat]),
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

    # Loads on the pipe and the gap behind the contact.
    for xy, xytext, label, off in (
        ((-R, yc), (-R - 0.33, yc), "$H$", (0.0, 0.07)),
        ((0.0, yc + R), (0.0, yc + R + 0.25), "$V$", (0.07, 0.0)),
    ):
        ax.annotate(
            "",
            xy=xy,
            xytext=xytext,
            arrowprops={"arrowstyle": "-|>", "color": LOAD, "lw": 2.2},
            zorder=9,
        )
        ax.text(
            xytext[0] + off[0],
            xytext[1] + off[1],
            label,
            color=LOAD,
            fontsize=FS + 1,
            ha="center",
            va="center",
        )
    k = len(t_gap) // 3
    ax.annotate(
        "ingen kontakt",
        xy=on_pipe(s, t_gap[k], 0.5 * GAP * (lam - t_gap[k]) / (lam + phi)),
        xytext=(x_l + 0.25, y_bot + 0.12),
        fontsize=FS - 3,
        color="0.25",
        ha="center",
        arrowprops={"arrowstyle": "-", "color": "0.35", "lw": 0.8},
        zorder=8,
    )
    if first:
        ax.legend(loc="upper right", fontsize=FS - 4, framealpha=0.9)
    ax.set_title(title, fontsize=FS + 1)
    ax.set_xlim(x_l - 0.02, x_r + 0.02)
    ax.set_ylim(y_bot - 0.02, yc + R + 0.33)
    ax.set_aspect("equal")
    ax.axis("off")
    print(
        f"{title[:3]} L/R code {s.L / R:.4f}, formula {(phi - lam) * np.sqrt(1 - ROUGHNESS):.4f}"
    )


# Same scale in both panels.
nets = [solve(lam) for lam in LAMBDAS]
x_l = -nets[0].B_s - GAP - 0.55
x_r = max(max(net_extent(s)[0], s.B_s + 0.5) for s in nets) + 0.14
y_bot = -(max(net_extent(s)[1] for s in nets) + 0.08)
fig, axes = plt.subplots(1, 2, figsize=(13, 5.0))
for ax, s, title in zip(
    axes,
    nets,
    (
        r"(a) $\lambda < 0$: forbi rørbunnen",
        r"(b) $\lambda > 0$: slutter før rørbunnen",
    ),
):
    panel(ax, s, title, ax is axes[0], x_l, x_r, y_bot)
fig.tight_layout()
out = PLOTS / "pipe_contact.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print(f"saved {out}")
