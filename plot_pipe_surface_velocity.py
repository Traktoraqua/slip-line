"""Velocity of the free surface over Region III (plots/pipe_surface_velocity.png).

v_n(s) is the velocity normal to the surface (out of the ground) at the
distance s from the contact edge E, for a pipe moving at unit speed.  Lines:
the closed form in pipe_formulas.tex.  Circles: InclinedDiskVelocityField on
the OneSidedDiskSolver net.  (a) lateral push, (b) vertical penetration.
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
from slip_line.lateral import InclinedDiskVelocityField, OneSidedDiskSolver

PLOTS = Path(__file__).parent / "plots"
PLOTS.mkdir(exist_ok=True)

D, P_OVER_D, THETA = 1.0, 0.2, np.radians(15.0)
CURVES = ((0.0, "black"), (0.5, "tab:blue"), (0.9, "tab:red"))
FS = 19


def comma(v: float) -> str:
    return f"{v:g}".replace(".", "{,}")


def surface_velocity(phi: float, omega: float, theta: float, beta: float):
    """Closed form: (s/R, v_n/v0) from E (s = 0) to F (s = L), and s0/R.

    omega > 0: lateral push (a), omega < 0: vertical penetration (b).  The
    pipe moves along (cos beta, sin beta), y down, so w(t) = sin(t + beta).
    s0 is the end of the parallel zone (0 without one).
    """
    k = np.cos(omega) - np.sin(omega)  # sqrt(1 -+ r)
    t0 = 0.5 * np.pi - omega - theta
    s, v = np.empty(0), np.empty(0)
    if phi > t0:  # parallel zone next to E, t from phi down to t0
        t = np.linspace(phi, t0, 200)
        s = (
            np.cos(t + theta)
            + np.sin(t + theta)
            - np.cos(phi + theta)
            - np.sin(phi + theta)
        )
        v = np.sin(t + beta) / (np.sin(t + theta) - np.cos(t + theta))
    s0 = s[-1] if s.size else 0.0
    t_end = min(phi, t0)
    t = np.linspace(t_end, 0.0, 400)
    s = np.concatenate((s, s0 + (t_end - t) * k))
    v = np.concatenate((v, np.sin(t + beta) / k))
    return s, v, s0


def code_velocity(r: float, omega: float, beta: float):
    """(s/R, v_n/v0) at the surface nodes of the code's field, and phi."""
    disk = DiskFoundationBoundary(cu=1.0, D=D, p=P_OVER_D, r=r)
    sol = OneSidedDiskSolver(disk, 0.0, omega, alpha=THETA, n_fan=12, n_disk=60)
    sol.solve()
    field = InclinedDiskVelocityField(sol, beta)
    field.solve()
    x, vx, up = field.surface_velocity()
    s = (x - sol.B_s) / np.cos(THETA) / sol.R
    return s, vx * np.sin(THETA) + up * np.cos(THETA), sol.theta_max


def panel(ax, sign: float, beta: float, title: str, label_at) -> None:
    for r, color in CURVES:
        omega = sign * 0.5 * np.arcsin(r)
        s_c, v_c, phi = code_velocity(r, omega, beta)
        s, v, s0 = surface_velocity(phi, omega, THETA, beta)
        # Formula, with the drop to rest at F.
        ax.plot(np.append(s, s[-1]), np.append(v, 0.0), color=color, lw=2.0)
        idx = np.unique(np.r_[np.arange(0, len(s_c), max(1, len(s_c) // 8)), -1])
        ax.plot(
            s_c[idx], v_c[idx], "o", ms=8, mfc="white", mec=color, mew=1.6, zorder=3
        )
        if s0 > 0.0:
            m = s <= s0 + 1e-12
            ax.fill_between(s[m], 0.0, v[m], color=color, alpha=0.15, lw=0)
        print(
            f"{title[:3]} r={r}: largest |code - formula| = "
            f"{np.max(np.abs(v_c - np.interp(s_c, s, v))):.1e}"
        )
        xy, ha = label_at[r]
        ax.text(*xy, f"$r = {comma(r)}$", fontsize=FS - 3, ha=ha, va="center")
    fmt = FuncFormatter(lambda v, _: f"{v:g}".replace(".", ","))
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)
    ax.set_xlabel("$s/R$", fontsize=FS)
    ax.set_ylabel("$v_n/v_0$", fontsize=FS)
    ax.tick_params(labelsize=FS - 4)
    ax.set_xlim(0.0, 1.45)
    ax.set_ylim(0.0, None)
    ax.grid(alpha=0.3)
    ax.set_title(title, fontsize=FS + 1)


fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 5.4))
panel(
    ax_a,
    1.0,
    0.0,
    "(a) Sideveis skyv",
    {
        0.0: ((0.62, 0.42), "left"),
        0.5: ((0.32, 0.82), "left"),
        0.9: ((0.2, 2.0), "left"),
    },
)
panel(
    ax_b,
    -1.0,
    0.5 * np.pi,
    "(b) Vertikal nedtrykking",
    {
        0.0: ((0.96, 1.0), "left"),
        0.5: ((1.17, 0.82), "left"),
        0.9: ((1.31, 0.72), "left"),
    },
)
ax_a.legend(
    handles=[
        Line2D([], [], color="0.3", lw=2.0, label="formel"),
        Line2D(
            [],
            [],
            color="0.3",
            marker="o",
            ms=8,
            mfc="white",
            mew=1.6,
            ls="none",
            label="kode",
        ),
        mpatches.Patch(color="tab:red", alpha=0.15, label="parallellsone"),
        *[
            Line2D([], [], color=c, lw=2.0, label=f"$r = {comma(r)}$")
            for r, c in CURVES
        ],
    ],
    loc="upper right",
    fontsize=FS - 5,
    framealpha=0.9,
)
fig.tight_layout()
out = PLOTS / "pipe_surface_velocity.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print(f"saved {out}")
