"""The time step of the dynamic model (plots/projection.png).

(a) The closest-point projection (slip_line/dynamics.py) on a schematic
yield envelope, V = W = 2 and M_x = M_z, so the distance is Euclidean: the
trial force F* = F_ext + M·v/Δt and the soil force F, the point of the
envelope closest to F*.  The pipe moves along F* − F, the outward normal
(normality): at rest inside, along the edge normal, between the normals at
a corner, and straight up at the origin (lift-off).  (b) Newmark's sliding
block as a check: the rectangle |H| ≤ H0 = 1 kN/m, 0 ≤ V ≤ 5 kN/m, W = 2,
M = 1000 kg/m, Fx = 1.5·H0·sin(2πt/T), T = 1 s, Δt = T/2000, against the
closed form of manual.tex.  Prints the numbers used there.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from slip_line.dynamics import PipeDynamics, closest_point
from slip_line.envelope import YieldEnvelope

PLOTS = Path(__file__).parent / "plots"
PLOTS.mkdir(exist_ok=True)

BLUE, LIGHT = "#2a6fc9", "#7fb0e8"
INK, INK2, GRID, GREEN = "#0b0b0b", "#52514e", "#e4e3df", "#1a8a5c"
FS = 14


def style(ax, xlabel, ylabel, title=None):
    ax.set_xlabel(xlabel, color=INK2, fontsize=FS)
    ax.set_ylabel(ylabel, color=INK2, fontsize=FS)
    if title:
        ax.set_title(title, fontsize=FS, loc="left", color=INK)
    ax.grid(color=GRID, lw=0.6)
    ax.tick_params(colors=INK2, labelsize=FS - 2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


# ----------------------------------------------------------------------
# (a) Projection
# ----------------------------------------------------------------------

W = 2.0
M = (1.0, 1.0)
env = YieldEnvelope(
    [
        [0.0, 0.0],
        [0.55, 0.45],
        [0.85, 1.3],
        [0.8, 2.3],
        [0.45, 3.0],
        [-0.45, 3.0],
        [-0.8, 2.3],
        [-0.85, 1.3],
        [-0.55, 0.45],
    ]
)
k = int(np.argmin(np.hypot(*(env.vertices - [0.85, 1.3]).T)))  # corner
n1, n2 = env.edge_delta[k - 1], env.edge_delta[k]
mid = np.arctan2(np.sin(n1) + np.sin(n2), np.cos(n1) + np.cos(n2))
cases = [
    ("rest", [0.3, W], (8, -18)),
    ("edge", [1.3, W], (8, 6)),
    ("corner", env.vertices[k] + 0.5 * np.array([np.cos(mid), np.sin(mid)]), (-45, 10)),
    ("lift-off", [0.15, -0.5], (-40, 9)),
]

fig = plt.figure(figsize=(14, 6.2))
gs = fig.add_gridspec(2, 2, width_ratios=(1.0, 1.15), hspace=0.12)
ax = fig.add_subplot(gs[:, 0])
P = np.vstack((env.vertices, env.vertices[:1]))
ax.fill(P[:, 1], P[:, 0], color=LIGHT, alpha=0.18, lw=0, zorder=1)
ax.plot(P[:, 1], P[:, 0], color=BLUE, lw=2.0, zorder=2, label="envelope")
print("(a) projection, M_x = M_z:")
for name, F_star, offset in cases:
    F_star = np.asarray(F_star, dtype=float)
    F = closest_point(env, F_star, M)
    ax.plot(F[1], F[0], "o", ms=8, color=INK, zorder=5)
    if np.allclose(F, F_star):
        print(f"  {name:8s} F* = ({F_star[0]:5.2f}, {F_star[1]:5.2f}) inside: v = 0")
        ax.annotate(
            "F = F*  (rest)",
            (F[1], F[0]),
            xytext=offset,
            textcoords="offset points",
            fontsize=FS - 2,
            color=INK,
        )
        continue
    ax.plot(F_star[1], F_star[0], "o", ms=8, mfc="white", mew=1.6, color=INK, zorder=5)
    ax.annotate(
        "",
        xy=(F_star[1], F_star[0]),
        xytext=(F[1], F[0]),
        arrowprops={"arrowstyle": "-|>", "color": INK, "lw": 1.4, "mutation_scale": 14},
        zorder=4,
    )
    ax.annotate(
        f"F* ({name})",
        (F_star[1], F_star[0]),
        xytext=offset,
        textcoords="offset points",
        fontsize=FS - 2,
        color=INK,
    )
    lo, hi = env.directions(F)
    d = np.degrees(np.arctan2(*(F_star - F)[::-1]))
    print(
        f"  {name:8s} F* = ({F_star[0]:5.2f}, {F_star[1]:5.2f})"
        f" -> F = ({F[0]:5.2f}, {F[1]:5.2f}), motion {d:6.1f} deg"
        f" (normals {np.degrees(lo):6.1f} to {np.degrees(hi):6.1f})"
    )
ax.axvline(W, color=GREEN, lw=1.2, ls=":", zorder=0)
ax.axvline(0.0, color=INK2, lw=0.7, zorder=0)
ax.axhline(0.0, color=INK2, lw=0.7, zorder=0)
ax.text(W + 0.05, -1.2, "V = W", color=GREEN, fontsize=FS - 2)
ax.plot([], [], "o", ms=8, mfc="white", mew=1.6, color=INK, label="trial force F*")
ax.plot([], [], "o", ms=8, color=INK, label="soil force F")
ax.plot([], [], color=INK, lw=1.4, label="motion along F* - F")
ax.set_aspect("equal")
ax.set_anchor("N")
ax.set_xlim(-0.8, 3.4)
ax.set_ylim(-1.3, 1.9)
ax.legend(loc="upper left", fontsize=FS - 3, frameon=False)
style(ax, "V [kN/m], downward", "H [kN/m], toward +x", "(a) Closest-point projection")

# ----------------------------------------------------------------------
# (b) Sliding block
# ----------------------------------------------------------------------

H0, F0, T, Mx = 1.0, 1.5, 1.0, 1.0  # kN/m, kN/m, s, t/m
omega = 2 * np.pi / T
box = YieldEnvelope([[-H0, 0.0], [H0, 0.0], [H0, 5.0], [-H0, 5.0]])
dyn = PipeDynamics(
    box, W=W, M=(1e3 * Mx, 1e3 * Mx), Fx=lambda t: F0 * np.sin(omega * t)
)
states = dyn.run(2.0 * T, T / 2000)
t = np.array([s.t for s in states])
u = np.array([s.u for s in states])

# Closed form (manual.tex): start t1, stop t2 (bisection), u while sliding.
t1 = np.arcsin(H0 / F0) / omega
c1, s1 = np.cos(omega * t1), np.sin(omega * t1)


def v_exact(t):
    return (F0 / omega * (c1 - np.cos(omega * t)) - H0 * (t - t1)) / Mx


def u_exact(t):
    tau = t - t1
    return (
        F0 / omega * (tau * c1 - (np.sin(omega * t) - s1) / omega) - 0.5 * H0 * tau**2
    ) / Mx


a, b = (np.pi - omega * t1) / omega, t1 + T
for _ in range(100):
    c = 0.5 * (a + b)
    a, b = (c, b) if v_exact(c) > 0.0 else (a, c)
t2 = 0.5 * (a + b)
du = u_exact(t2)
u_num = u[np.searchsorted(t, 0.5 * T + t1) - 1]
print(f"(b) sliding block, F0/H0 = {F0 / H0:g}:")
print(f"  t1 = {t1 / T:.4f} T, t2 = {t2 / T:.4f} T")
print(f"  du closed form = {1e3 * du:.3f} mm, numerical = {1e3 * u_num:.3f} mm")
print(f"  relative error = {(u_num - du) / du:.1e}")

ax_f = fig.add_subplot(gs[0, 1])
ax_f.plot(t / T, dyn.Fx(t) / H0, color=INK2, lw=2.0)
for y in (-1.0, 1.0):
    ax_f.axhline(y, color=GREEN, lw=1.2, ls=":", zorder=0)
ax_f.text(
    2.0, 1.08, r"$\pm H_0$", color=GREEN, fontsize=FS - 2, ha="right", va="bottom"
)
ax_f.set_ylim(-1.7, 1.7)
ax_f.tick_params(labelbottom=False)
style(ax_f, "", r"$F_x / H_0$", r"(b) Sliding block, $F_0 = 1.5\,H_0$")

ax_u = fig.add_subplot(gs[1, 1], sharex=ax_f)
ax_u.plot(t / T, 1e3 * u, color=BLUE, lw=2.0, label="PipeDynamics, Δt = T/2000")
for start, sign in ((0.0, 1.0), (0.5, -1.0), (1.0, 1.0), (1.5, -1.0)):
    tt = np.linspace(t1, t2, 60)
    base = 0.0 if sign > 0 else du
    ax_u.plot(
        start + tt / T,
        1e3 * (base + sign * u_exact(tt)),
        "o",
        ms=3.5,
        color=INK,
        markevery=4,
        label="closed form (Newmark)" if start == 0.0 else None,
    )
ax_u.legend(loc="upper center", ncol=2, fontsize=FS - 3, frameon=False)
ax_u.set_ylim(-1.0, 29.0)
ax_u.set_xlim(0.0, 2.0)
style(ax_u, "t / T", "u [mm]")

fig.subplots_adjust(left=0.06, right=0.98, wspace=0.22)
out = PLOTS / "projection.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print("wrote", out)
