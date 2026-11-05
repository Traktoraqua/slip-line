"""Pipe pushed down from P = 0 to P = PR_END*R with heave, using the closed
forms (PipePenetrationClosedForm).  Writes plots/penetration.gif (the pipe and
the soil surface, and the vertical force against penetration) and its last
frame as plots/penetration.png.

P is the settlement of the invert below the original surface; the load
depends on the effective embedment d below the contact level, which heave
makes larger than P.  The closed forms need d <= R, so PR_END must stay well
below 1.
"""

import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.ticker import FuncFormatter

from slip_line.penetration import PenetrationState, PipePenetrationClosedForm

D = 1.0  # pipe diameter [m]
CU = 1.0  # undrained shear strength [kPa]
ROUGHNESS = 0.0  # r = tau/cu at the pipe
PR_END = 0.5  # final P/R

PLOTS = Path(__file__).parent / "plots"
PLOTS.mkdir(exist_ok=True)
SOIL, PIPE, LINE = "#f1e6cc", "#d9d9d9", "saddlebrown"
FS = 15

# ---- analysis (displacement controlled: step until P = w reaches PR_END*R)
sim = PipePenetrationClosedForm(D=D, cu=CU, r=ROUGHNESS, W=np.inf, x_max=2.0 * D)
R = sim.R
state = PenetrationState(0.0, 0.0, 0.0, 0.0, 0.0, np.zeros_like(sim.x), 1.0)
states = [state]
t = time.perf_counter()
for _ in range(round(PR_END * R / sim.dw)):
    state = sim.advance(state, sim.dw)
    states.append(state)
print(f"{len(states) - 1} steps in {time.perf_counter() - t:.2f} s")

w = np.array([s.w for s in states])
V = np.array([s.V for s in states]) / (CU * D)
V_flat = np.array([sim.load(x) for x in w]) / (CU * D)  # no heave: d = w
end = states[-1]
print(
    f"P/R = {end.w / R:.3f}: d/R = {end.d / R:.3f}, "
    f"heave at the pipe = {1e3 * (end.d - end.w):.1f} mm, "
    f"largest heave = {1e3 * end.heave.max():.1f} mm, "
    f"plastic zone to x = {1e3 * end.x_end:.0f} mm"
)
print(f"V/(cu D) = {V[-1]:.3f} (no heave: {V_flat[-1]:.3f})")

# ---- animation
comma = FuncFormatter(lambda v, _: f"{v:g}".replace(".", ","))
X = max(1.3 * D, 1.05 * end.x_end)
Y_BOT, Y_TOP = -0.62 * D, 1.05 * D
last = len(states) - 1
frames = list(np.linspace(0, last, 121).astype(int)) + [last] * 25

fig, (ax_s, ax_v) = plt.subplots(
    1, 2, figsize=(13, 5.0), gridspec_kw={"width_ratios": [1.45, 1]}
)

# Force panel: the curve without heave is static, the rest moves.
ax_v.plot(w / R, V_flat, "--", color="0.55", lw=1.6, label="uten heving ($d = P$)")
(trace,) = ax_v.plot([], [], color="tab:blue", lw=2.2, label="med heving")
(dot,) = ax_v.plot([], [], "o", color="tab:blue", ms=8, mec="white", mew=1.5)
ax_v.set_xlim(0.0, PR_END)
ax_v.set_ylim(0.0, 1.08 * max(V.max(), V_flat.max()))
ax_v.set_xlabel("$P/R$", fontsize=FS)
ax_v.set_ylabel("$V/(c_u D)$", fontsize=FS)
ax_v.xaxis.set_major_formatter(comma)
ax_v.yaxis.set_major_formatter(comma)
ax_v.tick_params(labelsize=FS - 3)
ax_v.grid(alpha=0.3)
ax_v.legend(loc="lower right", fontsize=FS - 3, framealpha=0.9)
ax_v.set_title("Vertikalkraft mot nedtrykking", fontsize=FS + 1)


def surface(s: PenetrationState) -> tuple[np.ndarray, np.ndarray]:
    """Soil top (x, y up) across the full width: the pipe bottom inside the
    contact, the heave outside."""
    yc = R - s.w
    xr = sim.x[sim.x <= X]
    yr = np.where(
        xr < s.x_c, yc - np.sqrt(np.maximum(R * R - xr * xr, 0.0)), s.heave[: len(xr)]
    )
    return np.concatenate((-xr[::-1], xr[1:])), np.concatenate((yr[::-1], yr[1:]))


def draw(k: int) -> list:
    s = states[k]
    ax_s.clear()
    xs, ys = surface(s)
    ax_s.fill(np.r_[-X, xs, X], np.r_[Y_BOT, ys, Y_BOT], fc=SOIL, ec="none", zorder=1)
    ax_s.plot(xs, ys, color=LINE, lw=1.6, zorder=2)
    ax_s.axhline(0.0, color="0.45", ls="--", lw=1.0, zorder=3)
    ax_s.text(
        -X + 0.03 * D, 0.02 * D, "opprinnelig overflate", fontsize=FS - 5, color="0.35"
    )
    ax_s.add_patch(
        mpatches.Circle((0.0, R - s.w), R, fc=PIPE, ec="black", lw=1.6, zorder=4)
    )
    ax_s.plot(0.0, R - s.w, "k.", ms=4, zorder=5)
    ax_s.text(
        -X + 0.03 * D,
        Y_BOT + 0.04 * D,
        f"$P/R$ = {s.w / R:.3f}\n"
        f"$d/R$ = {s.d / R:.3f}\n"
        f"heving ved røret = {1e3 * (s.d - s.w):.0f} mm\n"
        f"$V/(c_u D)$ = {s.V / (CU * D):.2f}".replace(".", ","),
        fontsize=FS - 3,
        va="bottom",
        bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "0.7", "alpha": 0.9},
        zorder=6,
    )
    ax_s.set_xlim(-X, X)
    ax_s.set_ylim(Y_BOT, Y_TOP)
    ax_s.set_aspect("equal")
    ax_s.axis("off")
    kind = "Glatt rør" if ROUGHNESS == 0.0 else "Ru rør"
    ax_s.set_title(
        f"{kind}, $D$ = {D:g} m, $r$ = {ROUGHNESS:g}".replace(".", ","),
        fontsize=FS + 1,
    )
    trace.set_data(w[: k + 1] / R, V[: k + 1])
    dot.set_data([w[k] / R], [V[k]])
    return []


fig.tight_layout()
out = PLOTS / "penetration.gif"
FuncAnimation(fig, draw, frames=frames, blit=False).save(
    out, writer=PillowWriter(fps=20), dpi=80
)
draw(last)
fig.savefig(out.with_suffix(".png"), dpi=110, bbox_inches="tight")
print(f"saved {out} and {out.with_suffix('.png').name}")
