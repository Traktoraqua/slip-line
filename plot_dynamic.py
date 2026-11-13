"""The pipe pushed by a load in time with the geometry following it
(PipeMotion, step 3 of PLAN_DYNAMIKK.md; plots/dynamic_push.png).

D = 1 m, r = 0.8, W = 2.8 kN/m, cu = 1 kPa, mass from W with C_a = 1.  The
load Fx = 0.3·(1 − cos(πt/T_rise)) rises to 0.6 kN/m at T_rise, slowly
(T_rise = 50 s) and fast (T_rise = 1 s).  (a) u and w against t/T_rise, and
the load.  (b) The soil force (H, V) on its way through the envelopes of the
fast run, drawn at the start and after about 8, 16 and 25 mm.  Prints the
numbers used in manual.tex.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from slip_line.dynamics import PipeMotion, pipe_mass
from slip_line.lateral import PipeLateralClosedForm
from slip_line.penetration import PipePenetrationClosedForm

PLOTS = Path(__file__).parent / "plots"
PLOTS.mkdir(exist_ok=True)

BLUES = ("#c2d8f2", "#7fb0e8", "#2a6fc9", "#0e3470")  # light to dark: later
BLUE, GREEN = "#2a6fc9", "#1a8a5c"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
FS = 14
W, D, F_MAX = 2.8, 1.0, 0.6


def style(ax, xlabel, ylabel, title=None):
    ax.set_xlabel(xlabel, color=INK2, fontsize=FS)
    ax.set_ylabel(ylabel, color=INK2, fontsize=FS)
    if title:
        ax.set_title(title, fontsize=FS, loc="left", color=INK)
    ax.grid(color=GRID, lw=0.6)
    ax.tick_params(colors=INK2, labelsize=FS - 2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


sim = PipePenetrationClosedForm(D=D, cu=1.0, r=0.8, W=W)
lat, start = PipeLateralClosedForm.after_penetration(sim, sim.run()[-1])
runs = {}
for T_rise, dt in ((50.0, 0.02), (1.0, 2e-3)):
    motion = PipeMotion(
        lat,
        start,
        pipe_mass(W, D),
        Fx=lambda t, T=T_rise: 0.5 * F_MAX * (1.0 - np.cos(np.pi * t / T)),
    )
    states = motion.run(2.0 * T_rise, dt)
    runs[T_rise] = (motion, states)
    moving = [s for s in states if s.u_dot > 0.0]
    fastest = max(states, key=lambda s: s.u_dot)
    print(
        f"T_rise = {T_rise:g} s: moves from t = {moving[0].t / T_rise:.3f} T_rise"
        f" (Fx = {motion.dyn.Fx(moving[0].t):.3f}),"
        f" fastest {1e3 * fastest.u_dot:.1f} mm/s at {fastest.t / T_rise:.3f} T_rise,"
        f" stops moving at {moving[-1].t / T_rise:.3f} T_rise"
        f" (Fx = {motion.dyn.Fx(moving[-1].t):.3f}),"
        f" u = {1e3 * states[-1].u:.2f} mm, w = {1e3 * states[-1].w:.2f} mm,"
        f" V = {min(s.V for s in moving):.3f} to {max(s.V for s in moving):.3f},"
        f" run stops ({motion.stop_reason}) at {states[-1].t / T_rise:.3f} T_rise"
    )

fig = plt.figure(figsize=(14, 6.2))
gs = fig.add_gridspec(2, 2, height_ratios=(1.0, 2.3), hspace=0.12)
ax_f = fig.add_subplot(gs[0, 0])
ax_a = fig.add_subplot(gs[1, 0], sharex=ax_f)
ax_b = fig.add_subplot(gs[:, 1])

# (a) the load, and u(t), w(t)
tt = np.linspace(0.0, 2.0, 200)
ax_f.plot(tt, 0.5 * F_MAX * (1.0 - np.cos(np.pi * tt)), color=INK2, lw=2.0)
ax_f.set_ylim(0.0, 0.7)
ax_f.tick_params(labelbottom=False)
style(ax_f, "", "Fx [kN/m]", "(a) Load and displacement, cu = 1 kPa")
for T_rise, ls, name in ((50.0, "-", "slow, 50 s"), (1.0, "--", "fast, 1 s")):
    _, states = runs[T_rise]
    t = np.array([s.t for s in states]) / T_rise
    ax_a.plot(
        t, [1e3 * s.u for s in states], ls, color=BLUE, lw=2.0, label=f"u, {name}"
    )
    ax_a.plot(
        t, [1e3 * s.w for s in states], ls, color=GREEN, lw=2.0, label=f"w, {name}"
    )
ax_a.set_xlim(0.0, 2.0)
ax_a.set_ylim(0.0, 27.0)
ax_a.legend(loc="upper left", fontsize=FS - 3, frameon=False)
style(ax_a, "t / T_rise", "u, w [mm]")

# (b) the soil force through the envelopes of the fast run
motion, states = runs[1.0]
geometries = motion.geometries
for target, color in zip((0.0, 8e-3, 16e-3, 25e-3), BLUES):
    t, geo, env = min(geometries, key=lambda g: abs(g[1].u - target))
    P = np.vstack((env.vertices, env.vertices[:1]))
    ax_b.plot(P[:, 1], P[:, 0], color=color, lw=1.8, label=f"u = {1e3 * geo.u:.0f} mm")
for T_rise, ls, name in ((50.0, "-", "slow"), (1.0, "--", "fast")):
    _, sts = runs[T_rise]
    ax_b.plot(
        [s.V for s in sts],
        [s.H for s in sts],
        ls,
        color=INK,
        lw=1.6,
        label=f"soil force, {name}",
    )
ax_b.axvline(W, color=GREEN, lw=1.2, ls=":", zorder=0)
ax_b.axhline(0.0, color=INK2, lw=0.7, zorder=0)
ax_b.text(W + 0.03, -0.55, "V = W", color=GREEN, fontsize=FS - 2)
ax_b.set_xlim(0.0, 4.2)
ax_b.set_ylim(-0.6, 0.95)
ax_b.legend(loc="upper left", fontsize=FS - 3, frameon=False, ncol=2)
zoom = ax_b.inset_axes([0.75, 0.43, 0.24, 0.55], xlim=(2.72, 2.92), ylim=(0.33, 0.63))
for target, color in zip((0.0, 8e-3, 16e-3, 25e-3), BLUES):
    t, geo, env = min(geometries, key=lambda g: abs(g[1].u - target))
    P = np.vstack((env.vertices, env.vertices[:1]))
    zoom.plot(P[:, 1], P[:, 0], color=color, lw=1.8)
for T_rise, ls in ((50.0, "-"), (1.0, "--")):
    _, sts = runs[T_rise]
    zoom.plot([s.V for s in sts], [s.H for s in sts], ls, color=INK, lw=1.6)
zoom.axvline(W, color=GREEN, lw=1.2, ls=":", zorder=0)
zoom.tick_params(colors=INK2, labelsize=FS - 4)
ax_b.indicate_inset_zoom(zoom, edgecolor=INK2)
style(ax_b, "V [kN/m], downward", "H [kN/m], toward +x", "(b) Through the envelope")

fig.subplots_adjust(left=0.06, right=0.98, wspace=0.18)
out = PLOTS / "dynamic_push.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print("wrote", out)
