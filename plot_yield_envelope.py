"""Yield envelopes of the soil around the pipe (plots/yield_envelope.png).

YieldEnvelope.from_state on states of the quasi-static lateral push with
uniform friction (PipeLateralClosedForm): D = 1 m, r = 0.8, W = 2.8 kN/m.
(a) cu = 5 kPa, the pipe of PLAN_DYNAMIKK.md, at the start and after 4 and
10 mm.  (b) cu = 1 kPa, run_lateral_closed_form.py, at the start and after
50 and 200 mm.  H is the force from the pipe on the soil toward +x, V
downward.  Dots: the push force at V = W (push_force, the H of a fresh
mechanism search), filled where it lies on the envelope.  Crosses: the H of
the quasi-static run where its local search has not found the lowest
minimum.  Prints the numbers used in manual.tex.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from slip_line.envelope import YieldEnvelope
from slip_line.lateral import PipeLateralClosedForm
from slip_line.penetration import PipePenetrationClosedForm

PLOTS = Path(__file__).parent / "plots"
PLOTS.mkdir(exist_ok=True)

W = 2.8
BLUES = ("#7fb0e8", "#2a6fc9", "#0e3470")  # light to dark: later states
INK, INK2, GRID, GREEN = "#0b0b0b", "#52514e", "#e4e3df", "#1a8a5c"
FS = 14


def push(cu: float, marks: tuple[float, ...]):
    """Quasi-static push (0.25 mm steps for the first 50 mm, then 1 mm) and
    the states nearest to the displacements in marks (m)."""
    sim = PipePenetrationClosedForm(D=1.0, cu=cu, r=0.8, W=W)
    lat, start = PipeLateralClosedForm.after_penetration(sim, sim.run()[-1])
    states = lat.run(start, min(0.05, max(marks)), du=0.25e-3)
    if max(marks) > 0.05:
        states += lat.run(states[-1], max(marks) - 0.05)[1:]
    return lat, [min(states, key=lambda s: abs(s.u - u)) for u in marks]


def panel(ax, cu: float, marks: tuple[float, ...], title: str):
    lat, states = push(cu, marks)
    print(f"cu = {cu:g} kPa")
    print(
        "  u [mm]  w [mm]  d [mm]  run H  beta   push H  delta  H at V = W (back, fwd)"
        "   V at H = 0  V max   H max at V"
    )
    for s, color in zip(states, BLUES):
        env = YieldEnvelope.from_state(lat, s)
        P = np.vstack((env.vertices, env.vertices[:1]))
        label = f"u = {1e3 * s.u:.0f} mm"
        ax.fill(P[:, 1], P[:, 0], color=color, alpha=0.07, lw=0, zorder=1)
        ax.plot(P[:, 1], P[:, 0], color=color, lw=2.0, zorder=2, label=label)
        H, delta = env.push_force(W)
        on = env.contains([H, W], tol=1e-6)
        ax.plot(
            W,
            H,
            "o",
            ms=7,
            color=color,
            mfc=color if on else "white",
            mew=1.6,
            zorder=5,
        )
        if s.mech.H > H + 1e-3:
            ax.plot(W, s.mech.H, "x", ms=8, mew=2.0, color=color, zorder=5)
        hr = env.H_range(W)
        at_w = "outside" if hr is None else f"{hr[0]:7.3f} {hr[1]:6.3f}"
        j, k = np.argmax(env.vertices[:, 1]), np.argmax(env.vertices[:, 0])
        print(
            f"  {1e3 * s.u:6.1f} {1e3 * s.w:7.2f} {1e3 * s.d(lat.R):7.2f}"
            f" {s.mech.H:6.3f} {np.degrees(s.mech.beta):5.1f}"
            f" {H:8.3f} {np.degrees(delta):6.1f}  {at_w:>22}"
            f" {env.V_range(0.0)[1]:12.3f} {env.vertices[j, 1]:6.3f}"
            f" {env.vertices[k, 0]:7.3f} {env.vertices[k, 1]:5.3f}"
        )
    ax.axvline(W, color=GREEN, lw=1.2, ls=":", zorder=0)
    ax.axhline(0.0, color=INK2, lw=0.7, zorder=0)
    ax.set_xlabel("V [kN/m], downward", color=INK2, fontsize=FS)
    ax.set_ylabel("H [kN/m], toward +x", color=INK2, fontsize=FS)
    ax.set_title(title, fontsize=FS, loc="left", color=INK)
    ax.grid(color=GRID, lw=0.6)
    ax.tick_params(colors=INK2, labelsize=FS - 2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_xlim(left=0.0)


fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 5.8))
panel(ax_a, 5.0, (0.0, 0.004, 0.010), r"(a) $c_u = 5$ kPa")
panel(ax_b, 1.0, (0.0, 0.05, 0.20), r"(b) $c_u = 1$ kPa")
ax_a.set_ylim(-0.35, 0.78)
ax_b.set_ylim(-0.6, 1.2)
for ax in (ax_a, ax_b):
    lo, hi = ax.get_ylim()
    ax.text(
        W + 0.04, hi - 0.02 * (hi - lo), "V = W", color=GREEN, fontsize=FS - 2, va="top"
    )
ax_a.plot([], [], "o", color=INK2, ms=7, label="push force at V = W")
ax_a.plot(
    [],
    [],
    "o",
    color=INK2,
    mfc="white",
    ms=7,
    mew=1.6,
    label="  ... outside the envelope",
)
ax_a.plot(
    [], [], "x", color=INK2, ms=8, mew=2.0, label="quasi-static run (local search)"
)
ax_a.legend(loc="upper left", fontsize=FS - 3, frameon=False)
ax_b.legend(loc="upper left", fontsize=FS - 3, frameon=False)
fig.tight_layout()
out = PLOTS / "yield_envelope.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print("wrote", out)
