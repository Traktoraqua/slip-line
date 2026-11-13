"""Pipe (D = 1 m, roughness r = 0.8, weight W = 2.8 kN/m, mass from W with
C_a = 1) in weightless Tresca soil, pushed sideways by a load given in time,
with the geometry following the motion (PipeMotion, step 3 of
PLAN_DYNAMIKK.md).  It first sinks under its own weight
(PipePenetrationClosedForm).  Then Fx = 0.3·(1 − cos(πt/T_rise)) rises to
0.6 kN/m at T_rise and falls again:

  1. cu = 1 kPa, T_rise = 50 s (slow, T_n ≈ 2–3 s): the quasi-static push,
     force controlled.
  2. cu = 1 kPa, T_rise = 1 s: the same load, fast.
  3. cu = 5 kPa, T_rise = 50 s: the light pipe rises and leaves the trench
     floor behind it.

Each run stops when the motion would need a mechanism of the mirror image
(moving back, or sinking back into the gap behind), which comes in step 4.
Prints the steps, the energy balance and the time taken."""

import time

import numpy as np

from slip_line.dynamics import PipeMotion, pipe_mass
from slip_line.lateral import PipeLateralClosedForm
from slip_line.penetration import PipePenetrationClosedForm

W, D, F_MAX = 2.8, 1.0, 0.6
M = pipe_mass(W, D)  # C_a = 1


def rising(T_rise):
    return lambda t: 0.5 * F_MAX * (1.0 - np.cos(np.pi * t / T_rise))


def report(title, cu, T_rise, dt):
    sim = PipePenetrationClosedForm(D=D, cu=cu, r=0.8, W=W)
    lat, start = PipeLateralClosedForm.after_penetration(sim, sim.run()[-1])
    motion = PipeMotion(lat, start, M, Fx=rising(T_rise))
    states = motion.run(2.0 * T_rise, dt)
    Fx = motion.dyn.Fx
    print(f"\n{title}")
    print(
        "   t [s]  Fx [kN/m]  H [kN/m]  V [kN/m]  u [mm]  w [mm]  du/dt [mm/s]  delta"
    )
    moving = [s for s in states if s.u_dot > 0.0]
    shown = [states[0], *moving[:: max(1, len(moving) // 12)], states[-1]]
    for s in shown:
        print(
            f"{s.t:8.3f} {Fx(s.t):10.3f} {s.H:9.3f} {s.V:9.3f}"
            f" {1e3 * s.u:7.2f} {1e3 * s.w:7.3f} {1e3 * s.u_dot:13.2f}"
            f" {np.degrees(s.delta):6.1f}"
        )
    e = motion.dyn.energy(states)
    print(
        f"  stop: {motion.stop_reason} at t = {states[-1].t:.3f} s,"
        f" {len(motion.geometries)} envelopes"
    )
    print(
        f"  energy [J/m]: work {1e3 * e['work']:.2f}, dissipation"
        f" {1e3 * e['dissipation']:.2f}, kinetic {1e3 * e['kinetic']:.3f},"
        f" numerical {1e3 * e['numerical']:.3f}"
    )


t_start = time.perf_counter()
print(f"M = {M[0]:.0f} kg/m, load Fx = {F_MAX / 2:g}*(1 - cos(pi t / T_rise))")
report("1. cu = 1 kPa, T_rise = 50 s (slow):", 1.0, 50.0, 0.02)
report("2. cu = 1 kPa, T_rise = 1 s (fast):", 1.0, 1.0, 2e-3)
report("3. cu = 5 kPa, T_rise = 50 s (light pipe):", 5.0, 50.0, 0.02)
print(f"\n{time.perf_counter() - t_start:.1f} s in all")
