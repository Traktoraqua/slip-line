"""Pipe (D = 1 m, roughness r = 0.8, weight W = 2.8 kN/m) in weightless Tresca
soil (cu = 1 kPa), with the closed forms of pipe_formulas.tex instead of a
slip-line net.  Step 1: it sinks under its own weight from the surface
(PipePenetrationClosedForm).  Step 2: it is pushed 0.2 m sideways in small
steps while it is free to move vertically (PipeLateralClosedForm, V = W, with
the interface friction +r*cu toward the front edge on the whole contact).
Step 3: it is pushed 40 mm back toward -x, still with V = W.  Prints all
three and the time taken; H is positive toward +x.  The same model on the net
is PipeLateral(uniform=True)."""

import time

import numpy as np

from slip_line.lateral import PipeLateralClosedForm
from slip_line.penetration import PipePenetrationClosedForm

t_start = time.perf_counter()
sim = PipePenetrationClosedForm(D=1.0, cu=1.0, r=0.8, W=2.8)
s1 = sim.run()[-1]
print("Step 1, sinking under the pipe's weight:")
print(f"  settlement w       = {1e3 * s1.w:.2f} mm")
print(f"  heave at the pipe  = {1e3 * (s1.d - s1.w):.2f} mm")
print(f"  effective depth d  = {1e3 * s1.d:.2f} mm")


def on_grid(u, step_mm):
    """True when the displacement u (m) is a multiple of step_mm."""
    n = 1e3 * u / step_mm
    return abs(n - round(n)) < 1e-6


# Steps of 0.25 mm at the start of each push, where the direction changes fast.
lat, start = PipeLateralClosedForm.after_penetration(sim, s1, du=1e-3)
states = lat.run(start, 0.05, du=0.25e-3)
states += lat.run(states[-1], 0.15)[1:]
R = lat.R
print("\nStep 2, lateral push with V = W (uniform friction):")
print("  u [mm]   w [mm]   d [mm]  berm [mm]  beta [deg]  alpha [deg]  H [kN/m]")
for s in states:
    if on_grid(s.u, 10):
        m = s.mech
        print(
            f"{1e3 * s.u:8.1f} {1e3 * s.w:8.2f} {1e3 * s.d(R):8.2f}"
            f" {1e3 * (s.d(R) - s.w):10.2f} {np.degrees(m.beta):11.2f}"
            f" {np.degrees(s.alpha):12.2f} {m.H:9.3f}"
        )

rev, back = lat.reversed(states[-1])
states3 = rev.run(back, 0.01, du=0.25e-3)
states3 += rev.run(states3[-1], 0.03)[1:]
print("\nStep 3, pushed back toward -x with V = W (front = soil level at the")
print("contact on the -x side, above the original surface):")
print("  u [mm]   w [mm]   d [mm] front [mm]  beta [deg]  alpha [deg]  H [kN/m]")
for s in states3:
    if on_grid(s.u, 5):
        m = s.mech
        print(
            f"{-1e3 * s.u:8.1f} {1e3 * s.w:8.2f} {1e3 * s.d(R):8.2f}"
            f" {1e3 * (s.d(R) - s.w):10.2f} {np.degrees(m.beta):11.2f}"
            f" {np.degrees(s.alpha):12.2f} {-m.H:9.3f}"
        )
print(f"\n{time.perf_counter() - t_start:.1f} s in all")
