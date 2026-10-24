"""Pipe (D = 1 m, roughness r = 0.8, weight W = 2.8 kN/m) in weightless Tresca
soil (cu = 1 kPa).  Step 1: it sinks under its own weight from the surface
(PipePenetration).  Step 2: it is pushed 0.2 m sideways in 1 mm steps while it
is free to move vertically (PipeLateral, V = W).  Prints both."""

import numpy as np

from slip_line.lateral import PipeLateral
from slip_line.penetration import PipePenetration

sim = PipePenetration(D=1.0, cu=1.0, r=0.8, W=2.8)
s1 = sim.run()[-1]
print("Step 1, sinking under the pipe's weight:")
print(f"  settlement w       = {1e3 * s1.w:.2f} mm")
print(f"  heave at the pipe  = {1e3 * (s1.d - s1.w):.2f} mm")
print(f"  effective depth d  = {1e3 * s1.d:.2f} mm")

lat, start = PipeLateral.after_penetration(sim, s1, du=1e-3)
states = lat.run(start, 0.2)
R = lat.R
print("\nStep 2, lateral push with V = W:")
print("  u [mm]   w [mm]   d [mm]  berm [mm]  beta [deg]  alpha [deg]  H [kN/m]")
for k, s in enumerate(states):
    if k % 10 == 0:
        m = s.mech
        print(
            f"{1e3 * s.u:8.1f} {1e3 * s.w:8.2f} {1e3 * s.d(R):8.2f}"
            f" {1e3 * (s.d(R) - s.w):10.2f} {np.degrees(m.beta):11.2f}"
            f" {np.degrees(s.alpha):12.2f} {m.H:9.3f}"
        )
