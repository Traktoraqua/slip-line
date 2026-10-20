"""Velocity field (Geiringer) for a smooth strip footing pushed vertically
into a weightless Tresca soil; saves plots/velocity_strip_smooth.png."""

import itertools

import matplotlib

matplotlib.use("Agg")
from pathlib import Path

import numpy as np

from slip_line.boundary import FoundationBoundary
from slip_line.plot import plot_prandtl_velocity
from slip_line.solver import PrandtlSolver
from slip_line.velocity import PrandtlVelocityField

PLOTS = Path(__file__).parent / "plots"
PLOTS.mkdir(exist_ok=True)

cu, B, v0 = 1.0, 2.0, 1.0
solver = PrandtlSolver(FoundationBoundary(cu=cu, r=0.0, B=B), n_fan=12)
q = solver.solve()

vf = PrandtlVelocityField(solver, v0=v0)
U = vf.solve()

# Energy check (per half footing): q·B·V0 = dissipation.
#   jump along outer α-line:  cu·U·L_outer
#   shear in the fan (u = U, v = 0):  cu·U·R·Δψ   (R = fan radius)
outer = vf.outer_alpha_line()
L_outer = sum(np.hypot(b.x - a.x, b.y - a.y) for a, b in itertools.pairwise(outer))
R_fan = solver.n_radial * solver.dr
D = cu * U * (L_outer + R_fan * solver.delta_psi)

print(f"Velocity jump U / V0     = {U / v0:.6f}  (sqrt 2 = {np.sqrt(2):.6f})")
print(f"q / cu (stress field)    = {q / cu:.6f}")
print(f"q / cu (work = dissip.)  = {D / (B * v0) / cu:.6f}")

fig = plot_prandtl_velocity(vf, show=False)
fig.savefig(PLOTS / "velocity_strip_smooth.png", dpi=130, bbox_inches="tight")
print("Saved plots/velocity_strip_smooth.png")
