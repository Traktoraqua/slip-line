"""Velocity field (Geiringer) for a smooth disk pushed vertically into a
weightless Tresca soil; saves plots/velocity_disk_p025.png and prints the
work balance for a few penetrations."""

import matplotlib

matplotlib.use("Agg")
from pathlib import Path

from slip_line.boundary import DiskFoundationBoundary
from slip_line.plot import plot_disk_velocity
from slip_line.solver import DiskPrandtlSolver
from slip_line.velocity import DiskVelocityField

PLOTS = Path(__file__).parent / "plots"
PLOTS.mkdir(exist_ok=True)

cu, D, v0 = 1.0, 1.0, 1.0


def velocity_field(p: float) -> DiskVelocityField:
    disk = DiskFoundationBoundary(cu=cu, D=D, p=p, r=0.0)
    solver = DiskPrandtlSolver(disk=disk, n_fan=18, n_disk=30)
    solver.solve()
    vf = DiskVelocityField(solver, v0=v0)
    vf.solve()
    return vf


print("   p   Q/(D·cu): stress  virtual work  upper bound   negative-work area")
for p in (0.1, 0.25, 0.4, 0.5):
    e = velocity_field(p).energy_balance()
    print(
        f"{p:5.2f}            {e['F_stress'] / (D * cu):6.3f}  "
        f"{e['F_virtual_work'] / (D * cu):12.3f}  {e['F_upper_bound'] / (D * cu):11.3f}"
        f"   {100 * e['negative_area_fraction']:5.0f} %"
    )

fig = plot_disk_velocity(velocity_field(0.25), show=False)
fig.savefig(PLOTS / "velocity_disk_p025.png", dpi=130, bbox_inches="tight")
print("Saved plots/velocity_disk_p025.png")
