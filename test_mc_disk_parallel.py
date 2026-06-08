"""Quick sanity test for MC disk parallel-\u03b2 mechanism."""
import numpy as np
from slip_line.boundary import DiskFoundationBoundary
from slip_line.solver import MohrCoulombDiskPrandtlSolver

c = 1.0
print(f"{'phi':>5s} {'p':>5s} {'r':>4s}  {'parallel':>9s}  "
      f"{'th0_deg':>8s}  {'th_max_deg':>10s}  {'fan_deg':>8s}  "
      f"{'Q/(Dc)':>9s}  {'sn_bot':>8s}  {'sn_edge':>8s}")
cases = [
    (15.0, 0.10, 0.0),
    (15.0, 0.30, 0.5),
    (15.0, 0.40, 1.0),   # would previously give negative fan
    (15.0, 0.50, 1.0),
    (30.0, 0.20, 0.5),
    (30.0, 0.40, 1.0),   # would previously give negative fan
    (30.0, 0.50, 0.8),
    (30.0, 0.50, 1.0),
]
for phi_deg, p, r in cases:
    phi = np.radians(phi_deg)
    disk = DiskFoundationBoundary(cu=c, D=1.0, p=p, r=r)
    s = MohrCoulombDiskPrandtlSolver(disk=disk, c=c, phi=phi, n_fan=20, n_disk=40)
    out = s.solve()
    print(f"{phi_deg:5.1f} {p:5.2f} {r:4.1f}  "
          f"{str(out['parallel_beta_zone']):>9s}  "
          f"{out['theta_0_deg']:8.2f}  {out['theta_max_deg']:10.2f}  "
          f"{out['fan_angle_deg']:8.2f}  "
          f"{out['Q_over_Dc']:9.4f}  "
          f"{out['sigma_n_at_bottom']:8.3f}  {out['sigma_n_at_edge']:8.3f}")
