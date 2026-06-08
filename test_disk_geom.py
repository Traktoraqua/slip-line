"""Inspect disk Region I corner β-line endpoint vs fan_I expected angle."""
import matplotlib
matplotlib.use("Agg")
import numpy as np
from slip_line.boundary import DiskFoundationBoundary
from slip_line.solver import DiskPrandtlSolver

for p, r in [(0.1, 0.0), (0.25, 0.0), (0.5, 0.0), (0.4, 0.5)]:
    disk = DiskFoundationBoundary(cu=1.0, D=1.0, p=p, r=r)
    sv = DiskPrandtlSolver(disk, n_fan=20, n_disk=20)
    sv.solve()
    mesh = sv.region_I_mesh
    corner = mesh[0][-1]              # mesh[0][n_disk] = corner
    deep_end = mesh[-1][0]            # last row, first column
    dx = deep_end.x - corner.x
    dy = deep_end.y - corner.y
    dist = (dx*dx + dy*dy)**0.5
    angle_actual = np.degrees(np.arctan2(dy, dx))
    angle_expected = np.degrees(sv.psi_edge + np.pi/4)  # theta_beta
    print(f"p={p:.2f} r={r:.2f}: psi_edge={np.degrees(sv.psi_edge):.2f}deg, "
          f"theta_beta_exp={angle_expected:.2f}deg, actual={angle_actual:.2f}deg, "
          f"dist={dist:.4f}, n_radial·dr={sv.n_radial*sv.dr:.4f}, "
          f"B_s/|cos(theta_b)|={sv.B_s/abs(np.cos(sv.psi_edge+np.pi/4)):.4f}")
