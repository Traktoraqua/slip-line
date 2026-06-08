"""Test the piecewise-roughness (parallel-beta zone) regime."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from slip_line.boundary import DiskFoundationBoundary
from slip_line.solver import DiskPrandtlSolver
from slip_line.plot import plot_disk_slip_lines

cases = [
    ('plot_disk_p050_r050.png',  dict(p=0.50, r=0.50)),   # user's case
    ('plot_disk_p050_r020.png',  dict(p=0.50, r=0.20)),
    ('plot_disk_p045_r050.png',  dict(p=0.45, r=0.50)),
    ('plot_disk_p040_r050.png',  dict(p=0.40, r=0.50)),   # right at the threshold
    ('plot_disk_p035_r050.png',  dict(p=0.35, r=0.50)),   # below threshold (no parallel zone)
    ('plot_disk_p050_r100.png',  dict(p=0.50, r=1.00)),
]

print(f"{'case':35s} {'p':>5s} {'r':>5s} {'theta_max':>9s} {'theta_0':>8s} "
      f"{'parallel?':>9s} {'fan_deg':>8s} {'Q/(D*cu)':>10s}")
for name, kw in cases:
    d = DiskFoundationBoundary(cu=1.0, D=1.0, **kw)
    s = DiskPrandtlSolver(disk=d, n_fan=18, n_disk=30)
    res = s.solve()
    print(f"{name:35s} {kw['p']:5.2f} {kw['r']:5.2f} "
          f"{res['theta_max_deg']:9.2f} {res['theta_0_deg']:8.2f} "
          f"{str(res['parallel_beta_zone']):>9s} "
          f"{res['fan_angle_deg']:8.2f} {res['Q_over_D_over_cu']:10.4f}")
    fig = plot_disk_slip_lines(s, show=False)
    fig.savefig(rf'c:\Appl\Slip_line\{name}', dpi=120, bbox_inches='tight')
    plt.close(fig)
