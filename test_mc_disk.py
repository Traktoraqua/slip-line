"""Mohr-Coulomb disk foundation slip-line plots."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from slip_line.boundary import DiskFoundationBoundary
from slip_line.solver import MohrCoulombDiskPrandtlSolver
from slip_line.plot import plot_mc_disk_slip_lines

cases = [
    ('plot_mc_disk_phi20_p025_r0.png',  dict(phi_deg=20, p=0.25, r=0.0)),
    ('plot_mc_disk_phi30_p025_r0.png',  dict(phi_deg=30, p=0.25, r=0.0)),
    ('plot_mc_disk_phi30_p050_r0.png',  dict(phi_deg=30, p=0.50, r=0.0)),
    ('plot_mc_disk_phi30_p040_r05.png', dict(phi_deg=30, p=0.40, r=0.5)),
    ('plot_mc_disk_phi20_p010_r0.png',  dict(phi_deg=20, p=0.10, r=0.0)),
]

print(f"{'case':40s}  {'Q/(D*c)':>10s}  {'Nc_flat':>9s}  "
      f"{'sn_bot/c':>9s}  {'theta_max':>9s}  {'rows':>5s}")
for name, kw in cases:
    phi = np.radians(kw['phi_deg'])
    disk = DiskFoundationBoundary(cu=1.0, D=1.0, p=kw['p'], r=kw['r'])
    s = MohrCoulombDiskPrandtlSolver(disk=disk, c=1.0, phi=phi,
                                     n_fan=18, n_disk=30)
    res = s.solve()
    print(f"{name:40s}  {res['Q_over_Dc']:10.4f}  "
          f"{res['Nc_flat_smooth']:9.4f}  "
          f"{res['sigma_n_at_bottom']:9.4f}  "
          f"{res['theta_max_deg']:9.2f}  {len(s.region_I_mesh):5d}")
    fig = plot_mc_disk_slip_lines(s, show=False)
    fig.savefig(rf'c:\Appl\Slip_line\{name}', dpi=120, bbox_inches='tight')
    plt.close(fig)
