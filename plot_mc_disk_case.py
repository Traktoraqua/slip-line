import matplotlib
matplotlib.use("Agg")
import numpy as np
from slip_line.boundary import DiskFoundationBoundary
from slip_line.solver import MohrCoulombDiskPrandtlSolver
from slip_line.plot import plot_mc_disk_slip_lines

phi_deg = 30.0
p = 0.5
r = 0.5

disk = DiskFoundationBoundary(cu=1.0, D=1.0, p=p, r=r)
s = MohrCoulombDiskPrandtlSolver(disk=disk, c=1.0, phi=np.radians(phi_deg),
                                  n_fan=24, n_disk=40, n_radial=24)
out = s.solve()
print(out)

fig = plot_mc_disk_slip_lines(s, show=False)
fig.savefig(r'c:\Appl\Slip_line\mc_disk_phi30_p05_r05.png', dpi=130, bbox_inches='tight')
print("saved mc_disk_phi30_p05_r05.png")
