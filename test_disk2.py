"""Re-render disk plots with curved Region I characteristics."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from slip_line.boundary import DiskFoundationBoundary
from slip_line.solver import DiskPrandtlSolver
from slip_line.plot import plot_disk_slip_lines

# Debug: dump Region I mesh for p=0.25, n_disk=8
d = DiskFoundationBoundary(cu=1, D=1, p=0.10, r=0.8)
s_dbg = DiskPrandtlSolver(disk=d, n_fan=18, n_disk=8)
s_dbg.solve()
with open('dbg.out', 'w') as f:
    f.write(f'B_s={s_dbg.B_s:.4f}  y_c={s_dbg.y_c:.4f}  R={s_dbg.R:.4f}\n')
    f.write(f'Rows: {len(s_dbg.region_I_mesh)}\n')
    for k, row in enumerate(s_dbg.region_I_mesh):
        f.write(f'--- Row {k} ({len(row)} nodes) ---\n')
        for n in row:
            f.write(f'  ({n.x:+.4f}, {n.y:+.4f})  psi={n.psi*180/3.14159:6.2f} deg\n')

cases = [
    ('plot_disk_p025.png',    dict(p=0.25, r=0.0)),
    ('plot_disk_p050.png',    dict(p=0.50, r=0.0)),
    ('plot_disk_p010.png',    dict(p=0.10, r=0.0)),
    ('plot_disk_p04_r05.png', dict(p=0.40, r=0.5)),
]

for name, kw in cases:
    disk = DiskFoundationBoundary(cu=1.0, D=1.0, **kw)
    s = DiskPrandtlSolver(disk=disk, n_fan=18, n_disk=30)
    res = s.solve()
    print(f"{name}: p={kw['p']:.2f}, r={kw['r']:.2f}  "
          f"Q/(D*cu)={res['Q_over_D_over_cu']:.4f}  "
          f"theta_max={res['theta_max_deg']:.1f} deg  "
          f"mesh_rows={len(s.region_I_mesh)}")
    fig = plot_disk_slip_lines(s, show=False)
    fig.savefig(rf'c:\Appl\Slip_line\{name}', dpi=120)
    plt.close(fig)
