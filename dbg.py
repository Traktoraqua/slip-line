import sys
from slip_line.boundary import DiskFoundationBoundary
from slip_line.solver import DiskPrandtlSolver

with open('dbg.out', 'w') as f:
    d = DiskFoundationBoundary(cu=1, D=1, p=0.25, r=0)
    s = DiskPrandtlSolver(disk=d, n_fan=18, n_disk=8)
    s.solve()
    f.write(f'B_s={s.B_s:.4f}  y_c={s.y_c:.4f}  R={s.R:.4f}\n')
    f.write(f'Rows: {len(s.region_I_mesh)}\n')
    for k, row in enumerate(s.region_I_mesh):
        f.write(f'--- Row {k} ({len(row)} nodes) ---\n')
        for n in row:
            f.write(f'  ({n.x:+.4f}, {n.y:+.4f})  psi={n.psi*180/3.14159:6.2f} deg  sig={n.sigma:6.3f}\n')
print("done")
