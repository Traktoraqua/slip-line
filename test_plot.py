"""Regenerate flat foundation plots after Region I fix."""
import matplotlib
matplotlib.use("Agg")

from slip_line.boundary import FoundationBoundary
from slip_line.solver import PrandtlSolver
from slip_line.plot import plot_slip_lines

# Smooth (r=0)
fdn = FoundationBoundary(cu=1.0, r=0.0, B=2.0)
solver = PrandtlSolver(fdn, n_fan=18)
q = solver.solve()
print(f"Smooth: q/cu = {q:.6f}  (pi+2 = {3.141592653589793+2:.6f})")
fig = plot_slip_lines(solver, q, show=False)
fig.savefig("plot_smooth.png", dpi=120, bbox_inches="tight")

# Rough (r=0.5)
fdn = FoundationBoundary(cu=1.0, r=0.5, B=2.0)
solver = PrandtlSolver(fdn, n_fan=18)
q = solver.solve()
print(f"r=0.5: q/cu = {q:.6f}")
fig = plot_slip_lines(solver, q, show=False)
fig.savefig("plot_r05.png", dpi=120, bbox_inches="tight")

print("Done.")
