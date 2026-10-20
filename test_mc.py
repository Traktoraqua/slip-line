"""Generate Mohr-Coulomb Prandtl slip-line field plots."""
import matplotlib
matplotlib.use("Agg")
from pathlib import Path
PLOTS = Path(__file__).parent / 'plots'
PLOTS.mkdir(exist_ok=True)
import numpy as np

from slip_line.solver import MohrCoulombPrandtlSolver
from slip_line.plot import plot_mc_slip_lines

cases = [
    dict(phi_deg=20, label="phi20"),
    dict(phi_deg=30, label="phi30"),
    dict(phi_deg=40, label="phi40"),
]

for case in cases:
    phi = np.radians(case["phi_deg"])
    solver = MohrCoulombPrandtlSolver(c=1.0, phi=phi, B=2.0, n_fan=18)
    q = solver.solve()
    Nc = solver.analytical_Nc()
    print(
        f"phi={case['phi_deg']:2d}°  q/c = {q:.4f}  Nc_analyt = {Nc:.4f}"
        f"  diff = {abs(q/1.0 - Nc):.2e}"
    )
    fig = plot_mc_slip_lines(solver, show=False)
    fname = f"plot_mc_{case['label']}.png"
    fig.savefig(PLOTS / fname, dpi=120, bbox_inches="tight")
    print(f"  saved {fname}")
