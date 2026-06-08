"""F_v/(D*cu) vs penetration ratio p for a Tresca disk foundation.

Left panel  : penetration phase, p ∈ [0, 0.5] (real solver).
Right panel : extraction phase,  p ∈ [0.5, 1], by symmetry F_v(p) = F_v(1−p).

Three curves per panel: r = 0, 0.5, 1.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from slip_line.boundary import DiskFoundationBoundary
from slip_line.solver import DiskPrandtlSolver


def Fv_over_Dcu(p: float, r: float) -> float:
    """Half-disk (right side) vertical force, normalised by D·cu."""
    if p < 1e-9:
        return 0.0
    disk = DiskFoundationBoundary(cu=1.0, D=1.0, p=p, r=r)
    s = DiskPrandtlSolver(disk=disk, n_fan=18, n_disk=30)
    s.solve()
    return s.forces_on_half_disk()['F_v']      # cu=D=1 so already normalised


ps_pen = np.linspace(0.0, 0.5, 41)
rs     = [0.0, 0.5, 1.0]
colors = {0.0: 'tab:blue', 0.5: 'tab:orange', 1.0: 'tab:red'}

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)

for r in rs:
    Fv_pen = np.array([Fv_over_Dcu(p, r) for p in ps_pen])
    # Extraction: F_v(p) = F_v(1-p), 0.5 ≤ p ≤ 1
    ps_ext = 1.0 - ps_pen[::-1]            # [0.5 … 1.0]
    Fv_ext = Fv_pen[::-1]

    axL.plot(ps_pen, Fv_pen, color=colors[r], lw=2.0,
             marker='o', ms=4, label=f'r = {r:.1f}')
    axR.plot(ps_ext, Fv_ext, color=colors[r], lw=2.0,
             marker='o', ms=4, label=f'r = {r:.1f}')

for ax, title, xlim in [
    (axL, 'Penetration phase', (0.0, 0.5)),
    (axR, 'Extraction phase  (F_v(p) = F_v(1−p))', (0.5, 1.0)),
]:
    ax.set_xlim(xlim)
    ax.set_xlabel('penetration ratio  p = z / D', fontsize=11)
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.set_title(title, fontsize=12)
    ax.legend(loc='best', fontsize=10)
    ax.axvline(0.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)

axL.set_ylabel(r'$F_v / (D \cdot c_u)$  (half-disk, per unit length)', fontsize=11)
axL.set_ylim(bottom=0)
fig.suptitle('Tresca disk foundation — vertical force on right half-disk', fontsize=13)
fig.tight_layout()
fig.savefig(r'c:\Appl\Slip_line\plot_disk_Fv_vs_p.png',
            dpi=130, bbox_inches='tight')
print("saved plot_disk_Fv_vs_p.png")

# Print summary table
print(f"\n{'p':>5s}  " + "  ".join(f"r={r:.1f}".rjust(10) for r in rs))
for i, p in enumerate(ps_pen):
    row = [f"{Fv_over_Dcu(p, r):10.4f}" for r in rs]
    print(f"{p:5.3f}  " + "  ".join(row))
