"""F_v and F_h on the right half-disk vs penetration ratio p (Tresca disk).

Right panel : F_v/(D·cu) for p ∈ [0, 1].  Penetration 0..0.5 from solver,
              extraction 0.5..1 by symmetry F_v(p) = F_v(1−p).
Left  panel : F_h/(D·cu) for p ∈ [0, 0.5].

Three curves per panel: r = 0, 0.5, 1.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from slip_line.boundary import DiskFoundationBoundary
from slip_line.solver import DiskPrandtlSolver


def forces(p: float, r: float) -> tuple[float, float]:
    """(F_v, F_h) on the right half-disk, normalised by D·cu (cu=D=1)."""
    if p < 1e-9:
        return 0.0, 0.0
    disk = DiskFoundationBoundary(cu=1.0, D=1.0, p=p, r=r)
    s = DiskPrandtlSolver(disk=disk, n_fan=18, n_disk=30)
    s.solve()
    F = s.forces_on_half_disk()
    return F['F_v'], F['F_h']


ps_pen = np.linspace(0.0, 0.5, 41)
rs     = [0.0, 0.5, 1.0]
colors = {0.0: 'tab:blue', 0.5: 'tab:orange', 1.0: 'tab:red'}

fig, (axH, axV) = plt.subplots(1, 2, figsize=(13, 5.5))

for r in rs:
    Fv_pen = np.empty_like(ps_pen)
    Fh_pen = np.empty_like(ps_pen)
    for i, p in enumerate(ps_pen):
        Fv_pen[i], Fh_pen[i] = forces(p, r)

    # ---- Left panel: F_h for p ∈ [0, 0.5] ----
    axH.plot(ps_pen, Fh_pen, color=colors[r], lw=2.0,
             marker='o', ms=4, label=f'r = {r:.1f}')

    # ---- Right panel: F_v for p ∈ [0, 1] (mirror for extraction) ----
    ps_full = np.concatenate([ps_pen, 1.0 - ps_pen[::-1][1:]])
    Fv_full = np.concatenate([Fv_pen, Fv_pen[::-1][1:]])
    axV.plot(ps_full, Fv_full, color=colors[r], lw=2.0,
             marker='o', ms=3, label=f'r = {r:.1f}')

# Left panel formatting
axH.set_xlim(0.0, 0.5)
axH.set_xlabel('penetration ratio  p = z / D', fontsize=11)
axH.set_ylabel(r'$F_h / (D \cdot c_u)$  (half-disk, outward = +)', fontsize=11)
axH.set_title('Horizontal force on right half-disk', fontsize=12)
axH.grid(True, linestyle=':', alpha=0.5)
axH.legend(loc='best', fontsize=10)
axH.axhline(0.0, color='gray', linewidth=0.8, alpha=0.6)

# Right panel formatting
axV.set_xlim(0.0, 1.0)
axV.set_ylim(bottom=0)
axV.set_xlabel('penetration ratio  p = z / D', fontsize=11)
axV.set_ylabel(r'$F_v / (D \cdot c_u)$  (half-disk, upward)', fontsize=11)
axV.set_title('Vertical force — penetration (p<0.5) + extraction (p>0.5)', fontsize=12)
axV.grid(True, linestyle=':', alpha=0.5)
axV.legend(loc='best', fontsize=10)
axV.axvline(0.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)

fig.suptitle('Tresca disk foundation — half-disk forces vs penetration',
             fontsize=13)
fig.tight_layout()
fig.savefig(r'c:\Appl\Slip_line\plot_disk_Fv_Fh_vs_p.png',
            dpi=130, bbox_inches='tight')
print("saved plot_disk_Fv_Fh_vs_p.png")

# Summary table
print(f"\n{'p':>5s}  " +
      "  ".join(f"Fv(r={r:.1f})".rjust(12) for r in rs) + "  " +
      "  ".join(f"Fh(r={r:.1f})".rjust(12) for r in rs))
for i, p in enumerate(ps_pen[::4]):
    idx = i * 4
    Fv = [forces(p, r)[0] for r in rs]
    Fh = [forces(p, r)[1] for r in rs]
    print(f"{p:5.3f}  " +
          "  ".join(f"{v:12.4f}" for v in Fv) + "  " +
          "  ".join(f"{h:12.4f}" for h in Fh))
