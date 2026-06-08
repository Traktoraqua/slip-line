"""|F_v| and |F_h| on the right half-disk vs penetration depth z/R.

Both forces are normalised by cu·R where R = D/2 is the disk radius and
plotted as absolute values vs z/R = 2p, for z/R ∈ [0, 1] (i.e. p ∈ [0, 0.5]).

Three curves per panel: r = 0, 0.5, 1.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from slip_line.boundary import DiskFoundationBoundary
from slip_line.solver import DiskPrandtlSolver


def forces(p: float, r: float) -> tuple[float, float]:
    """(F_v, F_h) on the right half-disk, with cu = 1, D = 1 (R = 0.5)."""
    if p < 1e-9:
        return 0.0, 0.0
    disk = DiskFoundationBoundary(cu=1.0, D=1.0, p=p, r=r)
    s = DiskPrandtlSolver(disk=disk, n_fan=18, n_disk=30)
    s.solve()
    F = s.forces_on_half_disk()
    return F['F_v'], F['F_h']


ps   = np.linspace(0.0, 0.5, 51)
zR   = 2.0 * ps   # z/R = 2 p
rs     = [0.0, 0.5, 0.8, 1.0]
colors = {0.0: 'tab:blue', 0.5: 'tab:orange', 0.8: 'tab:green', 1.0: 'tab:red'}

fig, (axV, axH) = plt.subplots(1, 2, figsize=(13, 5.5))

# F/(cu·R) = F/(cu·D/2) = 2 · F/(cu·D)
norm = 2.0

for r in rs:
    Fv = np.empty_like(ps)
    Fh = np.empty_like(ps)
    for i, p in enumerate(ps):
        Fv[i], Fh[i] = forces(p, r)

    axV.plot(zR, np.abs(Fv) * norm, color=colors[r], lw=2.0,
             marker='o', ms=4, label=f'r = {r:.1f}')
    axH.plot(zR, np.abs(Fh) * norm, color=colors[r], lw=2.0,
             marker='o', ms=4, label=f'r = {r:.1f}')

for ax, ylabel, title in (
    (axV, r'$|F_v| / (c_u \cdot R)$', 'Vertical resistance on right half-disk'),
    (axH, r'$|F_h| / (c_u \cdot R)$', 'Horizontal resistance on right half-disk'),
):
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r'penetration depth  $z / R$', fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.legend(loc='best', fontsize=10)

fig.suptitle('Tresca disk — half-disk forces vs penetration depth (horizontal push)',
             fontsize=13)
fig.tight_layout()
fig.savefig(r'c:\Appl\Slip_line\plot_disk_forces_vs_z.png',
            dpi=130, bbox_inches='tight')
print("saved plot_disk_forces_vs_z.png")

# Summary table
print(f"\n{'z/R':>5s}  " +
      "  ".join(f"|Fv|(r={r:.1f})".rjust(14) for r in rs) + "  " +
      "  ".join(f"|Fh|(r={r:.1f})".rjust(14) for r in rs))
for i in range(0, len(ps), 5):
    p = ps[i]
    Fv = [abs(forces(p, r)[0]) * norm for r in rs]
    Fh = [abs(forces(p, r)[1]) * norm for r in rs]
    print(f"{zR[i]:5.3f}  " +
          "  ".join(f"{v:14.4f}" for v in Fv) + "  " +
          "  ".join(f"{h:14.4f}" for h in Fh))
