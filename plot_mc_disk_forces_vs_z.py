"""|F_v| and |F_h| on the right half-disk vs z/R — Mohr-Coulomb soil.

Both forces normalised by c·R vs z/R ∈ [0, 1] (p ∈ [0, 0.5]).
phi = 30°, four roughness values: r = 0, 0.5, 0.8, 1.0.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from slip_line.boundary import DiskFoundationBoundary
from slip_line.solver import MohrCoulombDiskPrandtlSolver

phi_deg = 30.0
phi     = np.radians(phi_deg)
c       = 1.0
D       = 1.0   # R = 0.5


def forces(p: float, r: float) -> tuple[float, float]:
    """(F_v, F_h) on right half-disk for given p and r."""
    if p < 1e-9:
        return 0.0, 0.0
    disk = DiskFoundationBoundary(cu=c, D=D, p=p, r=r)
    s = MohrCoulombDiskPrandtlSolver(disk=disk, c=c, phi=phi,
                                      n_fan=20, n_disk=40, n_radial=20)
    s.solve()
    F = s.forces_on_half_disk()
    return F['F_v'], F['F_h']


ps     = np.linspace(0.0, 0.5, 51)
zR     = 2.0 * ps          # z/R = p·D / R = 2p
norm   = 2.0               # F/(c·R) = F/(c·D/2) = 2·F/(c·D) = 2·F  (c=D=1)
rs     = [0.0, 0.5, 0.8, 1.0]
colors = {0.0: 'tab:blue', 0.5: 'tab:orange', 0.8: 'tab:green', 1.0: 'tab:red'}

fig, (axV, axH) = plt.subplots(1, 2, figsize=(13, 5.5))

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
    (axV, r'$|F_v| / (c \cdot R)$',
          f'Vertical resistance — Mohr-Coulomb  φ = {phi_deg:.0f}°'),
    (axH, r'$|F_h| / (c \cdot R)$',
          f'Horizontal resistance — Mohr-Coulomb  φ = {phi_deg:.0f}°'),
):
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r'penetration depth  $z / R$', fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.legend(loc='best', fontsize=10)

fig.suptitle(
    f'Mohr-Coulomb disk (φ = {phi_deg:.0f}°) — half-disk forces vs penetration depth (horizontal push)',
    fontsize=13,
)
fig.tight_layout()
outfile = rf'c:\Appl\Slip_line\plot_mc_disk_forces_vs_z_phi{int(phi_deg)}.png'
fig.savefig(outfile, dpi=130, bbox_inches='tight')
print(f"saved {outfile}")

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
