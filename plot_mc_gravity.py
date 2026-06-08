"""
Mohr-Coulomb strip footing with submerged soil weight (gravity).

Plots the slip-line field (weightless topology with gravity-corrected
stresses) and the contact-pressure distribution σ_yy(x) for several
values of the dimensionless weight γ·B/c.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from slip_line.solver import MohrCoulombGravityPrandtlSolver


def _draw_field(ax, solver, title):
    alpha_lines = solver.all_alpha_lines()
    beta_lines  = solver.all_beta_lines()
    for ln in alpha_lines:
        ax.plot([n.x for n in ln], [n.y for n in ln],
                color='royalblue', lw=0.7)
    for ln in beta_lines:
        ax.plot([n.x for n in ln], [n.y for n in ln],
                color='tomato', lw=0.7)
    B = solver.B
    ax.plot([0, B], [0, 0], color='black', lw=3.5, zorder=10)
    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.set_xlabel('x / B')
    ax.set_ylabel('depth y / B  (down)')
    ax.set_title(title)
    ax.grid(alpha=0.25)


def main():
    c = 1.0
    phi = np.deg2rad(30.0)
    B = 1.0
    n_fan = 24
    n_radial = 24

    gammas = [0.0, 1.0, 2.5, 5.0]   # γ·B/c values (since c = B = 1)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    print(f"Mohr-Coulomb strip footing, φ = 30°, c = {c}, B = {B}\n")
    print(f"{'γB/c':>6}  {'Q/(cB)':>10}  {'σ_yy(0)/c':>11}  "
          f"{'σ_yy(B)/c':>11}  {'q_avg/c':>10}")
    print('-' * 60)

    for ax, gamma in zip(axes.flat, gammas):
        s = MohrCoulombGravityPrandtlSolver(
            c=c, phi=phi, gamma=gamma, B=B,
            n_fan=n_fan, n_radial=n_radial,
        )
        res = s.solve()
        x, sn = s.contact_pressure_curve()
        Q = res['q_per_unit_length_half']
        print(f"{gamma:>6.2f}  {Q/(c*B):>10.4f}  {sn[0]/c:>11.4f}  "
              f"{sn[-1]/c:>11.4f}  {Q/(c*B):>10.4f}")

        _draw_field(
            ax, s,
            f'γB/c = {gamma:g}   '
            f'Q/(cB) = {Q/(c*B):.3f}',
        )

    fig.suptitle(
        'Mohr-Coulomb strip footing with submerged weight  '
        '(φ = 30°, smooth, rigorous MOC with gravity)'
    )
    fig.tight_layout()
    out = 'mc_gravity_phi30.png'
    fig.savefig(out, dpi=140)
    print(f"\nSaved: {out}")

    # Contact-pressure distributions (single panel)
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    for gamma in gammas:
        s = MohrCoulombGravityPrandtlSolver(
            c=c, phi=phi, gamma=gamma, B=B,
            n_fan=n_fan, n_radial=n_radial,
        )
        s.solve()
        x, sn = s.contact_pressure_curve()
        ax2.plot(x/B, sn/c, lw=1.8, label=f'γB/c = {gamma:g}')
    ax2.set_xlabel('x / B')
    ax2.set_ylabel(r'$\sigma_{yy}/c$')
    ax2.set_title('Contact pressure under smooth strip footing (MC, φ = 30°)')
    ax2.grid(alpha=0.3)
    ax2.legend()
    fig2.tight_layout()
    out2 = 'mc_gravity_contact_pressure.png'
    fig2.savefig(out2, dpi=140)
    print(f"Saved: {out2}")


if __name__ == '__main__':
    main()
