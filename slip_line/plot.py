"""
Matplotlib visualisation of the Prandtl slip line field.

Internal y is depth (positive down); we flip for display.
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from .solver import (
    PrandtlSolver, DiskPrandtlSolver, MohrCoulombPrandtlSolver,
    MohrCoulombDiskPrandtlSolver,
)


def plot_slip_lines(
    solver: PrandtlSolver,
    q: float,
    show: bool = True,
    title: str | None = None,
):
    """Draw α/β slip lines, foundation, free surface, and annotations."""
    cu = solver.cu
    B  = solver.B

    fig, ax = plt.subplots(figsize=(11, 6.5))

    alpha_lines = solver.all_alpha_lines()
    beta_lines  = solver.all_beta_lines()

    def _draw(lines, color, label):
        first = True
        for line in lines:
            xs = [n.x for n in line]
            ys = [n.y for n in line]
            ax.plot(xs, ys, color=color, linewidth=0.9,
                    label=label if first else None)
            first = False

    _draw(alpha_lines, 'royalblue', 'α-lines')
    _draw(beta_lines,  'tomato',    'β-lines')

    # --- Foundation: thick black bar (half, from symmetry axis to corner) ---
    ax.plot([0, B], [0, 0], color='black', linewidth=3.5, zorder=10)
    ax.fill_between([0, B], [0, 0], [-0.04 * B, -0.04 * B],
                    color='dimgray', alpha=0.7, zorder=9)

    # Hatched foundation
    for xi in np.linspace(0, B, 7):
        ax.plot([xi, xi - 0.04 * B], [-0.04 * B, -0.10 * B],
                color='dimgray', linewidth=1.0, zorder=9)

    # --- Pressure arrows on foundation (pointing down into soil) ---
    n_arrows = 5
    for xi in np.linspace(0.05 * B, 0.95 * B, n_arrows):
        ax.annotate(
            '', xy=(xi, 0.02 * B), xytext=(xi, -0.32 * B),
            arrowprops=dict(arrowstyle='->', color='darkgreen', lw=1.4),
            zorder=8,
        )

    # --- Free-surface line and ground markers ---
    x_max = _max_extent(alpha_lines, beta_lines)
    ax.plot([B, x_max], [0, 0], color='saddlebrown', linewidth=1.2)
    # --- Symmetry axis ---
    y_min_axis = max((n.y for ln in alpha_lines + beta_lines for n in ln),
                     default=B)
    ax.plot([0, 0], [0, y_min_axis], color='gray', linestyle='-.',
            linewidth=1.0, alpha=0.6)

    # --- Annotations ---
    fan_deg = np.degrees(solver.delta_psi)
    txt = (f'q / cᵤ = {q / cu:.4f}\n'
           f'(π + 2 = {np.pi + 2:.4f})\n'
           f'fan angle = {fan_deg:.2f}°\n'
           f'n_fan = {solver.n_fan},  r = {solver.r:.2f}')
    ax.text(0.02, 0.98, txt, transform=ax.transAxes,
            va='top', ha='left', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))

    ax.text(0.5 * B, -0.45 * B, f'q = {q / cu:.3f} · cᵤ',
            color='darkgreen', fontsize=12, ha='center', va='bottom')

    # --- Axes ---
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('depth y', fontsize=12)
    ax.invert_yaxis()  # display y downward
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.grid(True, linestyle=':', alpha=0.4)
    ax.set_title(
        title or f'Prandtl slip line field — Tresca, cᵤ = {cu:g}, B = {B:g}',
        fontsize=12,
    )
    ax.legend(loc='lower right', fontsize=10)
    plt.tight_layout()

    if show:
        plt.show()
    return fig


def plot_mc_slip_lines(
    solver: MohrCoulombPrandtlSolver,
    show: bool = True,
    title: str | None = None,
):
    """
    Draw the Mohr-Coulomb Prandtl slip-line field.

    Shows α/β characteristics, foundation, free surface, and a text box
    with q/c, Nc (analytical), φ, and mesh parameters.
    """
    c   = solver.c
    phi = solver.phi
    B   = solver.B
    q   = solver.collapse_load()
    Nc  = solver.analytical_Nc()

    fig, ax = plt.subplots(figsize=(11, 6.5))

    alpha_lines = solver.all_alpha_lines()
    beta_lines  = solver.all_beta_lines()

    def _draw(lines, color, label):
        first = True
        for line in lines:
            xs = [n.x for n in line]
            ys = [n.y for n in line]
            ax.plot(xs, ys, color=color, linewidth=0.9,
                    label=label if first else None)
            first = False

    _draw(alpha_lines, 'royalblue', 'α-lines')
    _draw(beta_lines,  'tomato',    'β-lines')

    # Foundation bar
    ax.plot([0, B], [0, 0], color='black', linewidth=3.5, zorder=10)
    ax.fill_between([0, B], [0, 0], [-0.04 * B, -0.04 * B],
                    color='dimgray', alpha=0.7, zorder=9)
    for xi in np.linspace(0, B, 7):
        ax.plot([xi, xi - 0.04 * B], [-0.04 * B, -0.10 * B],
                color='dimgray', linewidth=1.0, zorder=9)

    # Pressure arrows
    for xi in np.linspace(0.05 * B, 0.95 * B, 5):
        ax.annotate(
            '', xy=(xi, 0.02 * B), xytext=(xi, -0.32 * B),
            arrowprops=dict(arrowstyle='->', color='darkgreen', lw=1.4),
            zorder=8,
        )

    # Free surface
    x_max = _max_extent(alpha_lines, beta_lines)
    ax.plot([B, x_max], [0, 0], color='saddlebrown', linewidth=1.2)

    # Symmetry axis
    y_min_axis = max((n.y for ln in alpha_lines + beta_lines for n in ln),
                     default=B)
    ax.plot([0, 0], [0, y_min_axis], color='gray', linestyle='-.',
            linewidth=1.0, alpha=0.6)

    # Annotations
    phi_deg = np.degrees(phi)
    txt = (f'q / c  = {q / c:.4f}\n'
           f'Nc (analyt.) = {Nc:.4f}\n'
           f'φ = {phi_deg:.1f}°\n'
           f'n_fan = {solver.n_fan}')
    ax.text(0.02, 0.98, txt, transform=ax.transAxes,
            va='top', ha='left', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))

    ax.text(0.5 * B, -0.45 * B, f'q = {q / c:.3f} · c',
            color='darkgreen', fontsize=12, ha='center', va='bottom')

    # Axes
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('depth y', fontsize=12)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.grid(True, linestyle=':', alpha=0.4)
    ax.set_title(
        title or f'Prandtl slip-line field — Mohr-Coulomb, φ = {phi_deg:.1f}°, c = {c:g}, B = {B:g}',
        fontsize=12,
    )
    ax.legend(loc='lower right', fontsize=10)
    plt.tight_layout()

    if show:
        plt.show()
    return fig


def _max_extent(alpha_lines, beta_lines) -> float:
    vals = []
    for ln in alpha_lines + beta_lines:
        for n in ln:
            vals.append(n.x)
    return max(vals) if vals else 1.0


# ======================================================================
# Disk foundation plot
# ======================================================================


def plot_disk_slip_lines(
    solver: DiskPrandtlSolver,
    show: bool = True,
    title: str | None = None,
):
    """
    Two-panel figure for a partly-penetrated disk:

      Left  — slip line field with disk shown.
      Right — normal contact pressure σ_n(φ)/cu along the disk surface.
    """
    cu  = solver.cu
    D   = solver.D
    R   = solver.R
    y_c = solver.y_c
    p   = solver.p
    r   = solver.r
    tmax = solver.theta_max
    B_s = solver.B_s

    fig, (ax, axp) = plt.subplots(
        1, 2, figsize=(14, 6.5), gridspec_kw=dict(width_ratios=[2.2, 1.0])
    )

    # ---- slip lines (right half + mirror) ----
    alpha_lines = solver.all_alpha_lines()
    beta_lines  = solver.all_beta_lines()

    def _draw(lines, color, label):
        first = True
        for line in lines:
            xs = [n.x for n in line]
            ys = [n.y for n in line]
            ax.plot(xs, ys, color=color, linewidth=0.9,
                    label=label if first else None)
            first = False

    _draw(alpha_lines, 'royalblue', 'α-lines')
    _draw(beta_lines,  'tomato',    'β-lines')

    # ---- Disk (the buried portion only) — right half ----
    # Draw the full circle outline (right half), then highlight the contact arc.
    n_arc = 200
    phi_full = np.linspace(-np.pi / 2, np.pi / 2, n_arc)
    ax.plot(R * np.cos(phi_full), y_c + R * np.sin(phi_full),
            color='black', linewidth=2.0, zorder=10)
    # Light fill of the disk (right half via a wedge patch)
    disk_patch = mpatches.Wedge((0.0, y_c), R, -90, 90,
                                facecolor='none',
                                edgecolor='black', linewidth=2.0, zorder=9)
    ax.add_patch(disk_patch)

    # Contact arc (φ ∈ [0, θmax] from downward axis) — emphasize
    phi_contact = np.linspace(0.0, tmax, 100)
    ax.plot(R * np.sin(phi_contact), y_c + R * np.cos(phi_contact),
            color='darkred', linewidth=3.0, zorder=11)

    # ---- Free surface (right of the disk contact) ----
    x_max = _max_extent(alpha_lines, beta_lines) if (alpha_lines or beta_lines) else 1.5 * D
    ax.plot([B_s, x_max], [0, 0], color='saddlebrown', linewidth=1.2)
    # ---- Symmetry axis ----
    y_min_axis = max((n.y for ln in alpha_lines + beta_lines for n in ln),
                     default=R)
    ax.plot([0, 0], [y_c - R, max(y_min_axis, y_c + R)],
            color='gray', linestyle='-.', linewidth=1.0, alpha=0.6)

    # ---- Applied load arrow on disk ----
    ax.annotate('', xy=(0, y_c - 0.05 * D), xytext=(0, y_c - 0.55 * D),
                arrowprops=dict(arrowstyle='->', color='darkgreen', lw=2.5))
    ax.text(0.05 * D, y_c - 0.30 * D, 'Q',
            color='darkgreen', fontsize=14, va='center')

    # ---- Annotations ----
    Q = solver.collapse_load_Q()
    sn_bot = solver._sigma_n(0.0)
    sn_edge = solver._sigma_n(tmax)
    forces = solver.forces_on_half_disk()
    Fv = forces['F_v']
    Fh = forces['F_h']
    info = (f'p = {p:.3f}     D = {D:g}     r = {r:.2f}\n'
            f'θ_max = {np.degrees(tmax):.2f}°\n'
            f'fan angle = {np.degrees(solver.delta_psi_fan):.2f}°\n'
            f'B_s / D = {B_s / D:.4f}\n'
            f'\n'
            f'σ_n(bottom) / cᵤ = {sn_bot/cu:.4f}\n'
            f'σ_n(edge)   / cᵤ = {sn_edge/cu:.4f}\n'
            f'\n'
            f'Half-disk (right side, per unit length):\n'
            f'  F_v / (D·cᵤ) = {Fv / (D * cu):.4f}\n'
            f'  F_h / (D·cᵤ) = {Fh / (D * cu):.4f}\n'
            f'\n'
            f'Full disk:\n'
            f'  Q / D       = {Q / D:.4f}\n'
            f'  Q / (D·cᵤ)  = {Q / (D * cu):.4f}   (= 2·F_v)')
    ax.text(0.02, 0.98, info, transform=ax.transAxes,
            va='top', ha='left', fontsize=10, family='monospace',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('depth y', fontsize=12)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.grid(True, linestyle=':', alpha=0.4)
    ax.set_title(
        title or f'Disk foundation — Tresca, p = {p:.3f}, r = {r:.2f}',
        fontsize=12,
    )
    ax.legend(loc='lower right', fontsize=10)

    # ---- right panel: stress distribution σ_n(φ) and τ(φ) ----
    phis, sn, tau = solver.stress_distribution_curve(200)
    axp.plot(np.degrees(phis), sn / cu, color='darkred', linewidth=2.0,
             label='σ_n / cᵤ  (normal)')
    axp.fill_between(np.degrees(phis), 0, sn / cu, color='darkred', alpha=0.15)
    axp.plot(np.degrees(phis), tau / cu, color='steelblue', linewidth=2.0,
             linestyle='--', label=f'τ / cᵤ = r = {r:.2f}  (shear)')
    axp.axhline(np.pi + 2, color='gray', linestyle=':', linewidth=1,
                label=f'π + 2 ≈ {np.pi + 2:.3f}  (flat smooth)')
    axp.set_xlabel('φ  [deg, 0 = disk bottom]', fontsize=11)
    axp.set_ylabel('stress / cᵤ', fontsize=11)
    axp.set_title('Stress distribution', fontsize=12)
    axp.grid(True, linestyle=':', alpha=0.5)
    axp.set_xlim(0, max(1, np.degrees(tmax)))
    axp.set_ylim(0, max(np.pi + 2, sn.max() / cu) * 1.1)
    axp.legend(loc='upper right', fontsize=9)

    plt.tight_layout()
    if show:
        plt.show()
    return fig


# ======================================================================
# Mohr-Coulomb disk foundation plot
# ======================================================================


def plot_mc_disk_slip_lines(
    solver: MohrCoulombDiskPrandtlSolver,
    show: bool = True,
    title: str | None = None,
):
    """
    Two-panel figure for a partly-penetrated disk on Mohr-Coulomb soil.

      Left  — slip line field with disk shown.
      Right — normal σ_n(θ)/c and shear |τ(θ)|/c along the contact arc.
    """
    c    = solver.c
    phi  = solver.phi
    D    = solver.D
    R    = solver.R
    y_c  = solver.y_c
    p    = solver.p
    r    = solver.r
    tmax = solver.theta_max
    B_s  = solver.B_s

    fig, (ax, axp) = plt.subplots(
        1, 2, figsize=(14, 6.5), gridspec_kw=dict(width_ratios=[2.2, 1.0])
    )

    alpha_lines = solver.all_alpha_lines()
    beta_lines  = solver.all_beta_lines()

    def _draw(lines, color, label):
        first = True
        for line in lines:
            xs = [n.x for n in line]
            ys = [n.y for n in line]
            ax.plot(xs, ys, color=color, linewidth=0.9,
                    label=label if first else None)
            first = False

    _draw(alpha_lines, 'royalblue', 'α-lines')
    _draw(beta_lines,  'tomato',    'β-lines')

    # ---- Disk outline (right-half) ----
    n_arc = 200
    phi_full = np.linspace(-np.pi / 2, np.pi / 2, n_arc)
    ax.plot(R * np.cos(phi_full), y_c + R * np.sin(phi_full),
            color='black', linewidth=2.0, zorder=10)
    disk_patch = mpatches.Wedge((0.0, y_c), R, -90, 90,
                                facecolor='none',
                                edgecolor='black', linewidth=2.0, zorder=9)
    ax.add_patch(disk_patch)

    # Contact arc (θ ∈ [0, θmax])
    th_contact = np.linspace(0.0, tmax, 100)
    ax.plot(R * np.sin(th_contact), y_c + R * np.cos(th_contact),
            color='darkred', linewidth=3.0, zorder=11)

    # ---- Free surface ----
    x_max = _max_extent(alpha_lines, beta_lines) if (alpha_lines or beta_lines) else 1.5 * D
    ax.plot([B_s, x_max], [0, 0], color='saddlebrown', linewidth=1.2)
    # ---- Symmetry axis ----
    y_min_axis = max((n.y for ln in alpha_lines + beta_lines for n in ln),
                     default=R)
    ax.plot([0, 0], [y_c - R, max(y_min_axis, y_c + R)],
            color='gray', linestyle='-.', linewidth=1.0, alpha=0.6)

    # ---- Load arrow ----
    ax.annotate('', xy=(0, y_c - 0.05 * D), xytext=(0, y_c - 0.55 * D),
                arrowprops=dict(arrowstyle='->', color='darkgreen', lw=2.5))
    ax.text(0.05 * D, y_c - 0.30 * D, 'Q',
            color='darkgreen', fontsize=14, va='center')

    # ---- Annotations ----
    Q       = solver.collapse_load_Q()
    sn_bot  = solver._sigma_n(0.0)
    sn_edge = solver._sigma_n(tmax)
    forces  = solver.forces_on_half_disk()
    Fv = forces['F_v']
    Fh = forces['F_h']
    Nc_flat = solver.yield_fn.Nc()
    info = (f'φ = {np.degrees(phi):.1f}°    c = {c:g}    D = {D:g}\n'
            f'p = {p:.3f}    r = {r:.2f}    θ_max = {np.degrees(tmax):.2f}°\n'
            f'fan angle = {np.degrees(solver.delta_psi_fan):.2f}°\n'
            f'B_s / D = {B_s / D:.4f}\n'
            f'\n'
            f'σ_n(bottom) / c = {sn_bot/c:.4f}\n'
            f'σ_n(edge)   / c = {sn_edge/c:.4f}\n'
            f'Nc (flat smooth) = {Nc_flat:.4f}\n'
            f'\n'
            f'Half-disk (right side, per unit length):\n'
            f'  F_v / (D·c) = {Fv / (D * c):.4f}\n'
            f'  F_h / (D·c) = {Fh / (D * c):.4f}\n'
            f'\n'
            f'Full disk:\n'
            f'  Q / (D·c)   = {Q / (D * c):.4f}   (= 2·F_v)')
    ax.text(0.02, 0.98, info, transform=ax.transAxes,
            va='top', ha='left', fontsize=10, family='monospace',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('depth y', fontsize=12)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.grid(True, linestyle=':', alpha=0.4)
    ax.set_title(
        title or f'Disk foundation — Mohr-Coulomb, φ = {np.degrees(phi):.1f}°, '
                 f'c = {c:g}, p = {p:.3f}, r = {r:.2f}',
        fontsize=12,
    )
    ax.legend(loc='lower right', fontsize=10)

    # ---- Right panel: stress distribution ----
    ths, sn, tau = solver.stress_distribution_curve(300)
    axp.plot(np.degrees(ths), sn / c, color='darkred', linewidth=2.0,
             label='σ_n / c  (normal)')
    axp.fill_between(np.degrees(ths), 0, sn / c, color='darkred', alpha=0.15)
    axp.plot(np.degrees(ths), tau / c, color='steelblue', linewidth=2.0,
             linestyle='--', label=f'|τ| / c   (shear, r = {r:.2f})')
    axp.axhline(Nc_flat, color='gray', linestyle=':', linewidth=1,
                label=f'Nc = {Nc_flat:.3f}  (flat smooth)')
    axp.set_xlabel('θ  [deg, 0 = disk bottom]', fontsize=11)
    axp.set_ylabel('stress / c', fontsize=11)
    axp.set_title('Stress distribution', fontsize=12)
    axp.grid(True, linestyle=':', alpha=0.5)
    axp.set_xlim(0, max(1, np.degrees(tmax)))
    axp.set_ylim(0, max(Nc_flat, sn.max() / c) * 1.1)
    axp.legend(loc='upper right', fontsize=9)

    plt.tight_layout()
    if show:
        plt.show()
    return fig
