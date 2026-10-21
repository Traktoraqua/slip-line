"""
Matplotlib visualisation of the Prandtl slip line field.

Internal y is depth (positive down); we flip for display.
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PolyCollection

from .solver import (
    PrandtlSolver, DiskPrandtlSolver, MohrCoulombPrandtlSolver,
    MohrCoulombDiskPrandtlSolver,
)
from .velocity import DiskVelocityField, PrandtlVelocityField


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


# Shared colours for the velocity-field figures.
_C_ARROW = '#2a78d6'
_C_JUMP = '#eb6834'
_C_NET = '#c8c7c2'
_C_RIGID = '#efeeea'
_C_INK = '#52514e'
_C_NEG = '#f2b8b7'


def _velocity_background(ax, solver, outer, zones, x_lim, y_lim, jump):
    """Rigid soil, plastic zone(s), mirrored slip-line net and the
    velocity-discontinuity line `outer` (right half)."""
    ax.add_patch(mpatches.Rectangle(
        (-x_lim, 0), 2 * x_lim, y_lim, facecolor=_C_RIGID, hatch='///',
        edgecolor='#dddcd6', linewidth=0, zorder=0, label='rigid soil, V = 0',
    ))
    for zx, zy in zones:
        ax.fill(zx, zy, facecolor='white', edgecolor='none', zorder=1)
    for line in solver.all_alpha_lines() + solver.all_beta_lines():
        xs = np.array([nd.x for nd in line])
        ys = np.array([nd.y for nd in line])
        for sgn in (1, -1):
            ax.plot(sgn * xs, ys, color=_C_NET, linewidth=0.6, zorder=2)
    ox = np.array([nd.x for nd in outer])
    oy = np.array([nd.y for nd in outer])
    for sgn, lbl in ((1, f'velocity jump |ΔV| = {jump:.3f}·V₀'), (-1, None)):
        ax.plot(sgn * ox, oy, color=_C_JUMP, linewidth=2.2, zorder=4, label=lbl)
    ax.plot([-x_lim, x_lim], [0, 0], color='saddlebrown', linewidth=1.0,
            zorder=3)


def _fan_arrow_samples(add, fan_samples, n_fan, n_r, dpsi, step):
    """Fan nodes [k][j] for arrows: every `step`-th radius (corner skipped),
    thinning radials so the arc spacing stays ≈ step·dr."""
    it = iter(fan_samples)
    for k in range(n_fan + 1):
        for j in range(n_r + 1):
            nd, vx, vy = next(it)
            if j == 0 or j % step:
                continue
            if k % max(1, round(step / (j * dpsi))) == 0:
                add(nd, vx, vy)


def _velocity_arrows(ax, pts, arrow_len, v0):
    """Mirrored quiver of {key: (x, y, vx, vy)}; the fastest arrow has
    length `arrow_len`."""
    P = np.array(list(pts.values()))
    speed = np.hypot(P[:, 2], P[:, 3])
    speed_max = speed.max()
    q_kw = dict(angles='xy', scale_units='xy', scale=speed_max / arrow_len,
                color=_C_ARROW, width=0.0032, headwidth=4, headlength=4.5,
                zorder=5)
    if np.ptp(speed) < 1e-6 * speed_max:
        v_lbl = f'soil velocity, |V| = {speed_max / v0:.3f}·V₀'
    else:
        v_lbl = (f'soil velocity, |V| = {speed.min() / v0:.3f}'
                 f'–{speed_max / v0:.3f}·V₀')
    ax.quiver(P[:, 0], P[:, 1], P[:, 2], P[:, 3], label=v_lbl, **q_kw)
    ax.quiver(-P[:, 0], P[:, 1], -P[:, 2], P[:, 3], **q_kw)


def _region_labels(ax, samples, names, fan_shift=1.35):
    """Bold region names at the centroid of each region's nodes."""
    for name in names:
        nds = [smp[0] for smp in samples[name]]
        cx = np.mean([nd.x for nd in nds])
        cy = np.mean([nd.y for nd in nds])
        if name == 'II':
            cy *= fan_shift
        ax.text(cx, cy, name, fontsize=13, fontweight='bold', color=_C_INK,
                ha='center', va='center', zorder=8,
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                          edgecolor='none', alpha=0.85))


def _hodograph(axh, samples, v0, v_jump, lim, markersize=4):
    """Mirrored hodograph of all node velocities with O (rigid soil),
    F (foundation) and the jump vector O → v_jump (right half)."""
    H = np.array([(vx, vy) for grp in samples.values()
                  for _, vx, vy in grp]) / v0
    for sgn in (1, -1):
        axh.plot(sgn * H[:, 0], H[:, 1], 'o', color=_C_ARROW,
                 markersize=markersize, zorder=4)
        axh.annotate('', xy=(sgn * v_jump[0], v_jump[1]), xytext=(0, 0),
                     arrowprops=dict(arrowstyle='-|>', color=_C_JUMP, lw=1.8))
    axh.plot([0, 0], [0, 1], color='black', linewidth=2)
    axh.plot(0, 0, 'ks', markersize=7)
    axh.plot(0, 1, 'ks', markersize=7)
    axh.text(0.06, -0.06, 'O  rigid soil', ha='left', va='bottom',
             fontsize=10, color=_C_INK)
    axh.text(0.5 * v_jump[0], 0.5 * v_jump[1] - 0.08, 'jump', ha='right',
             va='bottom', fontsize=9, color=_C_JUMP)
    axh.set_xlim(-lim, lim)
    axh.set_ylim(lim, -lim)   # y down
    axh.set_aspect('equal')
    axh.axhline(0, color=_C_NET, linewidth=0.8, zorder=0)
    axh.axvline(0, color=_C_NET, linewidth=0.8, zorder=0)
    axh.grid(True, linestyle=':', alpha=0.3)
    axh.set_xlabel('vx / V₀', fontsize=11)
    axh.set_ylabel('vy / V₀  (down +)', fontsize=11)
    axh.set_title('Hodograph', fontsize=12)


def plot_prandtl_velocity(
    vf: PrandtlVelocityField,
    show: bool = True,
    title: str | None = None,
    step: int = 2,
):
    """
    Velocity field of a PrandtlVelocityField: arrows on the (mirrored)
    mechanism, velocity discontinuity and rigid soil (left), and the
    hodograph (right).  Arrows are drawn at every `step`-th mesh node.
    """
    solver = vf.solver
    B = solver.B
    v0 = vf.v0

    fig, (ax, axh) = plt.subplots(
        1, 2, figsize=(15, 6.2), gridspec_kw={'width_ratios': [2.3, 1]},
    )

    # ---- Rigid soil, plastic zone, net, velocity discontinuity ----
    outer = vf.outer_alpha_line()
    ox = np.array([nd.x for nd in outer])
    oy = np.array([nd.y for nd in outer])
    x_lim = 1.12 * ox.max()
    y_lim = 1.25 * oy.max()
    zone = (np.concatenate([-ox[::-1], ox[1:]]),
            np.concatenate([oy[::-1], oy[1:]]))
    _velocity_background(ax, solver, outer, [zone], x_lim, y_lim, vf.U / v0)

    # ---- Arrow samples: regular lattice on the net, skipping the corner ----
    pts: dict[tuple[float, float], tuple[float, float, float, float]] = {}

    def _add(nd, vx, vy):
        pts[(round(nd.x, 9), round(nd.y, 9))] = (nd.x, nd.y, vx, vy)

    samples = vf.region_samples()
    N_I = len(solver.region_I_nodes) - 1
    it = iter(samples['I'])
    for i in range(N_I + 1):
        for j in range(i + 1):
            nd, vx, vy = next(it)
            if i % step == 0 and j % step == 0:
                _add(nd, vx, vy)
    n_p = len(solver.region_III_nodes) - 1
    it = iter(samples['III'])
    for j in range(n_p + 1):
        for i in range(n_p + 1 - j):
            nd, vx, vy = next(it)
            if i % step == 0 and j % step == 0:
                _add(nd, vx, vy)
    _fan_arrow_samples(_add, samples['II'], solver.n_fan,
                       len(solver.region_II_nodes[0]) - 1,
                       solver.delta_psi / solver.n_fan, step)
    _velocity_arrows(ax, pts, 0.85 * step * solver.dr, v0)

    # ---- Footing and its velocity ----
    ax.fill_between([-B, B], [0, 0], [-0.06 * B, -0.06 * B],
                    color='dimgray', zorder=6)
    ax.plot([-B, B], [0, 0], color='black', linewidth=3, zorder=7)
    ax.annotate('', xy=(0, -0.08 * B), xytext=(0, -0.42 * B),
                arrowprops=dict(arrowstyle='-|>', color='black', lw=2),
                zorder=7)
    ax.text(0.05 * B, -0.27 * B, 'V₀', fontsize=13, va='center')

    _region_labels(ax, samples, ('I', 'II', 'III'))

    ax.set_xlim(-x_lim, x_lim)
    ax.set_ylim(y_lim, -0.55 * B)   # y down
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=11)
    ax.set_ylabel('depth y', fontsize=11)
    ax.grid(True, linestyle=':', alpha=0.3)
    ax.legend(loc='upper right', fontsize=9, framealpha=0.95)
    ax.set_title(
        title or f'Velocity field — smooth strip footing, Tresca, '
                 f'B = {B:g}, V₀ = {v0:g}',
        fontsize=12,
    )

    # ---- Hodograph (same y-down orientation as the physical plane) ----
    def _mean_v(name):
        V = np.array([(vx, vy) for _, vx, vy in samples[name]]) / v0
        return V.mean(axis=0)

    vI, vIII = _mean_v('I'), _mean_v('III')
    r_mid = vf.U / v0
    _hodograph(axh, samples, v0, vI, 1.25 * r_mid)
    for sgn in (1, -1):
        axh.plot([0, sgn * vI[0]], [1, vI[1]], color=_C_INK, linestyle='--',
                 linewidth=1)

    lbl = dict(fontsize=10, color=_C_INK)
    axh.text(-0.06, 1.06, 'F  footing', ha='right', va='top', **lbl)
    axh.text(vI[0] + 0.08, vI[1], 'I', ha='left', va='center', **lbl)
    axh.text(vIII[0] + 0.08, vIII[1], 'III', ha='left', va='center', **lbl)
    psi_mid = 0.5 * (solver.psi_A + solver.psi_B)
    axh.text((r_mid + 0.12) * np.cos(psi_mid - np.pi / 4),
             (r_mid + 0.12) * np.sin(psi_mid - np.pi / 4), 'II (fan)',
             ha='left', va='center', **lbl)
    axh.text(0.5 * vI[0], 1.06, 'slip under\nfooting', ha='center',
             va='top', fontsize=9, color=_C_INK)

    plt.tight_layout()
    if show:
        plt.show()
    return fig


def plot_disk_velocity(
    vf: DiskVelocityField,
    show: bool = True,
    title: str | None = None,
    step: int = 3,
):
    """
    Velocity field of a DiskVelocityField: arrows on the (mirrored)
    mechanism, velocity discontinuity, rigid soil and the cells where the
    plastic work against the slip-line stress field is negative (left),
    and the hodograph (right).  Arrows are drawn at every `step`-th node.
    """
    if vf.horizontal:
        raise NotImplementedError("plot_disk_velocity mirrors the field; vertical push only")
    solver = vf.solver
    D, R, y_c = solver.D, solver.R, solver.y_c
    v0 = vf.v0
    energy = vf.energy_balance()
    pw = vf.plastic_work()

    fig, (ax, axh) = plt.subplots(
        1, 2, figsize=(15, 6.6), gridspec_kw={'width_ratios': [2.0, 1]},
    )

    # ---- Rigid soil, plastic zone, net, velocity discontinuity ----
    outer = vf.outer_alpha_line()
    zone_nodes = outer + solver.disk_surface_nodes    # surface end → corner → bottom
    zx = np.array([nd.x for nd in zone_nodes])
    zy = np.array([nd.y for nd in zone_nodes])
    x_lim = 1.12 * zx.max()
    y_lim = 1.3 * zy.max()
    y_top = -0.4 * x_lim     # crop the disk above the surface
    _velocity_background(ax, solver, outer, [(zx, zy), (-zx, zy)],
                         x_lim, y_lim, vf.U / v0)

    # ---- Negative plastic work ----
    neg = pw['xy'][pw['d'] < 0]
    if len(neg):
        mirrored = neg * np.array([-1.0, 1.0])
        ax.add_collection(PolyCollection(
            list(neg) + list(mirrored), facecolor=_C_NEG, edgecolor=_C_NEG,
            linewidth=0.3, zorder=1.5, label='negative plastic work',
        ))

    # ---- Arrow samples: regular lattice on each region's net ----
    pts: dict[tuple[float, float], tuple[float, float, float, float]] = {}

    def _add(nd, vx, vy):
        pts[(round(nd.x, 9), round(nd.y, 9))] = (nd.x, nd.y, vx, vy)

    samples = vf.region_samples()
    it = iter(samples['I'])
    for k, row in enumerate(solver.region_I_mesh):
        for i in range(len(row)):
            nd, vx, vy = next(it)
            if i % step == 0 and (i + k) % step == 0:
                _add(nd, vx, vy)
    n_p = len(solver.region_III_nodes) - 1
    it = iter(samples['III'])
    for j in range(n_p + 1):
        for i in range(n_p + 1 - j):
            nd, vx, vy = next(it)
            if i % step == 0 and j % step == 0:
                _add(nd, vx, vy)
    fan = solver.region_II_nodes
    if len(fan) >= 2:
        _fan_arrow_samples(_add, samples['II'], len(fan) - 1, len(fan[0]) - 1,
                           solver.delta_psi_fan / (len(fan) - 1), step)
    _velocity_arrows(ax, pts, 0.85 * step * solver.dr, v0)

    # ---- Disk and its velocity ----
    ax.add_patch(mpatches.Circle((0.0, y_c), R, facecolor='#d9d8d3',
                                 edgecolor='black', linewidth=1.8, zorder=6))
    y_a = y_top + 0.06 * x_lim
    ax.annotate('', xy=(0, y_a + 0.22 * x_lim), xytext=(0, y_a),
                arrowprops=dict(arrowstyle='-|>', color='black', lw=2),
                zorder=7)
    ax.text(0.03 * x_lim, y_a + 0.11 * x_lim, 'V₀', fontsize=13, va='center',
            zorder=7)

    _region_labels(ax, samples, ('I', 'II', 'III') if len(fan) >= 2 else ('I', 'III'),
                   fan_shift=1.2)

    info = (f'p = {solver.p:.2f},  θ_max = {np.degrees(solver.theta_max):.1f}°\n'
            f'Q / (D·cᵤ):\n'
            f'  stress field      {energy["F_stress"] / (D * solver.cu):.3f}\n'
            f'  virtual work      {energy["F_virtual_work"] / (D * solver.cu):.3f}\n'
            f'  upper bound       {energy["F_upper_bound"] / (D * solver.cu):.3f}\n'
            f'negative-work area  {100 * energy["negative_area_fraction"]:.0f} %')
    ax.text(0.015, 0.985, info, transform=ax.transAxes, va='top', ha='left',
            fontsize=9, family='monospace', color=_C_INK, zorder=9,
            bbox=dict(boxstyle='round', facecolor='white', edgecolor=_C_NET))

    ax.set_xlim(-x_lim, x_lim)
    ax.set_ylim(y_lim, y_top)   # y down
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=11)
    ax.set_ylabel('depth y', fontsize=11)
    ax.grid(True, linestyle=':', alpha=0.3)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2,
              fontsize=9, frameon=False)
    ax.set_title(
        title or f'Velocity field — smooth disk, Tresca, D = {D:g}, '
                 f'p = {solver.p:g}, V₀ = {v0:g}',
        fontsize=12,
    )

    # ---- Hodograph ----
    bottom = solver.region_I_mesh[0][0]
    th = bottom.psi - np.pi / 4.0
    v_jump = (vf.U / v0 * np.cos(th), vf.U / v0 * np.sin(th))
    _hodograph(axh, samples, v0, v_jump, 1.25 * vf.U / v0, markersize=2)
    axh.text(-0.06, 1.06, 'F  disk', ha='right', va='top', fontsize=10,
             color=_C_INK)
    for name in samples:
        if name == 'II' and len(fan) < 2:
            continue
        V = np.array([(vx, vy) for _, vx, vy in samples[name]]) / v0
        ang = np.arctan2(V[:, 1].mean(), V[:, 0].mean())
        r_lbl = vf.U / v0 + 0.12
        axh.text(r_lbl * np.cos(ang), r_lbl * np.sin(ang), name, ha='left',
                 va='center', fontsize=10, color=_C_INK)

    plt.tight_layout()
    if show:
        plt.show()
    return fig
