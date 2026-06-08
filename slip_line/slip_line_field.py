"""
Core data structures for the slip line characteristic net.

Coordinate convention: y points DOWNWARD (depth).
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class Node:
    """A node in the slip line field."""
    x: float
    y: float
    sigma: float   # mean stress
    psi: float     # principal-stress angle [rad]


def alpha_slope(psi: float) -> float:
    """Slope dy/dx of α-line at orientation ψ (y-down convention)."""
    return np.tan(psi - np.pi / 4.0)


def beta_slope(psi: float) -> float:
    """Slope dy/dx of β-line at orientation ψ (y-down convention)."""
    return np.tan(psi + np.pi / 4.0)


def riemann_step(P: Node, Q: Node, cu: float) -> Node:
    """
    Riemann finite-difference step for the Hencky system.

      P : α-parent of R  (σ - 2cu·ψ const along α through P)
      Q : β-parent of R  (σ + 2cu·ψ const along β through Q)

    Stress update (exact for linear Hencky):
        σ_R - 2cu·ψ_R = σ_P - 2cu·ψ_P
        σ_R + 2cu·ψ_R = σ_Q + 2cu·ψ_Q
      ⇒ ψ_R = (σ_Q - σ_P)/(4cu) + (ψ_P + ψ_Q)/2
        σ_R = (σ_P + σ_Q)/2 + cu·(ψ_Q - ψ_P)

    Geometry: intersection of
        α-line through P with slope = tan(ψ_α - π/4)   (using mean ψ along α)
        β-line through Q with slope = tan(ψ_β + π/4)
    A one-step predictor-corrector is performed.
    """
    # --- Stress (exact) ---
    a = P.sigma - 2.0 * cu * P.psi
    b = Q.sigma + 2.0 * cu * Q.psi
    psi_R = (b - a) / (4.0 * cu)
    sigma_R = 0.5 * (a + b)

    # --- Geometry: predictor then corrector ---
    x_R, y_R = _intersect(P, Q, P.psi, Q.psi)
    psi_a = 0.5 * (P.psi + psi_R)
    psi_b = 0.5 * (Q.psi + psi_R)
    x_R, y_R = _intersect(P, Q, psi_a, psi_b)

    return Node(x=x_R, y=y_R, sigma=sigma_R, psi=psi_R)


def _intersect(P: Node, Q: Node, psi_a: float, psi_b: float) -> tuple[float, float]:
    """Intersection of α-line through P (slope tan(ψ_a-π/4)) and
    β-line through Q (slope tan(ψ_b+π/4))."""
    m_a = alpha_slope(psi_a)
    m_b = beta_slope(psi_b)
    denom = m_a - m_b
    if abs(denom) < 1e-14:
        return 0.5 * (P.x + Q.x), 0.5 * (P.y + Q.y)
    x = (Q.y - P.y + m_a * P.x - m_b * Q.x) / denom
    y = P.y + m_a * (x - P.x)
    return x, y


# ======================================================================
# Mohr-Coulomb characteristic directions and Riemann step
# ======================================================================

def alpha_slope_mc(psi: float, phi: float) -> float:
    """Slope dy/dx of α-line at ψ for Mohr-Coulomb with friction angle φ."""
    return np.tan(psi - np.pi / 4.0 + phi / 2.0)


def beta_slope_mc(psi: float, phi: float) -> float:
    """Slope dy/dx of β-line at ψ for Mohr-Coulomb with friction angle φ."""
    return np.tan(psi + np.pi / 4.0 - phi / 2.0)


def riemann_step_mc(P: Node, Q: Node, c: float, phi: float) -> Node:
    """
    Riemann finite-difference step for the Kötter system (Mohr-Coulomb).

    Node.sigma stores p̄ = p + c·cot(φ).

      P : α-parent  (p̄·exp(−2·tan φ·ψ) = const along α)
      Q : β-parent  (p̄·exp(+2·tan φ·ψ) = const along β)

    Stress update (exact):
        a = ln(p̄_P) − 2·tan φ·ψ_P
        b = ln(p̄_Q) + 2·tan φ·ψ_Q
        p̄_R = exp((a + b) / 2)
        ψ_R  = (b − a) / (4·tan φ)

    Geometry: predictor-corrector intersection of α through P and β through Q.
    """
    tanphi = np.tan(phi)
    a = np.log(P.sigma) - 2.0 * tanphi * P.psi
    b = np.log(Q.sigma) + 2.0 * tanphi * Q.psi
    pbar_R = np.exp(0.5 * (a + b))
    psi_R  = (b - a) / (4.0 * tanphi)

    # Predictor
    x_R, y_R = _intersect_mc(P, Q, P.psi, Q.psi, phi)
    # Corrector
    psi_a = 0.5 * (P.psi + psi_R)
    psi_b = 0.5 * (Q.psi + psi_R)
    x_R, y_R = _intersect_mc(P, Q, psi_a, psi_b, phi)

    return Node(x=x_R, y=y_R, sigma=pbar_R, psi=psi_R)


def _intersect_mc(
    P: Node, Q: Node, psi_a: float, psi_b: float, phi: float
) -> tuple[float, float]:
    """Intersection of α through P and β through Q (Mohr-Coulomb)."""
    m_a = alpha_slope_mc(psi_a, phi)
    m_b = beta_slope_mc(psi_b, phi)
    denom = m_a - m_b
    if abs(denom) < 1e-14:
        return 0.5 * (P.x + Q.x), 0.5 * (P.y + Q.y)
    x = (Q.y - P.y + m_a * P.x - m_b * Q.x) / denom
    y = P.y + m_a * (x - P.x)
    return x, y


# ======================================================================
# Mohr-Coulomb Riemann step WITH GRAVITY (body force +γ ĵ, y-down)
# ======================================================================

def riemann_step_mc_gravity(
    P: Node, Q: Node, c: float, phi: float, gamma: float, n_iter: int = 5
) -> Node:
    """
    Riemann finite-difference step for the Kötter system with weight.

    Kötter equations (compression positive, y-down, body force +γ ĵ):
        Along α:  d p̄ − 2 p̄ tan φ dψ = γ (dy − tan φ dx)
        Along β:  d p̄ + 2 p̄ tan φ dψ = γ (dy + tan φ dx)
    Dividing by p̄_avg:
        ln p̄_R − 2 tan φ ψ_R = ln p̄_P − 2 tan φ ψ_P
                              + (γ / p̄_avg_α)·[(y_R−y_P) − tan φ (x_R−x_P)]
        ln p̄_R + 2 tan φ ψ_R = ln p̄_Q + 2 tan φ ψ_Q
                              + (γ / p̄_avg_β)·[(y_R−y_Q) + tan φ (x_R−x_Q)]
    Solved by predictor–corrector iteration (geometry + Kötter coupled).
    Reduces to `riemann_step_mc` when γ = 0.
    """
    tanphi = np.tan(phi)

    # Weightless initial guess
    a0 = np.log(P.sigma) - 2.0 * tanphi * P.psi
    b0 = np.log(Q.sigma) + 2.0 * tanphi * Q.psi
    pbar_R = np.exp(0.5 * (a0 + b0))
    psi_R  = (b0 - a0) / (4.0 * tanphi)
    x_R, y_R = _intersect_mc(P, Q, P.psi, Q.psi, phi)

    for _ in range(int(n_iter)):
        pa = 0.5 * (P.sigma + pbar_R)
        pb = 0.5 * (Q.sigma + pbar_R)
        rhs_a = (gamma / pa) * ((y_R - P.y) - tanphi * (x_R - P.x))
        rhs_b = (gamma / pb) * ((y_R - Q.y) + tanphi * (x_R - Q.x))
        a = np.log(P.sigma) - 2.0 * tanphi * P.psi + rhs_a
        b = np.log(Q.sigma) + 2.0 * tanphi * Q.psi + rhs_b
        ln_pbar_R = 0.5 * (a + b)
        psi_R     = (b - a) / (4.0 * tanphi)
        pbar_R    = np.exp(ln_pbar_R)

        psi_a = 0.5 * (P.psi + psi_R)
        psi_b = 0.5 * (Q.psi + psi_R)
        x_R, y_R = _intersect_mc(P, Q, psi_a, psi_b, phi)

    return Node(x=x_R, y=y_R, sigma=pbar_R, psi=psi_R)
