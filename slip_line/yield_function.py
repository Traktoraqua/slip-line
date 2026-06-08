"""
Yield function classes for slip line analysis.

Sign convention (geotechnical, COMPRESSION POSITIVE):
  σ   : mean stress  (= (σ_1 + σ_3) / 2)
  ψ   : angle (CCW from +x axis) of the major principal stress σ_1
  y-axis points DOWNWARD (depth).  CCW in (x, y_down) is visually clockwise.

Mohr circle (Tresca, radius = cu):
    σ_xx = σ + cu·cos(2ψ)
    σ_yy = σ - cu·cos(2ψ)
    τ_xy = cu·sin(2ψ)

Slip line directions:
    α-line slope (dy/dx): tan(ψ - π/4)
    β-line slope (dy/dx): tan(ψ + π/4)

Hencky equations (for the convention chosen here):
    Along α-line:  σ - 2·cu·ψ = const
    Along β-line:  σ + 2·cu·ψ = const
"""

import numpy as np


class TrescaYield:
    """Tresca yield criterion: τ_max = cu."""

    def __init__(self, cu: float):
        if cu <= 0:
            raise ValueError("cu must be positive.")
        self.cu = cu

    # ------------------------------------------------------------------
    # Boundary stress states
    # ------------------------------------------------------------------

    def stress_state_free_surface(self) -> tuple[float, float]:
        """
        Boundary B: stress-free surface (σ_yy = τ_xy = 0).

        Passive state — major principal stress horizontal:
            ψ_B = 0      (from τ_xy = cu·sin(2ψ) = 0)
            σ_B = cu     (from σ_yy = σ - cu·cos(2ψ) = 0 with ψ=0)

        Returns
        -------
        (sigma_B, psi_B) = (cu, 0.0)
        """
        return self.cu, 0.0

    def stress_state_foundation(self, r: float) -> float:
        """
        Boundary A: foundation base.  Shear stress τ_xy = r·cu, r ∈ [0, 1].

        Active state — σ_1 close to vertical:
            ψ_A = π/2 − arcsin(r)/2   (from sin(2ψ_A) = r, active branch)

        Notes
        -----
        For r = 0 (smooth): ψ_A = π/2, σ_1 exactly vertical → exact Prandtl
        result q = (π + 2)·cu.

        For r > 0, this implementation uses a *uniform-stress* Region I,
        which is only fully consistent for the smooth case.  In reality the
        fully-rough strip-footing problem also yields q = (π+2)·cu (Hill
        mechanism with the foundation surface itself being a slip line);
        modelling that requires a non-uniform Region I and is not yet
        implemented here.

        Returns
        -------
        psi_A : float [radians]
        """
        if not 0.0 <= r <= 1.0:
            raise ValueError("Roughness ratio r must be in [0, 1].")
        return np.pi / 2.0 - np.arcsin(r) / 2.0

    # ------------------------------------------------------------------
    # Mohr-circle conversions
    # ------------------------------------------------------------------

    def sigma_yy(self, sigma: float, psi: float) -> float:
        return sigma - self.cu * np.cos(2.0 * psi)

    def sigma_xx(self, sigma: float, psi: float) -> float:
        return sigma + self.cu * np.cos(2.0 * psi)

    def tau_xy(self, sigma: float, psi: float) -> float:
        return self.cu * np.sin(2.0 * psi)


class MohrCoulombYield:
    """
    Mohr-Coulomb yield criterion: τ = c + σ_n·tan(φ).

    Internal stress variable: p̄ = p + c·cot(φ)  (modified mean stress).

    Kötter equations (weightless, along characteristics):
        Along α-line:  p̄ · exp(−2ψ·tan φ) = const
        Along β-line:  p̄ · exp(+2ψ·tan φ) = const

    Slip-line directions (y-down convention, ψ = angle of σ₁ from x-axis):
        α-slope (dy/dx) = tan(ψ − π/4 + φ/2)
        β-slope (dy/dx) = tan(ψ + π/4 − φ/2)

    The two families make angle (π/2 − φ) with each other.

    Stress conversions from p̄:
        R   = p̄ · sin φ
        p   = p̄ − c·cot φ
        σ_xx = p + R·cos(2ψ)
        σ_yy = p − R·cos(2ψ)
        τ_xy = R·sin(2ψ)
    """

    def __init__(self, c: float, phi: float):
        if c < 0.0:
            raise ValueError("Cohesion c must be non-negative.")
        if not (0.0 < phi < np.pi / 2.0):
            raise ValueError("Friction angle phi must be in (0, π/2) radians.")
        self.c = c
        self.phi = phi

    # ------------------------------------------------------------------
    # Boundary stress states  (returns p̄, ψ pairs)
    # ------------------------------------------------------------------

    def pbar_free_surface(self) -> float:
        """
        p̄ at a stress-free passive surface (σ_yy = τ_xy = 0, ψ = 0).

            p_B = c·cos φ / (1 − sin φ)
            p̄_B = p_B + c·cot φ = c·cos φ / (sin φ·(1 − sin φ))
        """
        sp = np.sin(self.phi)
        cp = np.cos(self.phi)
        return self.c * cp / (sp * (1.0 - sp))

    def stress_state_free_surface(self) -> tuple[float, float]:
        """Returns (p̄_B, ψ_B) = (pbar_free_surface(), 0.0)."""
        return self.pbar_free_surface(), 0.0

    def stress_state_foundation(self, r: float = 0.0) -> float:
        """
        ψ_A at a smooth (r = 0) foundation: τ_xy = 0 → ψ_A = π/2.

        For a rough interface (r > 0) the boundary condition couples
        ψ_A to φ and r in a non-trivial way; not yet implemented.

        Returns
        -------
        psi_A : float [radians]
        """
        if r != 0.0:
            raise NotImplementedError(
                "Rough Mohr-Coulomb foundation boundary (r > 0) not yet implemented."
            )
        return np.pi / 2.0

    # ------------------------------------------------------------------
    # Bearing-capacity factors (smooth strip footing, weightless soil)
    # ------------------------------------------------------------------

    def Nq(self) -> float:
        """Nq = tan²(π/4 + φ/2) · exp(π·tan φ)."""
        return np.tan(np.pi / 4.0 + self.phi / 2.0) ** 2 * np.exp(np.pi * np.tan(self.phi))

    def Nc(self) -> float:
        """Nc = (Nq − 1) / tan φ.  (Prandtl formula, weightless, no surcharge.)"""
        return (self.Nq() - 1.0) / np.tan(self.phi)

    # ------------------------------------------------------------------
    # Mohr-circle conversions from p̄
    # ------------------------------------------------------------------

    def sigma_yy(self, pbar: float, psi: float) -> float:
        """σ_yy = p − R·cos(2ψ),  R = p̄·sin φ,  p = p̄ − c·cot φ."""
        sp = np.sin(self.phi)
        R = pbar * sp
        p = pbar - self.c * np.cos(self.phi) / sp
        return p - R * np.cos(2.0 * psi)

    def sigma_xx(self, pbar: float, psi: float) -> float:
        """σ_xx = p + R·cos(2ψ)."""
        sp = np.sin(self.phi)
        R = pbar * sp
        p = pbar - self.c * np.cos(self.phi) / sp
        return p + R * np.cos(2.0 * psi)

    def tau_xy(self, pbar: float, psi: float) -> float:
        """τ_xy = R·sin(2ψ),  R = p̄·sin φ."""
        return pbar * np.sin(self.phi) * np.sin(2.0 * psi)
