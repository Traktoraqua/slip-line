"""
Boundary condition dataclasses for slip line analysis.
"""

from dataclasses import dataclass
import numpy as np


@dataclass
class FoundationBoundary:
    """
    Boundary A: foundation base.

    Attributes
    ----------
    cu : float
        Undrained shear strength.
    r : float
        Roughness ratio = τ / cu, in [0, 1].
        r=0  → smooth foundation (no shear traction)
        r=1  → fully rough foundation (shear = cu)
    B : float
        Foundation half-width [length units]; used for geometry scaling only.
        Stresses are independent of B in a weightless material.
    """
    cu: float
    r: float
    B: float = 1.0

    def __post_init__(self):
        if self.cu <= 0:
            raise ValueError("cu must be positive.")
        if not 0.0 <= self.r <= 1.0:
            raise ValueError("Roughness ratio r must be in [0, 1].")
        if self.B <= 0:
            raise ValueError("Foundation half-width B must be positive.")


@dataclass
class FreeSurfaceBoundary:
    """
    Boundary B: stress-free surface adjacent to the foundation.

    σ_yy = 0, τ_xy = 0 everywhere on this boundary.
    No parameters needed — state is fully determined by the Tresca yield criterion.
    """
    pass


@dataclass
class DiskFoundationBoundary:
    """
    Boundary A (disk variant): a rigid circular disk of diameter D partly
    penetrated into a flat soil surface to depth z = p·D.

    Geometry (y-axis pointing DOWN, free surface at y = 0):
      Disk radius:    R = D / 2
      Disk centre:    (0, y_c) with y_c = (p − 1/2)·D
                      (for p < 1/2 the centre is above the surface)
      Lowest point:   (0, p·D)
      Half-chord at the surface (disk-soil contact ends):
                      B_s = D · √(p·(1 − p))
      Contact half-angle (from the downward axis at the disk centre):
                      θ_max = arccos(1 − 2p)

    Attributes
    ----------
    cu : float
        Undrained shear strength.
    D  : float
        Disk diameter.
    p  : float
        Penetration ratio, p = z / D, in [0, 0.5].
        p = 0    →  no penetration (disk just touches the surface)
        p = 0.5  →  half disk (equator at the surface)
    r  : float
        Roughness ratio = |τ| / cu at the disk surface, in [0, 1].
    """
    cu: float
    D: float
    p: float
    r: float = 0.0

    def __post_init__(self):
        if self.cu <= 0:
            raise ValueError("cu must be positive.")
        if self.D <= 0:
            raise ValueError("Disk diameter D must be positive.")
        if not 0.0 <= self.p <= 0.5:
            raise ValueError("Penetration ratio p must be in [0, 0.5].")
        if not 0.0 <= self.r <= 1.0:
            raise ValueError("Roughness ratio r must be in [0, 1].")

    # ------------------------------------------------------------------
    # Derived geometric quantities
    # ------------------------------------------------------------------

    @property
    def R(self) -> float:
        """Disk radius."""
        return 0.5 * self.D

    @property
    def y_center(self) -> float:
        """y-coordinate of disk centre (y down, surface at y=0)."""
        return (self.p - 0.5) * self.D

    @property
    def B_s(self) -> float:
        """Half-width of the soil-disk contact at the free surface."""
        return self.D * np.sqrt(self.p * (1.0 - self.p))

    @property
    def theta_max(self) -> float:
        """Contact half-angle, measured from the downward axis at the centre."""
        return np.arccos(1.0 - 2.0 * self.p)

    def position(self, phi: float) -> tuple[float, float]:
        """(x, y) on the disk surface at angle φ from the downward axis."""
        return self.R * np.sin(phi), self.y_center + self.R * np.cos(phi)
