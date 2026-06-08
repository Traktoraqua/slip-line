"""
PrandtlSolver — slip line field for a strip footing on a weightless Tresca
half-space.  Right half (x ≥ 0) computed; symmetric mirror used for plotting.

Geometry (y axis pointing DOWN, foundation at y = 0, half-width B):

        |←———— B ————→|
        ===============  foundation (y=0, 0 ≤ x ≤ B)
       /   I   |  II  | III \\          free surface (y=0, x > B)
      /        |/fan\\|      \\
     /  active | (Δψ)| passive\\
    /          *—————*         \\
                ^corner (B,0)

For Tresca with γ = 0, the stress field is analytical:

  Region III (passive wedge, x > B):
      ψ = 0,  σ = cu                              (uniform)

  Region II (centred fan at the corner):
      ψ varies from 0 (on III side) to ψ_A (on I side)
      Along each radial (β-line), ψ = const and σ = const = cu·(1 + 2·ψ)

  Region I (active wedge under foundation):
      ψ = ψ_A,  σ = cu·(1 + 2·ψ_A)                (uniform)

Collapse load:
      q = σ_yy at foundation = σ_I − cu·cos(2·ψ_A)
        = cu·(1 + 2·ψ_A) − cu·cos(2·ψ_A)
      For ψ_A = π/2 (smooth):  q = cu·(2 + π)    ✓ Prandtl
"""

from __future__ import annotations
import numpy as np

from .yield_function import TrescaYield, MohrCoulombYield
from .boundary import FoundationBoundary, DiskFoundationBoundary
from .slip_line_field import (
    Node, alpha_slope, beta_slope, riemann_step,
    alpha_slope_mc, beta_slope_mc, riemann_step_mc, riemann_step_mc_gravity,
)


class PrandtlSolver:
    """
    Build the slip line field and compute the collapse load.

    Parameters
    ----------
    foundation : FoundationBoundary
    n_fan      : number of equal-angle increments in the Prandtl fan
    n_radial   : number of radial mesh steps in each region (default = n_fan)
    """

    def __init__(
        self,
        foundation: FoundationBoundary,
        n_fan: int = 20,
        n_radial: int | None = None,
    ):
        self.fdn = foundation
        self.cu = foundation.cu
        self.r = foundation.r
        self.B = foundation.B
        self.n_fan = int(n_fan)
        self.n_radial = int(n_radial) if n_radial is not None else self.n_fan

        self.yield_fn = TrescaYield(self.cu)
        self.sigma_B, self.psi_B = self.yield_fn.stress_state_free_surface()
        self.psi_A = self.yield_fn.stress_state_foundation(self.r)
        self.delta_psi = self.psi_A - self.psi_B  # fan rotation (positive)

        # Mesh radius step: choose so the passive wedge has horizontal extent B
        # along the free surface (i.e. apex at (B + B/2, B/2)).
        # Mesh sizing — chosen so the lower boundary (outermost α-line) is
        # CONTINUOUS through all three zones for any r ∈ [0, 1].
        #
        # Region I uses surface anchors at x = i·h_I, i = 0 … n_radial, with
        # h_I = dr / sin θ_α  (so that each α-line emanating from a foundation
        # anchor lands on a fan-I node spaced dr apart along fan_I).
        # Requiring N_I = n_radial AND x_{N_I} = B (corner) gives
        #     h_I = B / n_radial,
        #     dr  = h_I · sin θ_α = (B / n_radial) · sin θ_α.
        #
        # Region III's free-surface spacing h must give a corner-side β-line
        # that matches the k=0 fan radial (angle π/4, length n_radial·dr):
        #     h = dr · √2          (so n_radial·h/√2 = n_radial·dr matches fan).
        # The free-surface extent is L = n_radial · h.
        sin_alpha = np.sin(self.psi_A - np.pi / 4.0)
        self.dr = (self.B / self.n_radial) * sin_alpha
        self.h = self.dr * np.sqrt(2.0)
        self.L = self.n_radial * self.h

        # Storage
        self.region_III_nodes: list[list[Node]] = []   # [j][k] j=depth row, k=node along row
        self.region_II_nodes:  list[list[Node]] = []   # [k][j] k=fan radial index, j=radial step
        self.region_I_nodes:   list[list[Node]] = []   # [k][m] same layout as fan, m=step along α toward foundation

        self._built = False

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def solve(self) -> float:
        self._build_region_III()
        self._build_region_II()
        self._build_region_I()
        self._built = True
        return self.collapse_load()

    def collapse_load(self) -> float:
        """q = σ_yy at the foundation."""
        sigma_I = self.cu * (1.0 + 2.0 * self.psi_A)
        return sigma_I - self.cu * np.cos(2.0 * self.psi_A)

    def prandtl_analytical(self) -> float:
        """Reference value from the same closed-form expression (identical here
        because the field is analytical for weightless Tresca)."""
        return self.collapse_load()

    # ------------------------------------------------------------------
    # Region III — passive wedge (uniform stress)
    # ------------------------------------------------------------------

    def _build_region_III(self):
        """
        Mesh:
          j=0: free-surface nodes (B + i·h, 0), i = 0 … n_radial.
          j>0: Riemann grid using α-line slope -1 and β-line slope +1.
               Each new node uses P=(j-1, i+1) [α-parent] and Q=(j-1, i) [β-parent].
               Result: node (j, i) at (B + (i + j/2)·h, j·h/2)  for i = 0 … n_radial-j.

        Stress is uniform (ψ=0, σ=cu); we use the analytical positions directly.
        """
        n = self.n_radial
        h = self.h
        rows: list[list[Node]] = []
        for j in range(n + 1):
            row: list[Node] = []
            for i in range(n + 1 - j):
                x = self.B + (i + 0.5 * j) * h
                y = 0.5 * j * h
                row.append(Node(x=x, y=y, sigma=self.cu, psi=0.0))
            rows.append(row)
        self.region_III_nodes = rows

        # The k=0 column (i.e. row[0] in each j-row) is the β-line from the
        # corner; it is shared with the fan (Region II) at ψ = 0.

    # ------------------------------------------------------------------
    # Region II — centred fan (analytical)
    # ------------------------------------------------------------------

    def _build_region_II(self):
        """
        Centred fan at corner (B, 0).  Radials at angles
            θ_k = π/4 + ψ_k             (slope = tan(θ_k))
        where ψ_k = ψ_B + k·(Δψ / n_fan),  k = 0 … n_fan.

        Radial k=0 coincides with the corner-side boundary of Region III.
        Radial k=n_fan coincides with the corner-side boundary of Region I.

        Node (k, j) at radius r_j = j·dr from corner, along radial k:
            x = B + r_j·cos(θ_k),  y = r_j·sin(θ_k)
            ψ = ψ_k
            σ = cu·(1 + 2·ψ_k)         (constant along each radial)
        """
        n_p = self.n_radial
        cu = self.cu
        dpsi = self.delta_psi / self.n_fan

        radials: list[list[Node]] = []
        for k in range(self.n_fan + 1):
            psi_k = self.psi_B + k * dpsi
            theta_k = np.pi / 4.0 + psi_k          # β-line slope angle
            sigma_k = cu * (1.0 + 2.0 * psi_k)
            line: list[Node] = []
            for j in range(n_p + 1):
                r_j = j * self.dr
                x = self.B + r_j * np.cos(theta_k)
                y = r_j * np.sin(theta_k)
                line.append(Node(x=x, y=y, sigma=sigma_k, psi=psi_k))
            radials.append(line)
        self.region_II_nodes = radials

    # ------------------------------------------------------------------
    # Region I — active wedge under foundation (uniform stress)
    # ------------------------------------------------------------------

    def _build_region_I(self):
        """
        Uniform-stress region with ψ = ψ_A, σ = cu·(1 + 2·ψ_A).

        Region I is the triangle bounded by
          • the foundation surface  y = 0,  x ∈ [0, B]
          • the symmetry axis       x = 0
          • the fan-I β-line        from the corner (B, 0) along
                                    direction (cos θ_β, sin θ_β)

        Because the stress state is uniform, the α-lines (slope angle
        θ_α = ψ_A − π/4) and β-lines (slope angle θ_β = ψ_A + π/4) are
        each a family of straight, parallel lines.

        Mesh construction (parallelogram grid):
          Foundation anchor points  P_i = (i·h_I, 0),  i = 0 … N_I,
          where h_I = dr / sin θ_α  (so the corner-side α-line spacing
          matches the fan radial step dr along the fan-I β-line).
          Node (i, j) with i ≥ j is the intersection of the α-line
          through P_j and the β-line through P_i.

        Indexing: region_I_nodes[i][j],  0 ≤ j ≤ i ≤ N_I.
          • α-line through P_j  → column j: nodes (i, j) for i = j … N_I.
          • β-line through P_i  → row i:    nodes (i, j) for j = 0 … i.
          • Row i = N_I coincides with the fan-I β-line in reverse:
            node(N_I, j) == region_II_nodes[-1][N_I − j].
        """
        cu = self.cu
        psi_A = self.psi_A
        sigma_I = cu * (1.0 + 2.0 * psi_A)
        theta_alpha = psi_A - np.pi / 4.0
        theta_beta = psi_A + np.pi / 4.0
        sin_a = np.sin(theta_alpha)
        cos_b = np.cos(theta_beta)
        sin_b = np.sin(theta_beta)

        # Surface spacing chosen so each foundation anchor sends an
        # α-line that lands on a node of the fan-I β-line.
        if sin_a < 1e-12:
            # Degenerate (ψ_A = π/4) — shouldn't occur for r ∈ [0, 1].
            self.region_I_nodes = []
            return

        # By construction (see __init__): h_I = B / n_radial exactly,
        # so N_I = n_radial.
        h_I = self.B / self.n_radial
        N_I = self.n_radial

        # Foundation anchor x-coords:  x_0 = 0 (sym axis), x_{N_I} = corner side.
        x_anchor = [i * h_I for i in range(N_I + 1)]

        region_I: list[list[Node]] = []
        for i in range(N_I + 1):
            row: list[Node] = []
            for j in range(i + 1):
                if i == j:
                    # On the foundation surface
                    x = x_anchor[i]
                    y = 0.0
                else:
                    # Intersection of α-line from P_j and β-line from P_i.
                    # s_β = (x_i − x_j) · sin θ_α  (since sin(θ_β − θ_α)=1)
                    s_beta = (x_anchor[i] - x_anchor[j]) * sin_a
                    x = x_anchor[i] + s_beta * cos_b
                    y = s_beta * sin_b
                row.append(Node(x=x, y=y, sigma=sigma_I, psi=psi_A))
            region_I.append(row)
        self.region_I_nodes = region_I
        self.N_I = N_I

    # ------------------------------------------------------------------
    # Accessors for plotting
    # ------------------------------------------------------------------

    def all_alpha_lines(self) -> list[list[Node]]:
        """Return α-line polylines across all three regions."""
        lines: list[list[Node]] = []

        # ---- Region III: α-lines anchored on the free surface at (B + i·h, 0)
        # going down-LEFT (slope -1).  Node (j, k=k) sits on the α-line whose
        # surface anchor is at i_surface = j + k  (from the parametrisation).
        n = self.n_radial
        for i_surf in range(n + 1):
            line: list[Node] = []
            for j in range(min(i_surf, n) + 1):
                k = i_surf - j
                if k < 0 or k > n - j:
                    break
                line.append(self.region_III_nodes[j][k])
            if len(line) >= 2:
                lines.append(line)

        # ---- Region II: α-lines = arcs across the fan, joining nodes at the
        # same radius r_j across all fan radials.
        for j in range(1, n + 1):
            arc = [self.region_II_nodes[k][j]
                   for k in range(self.n_fan + 1)
                   if j < len(self.region_II_nodes[k])]
            if len(arc) >= 2:
                lines.append(arc)

        # ---- Region I: α-lines = columns (constant j) of region_I_nodes.
        # Each α-line goes from foundation point P_j at (j·h_I, 0) up-right
        # to fan_I.  Column j: nodes (i, j) for i = j … N_I.
        if self.region_I_nodes:
            N_I = len(self.region_I_nodes) - 1
            for j in range(N_I + 1):
                line = [self.region_I_nodes[i][j] for i in range(j, N_I + 1)]
                if len(line) >= 2:
                    lines.append(line)

        return lines

    def all_beta_lines(self) -> list[list[Node]]:
        """Return β-line polylines across all three regions."""
        lines: list[list[Node]] = []
        n = self.n_radial

        # ---- Region III: β-lines anchored on free surface at (B + i·h, 0)
        # going down-RIGHT (slope +1).  Node (j, k) lies on the β-line whose
        # surface anchor index is k (same k across all j).
        for k in range(n + 1):
            line = []
            for j in range(n + 1 - k):
                line.append(self.region_III_nodes[j][k])
            if len(line) >= 2:
                lines.append(line)

        # ---- Region II: β-lines = radials (each fan radial is a β-line).
        # Skip k=0 (already drawn as the corner-side β-line of Region III).
        for k in range(1, self.n_fan + 1):
            lines.append(self.region_II_nodes[k])

        # ---- Region I: β-lines = rows (constant i) of region_I_nodes.
        # Each β-line goes from foundation point P_i at (i·h_I, 0) down-left
        # to the symmetry axis.  Row i: nodes (i, j) for j = 0 … i.
        # Skip i = N_I (row coincides with fan_I, already drawn in Region II).
        if self.region_I_nodes:
            N_I = len(self.region_I_nodes) - 1
            for i in range(N_I):
                line = list(self.region_I_nodes[i])
                if len(line) >= 2:
                    lines.append(line)

        return lines


# ======================================================================
# Mohr-Coulomb flat foundation solver
# ======================================================================


class MohrCoulombPrandtlSolver:
    """
    Prandtl slip-line field for a smooth strip footing on weightless
    Mohr-Coulomb soil (c, φ).

    Geometry (y down, foundation at y = 0, half-width B):

        |←——— B ———→|
        =============  smooth foundation (τ = 0)
       /  I   | II  | III \\
      /       |fan  |      \\
     /        *————-*       \\
                ^corner (B, 0)

    Stress variable stored in Node.sigma: p̄ = p + c·cot(φ).

    Kötter equations:
      Along α:  p̄ · exp(−2ψ·tan φ) = const
      Along β:  p̄ · exp(+2ψ·tan φ) = const

    Characteristic directions (y-down):
      α-slope = tan(ψ − π/4 + φ/2)
      β-slope = tan(ψ + π/4 − φ/2)

    Boundary conditions:
      Free surface (Region III):  ψ = 0,    p̄ = c·cos φ / (sin φ·(1−sin φ))
      Foundation (smooth):        ψ = π/2,  p̄ = p̄_B · exp(π·tan φ)

    Analytical collapse load (Prandtl, smooth, weightless, no surcharge):
      q = c · Nc
      Nc = (Nq − 1) / tan φ
      Nq = tan²(π/4 + φ/2) · exp(π·tan φ)
    """

    def __init__(
        self,
        c: float,
        phi: float,
        B: float,
        n_fan: int = 20,
        n_radial: int | None = None,
    ):
        """
        Parameters
        ----------
        c       : cohesion
        phi     : friction angle [radians]
        B       : foundation half-width
        n_fan   : fan discretisation steps
        n_radial: mesh radial steps (default = n_fan)
        """
        self.c = c
        self.phi = phi
        self.B = B
        self.n_fan = int(n_fan)
        self.n_radial = int(n_radial) if n_radial is not None else self.n_fan

        self.yield_fn = MohrCoulombYield(c, phi)

        # Boundary stress states
        self.pbar_B, self.psi_B = self.yield_fn.stress_state_free_surface()
        self.psi_A = self.yield_fn.stress_state_foundation(r=0.0)   # = π/2
        self.delta_psi = self.psi_A - self.psi_B                    # = π/2
        self.pbar_A = self.pbar_B * np.exp(np.pi * np.tan(phi))

        # Slip-line angles in each region
        # Region I  (ψ_A = π/2):  α at π/4+φ/2,  β at 3π/4−φ/2
        # Region III (ψ_B = 0):   α at φ/2−π/4,  β at π/4−φ/2
        self.theta_alpha_I   =  np.pi / 4.0 + phi / 2.0   # α angle in Region I
        self.theta_beta_III  =  np.pi / 4.0 - phi / 2.0   # β angle in Region III

        # Mesh sizing (ensures continuous lower boundary across all zones)
        #   dr   : fan radial step (length of one β-line increment at Zone I boundary)
        #   h    : free-surface node spacing in Region III
        #   L    : total free-surface extent of Region III
        #
        # Derivation:
        #   In Zone I both families have slope ±tan(θ_α).  Intersection of the
        #   α-line from anchor x_j and β-line from anchor x_i = x_j + h_I:
        #       x = (x_i + x_j)/2,  y = tan(θ_α)·h_I/2
        #   The β-line runs from (x_i, 0) to the intersection, length:
        #       Δx = −h_I/2,  Δy = tan(θ_α)·h_I/2
        #       r  = h_I/2 · √(1 + tan²θ_α) = h_I / (2·cos θ_α)
        #   so  dr = (B/n_radial) / (2·cos(π/4 + φ/2))
        #
        #   The log-spiral growth of α-lines across the fan means the same α-line
        #   is at r·exp(+Δψ·tan φ) on the Zone-III side, giving:
        #       h = 2·dr·exp(+Δψ·tan φ)·cos(π/4 − φ/2)
        self.dr = B / (2.0 * self.n_radial * np.cos(self.theta_alpha_I))
        self.h  = 2.0 * self.dr * np.exp(+self.delta_psi * np.tan(phi)) * np.cos(self.theta_beta_III)
        self.L  = self.n_radial * self.h

        # Storage
        self.region_III_nodes: list[list[Node]] = []
        self.region_II_nodes:  list[list[Node]] = []
        self.region_I_nodes:   list[list[Node]] = []

        self._built = False

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def solve(self) -> float:
        self._build_region_III()
        self._build_region_II()
        self._build_region_I()
        self._built = True
        return self.collapse_load()

    def collapse_load(self) -> float:
        """q = σ_yy at foundation = c·Nc (analytical)."""
        return self.yield_fn.sigma_yy(self.pbar_A, self.psi_A)

    def analytical_Nc(self) -> float:
        return self.yield_fn.Nc()

    # ------------------------------------------------------------------
    # Region III — passive wedge (uniform ψ=0, p̄=p̄_B)
    # ------------------------------------------------------------------

    def _build_region_III(self):
        """
        Uniform-stress region: ψ = 0, p̄ = p̄_B.

        Characteristic directions:
          α-slope = tan(−π/4+φ/2) = −tan(π/4−φ/2)   [going right-upward]
          β-slope = tan(+π/4−φ/2) = +tan(π/4−φ/2)   [going right-downward]

        Node (j, i): intersection of
          β-line through free-surface point (B + i·h, 0)
          α-line through free-surface point (B + (i+j)·h, 0)

        Analytical position:
          x = B + (i + j/2)·h          [same formula as Tresca]
          y = j·h/2 · tan(π/4 − φ/2)   [Tresca: multiply by 1]
        """
        n = self.n_radial
        h = self.h
        m_beta = np.tan(self.theta_beta_III)
        rows: list[list[Node]] = []
        for j in range(n + 1):
            row: list[Node] = []
            for i in range(n + 1 - j):
                x = self.B + (i + 0.5 * j) * h
                y = 0.5 * j * h * m_beta
                row.append(Node(x=x, y=y, sigma=self.pbar_B, psi=0.0))
            rows.append(row)
        self.region_III_nodes = rows

    # ------------------------------------------------------------------
    # Region II — centred fan (analytical, logarithmic-spiral α-lines)
    # ------------------------------------------------------------------

    def _build_region_II(self):
        """
        Centred fan at corner (B, 0).

        β-lines: straight radials at angle θ_k = ψ_k + π/4 − φ/2.
        α-lines: logarithmic spirals (approximated by connecting equal-radius
                 nodes across radials for plotting).

        Node (k, j):
          ψ_k   = ψ_B + k·(Δψ/n_fan)
          p̄_k   = p̄_B · exp(2·tan φ·ψ_k)   [from α-invariant at Region III edge]
          θ_k   = ψ_k + π/4 − φ/2
          (x, y) = B + r_{k,j}·(cos θ_k, sin θ_k)
          r_{k,j} = j·dr·exp(−(Δψ − ψ_k)·tan φ)    [log-spiral distance]

        With this formula:
          k = n_fan (ψ = π/2):  r = j·dr              (matches Region I boundary)
          k = 0     (ψ = 0  ):  r = j·dr·exp(−Δψ·tan φ)  (matches Region III)
        α-lines (connecting nodes of equal j across radials) are exact log spirals.
        """
        n_p = self.n_radial
        tanphi = np.tan(self.phi)
        dpsi = self.delta_psi / self.n_fan

        radials: list[list[Node]] = []
        for k in range(self.n_fan + 1):
            psi_k  = self.psi_B + k * dpsi
            theta_k = psi_k + np.pi / 4.0 - self.phi / 2.0
            pbar_k  = self.pbar_B * np.exp(2.0 * tanphi * psi_k)
            # log-spiral scale: α-line at level j is at distance j·dr on the
            # last radial and r·exp(+(Δψ−ψ_k)·tanφ) on the k-th radial.
            spiral_scale = np.exp(+(self.delta_psi - psi_k) * tanphi)
            line: list[Node] = []
            for j in range(n_p + 1):
                r_j = j * self.dr * spiral_scale
                line.append(Node(
                    x=self.B + r_j * np.cos(theta_k),
                    y=r_j * np.sin(theta_k),
                    sigma=pbar_k,
                    psi=psi_k,
                ))
            radials.append(line)
        self.region_II_nodes = radials

    # ------------------------------------------------------------------
    # Region I — active wedge (uniform ψ=π/2, p̄=p̄_A)
    # ------------------------------------------------------------------

    def _build_region_I(self):
        """
        Uniform-stress region: ψ = π/2, p̄ = p̄_A.

        Characteristic directions:
          α-slope = tan(π/4 + φ/2)   [going right-downward steeply]
          β-slope = tan(3π/4 − φ/2)  [going left-downward steeply]

        Exact intersection of two straight characteristics:
          α-line from anchor j: y =  tan_α·(x − x_j)
          β-line from anchor i: y = −tan_α·(x − x_i)   [slope = −tan_α]
          → x = (x_i + x_j)/2,  y = tan_α·(x_i − x_j)/2

        Node (i, j) — i ≥ j:
          Foundation anchor P_i = (i·h_I, 0),  h_I = B/n_radial.
          x = (x_i + x_j) / 2
          y = tan(π/4+φ/2)·(x_i − x_j) / 2
        """
        n = self.n_radial
        h_I = self.B / n
        pbar_A = self.pbar_A
        tan_alpha = np.tan(self.theta_alpha_I)        # tan(π/4 + φ/2)

        x_anchor = [i * h_I for i in range(n + 1)]

        region_I: list[list[Node]] = []
        for i in range(n + 1):
            row: list[Node] = []
            for j in range(i + 1):
                x = 0.5 * (x_anchor[i] + x_anchor[j])
                y = 0.5 * tan_alpha * (x_anchor[i] - x_anchor[j])
                row.append(Node(x=x, y=y, sigma=pbar_A, psi=np.pi / 2.0))
            region_I.append(row)
        self.region_I_nodes = region_I
        self.N_I = n

    # ------------------------------------------------------------------
    # Accessors for plotting (identical structure to PrandtlSolver)
    # ------------------------------------------------------------------

    def all_alpha_lines(self) -> list[list[Node]]:
        """
        Return α-line polylines that are continuous across all three zones.

        Each full α-line (m = 1 … n_radial) is a single concatenated polyline:
          Region I column (foundation → fan boundary)
          + fan arc reversed (fan Zone-I side → fan Zone-III side)
          + Region III strip  (fan Zone-III side → free surface)

        At both junctions (I/II and II/III) the segments share exactly one
        endpoint (duplicates are removed), so the polyline is gap-free.
        Because the log-spiral tangent at the fan boundary exactly matches the
        characteristic slope in the adjacent uniform zone, the derivative is
        also continuous in the theoretical limit (discretisation error ∝ 1/n_fan).
        """
        lines: list[list[Node]] = []
        n       = self.n_radial
        n_fan   = self.n_fan

        for m in range(1, n + 1):
            j_I = n - m   # Region I column index for this α-line

            # 1. Region I: foundation anchor → fan boundary (Zone-I side of fan)
            seg_I = [self.region_I_nodes[i][j_I] for i in range(j_I, n + 1)]

            # 2. Fan arc reversed: Zone-I side (k=n_fan-1) → Zone-III side (k=0)
            #    Skip k=n_fan: already the last node of seg_I (same coordinates).
            seg_fan = [self.region_II_nodes[k][m]
                       for k in range(n_fan - 1, -1, -1)]

            # 3. Region III: fan Zone-III side → free surface
            #    i_surf = m; nodes are region_III_nodes[j][m-j] for j = 0..m.
            #    Reversed (fan side = j=m, surface = j=0); skip j=m (= last of
            #    seg_fan, same coordinates).
            seg_III = [self.region_III_nodes[j][m - j]
                       for j in range(m - 1, -1, -1)]

            line = seg_I + seg_fan + seg_III
            if len(line) >= 2:
                lines.append(line)

        return lines

    def all_beta_lines(self) -> list[list[Node]]:
        lines: list[list[Node]] = []
        n = self.n_radial

        # Region III (columns, constant i)
        for k in range(n + 1):
            line = [self.region_III_nodes[j][k] for j in range(n + 1 - k)]
            if len(line) >= 2:
                lines.append(line)

        # Region II (fan radials, skip k=0 — drawn in Region III)
        for k in range(1, self.n_fan + 1):
            lines.append(self.region_II_nodes[k])

        # Region I (rows, constant i; skip i=N — fan-I β-line already drawn)
        if self.region_I_nodes:
            N = len(self.region_I_nodes) - 1
            for i in range(N):
                line = list(self.region_I_nodes[i])
                if len(line) >= 2:
                    lines.append(line)

        return lines


# ======================================================================
# Mohr-Coulomb flat strip footing WITH SUBMERGED SOIL WEIGHT (gravity)
# ======================================================================


class MohrCoulombGravityPrandtlSolver:
    """
    Slip-line field for a smooth strip footing on weighted Mohr-Coulomb soil.

    Material:  c (cohesion), φ (friction angle), γ (submerged unit weight).
    Geometry:  half-width B, free surface at y = 0, foundation at y = 0,
               x ∈ [−B, B]; gravity body-force +γ ĵ (y-down).

    Strategy (this initial version):
      • Mesh topology is taken from the weightless Prandtl solution
        (`MohrCoulombPrandtlSolver`).  Region geometry is therefore
        identical to the γ = 0 case.
      • Stresses (p̄, ψ) are recomputed WITH gravity, by integrating the
        Kötter equations along characteristics:
          – Region III:  ψ ≡ 0, p̄(y) = p̄_B + γ y / (1 − sin φ).
          – Region II :  along each radial (β-line at ψ_k = const),
                         p̄(s) = p̄_B · exp(2 ψ_k tan φ)
                                + γ · s · [sin θ_k + tan φ · cos θ_k]
                         where θ_k = ψ_k + π/4 − φ/2, s = distance from corner.
          – Region I  :  ψ ≡ π/2 (smooth footing); along each α-line
                         from the inner-fan radial back to the foundation
                         surface, p̄ is updated by the α-Kötter equation.
                         (Straight-α assumption — exact only if ψ stays
                         constant, which holds in the limit γ → 0.)
      • The contact-pressure σ_yy(x) and the collapse load Q are evaluated
        from the foundation-row p̄ values via σ_yy = p̄(1 − sin φ) − c·cot φ.

    Limitations:
      • Slip-line geometry is the weightless topology — strictly valid only
        for moderate γ.  α-line consistency in Region II is not enforced;
        the field is therefore an approximation that becomes exact as γ → 0.
      • For larger γ a fully iterative method-of-characteristics scheme
        (`riemann_step_mc_gravity` used everywhere) is required.
    """

    def __init__(
        self,
        c: float,
        phi: float,
        gamma: float,
        B: float,
        n_fan: int = 20,
        n_radial: int | None = None,
    ):
        self.c = float(c)
        self.phi = float(phi)
        self.gamma = float(gamma)
        self.B = float(B)
        self.n_fan = int(n_fan)
        self.n_radial = int(n_radial) if n_radial is not None else self.n_fan

        self.yield_fn = MohrCoulombYield(self.c, self.phi)
        self.pbar_B, self.psi_B = self.yield_fn.stress_state_free_surface()
        self.psi_A = self.yield_fn.stress_state_foundation(r=0.0)   # = π/2
        self.delta_psi = self.psi_A - self.psi_B                    # = π/2
        # p̄ on the corner-fan side of the smooth foundation (γ = 0 limit at corner)
        self.pbar_A0 = self.pbar_B * np.exp(np.pi * np.tan(self.phi))

        # Mesh angles (same as weightless solver)
        self.theta_alpha_I   = np.pi / 4.0 + self.phi / 2.0
        self.theta_beta_III  = np.pi / 4.0 - self.phi / 2.0

        # Mesh sizing — chosen so the foundation chain (leftmost end of
        # Region I) exactly reaches the symmetry plane x = 0.  For γ = 0
        # this gives  dr_inner = B/n_p  (no iteration needed); for γ > 0
        # the field naturally extends further, so dr_inner must shrink.
        # solve() iterates `dr_inner` so that min_x(foundation) ≈ 0.
        self._set_dr_inner(B / self.n_radial)

        self.region_III_nodes: list[list[Node]] = []
        self.region_II_nodes:  list[list[Node]] = []
        self.region_I_nodes:   list[list[Node]] = []

        self._built = False

    def _set_dr_inner(self, dr_inner: float) -> None:
        """Update mesh-sizing attributes that depend on `dr_inner`."""
        self.dr_inner = float(dr_inner)
        self.dr_outer = self.dr_inner * np.exp(self.delta_psi * np.tan(self.phi))
        # Region III mesh size h chosen so its III/II diagonal coincides
        # one-to-one with the m = 0 row of the (k, m) grid.
        self.h  = 2.0 * self.dr_outer * np.cos(self.theta_beta_III)
        self.L  = self.n_radial * self.h
        # Backwards-compatible attribute
        self.dr = self.dr_outer

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def solve(self, fit_foundation: bool = True,
              fit_tol: float = 1e-4, fit_max_iter: int = 12) -> dict:
        """
        Build the slip-line field and return the bearing-capacity result.

        Parameters
        ----------
        fit_foundation : bool
            If True (default), iteratively rescale `dr_inner` so that the
            leftmost foundation node lands at x = 0 (the symmetry plane).
            This makes Region III adjust its size with γ so Region I
            exactly fits the foundation footprint [0, B].
        fit_tol : float
            Convergence tolerance on min(foundation_x) / B.
        fit_max_iter : int
            Maximum number of rescaling iterations.
        """
        if fit_foundation:
            for _ in range(int(fit_max_iter)):
                self._build_region_III()
                self._build_region_II()
                self._build_region_I()
                x_min = min(n.x for n in self.foundation_nodes)
                L_chain = self.B - x_min
                if abs(x_min) / self.B < fit_tol:
                    break
                # Scale dr so the new chain length equals B
                self._set_dr_inner(self.dr_inner * self.B / L_chain)
        else:
            self._build_region_III()
            self._build_region_II()
            self._build_region_I()
        self._built = True

        q_x, q_vals = self.contact_pressure_curve()
        Q = float(np.trapz(q_vals, q_x))   # half-footing load per unit length
        return {
            'q_per_unit_length_half': Q,
            'q_per_unit_length_full': 2.0 * Q,
            'q_over_2B': Q / self.B,                    # average pressure on half
            'q_over_cB': Q / (self.c * self.B),
            'q_over_gammaB2': Q / (self.gamma * self.B ** 2) if self.gamma > 0 else None,
            'Nc_weightless': self.yield_fn.Nc(),
            'sigma_yy_at_centre': q_vals[0],
            'sigma_yy_at_corner': q_vals[-1],
        }

    def contact_pressure_curve(self):
        """Return (x, σ_yy) along the foundation contact (smooth footing).

        ψ = π/2 along the smooth foundation, so cos(2ψ) = −1, hence
            σ_yy = (p̄ − c·cot φ) + p̄·sin φ = p̄·(1 + sin φ) − c·cot φ.

        Only the portion x ∈ [0, B] is returned (the actual half-footing
        footprint).  If the natural field extends to negative x (which
        happens for γ > 0 because the inner-fan-radial is shifted by
        gravity), the curve is linearly interpolated to a node at x = 0.
        """
        sp = np.sin(self.phi)
        cp = np.cos(self.phi)
        nodes = self.foundation_nodes   # already x ascending
        xs = np.array([n.x for n in nodes])
        sn = np.array([n.sigma * (1.0 + sp) - self.c * cp / sp for n in nodes])

        # Clip to [0, B] with linear interpolation at endpoints
        x_lo, x_hi = 0.0, self.B
        keep = (xs >= x_lo) & (xs <= x_hi)
        x_keep = xs[keep]
        s_keep = sn[keep]

        # Insert interpolated value at x = 0 if needed
        if xs[0] < 0.0 < (x_keep[0] if x_keep.size else x_hi):
            # find bracketing pair around x = 0
            i_hi = np.searchsorted(xs, 0.0)
            i_lo = i_hi - 1
            frac = (0.0 - xs[i_lo]) / (xs[i_hi] - xs[i_lo])
            s_at_0 = sn[i_lo] + frac * (sn[i_hi] - sn[i_lo])
            x_keep = np.concatenate(([0.0], x_keep))
            s_keep = np.concatenate(([s_at_0], s_keep))

        return x_keep, s_keep

    # ------------------------------------------------------------------
    # Region III — uniform-ψ wedge with linear p̄(y)
    # ------------------------------------------------------------------

    def _build_region_III(self):
        n = self.n_radial
        h = self.h
        m_beta = np.tan(self.theta_beta_III)
        denom = 1.0 - np.sin(self.phi)
        rows: list[list[Node]] = []
        for j in range(n + 1):
            row: list[Node] = []
            for i in range(n + 1 - j):
                x = self.B + (i + 0.5 * j) * h
                y = 0.5 * j * h * m_beta
                pbar = self.pbar_B + self.gamma * y / denom
                row.append(Node(x=x, y=y, sigma=pbar, psi=0.0))
            rows.append(row)
        self.region_III_nodes = rows

    # ------------------------------------------------------------------
    # Region II — RIGOROUS (k, m) MOC grid with gravity
    # ------------------------------------------------------------------

    def _build_region_II(self):
        """
        Build Region II as a rigorous (k, m) grid using `riemann_step_mc_gravity`.

        Index conventions match the existing weightless solver:
            region_II_nodes[k][m]
            k = 0 .. n_fan  →  "angular" index along the corner fan (ψ = k·dψ).
            m = 0 .. n_p    →  "radial" index along each β-line (m = 0 at corner / III/II β-line side).

        Boundary data:
            k = 0   : corner singularity — all nodes at (B, 0), ψ_k = k·dψ,
                      p̄_k = p̄_B · exp(2 ψ_k tan φ).         (γ contributes 0 at s = 0.)
            m = 0   : III/II β-line — nodes at distance k·dr from corner along
                      direction (cos μ, sin μ), ψ = 0, p̄ via β-Kötter with γ.

        Interior nodes (k ≥ 1, m ≥ 1):
            R = riemann_step_mc_gravity(P, Q, c, φ, γ)
            with  P = (k, m−1)  (α-parent on the same α-spiral)
                  Q = (k−1, m)  (β-parent on the same β-radial)

        With γ ≠ 0 the radials are no longer straight through the corner and
        the inner-fan-radial (k = n_fan) is shifted from the weightless case,
        which carries the γ-dependence into the foundation pressure.
        """
        n_p   = self.n_radial
        n_fan = self.n_fan
        dpsi  = self.delta_psi / n_fan
        tanphi = np.tan(self.phi)
        mu    = self.theta_beta_III          # π/4 − φ/2
        gamma = self.gamma

        # Allocate region_II_nodes[k][m]
        grid: list[list[Node | None]] = [
            [None] * (n_p + 1) for _ in range(n_fan + 1)
        ]

        # k = 0 : corner-fan column (all at (B, 0))
        for m in range(n_fan + 1):
            # Note: for m = 0 this is the corner of the III/II β-line too.
            psi_m = self.psi_B + m * dpsi
            pbar_m = self.pbar_B * np.exp(2.0 * tanphi * psi_m)
            grid[m][0] = Node(self.B, 0.0, pbar_m, psi_m)

        # m = 0 : III/II β-line (k > 0).  ψ = 0 constant, p̄ from β-Kötter.
        grad_along = gamma * (np.sin(mu) + tanphi * np.cos(mu))
        for k in range(1, n_p + 1):
            s = k * self.dr_outer
            x = self.B + s * np.cos(mu)
            y = s * np.sin(mu)
            pbar = self.pbar_B + grad_along * s
            grid[0][k] = Node(x, y, pbar, 0.0)

        # Interior nodes — Riemann sweep
        for m in range(1, n_fan + 1):
            for k in range(1, n_p + 1):
                P = grid[m - 1][k]      # α-parent: same β-radial (k), previous angular index
                Q = grid[m][k - 1]      # β-parent: same α-spiral (m), previous radial index
                grid[m][k] = riemann_step_mc_gravity(
                    P, Q, self.c, self.phi, gamma, n_iter=5,
                )

        # Cast back to list[list[Node]]
        self.region_II_nodes = [
            [n for n in row] for row in grid
        ]

    # ------------------------------------------------------------------
    # Region I — RIGOROUS triangular MOC mesh with gravity
    # ------------------------------------------------------------------

    def _build_region_I(self):
        """
        Build Region I as a triangular MOC mesh under the foundation.

        Index convention: `region_I_nodes[k][j]`
            k = 0 .. n_p   — α-line index (this α-line emanates from
                             `inner[k] = region_II_nodes[n_fan][k]`
                             and terminates at foundation point F[k]).
            j = 0 .. k     — step index along α-line k.
                             j = 0  → at inner-fan β-line  (= inner[k]).
                             j = k  → at foundation y = 0   (= F[k], ψ = π/2).

        β-line index j (for j ≥ 1) is the β-line emerging from F[j]; it
        contains nodes `region_I_nodes[k][j]` for k = j .. n_p.  β-line
        index j = 0 is the inner-fan β-line itself (handled by Region II).

        Recurrence (1 ≤ j ≤ k − 1, interior):
            region_I_nodes[k][j] = riemann_step_mc_gravity(
                P = region_I_nodes[k][j-1],   # α-parent
                Q = region_I_nodes[k-1][j],   # β-parent
            )

        Foundation node (j = k, on y = 0 with ψ = π/2 prescribed):
            α-shot from P = region_I_nodes[k][k-1].  Straight α-segment
            using ψ_avg = (ψ_P + π/2)/2 for slope; p̄_F solved from
            α-Kötter with gravity by Picard iteration on p̄_avg.
        """
        n_p    = self.n_radial
        n_fan  = self.n_fan
        tanphi = np.tan(self.phi)
        phi    = self.phi
        psi_F  = self.psi_A                  # π/2
        gamma  = self.gamma

        inner = [self.region_II_nodes[n_fan][k] for k in range(n_p + 1)]

        # region_I_nodes[k] is a list of length k+1
        nodes: list[list[Node]] = []

        # k = 0 — degenerate α-line at the corner, also the first foundation node
        F0 = Node(inner[0].x, 0.0, inner[0].sigma, psi_F)
        nodes.append([F0])
        foundation_nodes: list[Node] = [F0]

        # k = 1 .. n_p — build α-line k and foundation node F[k]
        for k in range(1, n_p + 1):
            chain: list[Node] = [inner[k]]   # j = 0
            # Interior: j = 1 .. k − 1
            for j in range(1, k):
                P = chain[j - 1]              # α-parent
                Q = nodes[k - 1][j]           # β-parent
                R = riemann_step_mc_gravity(
                    P, Q, self.c, phi, gamma, n_iter=5,
                )
                chain.append(R)
            # Foundation node F[k] (j = k) — α-shot to y = 0, ψ = π/2
            P = chain[k - 1]
            psi_avg = 0.5 * (P.psi + psi_F)
            slope_a = np.tan(psi_avg - np.pi / 4.0 + phi / 2.0)
            if abs(slope_a) < 1e-12 or P.y < 1e-12:
                F_k = Node(P.x, 0.0, P.sigma, psi_F)
            else:
                x_F  = P.x - P.y / slope_a
                dx_F = x_F - P.x
                dy_F = -P.y
                pbar_F = P.sigma
                a_lhs = (np.log(P.sigma) - 2.0 * tanphi * P.psi
                         + 2.0 * tanphi * psi_F)
                for _ in range(12):
                    pavg = 0.5 * (P.sigma + pbar_F)
                    ln_pbar = a_lhs + (gamma / pavg) * (dy_F - tanphi * dx_F)
                    pbar_new = np.exp(ln_pbar)
                    if abs(pbar_new - pbar_F) < 1e-12:
                        pbar_F = pbar_new
                        break
                    pbar_F = pbar_new
                F_k = Node(x_F, 0.0, pbar_F, psi_F)
            chain.append(F_k)
            nodes.append(chain)
            foundation_nodes.append(F_k)

        self.region_I_nodes = nodes
        self.region_I_inner = inner

        # Sort foundation nodes by x ASCENDING for the contact-pressure curve
        self.foundation_nodes = sorted(foundation_nodes, key=lambda n: n.x)

    # ------------------------------------------------------------------
    # Plot accessors
    # ------------------------------------------------------------------

    def all_alpha_lines(self) -> list[list[Node]]:
        """
        α-lines:
          Region III : diagonals (j, m−j) for m = 1..n.
          Region II  : α-spirals = rows of region_II_nodes (fixed m_angular,
                       varying k_radial).
          Region I   : α-line k = region_I_nodes[k][j] for j = 0..k.

        Each α-line at radial index k threads:
            Region III diagonal (k − 1..0)
          → Region II α-spiral at fixed k (m = 0..n_fan)
          → Region I α-chain (j = 0..k, from inner-fan to foundation F[k]).
        """
        lines: list[list[Node]] = []
        n     = self.n_radial
        n_fan = self.n_fan

        for k in range(1, n + 1):
            seg_III = [self.region_III_nodes[j][k - j]
                       for j in range(k - 1, -1, -1)]
            seg_II  = [self.region_II_nodes[m][k]
                       for m in range(n_fan + 1)]
            # Skip j = 0 of Region I — it equals inner[k] = seg_II[-1]
            seg_I   = [self.region_I_nodes[k][j] for j in range(1, k + 1)]
            line = list(reversed(seg_III)) + seg_II + seg_I
            if len(line) >= 2:
                lines.append(line)
        return lines

    def all_beta_lines(self) -> list[list[Node]]:
        """
        β-lines:
          Region III : columns (j, k) for j = 0..n−k.
          Region II  : β-radials = columns of region_II_nodes (fixed k_radial,
                       varying m_angular).
          Region I   : β-line j emerging from foundation point F[j],
                       = region_I_nodes[k][j] for k = j..n_p (j = 1..n_p).

        The inner-fan β-line (Region II β-radial at k = n_p) is already
        emitted by the Region II loop and naturally extends the j = 0
        β-boundary of Region I.
        """
        lines: list[list[Node]] = []
        n     = self.n_radial
        n_fan = self.n_fan

        # Region III β-lines
        for k in range(n + 1):
            line = [self.region_III_nodes[j][k] for j in range(n + 1 - k)]
            if len(line) >= 2:
                lines.append(line)

        # Region II β-radials (k > 0 to avoid duplicate at corner)
        for k in range(1, n_fan + 1):
            line = [self.region_II_nodes[k][m] for m in range(n + 1)]
            if len(line) >= 2:
                lines.append(line)

        # Region I β-lines emerging from foundation points F[j], j = 1..n_p
        for j in range(1, n + 1):
            line = [self.region_I_nodes[k][j] for k in range(j, n + 1)]
            if len(line) >= 2:
                lines.append(line)

        return lines


# ======================================================================
# Disk foundation solver
# ======================================================================


class DiskPrandtlSolver:
    """
    Slip line field for a partly-penetrated circular disk on Tresca soil.

    Mechanism (right half, by symmetry):
      Region III  — passive wedge at the free surface  (x > B_s),
                    ψ = 0, σ = cu, identical to the flat-footing case but
                    anchored at the contact edge (B_s, 0).
      Region II   — centred fan at the disk edge (B_s, 0).  ψ rotates from
                    ψ_B = 0 (free-surface side) to ψ_edge (disk-edge side).
      Region I    — soil under the disk.  The disk surface is the boundary
                    along which ψ varies smoothly:
                        ψ_disk(φ) = (π/2 − φ) − arcsin(r)/2
                    where φ ∈ [0, θ_max] is the angle from the disk's
                    downward axis at the centre.

    Analytical results (weightless Tresca):
      Along an α-line from the free surface to a disk-surface point at
      angle φ, σ − 2cu·ψ = const, giving
          σ(φ) = cu·(1 + 2·ψ_disk(φ))
      Normal contact stress on the disk surface:
          σ_n(φ) = cu·(1 + π − 2φ − arcsin(r) + √(1 − r²))
      Tangential contact stress:
          |τ(φ)| = r · cu
      Vertical reaction force per unit length of strip:
          Q = 2·R · ∫_0^{θ_max} (σ_n cos φ + |τ| sin φ) dφ
    """

    def __init__(
        self,
        disk: DiskFoundationBoundary,
        n_fan: int = 20,
        n_disk: int = 30,
        n_radial: int | None = None,
    ):
        self.disk = disk
        self.cu = disk.cu
        self.r = disk.r
        self.D = disk.D
        self.p = disk.p
        self.R = disk.R
        self.y_c = disk.y_center
        self.B_s = disk.B_s
        self.theta_max = disk.theta_max

        self.n_fan = int(n_fan)
        self.n_disk = int(n_disk)
        self.n_radial = int(n_radial) if n_radial is not None else self.n_fan

        self.yield_fn = TrescaYield(self.cu)
        self.sigma_B, self.psi_B = self.yield_fn.stress_state_free_surface()

        # Roughness angle ω (= arcsin(r)/2) corresponding to full mobilisation.
        self.omega_full = 0.5 * np.arcsin(self.r)

        # If π/2 − θ_max − ω < 0 the corner-fan mechanism is inadmissible.
        # Instead, the interface roughness is reduced on an outer arc so that
        # ψ_disk = 0 there (parallel β-lines, no fan).  Threshold:
        #     θ_0 = π/2 − ω_full
        #     ω(θ)  = ω_full           for θ ≤ θ_0      (inner zone)
        #     ω(θ)  = π/2 − θ         for θ ∈ (θ_0, θ_max]  (outer zone)
        # so ψ_disk(θ) is piecewise: (π/2 − θ − ω_full) inner, 0 outer.
        self.theta_0 = np.pi / 2.0 - self.omega_full
        self.use_parallel_zone = self.theta_max > self.theta_0 + 1e-12

        if self.use_parallel_zone:
            self.psi_edge = 0.0
            self.delta_psi_fan = 0.0
        else:
            self.psi_edge = (np.pi / 2.0 - self.theta_max) - self.omega_full
            self.delta_psi_fan = self.psi_edge - self.psi_B
        # ψ at θ = 0 (bottom) is always with full roughness (since θ_0 > 0 when r<1).
        self.psi_bot = (np.pi / 2.0) - self.omega_full

        # Mesh sizing — derived inside solve() from Region I geometry so
        # the fan_I β-line (zone-II boundary) exactly coincides with the
        # corner β-line of Region I, and Region III's corner-side β-line
        # exactly matches the k=0 fan radial.
        self.L = max(self.D, 1.0)            # placeholder; overwritten in solve()
        self.h = self.L / self.n_radial
        self.dr = self.h / np.sqrt(2.0)

        self.region_III_nodes: list[list[Node]] = []
        self.region_II_nodes:  list[list[Node]] = []
        self.disk_surface_nodes: list[Node] = []
        self.region_I_alpha_lines: list[list[Node]] = []

        self._built = False

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def solve(self) -> dict:
        """Build the field and return collapse-load summary."""
        # 1. Disk surface (Cauchy boundary for Region I).
        self._build_disk_surface()
        # 2. Region I — purely geometry-driven Riemann mesh from disk surface.
        self._build_region_I()
        # 3. Derive fan/Region-III sizing so the lower boundary is continuous
        #    through all three zones.  Region I's corner β-line is straight
        #    at angle θ_β = ψ_edge + π/4; set fan_I length = that β-line length.
        self._resize_fan_and_wedge_to_match_region_I()
        # 4. Fan (Region II) and passive wedge (Region III).
        self._build_region_II()
        self._build_region_III()
        self._built = True

        Q = self.collapse_load_Q()
        sn_sample = np.array([self._sigma_n(t)
                              for t in np.linspace(0.0, self.theta_max, 401)])
        return {
            'Q_per_unit_length': Q,
            'Q_over_D': Q / self.D,
            'Q_over_D_over_cu': Q / (self.D * self.cu),
            'sigma_n_max': float(sn_sample.max()),
            'sigma_n_at_bottom': self._sigma_n(0.0),
            'sigma_n_at_edge':   self._sigma_n(self.theta_max),
            'theta_max_deg': np.degrees(self.theta_max),
            'B_s': self.B_s,
            'fan_angle_deg': np.degrees(self.delta_psi_fan),
            'parallel_beta_zone': self.use_parallel_zone,
            'theta_0_deg': np.degrees(self.theta_0),
        }

    def collapse_load_Q(self) -> float:
        """
        Vertical reaction force per unit length of strip.

            Q = 2·R · ∫_0^{θ_max} [σ_n(θ) cos θ + |τ(θ)| sin θ] dθ

        Computed numerically (handles the piecewise interface law uniformly).
        """
        ths, sn, tau = self.stress_distribution_curve(2001)
        integrand = sn * np.cos(ths) + tau * np.sin(ths)
        return 2.0 * self.R * float(np.trapz(integrand, ths))

    # ------------------------------------------------------------------
    # Piecewise roughness, σ_n(θ), τ(θ)  (closed form)
    # ------------------------------------------------------------------

    def _omega_eff(self, theta: float) -> float:
        """Mobilised interface roughness angle ω at disk angle θ."""
        if self.use_parallel_zone and theta > self.theta_0:
            return 0.5 * np.pi - theta
        return self.omega_full

    def _psi_disk(self, theta: float) -> float:
        """ψ of the soil at the disk surface at angle θ."""
        return 0.5 * np.pi - theta - self._omega_eff(theta)

    def _r_eff(self, theta: float) -> float:
        """Mobilised interface shear ratio r_eff = sin(2ω_eff) at angle θ."""
        return np.sin(2.0 * self._omega_eff(theta))

    def _sigma_n(self, phi: float) -> float:
        """Normal contact stress at disk-surface angle φ (general piecewise form).

            σ_n = σ + cu·cos(2ω_eff) = cu·(1 + 2ψ_disk + √(1 − r_eff²))
        """
        cu = self.cu
        psi = self._psi_disk(phi)
        r_eff = self._r_eff(phi)
        return cu * (1.0 + 2.0 * psi + np.sqrt(max(0.0, 1.0 - r_eff * r_eff)))

    def _tau(self, phi: float) -> float:
        """Mobilised interface shear magnitude at angle φ: |τ| = r_eff·cu."""
        return self._r_eff(phi) * self.cu

    def _phi_of(self, node: Node) -> float:
        """Recover the disk-surface angle φ stored implicitly via node.psi.

        Valid only in the inner zone (uniform ω = ω_full).
        """
        return 0.5 * np.pi - node.psi - self.omega_full

    def contact_pressure_curve(self, n: int = 100):
        """Return arrays (phi, sigma_n/cu) for plotting."""
        phis = np.linspace(0.0, self.theta_max, n)
        sig  = np.array([self._sigma_n(p) for p in phis])
        return phis, sig

    def stress_distribution_curve(self, n: int = 200):
        """Return (phi, sigma_n, tau) arrays along the contact arc.

        sigma_n : normal stress (compression positive) at each phi
        tau     : shear stress |τ| = r_eff(θ)·cu  (piecewise; equals r·cu in inner zone)
        """
        phis    = np.linspace(0.0, self.theta_max, n)
        sigma_n = np.array([self._sigma_n(p) for p in phis])
        tau     = np.array([self._tau(p)     for p in phis])
        return phis, sigma_n, tau

    # ------------------------------------------------------------------
    # Resultant forces on the half-disk contact arc (per unit out-of-plane length)
    # ------------------------------------------------------------------

    def forces_on_half_disk(self, n: int = 400) -> dict:
        """
        Numerically integrate the contact tractions along the RIGHT-HALF
        disk-soil contact arc (φ ∈ [0, θ_max]) and return:

            F_v   — vertical force (upward, positive)
            F_h   — horizontal force on the half (outward = +x, positive)
            F_v_total = 2·F_v  — full-disk vertical force (= Q by symmetry)
            F_h_total = 0      — full-disk horizontal force (cancels by symmetry)

        Sign convention on the right half (φ measured from disk's downward
        axis toward the corner):
          outward normal to the disk (from disk into soil): n̂ = (sin φ, cos φ)
          tangent (from bottom toward corner):              t̂ = (cos φ, -sin φ)
          Soil pushes disk along -n̂ with magnitude σ_n.
          This is a HORIZONTAL-PUSH problem: the disk slides sideways, so the
          soil at the right-half contact slides along -t̂ relative to the disk
          and friction on the disk acts along -t̂ with magnitude |τ| = r·cu.
        Vertical (upward) component on the right half is therefore
            dF_v = (σ_n cos φ − r·cu · sin φ) · R dφ
        Horizontal (outward, +x) component:
            dF_h = (−σ_n sin φ − r·cu · cos φ) · R dφ
        """
        R = self.R
        tmax = self.theta_max

        phi = np.linspace(0.0, tmax, max(2, int(n)))
        sigma_n = np.array([self._sigma_n(p) for p in phi])
        tau     = np.array([self._tau(p)     for p in phi])

        dF_v = sigma_n * np.cos(phi) - tau * np.sin(phi)
        dF_h = -sigma_n * np.sin(phi) - tau * np.cos(phi)

        F_v = R * float(np.trapz(dF_v, phi))
        F_h = R * float(np.trapz(dF_h, phi))

        return {
            'F_v': F_v,
            'F_h': F_h,
            'F_v_total': 2.0 * F_v,
            'F_h_total': 0.0,
        }

    # ------------------------------------------------------------------
    # Sizing — set dr, h, L so the lower boundary is continuous across zones
    # ------------------------------------------------------------------

    def _resize_fan_and_wedge_to_match_region_I(self):
        """
        Set self.dr (fan radial step), self.h (Region III surface spacing) and
        self.L (Region III extent) so that:
          • fan_I's deep endpoint coincides with Region I's corner β-line
            deep endpoint (mesh[-1][0]);
          • Region III's corner-side β-line matches the k=0 fan radial.

        Region I's corner β-line is straight at angle θ_β = ψ_edge + π/4
        (verified analytically — ψ and σ are both constant along it).
        Hence we only need to match its LENGTH:  n_radial · dr = |mesh[-1][0] − corner|.
        Region III then uses h = dr·√2 so its 45° corner-side β-line (length
        n_radial·h/√2 = n_radial·dr) matches the k=0 fan radial.
        """
        mesh = getattr(self, 'region_I_mesh', None)
        if not mesh or len(mesh) < 2 or len(mesh[0]) < 2:
            # Nothing to match against — keep placeholder sizing.
            return
        corner = mesh[0][-1]          # disk corner (B_s, 0)
        deep_end = mesh[-1][0]        # tip of Region I corner β-line
        dx = deep_end.x - corner.x
        dy = deep_end.y - corner.y
        length = float(np.hypot(dx, dy))
        if length <= 1e-12 or self.n_radial < 1:
            return
        self.dr = length / self.n_radial
        self.h = self.dr * np.sqrt(2.0)
        self.L = self.n_radial * self.h

    # ------------------------------------------------------------------
    # Region III (identical layout to flat case, anchored at B_s)
    # ------------------------------------------------------------------

    def _build_region_III(self):
        n = self.n_radial
        h = self.h
        rows: list[list[Node]] = []
        for j in range(n + 1):
            row: list[Node] = []
            for i in range(n + 1 - j):
                x = self.B_s + (i + 0.5 * j) * h
                y = 0.5 * j * h
                row.append(Node(x=x, y=y, sigma=self.cu, psi=0.0))
            rows.append(row)
        self.region_III_nodes = rows

    # ------------------------------------------------------------------
    # Region II: centred fan at the disk edge (B_s, 0)
    # ------------------------------------------------------------------

    def _build_region_II(self):
        n_p = self.n_radial
        cu = self.cu

        if self.use_parallel_zone:
            # No Region II at all; Region III meets the parallel-β zone of
            # Region I directly along the corner β-line.
            self.region_II_nodes = []
            return

        if self.n_fan < 1 or abs(self.delta_psi_fan) < 1e-12:
            # Degenerate fan (e.g. p = 0.5 with r = 0 makes ψ_edge = ψ_B = 0).
            # Build a single radial at angle π/4 from the corner, of length
            # n_radial·dr — this coincides with both the Region III corner-side
            # β-line and the Region I corner β-line.
            theta0 = np.pi / 4.0 + self.psi_B
            sigma0 = cu * (1.0 + 2.0 * self.psi_B)
            radial = [
                Node(
                    x=self.B_s + (j * self.dr) * np.cos(theta0),
                    y=(j * self.dr) * np.sin(theta0),
                    sigma=sigma0,
                    psi=self.psi_B,
                )
                for j in range(n_p + 1)
            ]
            self.region_II_nodes = [radial]
            return

        dpsi = self.delta_psi_fan / self.n_fan

        radials: list[list[Node]] = []
        for k in range(self.n_fan + 1):
            psi_k = self.psi_B + k * dpsi
            theta_k = np.pi / 4.0 + psi_k
            sigma_k = cu * (1.0 + 2.0 * psi_k)
            line: list[Node] = []
            for j in range(n_p + 1):
                rj = j * self.dr
                x = self.B_s + rj * np.cos(theta_k)
                y = rj * np.sin(theta_k)
                line.append(Node(x=x, y=y, sigma=sigma_k, psi=psi_k))
            radials.append(line)
        self.region_II_nodes = radials

    # ------------------------------------------------------------------
    # Disk-surface nodes: ψ_disk(φ) varies; σ from Hencky
    # ------------------------------------------------------------------

    def _build_disk_surface(self):
        """Place n_disk + 1 nodes on the disk-soil contact (right half).

        When `use_parallel_zone` is active, the chain is densified at the
        transition θ_0 to capture the kink in dψ/dθ (slope discontinuity).
        """
        if self.use_parallel_zone:
            # Distribute n_disk + 1 nodes with one exactly on θ = θ_0.
            n_inner = max(1, int(round(self.n_disk * self.theta_0 / self.theta_max)))
            n_outer = max(1, self.n_disk - n_inner)
            inner = list(np.linspace(0.0, self.theta_0, n_inner + 1))
            outer = list(np.linspace(self.theta_0, self.theta_max, n_outer + 1))[1:]
            thetas = list(reversed(inner + outer))   # corner → bottom
        else:
            thetas = [self.theta_max * (1.0 - i / self.n_disk)
                      for i in range(self.n_disk + 1)]

        nodes: list[Node] = []
        for th in thetas:
            x, y = self.disk.position(th)
            psi = self._psi_disk(th)
            sigma = self.cu * (1.0 + 2.0 * psi)
            nodes.append(Node(x=x, y=y, sigma=sigma, psi=psi))
        self.disk_surface_nodes = nodes

    # ------------------------------------------------------------------
    # Region I (under disk): build a characteristic Riemann mesh emerging
    # from the disk surface (a non-characteristic stress boundary).
    # ψ varies along the disk, so α- and β-lines are curved.
    # ------------------------------------------------------------------

    def _build_region_I(self):
        """
        Goursat / Massau construction.

        Start with the chain of (n_disk + 1) disk-surface nodes ordered
        from the disk *bottom* (S_0) to the disk *edge / corner* (S_N).
        Each pair (S_i, S_{i+1}) → one new interior node by a Riemann
        step, where S_i is the α-parent and S_{i+1} the β-parent.
        Repeat until the triangular mesh is exhausted.

        Resulting indexing:
            mesh[k][i]  with  k = depth-row (0 = on disk), i = column
        α-line  :  fixed column i, k varies  →  mesh[0][i], mesh[1][i], …
        β-line  :  i + k = const,            →  mesh[0][j], mesh[1][j-1], …
        """
        cu = self.cu
        # disk_surface_nodes are stored corner→bottom; we want bottom→corner.
        chain = list(reversed(self.disk_surface_nodes))
        mesh: list[list[Node]] = [chain]

        for _ in range(len(chain) - 1):
            prev = mesh[-1]
            row: list[Node] = []
            for i in range(len(prev) - 1):
                try:
                    N = riemann_step(prev[i], prev[i + 1], cu)
                except Exception:
                    N = None
                if (N is None or not np.isfinite(N.x) or not np.isfinite(N.y)
                        or N.y < -1e-6):
                    break
                row.append(N)
            if len(row) < 1:
                break
            mesh.append(row)

        self.region_I_mesh = mesh
        # Keep legacy attribute name for any external references; unused now.
        self.region_I_alpha_lines = []

    # ------------------------------------------------------------------
    # Plot helpers
    # ------------------------------------------------------------------

    def all_alpha_lines(self) -> list[list[Node]]:
        lines: list[list[Node]] = []

        # ---- Region I (disk-under): α-line = fixed column across rows ----
        mesh = getattr(self, 'region_I_mesh', None)
        if mesh:
            n_cols = len(mesh[0])
            for i in range(n_cols):
                line: list[Node] = []
                for row in mesh:
                    if i < len(row):
                        line.append(row[i])
                    else:
                        break
                if len(line) >= 2:
                    lines.append(line)

        # ---- Region II (corner fan): α-lines = arcs at fixed radial index --
        n_p = self.n_radial
        if self.region_II_nodes and len(self.region_II_nodes) > 1:
            for j in range(1, n_p + 1):
                arc = [self.region_II_nodes[k][j]
                       for k in range(len(self.region_II_nodes))
                       if j < len(self.region_II_nodes[k])]
                if len(arc) >= 2:
                    lines.append(arc)

        # ---- Region III (passive wedge): diagonals ----
        n = self.n_radial
        for i_surf in range(n + 1):
            line = []
            for j in range(min(i_surf, n) + 1):
                k = i_surf - j
                if k < 0 or k > n - j:
                    break
                line.append(self.region_III_nodes[j][k])
            if len(line) >= 2:
                lines.append(line)

        return lines

    def all_beta_lines(self) -> list[list[Node]]:
        lines: list[list[Node]] = []

        # ---- Region I (disk-under): β-line = anti-diagonal (i+k = const) --
        mesh = getattr(self, 'region_I_mesh', None)
        if mesh:
            n_cols0 = len(mesh[0])
            for j_start in range(n_cols0):
                line: list[Node] = []
                for k in range(len(mesh)):
                    idx = j_start - k
                    if 0 <= idx < len(mesh[k]):
                        line.append(mesh[k][idx])
                    else:
                        break
                if len(line) >= 2:
                    lines.append(line)

        # ---- Region II (corner fan): β-lines = radials ----
        if self.region_II_nodes:
            for radial in self.region_II_nodes:
                if len(radial) >= 2:
                    lines.append(radial)

        # ---- Region III (passive wedge): vertical columns ----
        n = self.n_radial
        for k in range(n + 1):
            line = [self.region_III_nodes[j][k] for j in range(n + 1 - k)]
            if len(line) >= 2:
                lines.append(line)

        return lines


# ======================================================================
# Mohr-Coulomb disk foundation solver
# ======================================================================


class MohrCoulombDiskPrandtlSolver:
    """
    Slip-line field for a partly-penetrated circular disk on weightless
    Mohr-Coulomb soil (cohesion c, friction angle φ).

    Mechanism (right half, by symmetry):
      Region III  — passive wedge at the free surface (x > B_s),
                    ψ = 0, p̄ = p̄_B, anchored at the contact edge (B_s, 0).
      Region II   — centred fan at the disk edge (B_s, 0); ψ rotates from
                    ψ_B = 0 to ψ_edge, β-lines are straight radials, α-lines
                    are logarithmic spirals r = A·exp(+(Δψ−ψ_k)·tan φ).
      Region I    — soil under the disk; non-uniform ψ along the disk
                    surface drives a curved characteristic mesh built by
                    successive Riemann steps (`riemann_step_mc`).

    Boundary conditions on the disk surface (angle θ from downward axis):
        outward-normal angle from +x:  θ_n = π/2 − θ
        ψ_disk(θ) = θ_n − arcsin(r) / 2 = π/2 − θ − arcsin(r)/2
        p̄_disk(θ) = p̄_B · exp(2·ψ_disk(θ)·tan φ)
                  = p̄_B · exp((π − 2θ − arcsin(r))·tan φ)
        σ_n(θ)    = p̄·(1 + sin φ · √(1−r²)) − c·cot φ
        |τ(θ)|    = r · p̄ · sin φ

    Roughness ratio r ∈ [0, 1] is interpreted as the mobilisation
    fraction of the Mohr-circle radius (so r·p̄·sin φ is the actual
    interface shear — same convention as the Tresca disk solver).

    Vertical collapse load (full disk, per unit out-of-plane length):
        Q = 2·R · ∫_0^{θ_max} [σ_n(θ) cos θ + |τ(θ)| sin θ] dθ
    (computed numerically).
    """

    def __init__(
        self,
        disk: 'DiskFoundationBoundary',
        c: float,
        phi: float,
        n_fan: int = 20,
        n_disk: int = 30,
        n_radial: int | None = None,
    ):
        # Geometry from the existing boundary helper
        self.disk = disk
        self.D = disk.D
        self.p = disk.p
        self.R = disk.R
        self.y_c = disk.y_center
        self.B_s = disk.B_s
        self.theta_max = disk.theta_max
        self.r = disk.r

        # Material
        self.c = float(c)
        self.phi = float(phi)
        self.yield_fn = MohrCoulombYield(self.c, self.phi)

        self.n_fan = int(n_fan)
        self.n_disk = int(n_disk)
        self.n_radial = int(n_radial) if n_radial is not None else self.n_fan

        # Boundary stress states
        self.pbar_B, self.psi_B = self.yield_fn.stress_state_free_surface()

        # Roughness angle ω (= arcsin(r)/2) at full mobilisation.
        self.omega_full = 0.5 * np.arcsin(self.r)

        # When π/2 − θ_max − ω < 0, the corner fan would have negative
        # opening and is inadmissible. In that case the interface roughness
        # is reduced on an outer arc so that ψ_disk = 0 there (parallel
        # β-lines, no fan).  Same threshold as the Tresca disk:
        #     θ_0 = π/2 − ω_full
        #     ω(θ) = ω_full           for θ ≤ θ_0      (inner zone)
        #     ω(θ) = π/2 − θ         for θ ∈ (θ_0, θ_max]  (outer zone)
        # so ψ_disk(θ) is piecewise: (π/2 − θ − ω_full) inner, 0 outer.
        self.theta_0 = np.pi / 2.0 - self.omega_full
        self.use_parallel_zone = self.theta_max > self.theta_0 + 1e-12

        if self.use_parallel_zone:
            self.psi_edge = 0.0
            self.delta_psi_fan = 0.0
        else:
            self.psi_edge = (np.pi / 2.0 - self.theta_max) - self.omega_full
            self.delta_psi_fan = self.psi_edge - self.psi_B
        self.psi_bot = (np.pi / 2.0) - self.omega_full

        # Mesh sizing — placeholders; recomputed from Region I in solve().
        self.L = max(self.D, 1.0)
        self.h = self.L / self.n_radial
        self.dr = self.h / (2.0 * np.cos(np.pi / 4.0 - self.phi / 2.0)
                            * np.exp(self.delta_psi_fan * np.tan(self.phi)))

        self.region_III_nodes: list[list[Node]] = []
        self.region_II_nodes:  list[list[Node]] = []
        self.disk_surface_nodes: list[Node] = []
        self.region_I_mesh: list[list[Node]] = []

        self._built = False

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def solve(self) -> dict:
        self._build_disk_surface()
        self._build_region_I()
        self._resize_fan_and_wedge_to_match_region_I()
        self._build_region_II()
        self._build_region_III()
        self._built = True

        Q = self.collapse_load_Q()
        sn_sample = np.array([self._sigma_n(t)
                              for t in np.linspace(0.0, self.theta_max, 401)])
        return {
            'Q_per_unit_length': Q,
            'Q_over_D': Q / self.D,
            'Q_over_Dc': Q / (self.D * self.c),
            'sigma_n_max':       float(sn_sample.max()),
            'sigma_n_at_bottom': self._sigma_n(0.0),
            'sigma_n_at_edge':   self._sigma_n(self.theta_max),
            'theta_max_deg':     np.degrees(self.theta_max),
            'B_s':               self.B_s,
            'fan_angle_deg':     np.degrees(self.delta_psi_fan),
            'parallel_beta_zone': self.use_parallel_zone,
            'theta_0_deg':       np.degrees(self.theta_0),
            'Nc_flat_smooth':    self.yield_fn.Nc(),
        }

    # ------------------------------------------------------------------
    # Piecewise roughness, p̄(θ), σ_n(θ), τ(θ)
    # ------------------------------------------------------------------

    def _omega_eff(self, theta: float) -> float:
        """Mobilised interface roughness angle ω at disk angle θ."""
        if self.use_parallel_zone and theta > self.theta_0:
            return 0.5 * np.pi - theta
        return self.omega_full

    def _psi_disk(self, theta: float) -> float:
        """ψ of the soil at the disk surface at angle θ."""
        return 0.5 * np.pi - theta - self._omega_eff(theta)

    def _r_eff(self, theta: float) -> float:
        """Mobilised interface ratio r_eff = sin(2ω_eff) at angle θ.

        Interface shear is then τ = r_eff · p̄ · sinφ.
        """
        return np.sin(2.0 * self._omega_eff(theta))

    def _pbar_disk(self, theta: float) -> float:
        """p̄ on disk at angle θ from downward axis (piecewise)."""
        psi = self._psi_disk(theta)
        return self.pbar_B * np.exp(2.0 * psi * np.tan(self.phi))

    def _sigma_n(self, theta: float) -> float:
        """Normal contact stress on disk at angle θ (general piecewise form).

            σ_n = p̄·(1 + sinφ · √(1 − r_eff²)) − c·cotφ
        """
        sp = np.sin(self.phi)
        r_eff = self._r_eff(theta)
        rt = np.sqrt(max(0.0, 1.0 - r_eff * r_eff))
        pbar = self._pbar_disk(theta)
        return pbar * (1.0 + sp * rt) - self.c * np.cos(self.phi) / sp

    def _tau(self, theta: float) -> float:
        """Interface shear (magnitude) on disk at angle θ: τ = r_eff·p̄·sinφ."""
        return self._r_eff(theta) * self._pbar_disk(theta) * np.sin(self.phi)

    def contact_pressure_curve(self, n: int = 200):
        """Return (theta, sigma_n/c) for plotting."""
        ths = np.linspace(0.0, self.theta_max, n)
        sig = np.array([self._sigma_n(t) for t in ths])
        return ths, sig

    def stress_distribution_curve(self, n: int = 200):
        """Return (theta, sigma_n, tau) along the contact arc."""
        ths    = np.linspace(0.0, self.theta_max, n)
        sigma_n = np.array([self._sigma_n(t) for t in ths])
        tau     = np.array([self._tau(t)     for t in ths])
        return ths, sigma_n, tau

    def collapse_load_Q(self) -> float:
        """
        Vertical reaction Q (full disk, per unit length, upward positive).

        Q = 2·R · ∫_0^{θ_max} [σ_n(θ) cos θ + |τ(θ)| sin θ] dθ
        """
        ths, sn, tau = self.stress_distribution_curve(2001)
        integrand = sn * np.cos(ths) + tau * np.sin(ths)
        return 2.0 * self.R * float(np.trapz(integrand, ths))

    def forces_on_half_disk(self, n: int = 2001) -> dict:
        """
        Numerically integrate tractions on the RIGHT-HALF contact arc.

        Sign convention (right half, θ from downward axis):
          outward disk normal: n̂ = (sin θ, cos θ)
          tangent (bottom→corner): t̂ = (cos θ, −sin θ)
          soil pushes disk along −n̂ (magnitude σ_n).
          HORIZONTAL-PUSH: friction on disk along −t̂ (magnitude |τ|).
          dF_v = (σ_n cos θ − |τ| sin θ) · R dθ   (upward)
          dF_h = (−σ_n sin θ − |τ| cos θ) · R dθ  (resisting, −x)
        """
        ths    = np.linspace(0.0, self.theta_max, max(2, int(n)))
        sigma_n = np.array([self._sigma_n(t) for t in ths])
        tau     = np.array([self._tau(t)     for t in ths])

        dF_v = sigma_n * np.cos(ths) - tau * np.sin(ths)
        dF_h = -sigma_n * np.sin(ths) - tau * np.cos(ths)

        F_v = self.R * float(np.trapz(dF_v, ths))
        F_h = self.R * float(np.trapz(dF_h, ths))
        return {
            'F_v': F_v,
            'F_h': F_h,
            'F_v_total': 2.0 * F_v,
            'F_h_total': 0.0,
        }

    # ------------------------------------------------------------------
    # Sizing: match fan/wedge to Region I's corner β-line
    # ------------------------------------------------------------------

    def _resize_fan_and_wedge_to_match_region_I(self):
        """
        Region I's corner β-line is straight at angle θ_β = ψ_edge + π/4 − φ/2
        (ψ and p̄ are both constant along that β-line because the α-parent
        at the corner has ψ = ψ_edge).  We measure its length numerically
        from the mesh and set:

          dr = L_corner / n_radial            (fan_I radial = corner β-line)
          h  = 2·dr·exp(Δψ·tan φ)·cos(π/4 − φ/2)
                                              (Region III α-line spacing
                                              so that its corner-side β-line
                                              matches the k=0 fan radial).
        """
        mesh = self.region_I_mesh
        if not mesh or len(mesh) < 2 or len(mesh[0]) < 2:
            return
        corner   = mesh[0][-1]
        deep_end = mesh[-1][0]
        length = float(np.hypot(deep_end.x - corner.x, deep_end.y - corner.y))
        if length <= 1e-12 or self.n_radial < 1:
            return
        self.dr = length / self.n_radial
        self.h  = (2.0 * self.dr
                   * np.exp(self.delta_psi_fan * np.tan(self.phi))
                   * np.cos(np.pi / 4.0 - self.phi / 2.0))
        self.L  = self.n_radial * self.h

    # ------------------------------------------------------------------
    # Disk-surface boundary nodes
    # ------------------------------------------------------------------

    def _build_disk_surface(self):
        """Place (n_disk + 1) Cauchy-boundary nodes along the contact arc.

        When `use_parallel_zone` is active, the chain is densified at the
        transition θ_0 so that one node lies exactly on it (the disk-surface
        ψ has a kink there).
        """
        if self.use_parallel_zone:
            n_inner = max(1, int(round(self.n_disk * self.theta_0 / self.theta_max)))
            n_outer = max(1, self.n_disk - n_inner)
            inner = list(np.linspace(0.0, self.theta_0, n_inner + 1))
            outer = list(np.linspace(self.theta_0, self.theta_max, n_outer + 1))[1:]
            thetas = list(reversed(inner + outer))   # corner → bottom
        else:
            thetas = [self.theta_max * (1.0 - i / self.n_disk)
                      for i in range(self.n_disk + 1)]

        nodes: list[Node] = []
        for th in thetas:
            x, y = self.disk.position(th)
            psi = self._psi_disk(th)
            pbar = self.pbar_B * np.exp(2.0 * psi * np.tan(self.phi))
            nodes.append(Node(x=x, y=y, sigma=pbar, psi=psi))
        self.disk_surface_nodes = nodes

    # ------------------------------------------------------------------
    # Region I: characteristic Riemann mesh under the disk
    # ------------------------------------------------------------------

    def _build_region_I(self):
        """
        Goursat construction from the disk-surface chain (bottom → corner).
        Each adjacent pair drives one `riemann_step_mc` to produce one
        interior node; iterate until the triangular mesh is exhausted.

        Indexing identical to the Tresca disk solver:
            mesh[k][i]  with k = depth row (0 = on disk), i = column
            α-line  : fixed i, k varies   (parents are α-side)
            β-line  : i + k = const       (parents are β-side)
        """
        chain = list(reversed(self.disk_surface_nodes))     # bottom → corner
        mesh: list[list[Node]] = [chain]

        for _ in range(len(chain) - 1):
            prev = mesh[-1]
            row: list[Node] = []
            for i in range(len(prev) - 1):
                try:
                    N = riemann_step_mc(prev[i], prev[i + 1], self.c, self.phi)
                except Exception:
                    N = None
                if (N is None or not np.isfinite(N.x) or not np.isfinite(N.y)
                        or N.y < -1e-6):
                    break
                row.append(N)
            if len(row) < 1:
                break
            mesh.append(row)

        self.region_I_mesh = mesh

    # ------------------------------------------------------------------
    # Region II: log-spiral fan at the disk corner (B_s, 0)
    # ------------------------------------------------------------------

    def _build_region_II(self):
        """
        Centred fan at corner (B_s, 0).  Same construction as flat-MC fan
        but anchored at the disk corner and with Δψ = ψ_edge (not π/2).
            β-lines : straight radials at θ_k = ψ_k + π/4 − φ/2
            α-lines : log spirals — node (k, j) at distance
                      r_{k,j} = j·dr · exp((Δψ − ψ_k)·tan φ)
        """
        n_p   = self.n_radial
        phi   = self.phi
        tphi  = np.tan(phi)

        if self.use_parallel_zone:
            # No Region II at all; Region III meets the parallel-β zone of
            # Region I directly along the corner β-line.
            self.region_II_nodes = []
            return

        if self.n_fan < 1 or abs(self.delta_psi_fan) < 1e-12:
            # Degenerate fan (ψ_edge ≈ 0): one radial at angle π/4 − φ/2
            theta0 = self.psi_B + np.pi / 4.0 - phi / 2.0
            pbar0  = self.pbar_B * np.exp(2.0 * tphi * self.psi_B)
            radial = [
                Node(
                    x=self.B_s + (j * self.dr) * np.cos(theta0),
                    y=(j * self.dr) * np.sin(theta0),
                    sigma=pbar0,
                    psi=self.psi_B,
                )
                for j in range(n_p + 1)
            ]
            self.region_II_nodes = [radial]
            return

        dpsi = self.delta_psi_fan / self.n_fan
        radials: list[list[Node]] = []
        for k in range(self.n_fan + 1):
            psi_k = self.psi_B + k * dpsi
            theta_k = psi_k + np.pi / 4.0 - phi / 2.0
            pbar_k  = self.pbar_B * np.exp(2.0 * tphi * psi_k)
            spiral_scale = np.exp((self.delta_psi_fan - psi_k) * tphi)
            line: list[Node] = []
            for j in range(n_p + 1):
                rj = j * self.dr * spiral_scale
                line.append(Node(
                    x=self.B_s + rj * np.cos(theta_k),
                    y=rj * np.sin(theta_k),
                    sigma=pbar_k,
                    psi=psi_k,
                ))
            radials.append(line)
        self.region_II_nodes = radials

    # ------------------------------------------------------------------
    # Region III: passive wedge at the free surface (anchored at B_s)
    # ------------------------------------------------------------------

    def _build_region_III(self):
        """
        Uniform-stress region: ψ = 0, p̄ = p̄_B.
        Node (j, i): x = B_s + (i + j/2)·h,  y = j·h/2·tan(π/4 − φ/2).
        """
        n = self.n_radial
        h = self.h
        m_beta = np.tan(np.pi / 4.0 - self.phi / 2.0)
        rows: list[list[Node]] = []
        for j in range(n + 1):
            row: list[Node] = []
            for i in range(n + 1 - j):
                x = self.B_s + (i + 0.5 * j) * h
                y = 0.5 * j * h * m_beta
                row.append(Node(x=x, y=y, sigma=self.pbar_B, psi=0.0))
            rows.append(row)
        self.region_III_nodes = rows

    # ------------------------------------------------------------------
    # Plot accessors (same layout as Tresca disk solver)
    # ------------------------------------------------------------------

    def all_alpha_lines(self) -> list[list[Node]]:
        lines: list[list[Node]] = []

        # Region I: α-line = fixed column across mesh rows
        mesh = self.region_I_mesh
        if mesh:
            n_cols = len(mesh[0])
            for i in range(n_cols):
                line: list[Node] = []
                for row in mesh:
                    if i < len(row):
                        line.append(row[i])
                    else:
                        break
                if len(line) >= 2:
                    lines.append(line)

        # Region II: α-line = log spiral (one node per radial at fixed j)
        n_p = self.n_radial
        if self.region_II_nodes and len(self.region_II_nodes) > 1:
            for j in range(1, n_p + 1):
                arc = [self.region_II_nodes[k][j]
                       for k in range(len(self.region_II_nodes))
                       if j < len(self.region_II_nodes[k])]
                if len(arc) >= 2:
                    lines.append(arc)

        # Region III: α-line = diagonal (i + j = const)
        n = self.n_radial
        for i_surf in range(n + 1):
            line = []
            for j in range(min(i_surf, n) + 1):
                k = i_surf - j
                if k < 0 or k > n - j:
                    break
                line.append(self.region_III_nodes[j][k])
            if len(line) >= 2:
                lines.append(line)

        return lines

    def all_beta_lines(self) -> list[list[Node]]:
        lines: list[list[Node]] = []

        # Region I: β-line = anti-diagonal (i + k = const)
        mesh = self.region_I_mesh
        if mesh:
            n_cols0 = len(mesh[0])
            for j_start in range(n_cols0):
                line: list[Node] = []
                for k in range(len(mesh)):
                    idx = j_start - k
                    if 0 <= idx < len(mesh[k]):
                        line.append(mesh[k][idx])
                    else:
                        break
                if len(line) >= 2:
                    lines.append(line)

        # Region II: β-line = radial
        if self.region_II_nodes:
            for radial in self.region_II_nodes:
                if len(radial) >= 2:
                    lines.append(radial)

        # Region III: β-line = column (constant i)
        n = self.n_radial
        for k in range(n + 1):
            line = [self.region_III_nodes[j][k] for j in range(n + 1 - k)]
            if len(line) >= 2:
                lines.append(line)

        return lines
