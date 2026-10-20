"""
Velocity field on a Tresca slip-line net (Geiringer equations).

The velocity is resolved along the slip lines,  V = u·e_α + v·e_β,  with
(y-down, angles CCW in (x, y_down)):

    e_α = ( cos θ,  sin θ),   e_β = (−sin θ,  cos θ),   θ = ψ − π/4

Zero extension rate along both slip-line families (Tresca, associated flow):

    along α:  du − v·dψ = 0
    along β:  dv + u·dψ = 0

Velocity discontinuities can only occur across slip lines, are purely
tangential, and have constant magnitude along the line.
"""

from __future__ import annotations

import numpy as np

from .slip_line_field import Node
from .solver import PrandtlSolver


def alpha_beta_to_xy(psi: float, u: float, v: float) -> tuple[float, float]:
    """Convert slip-line components (u along α, v along β) to (vx, vy)."""
    theta = psi - np.pi / 4.0
    c, s = np.cos(theta), np.sin(theta)
    return u * c - v * s, u * s + v * c


def geiringer_step(
    P: Node,
    Q: Node,
    R: Node,
    uv_P: tuple[float, float],
    uv_Q: tuple[float, float],
) -> tuple[float, float]:
    """
    Velocity (u_R, v_R) at node R from its two parents (trapezoidal rule).

      P : α-neighbour of R  (same α-line)
      Q : β-neighbour of R  (same β-line)

        u_R − u_P = +½·(v_P + v_R)·(ψ_R − ψ_P)
        v_R − v_Q = −½·(u_Q + u_R)·(ψ_R − ψ_Q)
    """
    u_P, v_P = uv_P
    u_Q, v_Q = uv_Q
    a = 0.5 * (R.psi - P.psi)
    b = 0.5 * (R.psi - Q.psi)
    # u_R − a·v_R = u_P + a·v_P
    # b·u_R + v_R = v_Q − b·u_Q
    rhs_a = u_P + a * v_P
    rhs_b = v_Q - b * u_Q
    det = 1.0 + a * b
    return (rhs_a + a * rhs_b) / det, (rhs_b - b * rhs_a) / det


def geiringer_step_vertical_bc(
    Q: Node,
    R: Node,
    uv_Q: tuple[float, float],
    v0: float,
) -> tuple[float, float]:
    """
    Velocity at a foundation node R from its β-neighbour Q and the
    smooth-contact condition vy = v0 (normal velocity of the footing):

        u_R·sin θ_R + v_R·cos θ_R = v0
        v_R − v_Q = −½·(u_Q + u_R)·(ψ_R − ψ_Q)
    """
    u_Q, v_Q = uv_Q
    theta = R.psi - np.pi / 4.0
    s, c = np.sin(theta), np.cos(theta)
    b = 0.5 * (R.psi - Q.psi)
    rhs_b = v_Q - b * u_Q
    # s·u_R + c·v_R = v0  and  b·u_R + v_R = rhs_b
    u_R = (v0 - c * rhs_b) / (s - c * b)
    return u_R, rhs_b - b * u_R


class PrandtlVelocityField:
    """
    Velocity field for a PrandtlSolver net: the footing is pushed vertically
    into the soil with speed v0 (y-down, so v0 > 0 is downward).

    Boundary conditions (right half, x ≥ 0):
      • Outer α-line (Region I column 0 → fan outer arc → Region III
        boundary) borders soil at rest: v = 0 on the plastic side, and
        du = v·dψ = 0 gives u = U (the velocity jump) along the whole line.
      • Foundation (smooth contact): vy = v0.  At the symmetry point (0, 0)
        this fixes U = v0 / sin(ψ_A − π/4).
      • Free surface: no kinematic condition.

    Marching order: Region I (mixed problem) → fan → Region III (Goursat).

    For r > 0 the PrandtlSolver stress field is itself only approximate
    (see TrescaYield.stress_state_foundation); this class is meant for r = 0.

    Velocities are stored as (u, v) per node, in the same layout as the
    solver's region_*_nodes lists.
    """

    def __init__(self, solver: PrandtlSolver, v0: float = 1.0):
        if not solver._built:
            solver.solve()
        self.solver = solver
        self.v0 = float(v0)
        self.U = 0.0
        self.region_I_uv: list[list[tuple[float, float]]] = []
        self.region_II_uv: list[list[tuple[float, float]]] = []
        self.region_III_uv: list[list[tuple[float, float]]] = []

    def solve(self) -> float:
        """Build the velocity field; return U, the jump across the outer α-line."""
        self._solve_region_I()
        self._solve_region_II()
        self._solve_region_III()
        return self.U

    # ------------------------------------------------------------------
    # Region I — mixed problem: outer α-line (column 0) + foundation
    # ------------------------------------------------------------------

    def _solve_region_I(self):
        nodes = self.solver.region_I_nodes
        N = len(nodes) - 1

        # U from the symmetry point (0, 0): v = 0 and vy = v0.
        origin = nodes[0][0]
        self.U = self.v0 / np.sin(origin.psi - np.pi / 4.0)

        uv: list[list[tuple[float, float]]] = [[None] * (i + 1) for i in range(N + 1)]
        for i in range(N + 1):
            uv[i][0] = (self.U, 0.0)  # column 0 = outer α-line

        # Column j (α-line from foundation anchor j), marching away from it.
        # α-neighbour of (i, j) is (i−1, j); β-neighbour is (i, j−1).
        for j in range(1, N + 1):
            uv[j][j] = geiringer_step_vertical_bc(
                nodes[j][j - 1], nodes[j][j], uv[j][j - 1], self.v0
            )
            for i in range(j + 1, N + 1):
                uv[i][j] = geiringer_step(
                    nodes[i - 1][j],
                    nodes[i][j - 1],
                    nodes[i][j],
                    uv[i - 1][j],
                    uv[i][j - 1],
                )
        self.region_I_uv = uv

    # ------------------------------------------------------------------
    # Region II — centred fan, Goursat problem from fan-I radial + outer arc
    # ------------------------------------------------------------------

    def _solve_region_II(self):
        fan = self.solver.region_II_nodes  # [k][j]
        n_fan = len(fan) - 1
        n = len(fan[0]) - 1
        N_I = len(self.region_I_uv) - 1

        uv: list[list[tuple[float, float]]] = [
            [None] * (n + 1) for _ in range(n_fan + 1)
        ]
        # Radial k = n_fan coincides with Region I row N_I (reversed).
        for j in range(n + 1):
            uv[n_fan][j] = self.region_I_uv[N_I][N_I - j]

        # α-neighbour of (k, j) is (k+1, j) on the same arc;
        # β-neighbour is (k, j+1) on the same radial.
        for k in range(n_fan - 1, -1, -1):
            uv[k][n] = (self.U, 0.0)  # outer arc = outer α-line
            for j in range(n - 1, -1, -1):
                uv[k][j] = geiringer_step(
                    fan[k + 1][j],
                    fan[k][j + 1],
                    fan[k][j],
                    uv[k + 1][j],
                    uv[k][j + 1],
                )
        self.region_II_uv = uv

    # ------------------------------------------------------------------
    # Region III — Goursat problem from fan radial 0 + outer α-line
    # ------------------------------------------------------------------

    def _solve_region_III(self):
        rows = self.solver.region_III_nodes  # [j][i], i = 0 … n − j
        n = len(rows) - 1

        uv: list[list[tuple[float, float]]] = [[None] * len(row) for row in rows]
        for j in range(n + 1):
            uv[j][0] = self.region_II_uv[0][j]  # corner β-line = fan radial 0
            uv[j][n - j] = (self.U, 0.0)  # outer α-line

        # α-neighbour of (j, i) is (j+1, i−1); β-neighbour is (j+1, i).
        for j in range(n - 1, -1, -1):
            for i in range(1, n - j):
                uv[j][i] = geiringer_step(
                    rows[j + 1][i - 1],
                    rows[j + 1][i],
                    rows[j][i],
                    uv[j + 1][i - 1],
                    uv[j + 1][i],
                )
        self.region_III_uv = uv

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def region_samples(self) -> dict[str, list[tuple[Node, float, float]]]:
        """(node, vx, vy) for every node, grouped by region ('I', 'II', 'III')."""
        groups = {
            "I": (self.solver.region_I_nodes, self.region_I_uv),
            "II": (self.solver.region_II_nodes, self.region_II_uv),
            "III": (self.solver.region_III_nodes, self.region_III_uv),
        }
        out: dict[str, list[tuple[Node, float, float]]] = {}
        for name, (nodes, uv) in groups.items():
            samples = []
            for node_row, uv_row in zip(nodes, uv):
                for node, (u, v) in zip(node_row, uv_row):
                    samples.append((node, *alpha_beta_to_xy(node.psi, u, v)))
            out[name] = samples
        return out

    def outer_alpha_line(self) -> list[Node]:
        """The velocity-discontinuity line from (0, 0) to the free surface."""
        s = self.solver
        N_I = len(s.region_I_nodes) - 1
        n = len(s.region_III_nodes) - 1
        line = [s.region_I_nodes[i][0] for i in range(N_I + 1)]
        line += [radial[-1] for radial in reversed(s.region_II_nodes)][1:]
        line += [s.region_III_nodes[j][n - j] for j in range(n - 1, -1, -1)]
        return line
