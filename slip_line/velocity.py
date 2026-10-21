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

import itertools

import numpy as np

from .slip_line_field import Node
from .solver import DiskPrandtlSolver, PrandtlSolver


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


def geiringer_step_contact_bc(
    Q: Node,
    R: Node,
    uv_Q: tuple[float, float],
    normal: tuple[float, float],
    v_n: float,
) -> tuple[float, float]:
    """
    Velocity at a contact node R from its β-neighbour Q and the
    smooth-contact condition V·n = v_n (normal velocity of the foundation):

        u_R·(e_α·n) + v_R·(e_β·n) = v_n
        v_R − v_Q = −½·(u_Q + u_R)·(ψ_R − ψ_Q)
    """
    u_Q, v_Q = uv_Q
    a_n = np.dot(alpha_beta_to_xy(R.psi, 1.0, 0.0), normal)
    b_n = np.dot(alpha_beta_to_xy(R.psi, 0.0, 1.0), normal)
    b = 0.5 * (R.psi - Q.psi)
    rhs_b = v_Q - b * u_Q
    # a_n·u_R + b_n·v_R = v_n  and  b·u_R + v_R = rhs_b
    u_R = (v_n - b_n * rhs_b) / (a_n - b_n * b)
    return u_R, rhs_b - b * u_R


def _march_fan(
    fan: list[list[Node]],
    uv_fan_I: list[tuple[float, float]],
    U: float,
) -> list[list[tuple[float, float]]]:
    """
    Goursat problem in a centred fan, fan[k][j] (radial k, radius index j).
    Data: radial k = −1 (Region I side) and the outer arc j = −1 (v = 0, u = U).
    α-neighbour of (k, j) is (k+1, j) on the same arc; β-neighbour is
    (k, j+1) on the same radial.
    """
    n_fan = len(fan) - 1
    n = len(fan[0]) - 1
    uv: list[list[tuple[float, float]]] = [[None] * (n + 1) for _ in range(n_fan + 1)]
    uv[n_fan] = list(uv_fan_I)
    for k in range(n_fan - 1, -1, -1):
        uv[k][n] = (U, 0.0)
        for j in range(n - 1, -1, -1):
            uv[k][j] = geiringer_step(
                fan[k + 1][j], fan[k][j + 1], fan[k][j], uv[k + 1][j], uv[k][j + 1]
            )
    return uv


def _march_region_III(
    rows: list[list[Node]],
    uv_corner: list[tuple[float, float]],
    U: float,
) -> list[list[tuple[float, float]]]:
    """
    Goursat problem in the passive wedge, rows[j][i] (i = 0 … n − j).
    Data: corner β-line i = 0 and the outer α-line i = n − j (v = 0, u = U).
    α-neighbour of (j, i) is (j+1, i−1); β-neighbour is (j+1, i).
    """
    n = len(rows) - 1
    uv: list[list[tuple[float, float]]] = [[None] * len(row) for row in rows]
    for j in range(n + 1):
        uv[j][0] = uv_corner[j]
        uv[j][n - j] = (U, 0.0)
    for j in range(n - 1, -1, -1):
        for i in range(1, n - j):
            uv[j][i] = geiringer_step(
                rows[j + 1][i - 1],
                rows[j + 1][i],
                rows[j][i],
                uv[j + 1][i - 1],
                uv[j + 1][i],
            )
    return uv


def _samples(nodes, uv) -> list[tuple[Node, float, float]]:
    """Flatten matching node / (u, v) layouts into (node, vx, vy)."""
    out = []
    for node_row, uv_row in zip(nodes, uv):
        for node, (u, v) in zip(node_row, uv_row):
            out.append((node, *alpha_beta_to_xy(node.psi, u, v)))
    return out


def _row_mesh_triangles(rows, uv) -> list:
    """Triangles of a mesh where node [k+1][i] lies between [k][i] and [k][i+1]."""
    tris = []
    for k in range(len(rows) - 1):
        for i in range(len(rows[k + 1])):
            a, b, c = (k, i), (k, i + 1), (k + 1, i)
            tris.append([(rows[r][q], uv[r][q]) for r, q in (a, b, c)])
            if i + 1 < len(rows[k + 1]):
                d = (k + 1, i + 1)
                tris.append([(rows[r][q], uv[r][q]) for r, q in (b, c, d)])
    return tris


def _fan_triangles(fan, uv) -> list:
    """Two triangles per fan cell (k, j)–(k+1, j+1)."""
    tris = []
    for k in range(len(fan) - 1):
        for j in range(len(fan[k]) - 1):
            a, b, c, d = (k, j), (k + 1, j), (k + 1, j + 1), (k, j + 1)
            for tri in ((a, b, c), (a, c, d)):
                tris.append([(fan[r][q], uv[r][q]) for r, q in tri])
    return tris


def plastic_work_rate(triangles, cu: float) -> dict[str, np.ndarray]:
    """
    Plastic work rate per unit area in each triangle (velocity linear over it).

    With ε̇_αα = ε̇_ββ = 0 the strain rate is pure shear γ̇_αβ, and the work
    of the slip-line stress (compression positive, shear +cu on α-lines) is

        d = −cu·γ̇_αβ,   γ̇_αβ = e_α·(L + Lᵀ)·e_β,   L = ∇V

    d < 0 means the velocity field does negative work against the stress
    field there, i.e. the two are not compatible (no complete solution).
    Degenerate (zero-area) triangles are skipped.

    Returns 'xy' (n, 3, 2) vertices, 'area' (n,), 'd' (n,).
    """
    xy, area, d = [], [], []
    for tri in triangles:
        pts = np.array([(nd.x, nd.y) for nd, _ in tri])
        vel = np.array([alpha_beta_to_xy(nd.psi, u, v) for nd, (u, v) in tri])
        J = pts[1:] - pts[0]
        det = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
        if abs(det) < 1e-14:
            continue
        L = np.linalg.solve(J, vel[1:] - vel[0]).T  # L[i, j] = ∂v_i/∂x_j
        theta = np.mean([nd.psi for nd, _ in tri]) - np.pi / 4.0
        e_a = np.array([np.cos(theta), np.sin(theta)])
        e_b = np.array([-np.sin(theta), np.cos(theta)])
        gamma = e_a @ (L + L.T) @ e_b
        xy.append(pts)
        area.append(0.5 * abs(det))
        d.append(-cu * gamma)
    return {"xy": np.array(xy), "area": np.array(area), "d": np.array(d)}


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
            uv[j][j] = geiringer_step_contact_bc(
                nodes[j][j - 1], nodes[j][j], uv[j][j - 1], (0.0, 1.0), self.v0
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
        # Radial k = n_fan coincides with Region I row N_I (reversed).
        row = self.region_I_uv[-1]
        self.region_II_uv = _march_fan(self.solver.region_II_nodes, row[::-1], self.U)

    # ------------------------------------------------------------------
    # Region III — Goursat problem from fan radial 0 + outer α-line
    # ------------------------------------------------------------------

    def _solve_region_III(self):
        # Corner β-line of Region III = fan radial 0.
        self.region_III_uv = _march_region_III(
            self.solver.region_III_nodes, self.region_II_uv[0], self.U
        )

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
        return {name: _samples(nodes, uv) for name, (nodes, uv) in groups.items()}

    def outer_alpha_line(self) -> list[Node]:
        """The velocity-discontinuity line from (0, 0) to the free surface."""
        s = self.solver
        N_I = len(s.region_I_nodes) - 1
        n = len(s.region_III_nodes) - 1
        line = [s.region_I_nodes[i][0] for i in range(N_I + 1)]
        line += [radial[-1] for radial in reversed(s.region_II_nodes)][1:]
        line += [s.region_III_nodes[j][n - j] for j in range(n - 1, -1, -1)]
        return line


class DiskVelocityField:
    """
    Velocity field for a DiskPrandtlSolver net.  The disk moves with speed
    v0 either vertically into the soil (y-down, so v0 > 0 is downward; the
    left half mirrors the right) or, with horizontal=True, horizontally
    toward +x.  A horizontal push only loads the right half: the left half
    of the disk moves away from the soil, which separates and stays rigid.

    Boundary conditions (right half, x ≥ 0):
      • Outer α-line (Region I column 0 from the disk bottom → fan outer arc
        → Region III boundary) borders soil at rest: v = 0, u = U.
      • Disk contact: V·n = V0·n with V0 the disk velocity and n the outward
        disk normal; tangential slip is free (smooth disk, r = 0).  At the
        disk bottom this fixes U (U = 0 for a horizontal push).
      • Free surface: no kinematic condition.

    Region I's corner β-line is straight and coincides with the fan-I radial,
    but the two meshes have different nodes, so velocities are interpolated
    linearly in distance from the corner.

    Velocities are stored as (u, v) per node, in the same layout as the
    solver's region_I_mesh / region_II_nodes / region_III_nodes.
    """

    def __init__(
        self, solver: DiskPrandtlSolver, v0: float = 1.0, horizontal: bool = False
    ):
        if not solver._built:
            solver.solve()
        self.solver = solver
        self.v0 = float(v0)
        self.horizontal = horizontal
        # Unit direction of the disk velocity V0 (y-down).
        self.direction = (1.0, 0.0) if horizontal else (0.0, 1.0)
        self.U = 0.0
        self.region_I_uv: list[list[tuple[float, float]]] = []
        self.region_II_uv: list[list[tuple[float, float]]] = []
        self.region_III_uv: list[list[tuple[float, float]]] = []

    def solve(self) -> float:
        """Build the velocity field; return U, the jump across the outer α-line."""
        self._solve_region_I()
        uv_corner = self._corner_beta_line_uv()
        fan = self.solver.region_II_nodes
        if len(fan) >= 2:
            self.region_II_uv = _march_fan(fan, uv_corner, self.U)
            uv_corner = self.region_II_uv[0]
        else:
            # Degenerate fan (one radial) or parallel-β zone (no fan).
            self.region_II_uv = [uv_corner] if fan else []
        self.region_III_uv = _march_region_III(
            self.solver.region_III_nodes, uv_corner, self.U
        )
        return self.U

    def _normal(self, node: Node) -> tuple[float, float]:
        """Outward disk normal (disk → soil) at a disk-surface node."""
        s = self.solver
        return node.x / s.R, (node.y - s.y_c) / s.R

    # ------------------------------------------------------------------
    # Region I — mixed problem: outer α-line (column 0) + disk surface
    # ------------------------------------------------------------------

    def _solve_region_I(self):
        """
        mesh[k][i]: α-line = column i (from disk node mesh[0][i]),
        β-line = i + k const.  α-neighbour of (k, i) is (k−1, i), toward
        the disk; β-neighbour is (k+1, i−1), toward column 0.
        """
        mesh = self.solver.region_I_mesh
        N = len(mesh[0]) - 1

        # U from the disk bottom: v = 0 and V·n = V0·n.
        bottom = mesh[0][0]
        n_b = self._normal(bottom)
        self.U = (
            self.v0
            * np.dot(self.direction, n_b)
            / np.dot(alpha_beta_to_xy(bottom.psi, 1.0, 0.0), n_b)
        )

        uv: list[list[tuple[float, float]]] = [[None] * len(row) for row in mesh]
        for k in range(len(mesh)):
            uv[k][0] = (self.U, 0.0)  # column 0 = outer α-line

        for i in range(1, N + 1):
            n_i = self._normal(mesh[0][i])
            uv[0][i] = geiringer_step_contact_bc(
                mesh[1][i - 1],
                mesh[0][i],
                uv[1][i - 1],
                n_i,
                self.v0 * np.dot(self.direction, n_i),
            )
            for k in range(1, N + 1 - i):
                uv[k][i] = geiringer_step(
                    mesh[k - 1][i],
                    mesh[k + 1][i - 1],
                    mesh[k][i],
                    uv[k - 1][i],
                    uv[k + 1][i - 1],
                )
        self.region_I_uv = uv

    def _corner_beta_line_uv(self) -> list[tuple[float, float]]:
        """(u, v) at r_j = j·dr along the corner β-line, from Region I."""
        mesh = self.solver.region_I_mesh
        N = len(mesh[0]) - 1
        corner = mesh[0][N]
        line = [(mesh[k][N - k], self.region_I_uv[k][N - k]) for k in range(N + 1)]
        dist = [np.hypot(nd.x - corner.x, nd.y - corner.y) for nd, _ in line]
        u = [uv[0] for _, uv in line]
        v = [uv[1] for _, uv in line]
        r = np.arange(self.solver.n_radial + 1) * self.solver.dr
        return list(zip(np.interp(r, dist, u), np.interp(r, dist, v)))

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def region_samples(self) -> dict[str, list[tuple[Node, float, float]]]:
        """(node, vx, vy) for every node, grouped by region ('I', 'II', 'III')."""
        s = self.solver
        return {
            "I": _samples(s.region_I_mesh, self.region_I_uv),
            "II": _samples(s.region_II_nodes, self.region_II_uv),
            "III": _samples(s.region_III_nodes, self.region_III_uv),
        }

    def outer_alpha_line(self) -> list[Node]:
        """The velocity-discontinuity line from the disk bottom to the free surface."""
        s = self.solver
        n = len(s.region_III_nodes) - 1
        line = [row[0] for row in s.region_I_mesh]
        line += [radial[-1] for radial in reversed(s.region_II_nodes)][1:]
        line += [s.region_III_nodes[j][n - j] for j in range(n - 1, -1, -1)]
        return line

    def plastic_work(self) -> dict[str, np.ndarray]:
        """plastic_work_rate() over all triangles of Regions I, II and III."""
        s = self.solver
        tris = _row_mesh_triangles(s.region_I_mesh, self.region_I_uv)
        if len(s.region_II_nodes) >= 2:
            tris += _fan_triangles(s.region_II_nodes, self.region_II_uv)
        tris += _row_mesh_triangles(s.region_III_nodes, self.region_III_uv)
        return plastic_work_rate(tris, s.cu)

    def energy_balance(self) -> dict[str, float]:
        """
        Full-disk load F per unit length along the push direction (vertical Q,
        or horizontal H), from the work equation  F·v0 = halves·(P + cu·U·L)
        with halves = 2 for a vertical push and 1 for a horizontal push:

          F_stress        — from the solver's contact stresses
          F_virtual_work  — P = Σ d·A (signed); equals F_stress when the
                            velocity field is correct (virtual work identity)
          F_upper_bound   — P = Σ |d|·A; the upper bound of this mechanism
          negative_area_fraction — share of the plastic area with d < 0
        """
        pw = self.plastic_work()
        outer = self.outer_alpha_line()
        L = sum(np.hypot(b.x - a.x, b.y - a.y) for a, b in itertools.pairwise(outer))
        P_jump = self.solver.cu * self.U * L
        P_signed = float(np.sum(pw["d"] * pw["area"]))
        P_abs = float(np.sum(np.abs(pw["d"]) * pw["area"]))
        if self.horizontal:
            halves = 1.0
            F_stress = -self.solver.forces_on_half_disk()["F_h"]
        else:
            halves = 2.0
            F_stress = self.solver.collapse_load_Q()
        return {
            "F_stress": F_stress,
            "F_virtual_work": halves * (P_signed + P_jump) / self.v0,
            "F_upper_bound": halves * (P_abs + P_jump) / self.v0,
            "negative_area_fraction": float(
                pw["area"][pw["d"] < 0].sum() / pw["area"].sum()
            ),
        }
