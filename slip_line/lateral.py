"""
Lateral displacement of a partly embedded pipe that is free to move
vertically: the vertical load W stays constant while the pipe is pushed toward
+x in small steps du.  Weightless Tresca soil (cu), interface roughness r.  The
soil surface (berm in front, trench behind) is tracked step by step; the start
is the end of a PipePenetration run.

Mechanism (kinematic, upper bound).  The pipe moves with unit speed along
(cos β, sin β) (y down, so β > 0 is downward).  Its trailing (left) side
separates where it moves away from the soil, beyond φ = −β from the invert
(φ measured from the downward axis, positive toward +x), so the soil touches
the pipe on φ ∈ [φ_start, θ_R]:

    φ_start = −min(β, φ_L)    (β ≥ 0)
    φ_start = −β > 0          (β < 0: the pipe rises and leaves the invert)

φ_L is the largest trailing contact angle.  Against the original surface it
is the current left contact angle (diving steeper than that cuts into the
soil, with a velocity jump U along the outer α-line).  Against the trench
wall the pipe has cut itself there is no limit: on a smooth path the contact
ends exactly where the pipe leaves the wall, at φ = −β.

The net is that of DiskPrandtlSolver on this chain, one-sided (all soil flows
toward +x through the fan and passive wedge at the right edge), with a
uniform net angle ω: ψ = π/2 − φ − ω.  ω only shapes the mechanism; the
interface shear enters the dissipation as r·cu·|slip|.  The free surface in
front is a straight line from the right contact point, falling away from the
pipe at the angle α of a least-squares fit to the berm over the passive
wedge's extent (α = 0 is the flat surface of DiskPrandtlSolver); σ₁ is
parallel to it there, ψ = α.  The Geiringer field on the net gives, per
unit pipe speed,

    D = Σ|d|·A + cu·U·L + r·cu·∫|slip| ds

(plastic shear, the jump U along the outer α-line, interface slip), and the
work equation at constant V = W gives the horizontal force

    H·cos β + W·sin β = D   →   H = (D − W·sin β) / cos β.

The pipe moves in the direction β (and with the net ω) that minimises H, the
normality rule of the upper-bound yield envelope.  Each step:

  1. Mechanism: β, ω minimising H at the current geometry.
  2. Move: du sideways and du·tan β down (path length ds = du / cos β).
  3. Surface: the free surface in front rises by the Geiringer velocity
     through the sloping surface × ds (v_up + v_x·tan α, per unit width), and
     the soil the pipe moves into is removed; the rise is scaled so the soil
     volume is conserved.  α is then refitted to the new surface.

The left soil (trench wall, or the original surface) is rigid.  x is global,
with the pipe centre at u; the surface elevation η is positive upward.

A push back toward −x (PipeLateral.reversed) is computed in the mirror image
x → −x, where it is again a push toward +x: the trench is then in front (its
rising wall gives α < 0) and the berm behind.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np

from .boundary import DiskFoundationBoundary
from .penetration import PenetrationState, PipePenetration
from .slip_line_field import Node
from .solver import DiskPrandtlSolver
from .velocity import DiskVelocityField, alpha_beta_to_xy

OMEGA_MAX = np.pi / 4.0 - 0.02  # |ω| bound: the chain must not be a slip line


class OneSidedDiskSolver(DiskPrandtlSolver):
    """
    DiskPrandtlSolver net on the contact chain φ ∈ [phi_start, θ_max]
    (phi_start < 0 lies left of the invert) with a uniform net angle ω:
    ψ = π/2 − φ − ω, σ = cu·(1 + 2(ψ − α)).  The fan and passive wedge sit at
    the right edge only, the wedge under a free surface falling at angle α
    (ψ = α); the outer α-line starts at the chain's first node.  Only the net
    is meant to be used (the contact-stress methods of the parent cover
    φ ∈ [0, θ_max] under a flat surface).
    """

    def __init__(
        self,
        disk: DiskFoundationBoundary,
        phi_start: float,
        omega: float,
        alpha: float = 0.0,
        n_fan: int = 18,
        n_disk: int = 30,
    ):
        super().__init__(disk, n_fan=n_fan, n_disk=n_disk)
        self.phi_start = float(phi_start)
        self.omega_signed = float(omega)
        self.alpha = self.psi_B = float(alpha)
        self.use_parallel_zone = False
        self.psi_edge = 0.5 * np.pi - self.theta_max - self.omega_signed
        self.delta_psi_fan = self.psi_edge - self.psi_B
        self.psi_bot = 0.5 * np.pi - self.phi_start - self.omega_signed
        if self.delta_psi_fan < 0.0:
            raise ValueError("ω too large: the corner fan would be negative.")

    def solve(self) -> None:
        """Build the net (without the parent's collapse-load summary)."""
        self._build_disk_surface()
        self._build_region_I()
        self._resize_fan_and_wedge_to_match_region_I()
        self._build_region_II()
        self._build_region_III()
        self._built = True

    def _build_disk_surface(self):
        nodes = []
        for th in np.linspace(self.theta_max, self.phi_start, self.n_disk + 1):
            x, y = self.disk.position(th)
            psi = self._psi_disk(th)
            sigma = self.cu * (1.0 + 2.0 * (psi - self.alpha))
            nodes.append(Node(x=x, y=y, sigma=sigma, psi=psi))
        self.disk_surface_nodes = nodes

    def _build_region_III(self):
        """Passive wedge under the sloping surface: the flat-surface layout
        rotated by α about the contact edge."""
        n, h = self.n_radial, self.h
        c, s = np.cos(self.alpha), np.sin(self.alpha)
        rows = []
        for j in range(n + 1):
            row = []
            for i in range(n + 1 - j):
                a, b = (i + 0.5 * j) * h, 0.5 * j * h
                row.append(
                    Node(
                        x=self.B_s + a * c - b * s,
                        y=a * s + b * c,
                        sigma=self.cu,
                        psi=self.alpha,
                    )
                )
            rows.append(row)
        self.region_III_nodes = rows

    def complete(self) -> bool:
        """True when Region I is the full triangular net (no node cut off)."""
        mesh = self.region_I_mesh
        return len(mesh) == self.n_disk + 1 and all(
            len(row) == self.n_disk + 1 - k for k, row in enumerate(mesh)
        )


class InclinedDiskVelocityField(DiskVelocityField):
    """Geiringer field on a OneSidedDiskSolver net for a pipe moving with
    speed v0 along (cos β, sin β)."""

    def __init__(self, solver: OneSidedDiskSolver, beta: float, v0: float = 1.0):
        super().__init__(solver, v0=v0)
        self.beta = float(beta)
        self.direction = (np.cos(beta), np.sin(beta))

    def dissipation(self) -> dict[str, float]:
        """Upper-bound dissipation rate (per unit length) and its parts."""
        s = self.solver
        pw = self.plastic_work()
        outer = self.outer_alpha_line()
        L = sum(np.hypot(b.x - a.x, b.y - a.y) for a, b in itertools.pairwise(outer))
        it = self.interface_slip()
        plastic = float(np.sum(np.abs(pw["d"]) * pw["area"]))
        jump = s.cu * abs(self.U) * L
        interface = (
            s.R * s.r * s.cu * float(np.trapezoid(np.abs(it["slip"]), it["phi"]))
        )
        return {
            "plastic": plastic,
            "jump": jump,
            "interface": interface,
            "total": plastic + jump + interface,
        }

    def surface_velocity(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Local x, horizontal and upward velocity of the free-surface nodes
        (Region III)."""
        row = self.solver.region_III_nodes[0]
        v = np.array(
            [
                alpha_beta_to_xy(nd.psi, u, v)
                for nd, (u, v) in zip(row, self.region_III_uv[0])
            ]
        )
        return np.array([nd.x for nd in row]), v[:, 0], -v[:, 1]


@dataclass
class Mechanism:
    """Best mechanism at a geometry."""

    beta: float  # direction of motion below the horizontal (rad)
    omega: float  # net angle (rad)
    H: float  # horizontal force, full pipe per unit length
    work: dict[str, float]  # dissipation parts per unit pipe speed
    phi_start: float  # first contact angle of the chain (rad)
    U: float  # velocity jump along the outer α-line per unit pipe speed


@dataclass
class LateralState:
    """Pipe, soil surface and the mechanism found at that geometry."""

    u: float  # horizontal displacement of the pipe centre (toward +x)
    w: float  # invert depth below the original surface
    eta: np.ndarray  # soil surface elevation on PipeLateral.x
    x_left: float  # left end of the contact (global x)
    x_right: float  # right end of the contact (global x)
    cut: np.ndarray  # surface made by the pipe (excavated), on PipeLateral.x
    mech: Mechanism | None = field(default=None, repr=False)
    scale: float = 1.0  # volume correction applied to the uplift that led here
    alpha: float = 0.0  # slope of the free surface in front (falling away, rad)

    def d(self, R: float) -> float:
        """Effective embedment: invert depth below the right contact level."""
        dx = self.x_right - self.u
        return R - np.sqrt(R * R - dx * dx)

    def phi_left(self, R: float) -> float:
        """Left contact angle from the invert."""
        return float(np.arcsin(np.clip((self.u - self.x_left) / R, 0.0, 1.0)))


class PipeLateral:
    """Pipe of diameter D, roughness r and weight W pushed sideways in Tresca soil."""

    def __init__(
        self,
        D: float,
        cu: float,
        r: float,
        W: float,
        x: np.ndarray,
        du: float = 1e-3,
        n_fan: int = 18,
        n_disk: int = 30,
    ):
        self.D, self.R, self.cu, self.r, self.W, self.du = D, D / 2.0, cu, r, W, du
        self.x = x
        self.n_fan, self.n_disk = n_fan, n_disk

    @classmethod
    def after_penetration(
        cls,
        sim: PipePenetration,
        state: PenetrationState,
        x_min: float = -1.5,
        x_max: float = 2.5,
        **kw,
    ) -> tuple[PipeLateral, LateralState]:
        """Lateral model and start state from the end of a penetration run
        (the right-half heave mirrored to the left)."""
        dx = sim.x[1] - sim.x[0]
        x = np.arange(x_min, x_max + 0.5 * dx, dx)
        lat = cls(sim.D, sim.cu, sim.r, sim.W, x, **kw)
        heave = np.interp(np.abs(x), sim.x, state.heave, right=0.0)
        eta = lat._excavate(heave, 0.0, state.w)
        start = LateralState(0.0, state.w, eta, -state.x_c, state.x_c, heave > eta)
        start.alpha = lat._slope(eta, state.x_c, state.x_end - state.x_c)
        return lat, start

    def reversed(self, state: LateralState) -> tuple[PipeLateral, LateralState]:
        """Model and start state (with its mechanism) for pushing the pipe
        back toward −x from state, in the mirror image x → −x, where the push
        is again toward +x.  α of the new front (the trench) is fitted over
        the passive wedge of the mechanism found there."""
        lat = PipeLateral(
            self.D,
            self.cu,
            self.r,
            self.W,
            -self.x[::-1],
            self.du,
            self.n_fan,
            self.n_disk,
        )
        start = LateralState(
            -state.u,
            state.w,
            state.eta[::-1].copy(),
            -state.x_right,
            -state.x_left,
            state.cut[::-1].copy(),
        )
        length, guess = self.R, None
        for _ in range(10):
            start.alpha = lat._slope(start.eta, start.x_right, length)
            start.mech = lat.mechanism(start, guess)
            guess = (start.mech.beta, start.mech.omega)
            wedge = lat.field(start, *guess).solver.L
            if abs(wedge - length) < 1e-3:
                break
            length = wedge
        return lat, start

    # ------------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------------

    def bottom(self, x, u: float, w: float):
        """Depth of the pipe's lower surface at x (−inf outside the pipe)."""
        dx = np.asarray(x, dtype=float) - u
        inside = np.abs(dx) < self.R
        root = np.sqrt(np.where(inside, self.R * self.R - dx * dx, 0.0))
        return np.where(inside, w - self.R + root, -np.inf)

    def _excavate(self, eta: np.ndarray, u: float, w: float) -> np.ndarray:
        """Remove the soil above the pipe's lower surface."""
        return np.minimum(eta, -self.bottom(self.x, u, w))

    def _slope(self, eta: np.ndarray, x0: float, length: float) -> float:
        """Angle α of the least-squares line through the surface on
        [x0, x0 + length], positive when it falls away from the pipe."""
        on = (self.x >= x0) & (self.x <= x0 + length)
        return float(-np.arctan(np.polyfit(self.x[on], eta[on], 1)[0]))

    def _contact_end(self, eta: np.ndarray, u: float, w: float, side: int) -> float:
        """Outer end of the contact on one side (side = ±1): the outermost
        crossing of the soil surface with the pipe's lower surface, the soil
        lying above it inside the contact.  Returns u if there is none."""

        def gap(x):
            return float(self.bottom(x, u, w)) + float(np.interp(x, self.x, eta))

        out = side * (self.x - u)
        idx = np.nonzero((out > 0.0) & (out < self.R))[0]
        idx = idx[::-1] if side > 0 else idx  # outermost first
        hit = np.nonzero(self.bottom(self.x[idx], u, w) + eta[idx] >= 0.0)[0]
        if len(hit) == 0:
            return u
        if hit[0] == 0:
            raise ValueError("Soil above the pipe's springline.")
        lo, hi = self.x[idx[hit[0]]], self.x[idx[hit[0] - 1]]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if gap(mid) >= 0.0 else (lo, mid)
        return 0.5 * (lo + hi)

    def left_limit(self, state: LateralState) -> float:
        """Largest trailing contact angle φ_L: none (π/2) on the trench wall
        the pipe has cut, else the current left contact angle."""
        if state.x_left >= state.u:
            return 0.0
        outside = max(int(np.searchsorted(self.x, state.x_left)) - 1, 0)
        return 0.5 * np.pi if state.cut[outside] else state.phi_left(self.R)

    # ------------------------------------------------------------------
    # Mechanism
    # ------------------------------------------------------------------

    def field(
        self, state: LateralState, beta: float, omega: float
    ) -> InclinedDiskVelocityField:
        """Net and Geiringer field for direction β and net angle ω."""
        R = self.R
        phi_start = -min(beta, self.left_limit(state)) if beta >= 0.0 else -beta
        disk = DiskFoundationBoundary(
            cu=self.cu, D=self.D, p=state.d(R) / self.D, r=self.r
        )
        s = OneSidedDiskSolver(
            disk, phi_start, omega, state.alpha, self.n_fan, self.n_disk
        )
        if phi_start >= s.theta_max:
            raise ValueError("No contact left.")
        s.solve()
        if not s.complete():
            raise ValueError("Incomplete net.")
        vf = InclinedDiskVelocityField(s, beta)
        vf.solve()
        return vf

    def horizontal_force(self, state: LateralState, beta: float, omega: float):
        """H for direction β and net angle ω (inf if the net is not valid)."""
        try:
            vf = self.field(state, beta, omega)
        except (ValueError, IndexError, FloatingPointError):
            return np.inf, None
        D = vf.dissipation()
        return (D["total"] - self.W * np.sin(beta)) / np.cos(beta), vf

    def _omega_max(self, state: LateralState) -> float:
        theta = np.arccos(1.0 - 2.0 * state.d(self.R) / self.D)
        return min(OMEGA_MAX, 0.5 * np.pi - theta - state.alpha - 1e-6)

    def _best_omega(self, state, force, beta, omega0, span):
        lo = max(-OMEGA_MAX, omega0 - span)
        hi = min(self._omega_max(state), omega0 + span)
        omega = _golden(lambda om: force(beta, om), lo, hi)
        return omega, force(beta, omega)

    def mechanism(
        self,
        state: LateralState,
        guess: tuple[float, float] | None = None,
    ) -> Mechanism:
        """Mechanism minimising H.  Without a guess (β, ω) the directions are
        scanned first; the local search is recentred until the minimum lies
        inside its bracket.  If the search ends on an invalid net (the
        minimum lies on the edge of the valid nets), the best valid
        mechanism it has seen is taken, or the scan is run when it has seen
        none."""
        best = [np.inf, 0.0, 0.0]  # lowest H seen, with its β and ω

        def force(b, om):
            H = self.horizontal_force(state, b, om)[0]
            if H < best[0]:
                best[:] = [H, b, om]
            return H

        if guess is None:
            for beta in np.radians(np.arange(-20.0, 86.0, 5.0)):
                for om in np.linspace(-OMEGA_MAX, self._omega_max(state), 12):
                    force(beta, om)
            beta, omega = best[1], best[2]
            span_b, span_o = np.radians(5.0), 0.15
        else:
            (beta, omega), span_b, span_o = guess, np.radians(1.5), 0.04
        for _ in range(20):
            lo, hi = beta - span_b, beta + span_b
            om0 = omega
            new_beta = _golden(
                lambda b, om0=om0: self._best_omega(state, force, b, om0, span_o)[1],
                lo,
                hi,
                tol=np.radians(0.01),
            )
            new_omega, _ = self._best_omega(state, force, new_beta, omega, span_o)
            moved = min(new_beta - lo, hi - new_beta) < 0.05 * span_b or (
                abs(new_omega - omega) > 0.95 * span_o
                and -OMEGA_MAX < new_omega < self._omega_max(state)
            )
            beta, omega = new_beta, new_omega
            if not moved:
                break
        H, vf = self.horizontal_force(state, beta, omega)
        if vf is None:
            if not np.isfinite(best[0]):
                if guess is None:
                    raise ValueError("No valid mechanism.")
                return self.mechanism(state)
            beta, omega = best[1], best[2]
            H, vf = self.horizontal_force(state, beta, omega)
        return Mechanism(beta, omega, H, vf.dissipation(), vf.solver.phi_start, vf.U)

    # ------------------------------------------------------------------
    # Steps
    # ------------------------------------------------------------------

    def advance(self, state: LateralState, du: float) -> LateralState:
        """Move the pipe du sideways along the state's mechanism and update
        the soil surface."""
        m = state.mech
        vf = self.field(state, m.beta, m.omega)
        ds = du / np.cos(m.beta)
        u, w = state.u + du, state.w + du * np.tan(m.beta)

        xs, vx, up = vf.surface_velocity()
        flux = up + vx * np.tan(state.alpha)
        order = np.argsort(xs)
        inc = ds * np.interp(
            self.x - state.u, xs[order], flux[order], left=0.0, right=0.0
        )
        if self.x[-1] - state.u <= xs.max():
            raise ValueError("Plastic zone reaches x_max; increase x_max.")
        cap = -self.bottom(self.x, u, w)
        base = np.trapezoid(state.eta, self.x)
        scale = 1.0
        for _ in range(30):
            err = np.trapezoid(np.minimum(state.eta + scale * inc, cap), self.x) - base
            if abs(err) < 1e-14:
                break
            free = (state.eta + scale * inc) < cap
            scale -= err / np.trapezoid(np.where(free, inc, 0.0), self.x)
        heaved = state.eta + scale * inc
        x_left = self._contact_end(heaved, u, w, -1)
        x_right = self._contact_end(heaved, u, w, +1)
        eta = np.minimum(heaved, cap)
        cut = state.cut | (heaved > cap)
        alpha = self._slope(eta, x_right, vf.solver.L)
        return LateralState(u, w, eta, x_left, x_right, cut, scale=scale, alpha=alpha)

    def run(self, state: LateralState, distance: float) -> list[LateralState]:
        """Push the pipe `distance` sideways from state; return every state
        (each with its mechanism)."""
        if state.mech is None:
            state.mech = self.mechanism(state)
        states = [state]
        n = round(distance / self.du)
        for _ in range(n):
            prev = states[-1]
            nxt = self.advance(prev, self.du)
            nxt.mech = self.mechanism(nxt, guess=(prev.mech.beta, prev.mech.omega))
            states.append(nxt)
        return states


def _golden(f, a: float, b: float, tol: float = 2e-4) -> float:
    """Minimum of a unimodal f on [a, b] (golden-section search)."""
    g = 0.5 * (np.sqrt(5.0) - 1.0)
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    while b - a > tol:
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
    return 0.5 * (a + b)
