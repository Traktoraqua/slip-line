"""
Equations of motion of the pipe on a fixed yield envelope (step 2 of
PLAN_DYNAMIKK.md; manual.tex, "Equations of motion").

The pipe is a rigid body with the mass M = diag(M_x, M_z) and two degrees
of freedom, u (toward +x) and w (down).  The soil is rigid-plastic: its
force on the pipe is −F, with F = (H, V) on the yield envelope E
(envelope.py), and the pipe moves only when F is on the boundary, along the
outward normal there (normality):

    M·v̇ = F_ext(t) − F,   F ∈ E,   v in the normal cone of E at F,

with F_ext = (Fx(t), W − F_L(t)).  Backward Euler turns each step into a
closest-point projection: with the trial force F* = F_ext(t_k+1) + M·v_k/Δt
the soil force F_k+1 is the point of E closest to F* in the metric of M⁻¹,
and v_k+1 = Δt·M⁻¹·(F* − F_k+1).  F* inside E means rest.

Units: forces in kN/m, masses in kg/m (t/m inside, so that kN = t·m/s²),
time in s, lengths in m, energies in kJ/m.  The envelope (the geometry) is
fixed here; it follows the motion from step 3 of the plan.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise

import numpy as np

from .envelope import YieldEnvelope

RHO_W = 1025.0  # sea water, kg/m³


def pipe_mass(W, D, C_a=1.0, rho_w=RHO_W, g=9.81) -> tuple[float, float]:
    """Masses (M_x, M_z) in kg/m of a pipe with the submerged weight W
    (kN/m) and outer diameter D (m): the pipe with contents and coating,
    m = W/g + ρ_w·πD²/4 (its weight in water plus the buoyancy), and the
    added mass C_a·ρ_w·πD²/4 of the water.  C_a is one value or
    (C_a,x, C_a,z)."""
    water = rho_w * np.pi * D**2 / 4.0
    m = 1e3 * W / g + water
    ca_x, ca_z = np.broadcast_to(np.asarray(C_a, dtype=float), 2)
    return float(m + ca_x * water), float(m + ca_z * water)


def closest_point(env: YieldEnvelope, F, M) -> np.ndarray:
    """Point of the envelope closest to F = (H, V) in the metric of M⁻¹,
    M = (M_x, M_z): F itself if it lies inside.  Scaling the forces by
    M^(−1/2) makes the distance Euclidean and keeps the polygon convex; the
    point is then the nearest of the closest points on all edges."""
    F = np.asarray(F, dtype=float)
    if env.contains(F):
        return F.copy()
    s = 1.0 / np.sqrt(np.asarray(M, dtype=float))
    P, f = env.vertices * s, F * s
    E = np.roll(P, -1, axis=0) - P
    t = np.clip(np.sum((f - P) * E, axis=1) / np.sum(E * E, axis=1), 0.0, 1.0)
    X = P + t[:, None] * E
    return X[np.argmin(np.sum((X - f) ** 2, axis=1))] / s


@dataclass
class DynamicState:
    """The pipe at time t: displacement, velocity and the force (H, V) from
    the pipe on the soil."""

    t: float  # time, s
    u: float  # horizontal displacement of the pipe centre (toward +x), m
    w: float  # vertical displacement of the pipe centre (down), m
    u_dot: float  # m/s
    w_dot: float  # m/s
    H: float  # kN/m, toward +x
    V: float  # kN/m, down

    @property
    def delta(self) -> float:
        """Direction of motion below the horizontal (rad); nan at rest."""
        if self.u_dot == 0.0 and self.w_dot == 0.0:
            return float("nan")
        return float(np.arctan2(self.w_dot, self.u_dot))


class PipeDynamics:
    """Rigid pipe with mass on rigid-plastic soil whose yield envelope env
    is fixed, loaded by its submerged weight W, the horizontal load Fx(t)
    and the lift F_L(t) (kN/m; functions of time, zero if not given).
    M = (M_x, M_z) in kg/m, e.g. pipe_mass(W, D)."""

    def __init__(self, env: YieldEnvelope, W: float, M, Fx=None, F_L=None):
        self.env = env
        self.W = W
        self.M = np.array(M, dtype=float)
        self.Fx = Fx if Fx is not None else (lambda t: 0.0)
        self.F_L = F_L if F_L is not None else (lambda t: 0.0)
        self._Mt = 1e-3 * self.M  # t/m

    def load(self, t: float) -> np.ndarray:
        """F_ext = (Fx, W − F_L) on the pipe at time t (kN/m, V down)."""
        return np.array([self.Fx(t), self.W - self.F_L(t)])

    def rest(self, t: float = 0.0) -> DynamicState:
        """The pipe at rest at u = w = 0, with the soil force closest to the
        load at t."""
        H, V = closest_point(self.env, self.load(t), self.M)
        return DynamicState(t, 0.0, 0.0, 0.0, 0.0, H, V)

    def step(self, s: DynamicState, dt: float) -> DynamicState:
        """One backward Euler step of dt from s: the soil force is the point
        of the envelope closest to F* = F_ext(t + dt) + M·v/dt (metric
        M⁻¹), and the new velocity v = dt·M⁻¹·(F* − F)."""
        t = s.t + dt
        v = np.array([s.u_dot, s.w_dot])
        F_star = self.load(t) + self._Mt * v / dt
        F = closest_point(self.env, F_star, self.M)
        v = dt * (F_star - F) / self._Mt
        return DynamicState(t, s.u + dt * v[0], s.w + dt * v[1], v[0], v[1], F[0], F[1])

    def run(self, t_end: float, dt: float, start=None) -> list[DynamicState]:
        """States from start (default: at rest at t = 0) to t_end, in steps
        of dt (the times are start.t + k·dt, without round-off drift)."""
        s = start if start is not None else self.rest()
        t0, n = s.t, round((t_end - s.t) / dt)
        states = [s]
        for k in range(1, n + 1):
            states.append(self.step(states[-1], t0 + k * dt - states[-1].t))
        return states

    def energy(self, states: list[DynamicState]) -> dict:
        """Energy balance over consecutive states of run() (kJ/m), as in
        manual.tex: work of the loads (with the weight) = dissipation +
        change of kinetic energy + numerical damping of backward Euler,
        exactly up to round-off."""
        work = dissipation = numerical = 0.0
        for a, b in pairwise(states):
            dt = b.t - a.t
            v0, v1 = np.array([a.u_dot, a.w_dot]), np.array([b.u_dot, b.w_dot])
            work += dt * float(self.load(b.t) @ v1)
            dissipation += dt * float(np.array([b.H, b.V]) @ v1)
            numerical += 0.5 * float((v1 - v0) @ (self._Mt * (v1 - v0)))
        v0 = np.array([states[0].u_dot, states[0].w_dot])
        v1 = np.array([states[-1].u_dot, states[-1].w_dot])
        kinetic = 0.5 * float(v1 @ (self._Mt * v1) - v0 @ (self._Mt * v0))
        return {
            "work": work,
            "dissipation": dissipation,
            "kinetic": kinetic,
            "numerical": numerical,
        }
