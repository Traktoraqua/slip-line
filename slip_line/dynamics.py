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
    F, M = np.asarray(F, dtype=float), np.asarray(M, dtype=float)
    if np.any(M <= 0.0):
        raise ValueError("The masses must be positive.")
    if env.contains(F):
        return F.copy()
    s = 1.0 / np.sqrt(M)
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
    M = (M_x, M_z) in kg/m, e.g. pipe_mass(W, D).  env = None means no
    contact: the soil carries nothing."""

    def __init__(self, env: YieldEnvelope, W: float, M, Fx=None, F_L=None):
        self.env = env
        self.W = W
        self.M = np.array(M, dtype=float)
        if np.any(self.M <= 0.0):
            raise ValueError("The masses must be positive.")
        self.Fx = Fx if Fx is not None else (lambda t: 0.0)
        self.F_L = F_L if F_L is not None else (lambda t: 0.0)
        self._Mt = 1e-3 * self.M  # t/m

    def load(self, t: float) -> np.ndarray:
        """F_ext = (Fx, W − F_L) on the pipe at time t (kN/m, V down)."""
        return np.array([self.Fx(t), self.W - self.F_L(t)])

    def rest(self, t: float = 0.0) -> DynamicState:
        """The pipe at rest at u = w = 0, with the soil force closest to the
        load at t."""
        H, V = self.soil_force(self.load(t))
        return DynamicState(t, 0.0, 0.0, 0.0, 0.0, H, V)

    def soil_force(self, F_star) -> np.ndarray:
        """The point of the envelope closest to F_star (metric M⁻¹); zero
        without contact."""
        if self.env is None:
            return np.zeros(2)
        return closest_point(self.env, F_star, self.M)

    def step(self, s: DynamicState, dt: float) -> DynamicState:
        """One backward Euler step of dt from s: the soil force is the point
        of the envelope closest to F* = F_ext(t + dt) + M·v/dt (metric
        M⁻¹), and the new velocity v = dt·M⁻¹·(F* − F)."""
        if dt <= 0.0:
            raise ValueError("The time step must be positive.")
        t = s.t + dt
        v = np.array([s.u_dot, s.w_dot])
        F_star = self.load(t) + self._Mt * v / dt
        F = self.soil_force(F_star)
        v = dt * (F_star - F) / self._Mt
        return DynamicState(t, s.u + dt * v[0], s.w + dt * v[1], v[0], v[1], F[0], F[1])

    def run(self, t_end: float, dt: float, start=None) -> list[DynamicState]:
        """States from start (default: at rest at t = 0) to t_end, in steps
        of dt at the times start.t + k·dt (without round-off drift).  The
        last step is shorter if t_end is not a multiple of dt."""
        s = start if start is not None else self.rest()
        states = [s]
        for t in _times(s.t, t_end, dt):
            states.append(self.step(states[-1], t - states[-1].t))
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


def _times(t0: float, t_end: float, dt: float) -> list[float]:
    """The times t0 + k·dt up to t_end (without round-off drift), and t_end
    itself if it is not one of them."""
    if dt <= 0.0:
        raise ValueError("The time step must be positive.")
    n = int(np.floor((t_end - t0) / dt + 1e-9))
    times = [t0 + k * dt for k in range(1, n + 1)]
    if t_end - (times[-1] if times else t0) > 1e-9 * dt:
        times.append(t_end)
    return times


class _Backward(Exception):
    """The motion needs a mechanism of the mirror image (step 4)."""


class PipeMotion:
    """The pipe with mass, loaded in time, while the geometry follows the
    motion forward (step 3 of PLAN_DYNAMIKK.md; manual.tex, "Geometry that
    follows the motion").

    lat and geo are a PipeLateral model with uniform friction (e.g.
    PipeLateralClosedForm) and its state, e.g. from after_penetration.  The
    time steps are those of PipeDynamics on the envelope of the current
    geometry (YieldEnvelope.from_state with n directions).  When the pipe
    has moved the path length ds, the geometry is updated with
    PipeLateral.move, the soil flowing with the mechanisms of the envelope
    edges the velocity has followed (two at a corner), and a new envelope
    is made.  Without contact the soil force is zero, and the geometry is
    updated as soon as the pipe touches the soil again.  The run stops
    (stop_reason = "backward") before a step whose motion needs a
    mechanism of the mirror image, which comes in step 4.

    run() returns DynamicStates with u and w from the start; geometries
    lists (t, LateralState, envelope or None) at the start and at every
    update."""

    def __init__(self, lat, geo, M, Fx=None, F_L=None, ds=0.25e-3, n=720):
        if not lat.uniform:
            raise ValueError("PipeMotion needs uniform friction (uniform=True).")
        self.lat, self.geo0, self.ds, self.n = lat, geo, ds, n
        self.dyn = PipeDynamics(None, lat.W, M, Fx, F_L)
        self.geometries = []
        self.stop_reason = None

    def envelope(self, geo) -> YieldEnvelope | None:
        """Envelope of a geometry; None if the pipe touches no soil."""
        if geo.x_right - geo.x_left <= 1e-12:
            return None
        return YieldEnvelope.from_state(self.lat, geo, self.n)

    def run(self, t_end: float, dt: float) -> list[DynamicState]:
        """Time steps of dt from rest at t = 0 to t_end (shorter last step),
        or until the run stops."""
        dyn, geo = self.dyn, self.geo0
        dyn.env = self.envelope(geo)
        states = [dyn.rest()]
        self.geometries = [(0.0, geo, dyn.env)]
        self.stop_reason = None
        paths, path = {}, 0.0  # path along each sampled mechanism, and in all
        for t in _times(0.0, t_end, dt):
            prev = states[-1]
            s = dyn.step(prev, t - prev.t)
            v = np.array([s.u_dot, s.w_dot])
            if v.any():
                try:
                    split = self._split(dyn.env, s, v)
                except _Backward:
                    self.stop_reason = "backward"
                    break
                for j, a in split:
                    paths[j] = paths.get(j, 0.0) + a * (s.t - prev.t)
                path += float(np.hypot(*v)) * (s.t - prev.t)
            states.append(s)
            landed = dyn.env is None and self._touches(geo, s)
            if path >= self.ds or landed:
                geo = self._update(geo, states[-1], paths)
                paths, path = {}, 0.0
        if path > 0.0:
            self._update(geo, states[-1], paths)
        return states

    @staticmethod
    def _split(env, s: DynamicState, v: np.ndarray) -> list[tuple[int, float]]:
        """The velocity v at the soil force of s as speeds a_k ≥ 0 along the
        directions of motion there: one edge, or the two edges of a corner.
        Returns (index of the sampled direction, a_k) for the edges whose
        mechanism moves soil; raises _Backward if one of them is a
        mechanism of the mirror image."""
        if env is None:
            return []
        F, speed = np.array([s.H, s.V]), float(np.hypot(*v))
        lo, hi = env.directions(F, tol=1e-7)
        if hi - lo < 1e-9:
            dirs, speeds = [lo], [speed]
        else:
            A = np.array([[np.cos(lo), np.cos(hi)], [np.sin(lo), np.sin(hi)]])
            dirs, speeds = [lo, hi], np.linalg.solve(A, v)
        out = []
        for d, a in zip(dirs, speeds):
            j = PipeMotion._sample(env, d, F)
            if a <= 1e-9 * speed or env.D[j] <= 0.0:
                continue  # no motion along it (round-off), or no soil
            if j >= len(env.deltas) // 2:
                raise _Backward
            out.append((j, float(a)))
        return out

    @staticmethod
    def _sample(env, d: float, F) -> int:
        """Index of the sampled direction whose half-plane forms the edge
        with normal d through F (forward first, if both halves give it)."""
        diff = np.abs((env.deltas - d + np.pi) % (2.0 * np.pi) - np.pi)
        near = np.nonzero((diff < 1e-4) & np.isfinite(env.D))[0]
        if len(near) == 0:
            raise ValueError("No sampled direction for an edge of the envelope.")
        work = np.cos(d) * F[0] + np.sin(d) * F[1]
        return int(near[np.argmin(np.abs(env.D[near] - work))])

    def _touches(self, geo, s: DynamicState) -> bool:
        """Whether the pipe at s reaches below the soil surface anywhere."""
        u, w = self.geo0.u + s.u, self.geo0.w + s.w
        return bool(np.any(self.lat.bottom(self.lat.x, u, w) + geo.eta >= 0.0))

    def _update(self, geo, s: DynamicState, paths: dict):
        """Move the geometry to the pipe at s, the soil flowing with the
        mechanisms and paths since the last update, and make its envelope."""
        lat, env = self.lat, self.dyn.env
        flow = [
            (lat.field(geo, float(env.deltas[j]), lat.omega_r), p)
            for j, p in sorted(paths.items())
        ]
        du, dw = self.geo0.u + s.u - geo.u, self.geo0.w + s.w - geo.w
        geo = lat.move(geo, du, dw, flow)
        self.dyn.env = self.envelope(geo)
        self.geometries.append((s.t, geo, self.dyn.env))
        return geo
