"""
Yield envelope of the soil around a pipe in the (H, V) plane, for one
geometry (the first step of the dynamic model, PLAN_DYNAMIKK.md).

(H, V) is the force from the pipe on the soil, H toward +x and V downward
(y down).  When the pipe moves with unit speed along n(δ) = (cos δ, sin δ),
δ below the horizontal, the mechanism of PipeLateral for that direction
dissipates D(δ) per unit speed.  By the upper-bound theorem the soil can
only carry forces with

    H·cos δ + V·sin δ ≤ D(δ)    for all δ,

so the envelope is the intersection of these half-planes: a convex polygon
when D is sampled in finitely many directions.  D comes in three parts:

  * forward (−90° ≤ δ ≤ 90°): the state's one-sided mechanism toward +x;
  * backward (90° ≤ δ ≤ 270°): the mirror image (PipeLateral.reversed),
    with δ = 180° − δ', δ' the direction in the mirror image;
  * D = 0 where the pipe leaves all the soil it touches (no tension).

At V = W the largest H on the envelope is the minimum over δ of
(D − W·sin δ)/cos δ, the H of PipeLateral.mechanism(), and the direction
of the edge there is β (normality).  At a corner the pipe can move in any
direction between the normals of the two edges that meet there.
"""

from __future__ import annotations

import numpy as np

from .lateral import LateralState, PipeLateral


class YieldEnvelope:
    """Convex polygon of the forces (H, V) the soil can carry.

    vertices are counter-clockwise in the (H, V) plane (H to the right, V
    up).  Edge k runs from vertex k to vertex k + 1; edge_delta[k] is the
    direction δ of its outward normal (the direction of motion on that
    edge) and edge_D[k] the dissipation n(δ)·F there.  deltas and D are the
    sampled directions and dissipations the envelope was made from (None
    when it was made from vertices).
    """

    def __init__(self, vertices, deltas=None, D=None):
        P = _clean(np.asarray(vertices, dtype=float))
        if len(P) < 3:
            raise ValueError("An envelope needs at least three vertices.")
        if _signed_area(P) < 0.0:
            P = P[::-1]
        E = np.roll(P, -1, axis=0) - P
        turn = (
            E[:, 0] * np.roll(E, -1, axis=0)[:, 1]
            - E[:, 1] * np.roll(E, -1, axis=0)[:, 0]
        )
        self.scale = float(np.max(np.abs(P)))
        if np.any(turn < -1e-9 * self.scale**2):
            raise ValueError("The vertices do not form a convex polygon.")
        self.vertices = P
        self.edge_delta = np.arctan2(-E[:, 0], E[:, 1])
        self.edge_D = P[:, 0] * np.cos(self.edge_delta) + P[:, 1] * np.sin(
            self.edge_delta
        )
        self.deltas = None if deltas is None else np.asarray(deltas, dtype=float)
        self.D = None if D is None else np.asarray(D, dtype=float)

    @classmethod
    def from_dissipation(cls, deltas, D) -> YieldEnvelope:
        """Envelope from D(δ) sampled in the directions deltas (rad).
        Directions with D = nan are skipped.  The directions must surround
        the envelope (no gap of 180° or more), or it is unbounded."""
        deltas, D = np.asarray(deltas, dtype=float), np.asarray(D, dtype=float)
        ok = np.isfinite(D)
        if np.any(D[ok] < 0.0):
            raise ValueError("The dissipation must not be negative.")
        big = 10.0 * (float(np.max(D[ok], initial=0.0)) + 1.0)
        P = big * np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])
        for d, h in zip(deltas[ok], D[ok]):
            P = _clip(P, np.array([np.cos(d), np.sin(d)]), h)
        if np.max(np.abs(P)) > 0.5 * big:
            raise ValueError("Unbounded envelope: the directions leave a gap.")
        return cls(P, deltas, D)

    @classmethod
    def from_state(
        cls, lat: PipeLateral, state: LateralState, n: int = 720
    ) -> YieldEnvelope:
        """Envelope of a PipeLateral state with uniform friction
        (PipeLateralClosedForm or PipeLateral(uniform=True)), from n + 2
        directions: n/2 + 1 forward from −90° to 90° and as many backward,
        in the mirror image (lat.reversed(state))."""
        if not lat.uniform:
            raise ValueError("The envelope needs uniform friction (uniform=True).")
        half = np.linspace(-0.5 * np.pi, 0.5 * np.pi, n // 2 + 1)
        forward = dissipation(lat, state, half)
        rev, back = lat.reversed(state)
        backward = dissipation(rev, back, half)
        deltas = np.concatenate((half, np.pi - half[::-1]))
        return cls.from_dissipation(deltas, np.concatenate((forward, backward[::-1])))

    # ------------------------------------------------------------------

    def support(self, delta: float) -> tuple[np.ndarray, float]:
        """Force F on the envelope for motion along δ (normality) and the
        dissipation n(δ)·F.  If δ is the normal of an edge, F is the
        middle of that edge."""
        n = np.array([np.cos(delta), np.sin(delta)])
        work = self.vertices @ n
        top = work.max()
        F = self.vertices[work >= top - 1e-12 * self.scale].mean(axis=0)
        return F, float(top)

    def contains(self, F, tol: float = 1e-9) -> bool:
        """True if F = (H, V) lies inside the envelope or on it (within
        tol·scale)."""
        return bool(np.max(self._excess(F)) <= tol * self.scale)

    def H_range(self, V: float) -> tuple[float, float] | None:
        """(H_min, H_max) on the envelope at V, or None if V is outside it.
        H_max is the force that moves the pipe forward at V (the H of
        PipeLateral.mechanism() at V = W), H_min the one that moves it back."""
        return self._chord(1, V)

    def V_range(self, H: float) -> tuple[float, float] | None:
        """(V_min, V_max) on the envelope at H, or None if H is outside it."""
        return self._chord(0, H)

    def push_force(self, V: float) -> tuple[float, float]:
        """H and δ for pushing the pipe forward at V with its horizontal
        motion imposed (u̇ > 0), as PipeLateral does at V = W: the smallest
        (D − V·sin δ)/cos δ over the sampled forward directions.  The
        backward half-planes do not limit it, since the pipe may not move
        back.  H is very negative if the forward mechanisms cannot carry V
        (the pipe sinks)."""
        if self.D is None:
            raise ValueError("The envelope was not made from sampled directions.")
        fwd = np.isfinite(self.D) & (np.cos(self.deltas) > 1e-9)
        d = self.deltas[fwd]
        H = (self.D[fwd] - V * np.sin(d)) / np.cos(d)
        i = int(np.argmin(H))
        return float(H[i]), float(d[i])

    def directions(self, F, tol: float = 1e-9) -> tuple[float, float]:
        """Range (δ_lo, δ_hi) of the directions of motion at F on the
        boundary: equal on an edge, the normals of the two edges at a
        corner."""
        on = np.nonzero(np.abs(self._excess(F)) <= tol * self.scale)[0]
        if len(on) == 0:
            raise ValueError("F is not on the envelope.")
        d = self.edge_delta[on]
        d = d[0] + (d - d[0] + np.pi) % (2.0 * np.pi) - np.pi
        return float(d.min()), float(d.max())

    def active(self, tol: float = 1e-9) -> np.ndarray:
        """For each sampled direction: True if its half-plane touches the
        envelope.  False inside a corner, where the direction's own
        mechanism dissipates more than the corner force does."""
        if self.D is None:
            raise ValueError("The envelope was not made from sampled directions.")
        work = np.array([self.support(d)[1] for d in self.deltas])
        return np.isfinite(self.D) & (work >= self.D - tol * self.scale)

    # ------------------------------------------------------------------

    def _excess(self, F) -> np.ndarray:
        """n_k·F − D_k for every edge (≤ 0 inside)."""
        H, V = np.asarray(F, dtype=float)
        return H * np.cos(self.edge_delta) + V * np.sin(self.edge_delta) - self.edge_D

    def _chord(self, axis: int, value: float):
        P = self.vertices
        if not P[:, axis].min() <= value <= P[:, axis].max():
            return None
        Q = np.roll(P, -1, axis=0)
        a, b = P[:, axis] - value, Q[:, axis] - value
        hit = (np.minimum(a, b) <= 0.0) & (np.maximum(a, b) >= 0.0)
        flat = a == b
        t = np.divide(a, a - b, out=np.zeros_like(a), where=~flat)
        X = P + t[:, None] * (Q - P)
        other = np.concatenate((X[hit, 1 - axis], Q[hit & flat, 1 - axis]))
        return float(other.min()), float(other.max())


def dissipation(lat: PipeLateral, state: LateralState, deltas) -> np.ndarray:
    """D(δ) of the state's one-sided mechanism toward +x with uniform
    friction (ω = ½·asin r), per unit pipe speed; 0 where the pipe has no
    contact left, nan where the net is not valid."""
    out = np.full(len(deltas), np.nan)
    for i, d in enumerate(deltas):
        try:
            out[i] = lat.field(state, d, lat.omega_r).dissipation()["total"]
        except ValueError as e:
            if "No contact" in str(e):
                out[i] = 0.0
        except (IndexError, FloatingPointError):
            pass
    return out


def _clip(P: np.ndarray, n: np.ndarray, h: float) -> np.ndarray:
    """Part of the convex polygon P with n·p ≤ h (Sutherland-Hodgman)."""
    f = P @ n - h
    if np.all(f <= 0.0):
        return P
    Q, fq = np.roll(P, -1, axis=0), np.roll(f, -1)
    inside, cross = f <= 0.0, f * fq < 0.0
    t = np.divide(f, f - fq, out=np.zeros_like(f), where=cross)
    X = P + t[:, None] * (Q - P)
    keep = np.stack((inside, cross), axis=1).reshape(-1)
    return np.stack((P, X), axis=1).reshape(-1, 2)[keep]


def _clean(P: np.ndarray) -> np.ndarray:
    """Drop vertices that repeat the previous one (within round-off)."""
    tol = 1e-10 * max(float(np.max(np.abs(P), initial=0.0)), 1.0)
    keep = np.hypot(*(P - np.roll(P, 1, axis=0)).T) > tol
    return P[keep] if keep.any() else P[:1]


def _signed_area(P: np.ndarray) -> float:
    Q = np.roll(P, -1, axis=0)
    return 0.5 * float(np.sum(P[:, 0] * Q[:, 1] - Q[:, 0] * P[:, 1]))
