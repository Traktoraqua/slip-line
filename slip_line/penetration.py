"""
Quasi-static self-weight penetration of a smooth pipe (disk) into weightless
Tresca soil, with the heave around the pipe tracked step by step.

The pipe is pushed down in small steps dw.  In each step:

  1. Contact: the pipe meets the current soil surface at x_c.  The slip-line
     field is that of DiskPrandtlSolver with a FLAT surface at the contact
     level, i.e. at the effective embedment d = invert depth below that level
     (heave raises the contact level, so d > w).
  2. Heave: the Geiringer field for a unit downward pipe velocity gives the
     upward velocity v_up(x) of the free surface (Region III).  The surface is
     raised by v_up·dw; the horizontal component is not used, since the field
     assumes a flat surface.
  3. Volume: the soil is incompressible, so the integral along the surface of
     the downward displacement of the soil's top boundary (pipe contact for
     x < x_c, free surface beyond) stays zero.  The Geiringer field already
     satisfies this to first order (∫ v_up dx = x_c); the heave increment is
     scaled to remove the discretisation error.
  4. Stop: the pipe sinks while the collapse load V(d) is below the pipe
     weight W; the last step is shortened so V = W.

Right half only (x ≥ 0).  y points down from the original surface; the heave
η is positive upward.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .boundary import DiskFoundationBoundary
from .solver import DiskPrandtlSolver
from .velocity import DiskVelocityField, alpha_beta_to_xy


@dataclass
class PenetrationState:
    """Pipe and soil surface after a step (right half)."""

    w: float  # invert depth below the original surface
    d: float  # effective embedment: invert depth below the contact level
    x_c: float  # contact half-width
    x_end: float  # outer end of the plastic zone at the surface (this step's field)
    V: float  # collapse load at d, full pipe per unit length
    heave: np.ndarray  # η on PipePenetration.x
    scale: float  # volume correction applied to the heave increment


class PipePenetration:
    """Smooth pipe of diameter D and weight W sinking into Tresca soil (cu)."""

    def __init__(
        self,
        D: float = 1.0,
        cu: float = 1.0,
        W: float = 1.0,
        dw: float = 2.5e-5,
        dx: float = 2.5e-4,
        x_max: float | None = None,
        n_fan: int = 18,
        n_disk: int = 30,
    ):
        self.D, self.R, self.cu, self.W, self.dw = D, D / 2.0, cu, W, dw
        x_max = 1.5 * D if x_max is None else x_max
        self.x = np.arange(0.0, x_max + 0.5 * dx, dx)
        self.n_fan, self.n_disk = n_fan, n_disk

    def solver(self, d: float) -> DiskPrandtlSolver:
        disk = DiskFoundationBoundary(cu=self.cu, D=self.D, p=d / self.D, r=0.0)
        return DiskPrandtlSolver(disk, n_fan=self.n_fan, n_disk=self.n_disk)

    def load(self, d: float) -> float:
        """Collapse load V at effective embedment d (closed form, no mesh)."""
        return self.solver(d).collapse_load_Q() if d > 0.0 else 0.0

    def contact(self, w: float, heave: np.ndarray) -> tuple[float, float]:
        """Contact half-width x_c and effective embedment d = w + η(x_c)."""
        R = self.R

        def gap(x):
            # Pipe bottom depth below the soil surface; > 0 where embedded.
            return w - R + np.sqrt(R * R - x * x) + np.interp(x, self.x, heave)

        if w <= 0.0:
            return 0.0, 0.0
        xs = self.x[self.x < R]
        i = int(np.argmax(gap(xs) <= 0.0))
        lo, hi = xs[i - 1], xs[i]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if gap(mid) > 0.0 else (lo, mid)
        x_c = 0.5 * (lo + hi)
        return x_c, w + float(np.interp(x_c, self.x, heave))

    def _integral(self, f: np.ndarray, xq: float) -> float:
        """∫_0^xq of the piecewise-linear interpolant of f on self.x."""
        x = self.x
        i = min(int(np.searchsorted(x, xq, side="right")) - 1, len(x) - 2)
        cum = np.sum(0.5 * (f[1 : i + 1] + f[:i]) * np.diff(x[: i + 1]))
        return float(cum + 0.5 * (f[i] + np.interp(xq, x, f)) * (xq - x[i]))

    def areas(self, w: float, heave: np.ndarray) -> dict[str, float]:
        """Volume balance of the right half.

        pipe  — ∫_0^x_c of the pipe bottom depth (the soil pushed down)
        heave — ∫_x_c^∞ η (the soil pushed up)
        error — pipe − heave, zero when the soil volume is conserved
        """
        x_c, _ = self.contact(w, heave)
        R, y_c = self.R, w - self.R
        pipe = y_c * x_c + 0.5 * (
            x_c * np.sqrt(R * R - x_c * x_c) + R * R * np.arcsin(x_c / R)
        )
        up = self._integral(heave, self.x[-1]) - self._integral(heave, x_c)
        return {"pipe": float(pipe), "heave": up, "error": float(pipe - up)}

    def surface_uplift_rate(self, d: float) -> tuple[np.ndarray, float]:
        """Upward free-surface velocity on self.x for a unit downward pipe
        velocity at effective embedment d (zero outside the plastic zone),
        and the outer end of the zone."""
        s = self.solver(d)
        s.solve()
        vf = DiskVelocityField(s, v0=1.0)
        vf.solve()
        row = s.region_III_nodes[0]
        xs = np.array([nd.x for nd in row])
        up = np.array(
            [
                -alpha_beta_to_xy(nd.psi, u, v)[1]
                for nd, (u, v) in zip(row, vf.region_III_uv[0])
            ]
        )
        if xs.max() >= self.x[-1]:
            raise ValueError("Plastic zone reaches x_max; increase x_max.")
        order = np.argsort(xs)
        rate = np.interp(self.x, xs[order], up[order], left=0.0, right=0.0)
        return rate, float(xs.max())

    def advance(self, state: PenetrationState, dw: float) -> PenetrationState:
        """Push the pipe down by dw from state and raise the surface."""
        w = state.w + dw
        _, d = self.contact(w, state.heave)
        rate, x_end = self.surface_uplift_rate(d)
        inc = rate * dw
        total = self._integral(inc, self.x[-1])
        scale = 1.0
        for _ in range(20):
            err = self.areas(w, state.heave + scale * inc)["error"]
            if abs(err) < 1e-15:
                break
            scale += err / total
        heave = state.heave + scale * inc
        x_c, d = self.contact(w, heave)
        return PenetrationState(w, d, x_c, x_end, self.load(d), heave, scale)

    def run(self, max_steps: int = 100_000) -> list[PenetrationState]:
        """Sink the pipe from the surface until V = W; return every state."""
        states = [PenetrationState(0.0, 0.0, 0.0, 0.0, 0.0, np.zeros_like(self.x), 1.0)]
        while len(states) <= max_steps:
            prev = states[-1]
            nxt = self.advance(prev, self.dw)
            if nxt.V >= self.W:
                # Shorten the last step so V = W (regula falsi on its length).
                lo, v_lo, hi, v_hi = 0.0, prev.V, self.dw, nxt.V
                for _ in range(50):
                    dw = lo + (hi - lo) * (self.W - v_lo) / (v_hi - v_lo)
                    nxt = self.advance(prev, dw)
                    if abs(nxt.V - self.W) < 1e-8 * self.W:
                        break
                    if nxt.V < self.W:
                        lo, v_lo = dw, nxt.V
                    else:
                        hi, v_hi = dw, nxt.V
                states.append(nxt)
                return states
            states.append(nxt)
        raise RuntimeError("No equilibrium within max_steps.")
