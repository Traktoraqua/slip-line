"""Checks of the equations of motion on a fixed yield envelope
(slip_line/dynamics.py), step 2 of PLAN_DYNAMIKK.md.

Run with pytest (uv run pytest) or as a script:
    uv run tests/test_dynamics.py
"""

import sys
from functools import cache
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from slip_line.dynamics import PipeDynamics, closest_point, pipe_mass
from slip_line.envelope import YieldEnvelope
from slip_line.lateral import PipeLateralClosedForm
from slip_line.penetration import PipePenetrationClosedForm

W = 2.8  # kN/m
BOX = YieldEnvelope([[-1.0, 0.0], [1.0, 0.0], [1.0, 5.0], [-1.0, 5.0]])


@cache
def pipe_envelope():
    """The pipe of PLAN_DYNAMIKK.md (D = 1 m, cu = 5 kPa, r = 0.8) at the
    start of the push: H from −0.1556 to 0.1550 kN/m at V = W."""
    sim = PipePenetrationClosedForm(D=1.0, cu=5.0, r=0.8, W=W)
    lat, start = PipeLateralClosedForm.after_penetration(sim, sim.run()[-1])
    return YieldEnvelope.from_state(lat, start)


def newmark(F0, H0, omega, Mx):
    """Newmark's sliding block under F0·sin ωt from rest (manual.tex,
    Equations of motion): start t1, stop t2 and the displacement u(t2),
    with Mx in t/m.  Also returns u(t) while it slides."""
    t1 = np.arcsin(H0 / F0) / omega
    c1 = np.cos(omega * t1)

    def v(t):
        return (F0 / omega * (c1 - np.cos(omega * t)) - H0 * (t - t1)) / Mx

    def u(t):
        tau, s1 = t - t1, np.sin(omega * t1)
        sliding = F0 / omega * (tau * c1 - (np.sin(omega * t) - s1) / omega)
        return (sliding - 0.5 * H0 * tau**2) / Mx

    a, b = (np.pi - omega * t1) / omega, t1 + 2.0 * np.pi / omega  # v(a) > 0
    for _ in range(100):
        c = 0.5 * (a + b)
        a, b = (c, b) if v(c) > 0.0 else (a, c)
    t2 = 0.5 * (a + b)
    return t1, t2, u(t2), u


def block(F0, dt, T=1.0, periods=1.0):
    """Sliding block on BOX (H0 = 1, V0 = 5) with W = 2 and M = 1 t/m."""
    dyn = PipeDynamics(
        BOX, W=2.0, M=(1000.0, 1000.0), Fx=lambda t: F0 * np.sin(2 * np.pi * t / T)
    )
    return dyn, dyn.run(periods * T, dt)


# ----------------------------------------------------------------------
# Mass and projection
# ----------------------------------------------------------------------


def test_pipe_mass():
    """m = W/g + ρ_w·πD²/4 = 285.42 + 805.03 = 1090.46 kg/m, plus C_a·805.03."""
    Mx, Mz = pipe_mass(W, 1.0, C_a=0.0)
    assert abs(Mx - 1090.46) < 0.01 and Mz == Mx
    Mx, Mz = pipe_mass(W, 1.0)
    assert abs(Mx - 1895.49) < 0.01 and Mz == Mx
    Mx, Mz = pipe_mass(W, 1.0, C_a=(1.0, 2.0))
    assert abs(Mz - Mx - 805.03) < 0.01


def test_closest_point_box():
    """For a box and a diagonal M the closest point clamps each
    coordinate, whatever the masses; a point inside stays."""
    M = (1000.0, 4000.0)
    for F, expect in [
        ([3.0, 1.0], [1.0, 1.0]),
        ([3.0, 7.0], [1.0, 5.0]),
        ([0.5, -2.0], [0.5, 0.0]),
        ([-4.0, -1.0], [-1.0, 0.0]),
        ([0.2, 1.0], [0.2, 1.0]),
    ]:
        assert np.allclose(closest_point(BOX, F, M), expect), F


def test_closest_point_normality():
    """Outside the envelope the closest point lies on it, and
    M⁻¹(F* − F) points along a direction of motion there (normality),
    which makes it the closest point.  A triangle with unequal masses,
    and the pipe envelope."""
    rng = np.random.default_rng(1)
    tri = YieldEnvelope([[0.0, 0.0], [3.0, 3.0], [-3.0, 3.0]])
    for env, M, size in [
        (tri, (1000.0, 4000.0), 6.0),
        (pipe_envelope(), (1895.0, 1895.0), 4.0),
    ]:
        n_out = 0
        for F_star in rng.uniform(-size, size, (200, 2)):
            F = closest_point(env, F_star, M)
            if env.contains(F_star):
                assert np.allclose(F, F_star)
                continue
            n_out += 1
            assert np.max(env._excess(F)) < 1e-9 * env.scale
            g = (F_star - F) / np.asarray(M)
            lo, hi = env.directions(F, tol=1e-7)
            d = np.arctan2(g[1], g[0])
            d = lo + (d - lo + np.pi) % (2.0 * np.pi) - np.pi
            assert lo - 1e-6 <= d <= hi + 1e-6, np.degrees([lo, d, hi])
        assert n_out > 50


# ----------------------------------------------------------------------
# Motion
# ----------------------------------------------------------------------


def test_rest_inside_envelope():
    """|Fx| = 0.1 < 0.155 kN/m at V = W: the pipe stays put, and the soil
    carries the load."""
    dyn = PipeDynamics(
        pipe_envelope(),
        W=W,
        M=pipe_mass(W, 1.0),
        Fx=lambda t: 0.1 * np.sin(2 * np.pi * t / 10.0),
    )
    for s in dyn.run(10.0, 0.01):
        assert s.u == 0.0 and s.w == 0.0 and s.u_dot == 0.0 and s.w_dot == 0.0
        assert abs(s.H - dyn.Fx(s.t)) < 1e-12 and abs(s.V - W) < 1e-12
        assert np.isnan(s.delta)


def test_sliding_block():
    """F0 = 1.5·H0: the block slides forward Δu, rests, slides the same
    back half a period later, and forward again (Newmark), within 0.1 %.
    It never moves vertically, and V = W throughout."""
    F0, T = 1.5, 1.0
    t1, t2, du, _ = newmark(F0, 1.0, 2 * np.pi / T, 1.0)
    _, states = block(F0, T / 2000, T, periods=1.6)
    t = np.array([s.t for s in states])
    u = np.array([s.u for s in states])
    assert abs(u[np.searchsorted(t, 0.5 * T + t1) - 1] - du) < 1e-3 * du
    assert abs(u[np.searchsorted(t, T + t1) - 1]) < 1e-3 * du
    assert abs(u[-1] - du) < 1e-3 * du  # 1.6 T > T + t2
    for s in states:
        assert abs(s.w) < 1e-12 and abs(s.w_dot) < 1e-12 and abs(s.V - 2.0) < 1e-12
    moving = np.array([abs(s.u_dot) > 1e-12 for s in states])
    assert not moving[t < t1 - T / 2000].any()
    assert not moving[(t > t2 + T / 2000) & (t < 0.5 * T + t1 - T / 2000)].any()


def test_sliding_block_converges():
    """Backward Euler is first order: halving Δt halves the error of u
    in the middle of the slide (t = 0.3 T)."""
    F0, T = 1.5, 1.0
    *_, u = newmark(F0, 1.0, 2 * np.pi / T, 1.0)
    err = []
    for n in (1000, 2000, 4000):
        _, states = block(F0, T / n, T, periods=0.3)
        assert abs(states[-1].t - 0.3 * T) < 1e-12
        err.append(states[-1].u - u(0.3 * T))
    ratio = np.array(err[:-1]) / np.array(err[1:])
    assert np.all(np.abs(ratio - 2.0) < 0.05), ratio
    assert abs(err[-1]) < 5e-3 * u(0.3 * T)


def test_lift_off():
    """Lift F_L = 2W: the soil force is zero, and the pipe accelerates
    upward at (F_L − W)/M_z (exact in velocity for a constant force)."""
    M = pipe_mass(W, 1.0)
    dyn = PipeDynamics(pipe_envelope(), W=W, M=M, F_L=lambda t: 2.0 * W)
    for s in dyn.run(0.5, 0.01)[1:]:
        assert abs(s.H) < 1e-12 and abs(s.V) < 1e-12 and abs(s.u) < 1e-12
        assert abs(s.w_dot + W / (M[1] / 1000.0) * s.t) < 1e-12
        assert abs(np.degrees(s.delta) + 90.0) < 1e-9


def test_energy_balance():
    """Work of the loads = dissipation + kinetic energy + numerical
    damping, exactly; the numerical damping halves with Δt."""
    T = 2.0
    num = []
    for n in (1000, 2000):
        dyn = PipeDynamics(
            pipe_envelope(),
            W=W,
            M=pipe_mass(W, 1.0),
            Fx=lambda t: 0.3 * np.sin(2 * np.pi * t / T),
        )
        states = dyn.run(T, T / n)
        e = dyn.energy(states)
        assert e["dissipation"] > 0.0 and e["work"] > 0.0
        total = e["dissipation"] + e["kinetic"] + e["numerical"]
        assert abs(e["work"] - total) < 1e-9 * e["work"], e
        num.append(e["numerical"])
    assert 1.8 < num[0] / num[1] < 2.2, num
    assert num[1] < 1e-2 * e["work"]


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_") and callable(f):
            f()
            print("ok ", name)
