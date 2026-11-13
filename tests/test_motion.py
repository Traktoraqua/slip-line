"""Checks of the geometry that follows the motion (PipeLateral.move and
slip_line/dynamics.py PipeMotion), step 3 of PLAN_DYNAMIKK.md.

Run with pytest (uv run pytest) or as a script:
    uv run tests/test_motion.py
"""

import sys
from functools import cache
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from slip_line.dynamics import PipeMotion, pipe_mass
from slip_line.lateral import PipeLateralClosedForm
from slip_line.penetration import PipePenetrationClosedForm

W, D = 2.8, 1.0  # kN/m, m
M = pipe_mass(W, D)  # C_a = 1


@cache
def start(cu):
    """Lateral model and start state after penetration (r = 0.8)."""
    sim = PipePenetrationClosedForm(D=D, cu=cu, r=0.8, W=W)
    return PipeLateralClosedForm.after_penetration(sim, sim.run()[-1])


@cache
def quasi_static_cu1():
    """The quasi-static push of run_lateral_closed_form.py, 30 mm:
    arrays u, w (m) and H (kN/m)."""
    lat, s0 = start(1.0)
    states = lat.run(s0, 0.03, du=0.25e-3)
    return (
        np.array([s.u for s in states]),
        np.array([s.w for s in states]),
        np.array([s.mech.H for s in states]),
    )


def rising(F_max, T_rise):
    """Fx = F_max·(1 − cos(πt/T_rise))/2: rises from 0 to F_max at
    T_rise, then falls back."""
    return lambda t: 0.5 * F_max * (1.0 - np.cos(np.pi * t / T_rise))


@cache
def slow_cu1(ds=0.25e-3):
    """cu = 1 kPa, Fx rising to 0.6 kN/m over 50 s and falling again
    (T ≫ T_n ≈ 2–3 s), run until it stops or t = 100 s."""
    lat, s0 = start(1.0)
    motion = PipeMotion(lat, s0, M, Fx=rising(0.6, 50.0), ds=ds)
    return motion, motion.run(100.0, 0.02)


@cache
def fast_cu1(ds=0.25e-3, dt=2e-3):
    """cu = 1 kPa, the same load rising over 1 s, run to t = 3 s."""
    lat, s0 = start(1.0)
    motion = PipeMotion(lat, s0, M, Fx=rising(0.6, 1.0), ds=ds)
    return motion, motion.run(3.0, dt)


# ----------------------------------------------------------------------
# PipeLateral.move
# ----------------------------------------------------------------------


def test_advance_is_move_along_the_mechanism():
    """advance(du) is move() with the state's mechanism: same position,
    surface and contact."""
    lat, s0 = start(1.0)
    s0.mech = lat.mechanism(s0)
    m = s0.mech
    a = lat.advance(s0, 1e-3)
    vf = lat.field(s0, m.beta, m.omega)
    b = lat.move(s0, 1e-3, 1e-3 * np.tan(m.beta), [(vf, 1e-3 / np.cos(m.beta))])
    assert (a.u, a.w, a.x_left, a.x_right, a.alpha) == (
        b.u,
        b.w,
        b.x_left,
        b.x_right,
        b.alpha,
    )
    assert np.array_equal(a.eta, b.eta) and np.array_equal(a.cut, b.cut)


def test_move_without_soil_flow():
    """Lifting the pipe 1 mm straight up moves no soil: the surface stays,
    and the volume with it."""
    lat, s0 = start(5.0)
    s1 = lat.move(s0, 0.0, -1e-3, [])
    assert s1.u == s0.u and abs(s1.w - (s0.w - 1e-3)) < 1e-15
    assert np.array_equal(s1.eta, s0.eta) and s1.alpha == s0.alpha


def test_move_mix_conserves_volume():
    """A corner mix of two mechanisms (paths 0.3 and 0.2 mm) moves the
    pipe by their sum and keeps the soil volume."""
    lat, s0 = start(1.0)
    d1, d2 = np.radians(10.0), np.radians(20.0)
    f1, f2 = lat.field(s0, d1, lat.omega_r), lat.field(s0, d2, lat.omega_r)
    du = 0.3e-3 * np.cos(d1) + 0.2e-3 * np.cos(d2)
    dw = 0.3e-3 * np.sin(d1) + 0.2e-3 * np.sin(d2)
    s1 = lat.move(s0, du, dw, [(f1, 0.3e-3), (f2, 0.2e-3)])
    assert abs(s1.u - du) < 1e-15 and abs(s1.w - s0.w - dw) < 1e-15
    v0, v1 = np.trapezoid(s0.eta, lat.x), np.trapezoid(s1.eta, lat.x)
    assert abs(v1 - v0) < 1e-12
    assert np.max(s1.eta - s0.eta) > 0.0  # the berm has grown


# ----------------------------------------------------------------------
# PipeMotion
# ----------------------------------------------------------------------


def test_slow_load_follows_quasi_static_push():
    """T ≫ T_n: while the pipe moves, the load equals the quasi-static
    push force at the same u, within 1 %, and the path w(u) is the
    quasi-static one, within 0.05 mm."""
    motion, states = slow_cu1()
    u_qs, w_qs, H_qs = quasi_static_cu1()
    _, s0 = start(1.0)
    moving = [s for s in states if s.u_dot > 0.0]
    assert len(moving) > 100
    for s in moving:
        H = np.interp(s.u, u_qs, H_qs)
        assert abs(motion.dyn.Fx(s.t) - H) < 0.01 * H, (s.t, s.u)
        assert abs(s0.w + s.w - np.interp(s.u, u_qs, w_qs)) < 5e-5, (s.t, s.u)
    peak = max(states, key=lambda s: motion.dyn.Fx(s.t))
    u_peak = np.interp(0.6, H_qs, u_qs)
    assert abs(peak.u - u_peak) < 0.5e-3, (peak.u, u_peak)


def test_falling_load_stops_before_moving_back():
    """After the peak the pipe rests until the load falls to the holding
    force H_min at V = W, and the run stops there: moving back needs the
    mirror image (step 4).  The pipe never moves back."""
    motion, states = slow_cu1()
    assert motion.stop_reason == "backward"
    assert all(s.u_dot >= 0.0 for s in states)
    env = motion.geometries[-1][2]
    H_min = env.H_range(W)[0]
    assert H_min > 0.3
    assert abs(motion.dyn.Fx(states[-1].t) - H_min) < 0.01, H_min
    assert 50.0 < states[-1].t < 100.0


def test_volume_is_conserved():
    """The soil volume is the same in every geometry of the slow run."""
    motion, _ = slow_cu1()
    lat = motion.lat
    v0 = np.trapezoid(motion.geometries[0][1].eta, lat.x)
    for _, geo, _ in motion.geometries:
        assert abs(np.trapezoid(geo.eta, lat.x) - v0) < 1e-12


def test_energy_balance_with_changing_envelope():
    """Work of the loads = dissipation + kinetic + numerical damping,
    exactly, although the envelope changes between steps."""
    motion, states = fast_cu1()
    e = motion.dyn.energy(states)
    assert e["dissipation"] > 0.0
    total = e["dissipation"] + e["kinetic"] + e["numerical"]
    assert abs(e["work"] - total) < 1e-9 * e["work"], e


def test_fast_load_lags_and_moves_on():
    """The same load rising over 1 s (T_n ≈ 2–3 s).  The pipe starts when
    the load reaches the push force, as for the slow load, but its mass
    makes it lag: it is fastest after the load has peaked, and it moves on
    while the falling load is below the soil's resistance, until it
    stops.  It ends within 10 % of the slow displacement."""
    motion, fast = fast_cu1()
    slow_motion, slow = slow_cu1()
    Fx, Fx_slow = motion.dyn.Fx, slow_motion.dyn.Fx
    first = next(s for s in fast if s.u_dot > 0.0)
    first_slow = next(s for s in slow if s.u_dot > 0.0)
    assert abs(Fx(first.t) - Fx_slow(first_slow.t)) < 0.02
    fastest = max(fast, key=lambda s: s.u_dot)
    assert fastest.t > 1.0  # after the load peak
    moving = [s for s in fast if s.u_dot > 0.0]
    assert Fx(moving[-1].t) < moving[-1].H  # braked by the soil
    u_fast, u_slow = max(s.u for s in fast), max(s.u for s in slow)
    assert abs(u_fast / u_slow - 1.0) < 0.1, (u_fast, u_slow)


def test_converges_with_ds_and_dt():
    """Halving Δs, and halving Δt, changes the largest displacement of the
    fast run less each time; at the finest, by less than 1 %."""
    u = {
        k: max(s.u for s in fast_cu1(*k)[1])
        for k in [(0.5e-3, 2e-3), (0.25e-3, 2e-3), (0.125e-3, 2e-3), (0.125e-3, 1e-3)]
    }
    d1 = abs(u[(0.5e-3, 2e-3)] - u[(0.25e-3, 2e-3)])
    d2 = abs(u[(0.25e-3, 2e-3)] - u[(0.125e-3, 2e-3)])
    d3 = abs(u[(0.125e-3, 2e-3)] - u[(0.125e-3, 1e-3)])
    assert d2 < d1, u
    assert d2 < 0.01 * u[(0.125e-3, 2e-3)] and d3 < 0.01 * u[(0.125e-3, 1e-3)], u


def test_cu5_lifts_loses_trench_floor_and_stops():
    """cu = 5 kPa (light pipe), load rising slowly to 0.5 kN/m.  From
    u ≈ 1.25 mm the pipe rises (δ < 0, the lift branch the quasi-static run
    misses).  Rising, it leaves the trench floor behind it, by microns: the
    contact then starts at the invert, and V = W lies outside the envelope,
    so the pipe would sink back into the gap, which mixes in a mechanism of
    the mirror image.  The run stops there (step 4), at the same place
    with half the Δs, within the Δs."""
    lat, s0 = start(5.0)
    stops = []
    for ds in (0.25e-3, 0.125e-3):
        motion = PipeMotion(lat, s0, M, Fx=rising(0.5, 50.0), ds=ds)
        states = motion.run(50.0, 0.02)
        assert motion.stop_reason == "backward"
        assert any(s.delta < 0.0 for s in states if not np.isnan(s.delta))
        _, geo, env = motion.geometries[-1]
        assert geo.x_left >= geo.u - 1e-9 and env.H_range(W) is None
        stops.append(states[-1].u)
    assert all(1.0e-3 < u < 3.0e-3 for u in stops), stops
    assert abs(stops[0] - stops[1]) < 0.25e-3, stops


def test_free_flight_and_landing():
    """Lift 2W for 0.05 s (cu = 5 kPa): the pipe leaves the soil, flies
    with no soil force, and lands straight down in its own hollow, which
    needs the mechanisms on both sides, so the run stops there."""
    lat, s0 = start(5.0)
    motion = PipeMotion(lat, s0, M, F_L=lambda t: 2.0 * W if t < 0.05 else 0.0)
    states = motion.run(2.0, 1e-3)
    times = [t for t, _, env in motion.geometries]
    off = [i for i, (_, _, env) in enumerate(motion.geometries) if env is None]
    assert off, "the pipe never left the soil"
    t_off = times[off[0]]
    t_on = next(times[i] for i in range(off[0], len(times)) if i not in off)
    in_air = [s for s in states if t_off < s.t < t_on]
    assert in_air and all(s.H == 0.0 and s.V == 0.0 for s in in_air)
    assert motion.stop_reason == "backward"
    assert motion.geometries[-1][2] is not None  # it touched the soil again
    assert abs(states[-1].u) < 1e-12


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_") and callable(f):
            f()
            print("ok ", name)
