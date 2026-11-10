"""Checks of YieldEnvelope (slip_line/envelope.py), step 1 of PLAN_DYNAMIKK.md.

Written as pytest tests, but runnable without pytest:
    uv run tests/test_yield_envelope.py
"""

import sys
from functools import cache
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from slip_line.envelope import YieldEnvelope
from slip_line.lateral import PipeLateralClosedForm
from slip_line.penetration import PipePenetrationClosedForm

W = 2.8  # kN/m, both cases
STEP = np.radians(0.5)  # direction step of the default envelope


@cache
def cu5():
    """The pipe of PLAN_DYNAMIKK.md: D = 1 m, cu = 5 kPa, r = 0.8, W = 2.8,
    pushed 8 mm with the quasi-static model."""
    sim = PipePenetrationClosedForm(D=1.0, cu=5.0, r=0.8, W=W)
    lat, start = PipeLateralClosedForm.after_penetration(sim, sim.run()[-1])
    states = lat.run(start, 0.008, du=0.25e-3)
    return lat, states


@cache
def cu1():
    """run_lateral_closed_form.py (cu = 1 kPa): 50 mm, 200 mm, 40 mm back."""
    sim = PipePenetrationClosedForm(D=1.0, cu=1.0, r=0.8, W=W)
    lat, start = PipeLateralClosedForm.after_penetration(sim, sim.run()[-1])
    states = lat.run(start, 0.05, du=0.25e-3)
    s50 = states[-1]
    states += lat.run(s50, 0.15)[1:]
    rev, back = lat.reversed(states[-1])
    states3 = rev.run(back, 0.01, du=0.25e-3)
    states3 += rev.run(states3[-1], 0.03)[1:]
    return [(lat, s50), (lat, states[-1]), (rev, states3[-1])]


def matches_mechanism(lat, state):
    """push_force(W) against a fresh mechanism search (a scan of β and a
    golden-section search): same H, and δ within one direction step of β.
    Returns the envelope and the mechanism."""
    env = YieldEnvelope.from_state(lat, state)
    m = lat.mechanism(state)
    H, delta = env.push_force(W)
    assert abs(H - m.H) < 2e-4 * max(1.0, abs(m.H)), (H, m.H)
    assert abs(delta - m.beta) <= STEP, np.degrees([delta, m.beta])
    return env, m


def forward_edge_at_W(env, m):
    """The envelope reaches V = W, and its largest H there is the push
    force, with β among the directions of motion at that point."""
    H = env.H_range(W)[1]
    lo, hi = env.directions([H, W])
    assert abs(H - m.H) < 2e-4 * max(1.0, abs(m.H)), (H, m.H)
    assert lo - STEP <= m.beta <= hi + STEP, np.degrees([lo, m.beta, hi])


# ----------------------------------------------------------------------
# The polygon itself
# ----------------------------------------------------------------------


def test_rectangle():
    env = YieldEnvelope([[-1.0, 0.0], [1.0, 0.0], [1.0, 2.0], [-1.0, 2.0]])
    F, D = env.support(0.0)  # normal of the right edge: its middle
    assert np.allclose(F, [1.0, 1.0]) and np.isclose(D, 1.0)
    F, D = env.support(np.pi / 4)  # inside the corner (1, 2)
    assert np.allclose(F, [1.0, 2.0]) and np.isclose(D, 3.0 / np.sqrt(2.0))
    assert env.H_range(1.0) == (-1.0, 1.0) and env.V_range(0.0) == (0.0, 2.0)
    assert env.H_range(2.5) is None
    assert np.allclose(env.directions([1.0, 2.0]), [0.0, np.pi / 2])
    assert np.allclose(env.directions([1.0, 1.0]), [0.0, 0.0])
    assert env.contains([0.5, 1.0]) and env.contains([1.0, 2.0])
    assert not env.contains([1.01, 1.0])


def test_from_dissipation_recovers_polygon():
    """Coulomb-type triangle |H| ≤ V ≤ 3: its support function sampled in
    720 directions (which include its edge normals, −135°, −45° and 90°)
    gives back the triangle."""
    tri = np.array([[0.0, 0.0], [3.0, 3.0], [-3.0, 3.0]])
    deltas = np.linspace(-np.pi, np.pi, 721)[:-1]
    D = np.max(tri @ np.array([np.cos(deltas), np.sin(deltas)]), axis=0)
    env = YieldEnvelope.from_dissipation(deltas, np.maximum(D, 0.0))
    assert len(env.vertices) == 3
    assert all(np.min(np.hypot(*(env.vertices - p).T)) < 1e-9 for p in tri)


def test_circle_converges():
    """D = 1 in every direction: the polygon lies just outside the unit
    circle, by 1/cos(Δδ/2) − 1."""
    deltas = np.linspace(-np.pi, np.pi, 721)[:-1]
    env = YieldEnvelope.from_dissipation(deltas, np.ones_like(deltas))
    r = np.hypot(*env.vertices.T)
    assert np.all(r >= 1.0 - 1e-12)
    assert np.max(r) < 1.0 / np.cos(0.5 * (deltas[1] - deltas[0])) + 1e-12
    assert env.active().all()


# ----------------------------------------------------------------------
# Envelopes of pipe states
# ----------------------------------------------------------------------


def test_start_matches_mechanism():
    """cu = 5 kPa at the start: H = 0.155 kN/m at β = 17°, as the
    quasi-static push."""
    lat, states = cu5()
    env, m = matches_mechanism(lat, states[0])
    forward_edge_at_W(env, m)
    assert abs(m.H - 0.155) < 5e-4 and abs(np.degrees(m.beta) - 17.0) < 0.1


def test_start_forward_and_back():
    """At the start the pipe sits in its own symmetric hollow, so pushing
    back takes nearly the same force (0.1555 against 0.1550 kN/m), and
    that matches the mechanism of the mirror image."""
    lat, states = cu5()
    env = YieldEnvelope.from_state(lat, states[0])
    H_back, H_fwd = env.H_range(W)
    _, back = lat.reversed(states[0])
    assert abs(-H_back - back.mech.H) < 2e-4, (H_back, back.mech.H)
    assert abs(H_fwd + H_back) < 2e-3


def test_no_tension_and_convex():
    """The origin is a corner (lift-off), nothing below V = 0, and the
    polygon is convex (checked when it is made)."""
    lat, states = cu5()
    env = YieldEnvelope.from_state(lat, states[0])
    assert env.H_range(0.0) == (0.0, 0.0)
    assert env.V_range(0.0)[0] == 0.0
    assert not env.contains([0.0, -1e-3])
    lo, hi = env.directions([0.0, 0.0])
    assert lo < -0.5 * np.pi < hi
    YieldEnvelope(env.vertices)  # raises if not convex


def test_local_minimum_cu5():
    """At u = 8 mm the quasi-static run (local search from the previous
    step) is still on the downward branch, H = 0.68 kN/m.  The envelope,
    like a fresh search, finds the lift branch at δ ≈ −4.8°, H = 0.47."""
    lat, states = cu5()
    s = states[-1]
    env, m = matches_mechanism(lat, s)
    forward_edge_at_W(env, m)
    assert s.mech.H > 0.67 and np.degrees(s.mech.beta) > 8.0
    assert abs(m.H - 0.474) < 2e-3 and np.degrees(m.beta) < -4.0


def test_states_cu1():
    """run_lateral_closed_form.py.  After 50 mm the pipe needs a holding
    force not to slide back (H_min > 0 at V = W).  After 200 mm the
    contact behind has gone, and V = W lies outside the envelope: with the
    force given, the pipe would sink back.  After 40 mm back (mirror
    image, trench in front) the forward edge is at V = W again."""
    (lat, s50), (_, s200), (rev, s_back) = cu1()
    env, m = matches_mechanism(lat, s50)
    forward_edge_at_W(env, m)
    assert env.H_range(W)[0] > 0.5
    env, m = matches_mechanism(lat, s200)
    assert env.H_range(W) is None
    env, m = matches_mechanism(rev, s_back)
    forward_edge_at_W(env, m)


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_") and callable(f):
            f()
            print("ok ", name)
