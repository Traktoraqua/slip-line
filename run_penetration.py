"""Smooth pipe (D = 1 m, W = 1 kN/m) sinking under its own weight into
weightless Tresca soil (cu = 1 kPa), with the heave updated in every step;
prints the steps and the final state with and without heave."""

from slip_line.penetration import PipePenetration

sim = PipePenetration(D=1.0, cu=1.0, W=1.0)
states = sim.run()

print(" step   w [mm]   d [mm]  heave at pipe [mm]  x_c [mm]  x_end [mm]  V [kN/m]")
for k, s in enumerate(states):
    if k % 20 == 0 or k == len(states) - 1:
        print(
            f"{k:5d} {1e3 * s.w:8.3f} {1e3 * s.d:8.3f} {1e3 * (s.d - s.w):19.3f}"
            f" {1e3 * s.x_c:9.1f} {1e3 * s.x_end:11.1f} {s.V:9.4f}"
        )

# Without heave the contact level stays at the original surface (d = w).
lo, hi = 0.0, sim.D
for _ in range(60):
    mid = 0.5 * (lo + hi)
    lo, hi = (mid, hi) if sim.load(mid) < sim.W else (lo, mid)

final = states[-1]
a = sim.areas(final.w, final.heave)
print(f"\nEquilibrium V = W = {sim.W} kN/m after {len(states) - 1} steps:")
print(f"  settlement w       = {1e3 * final.w:.2f} mm  (no heave: {1e3 * lo:.2f} mm)")
print(f"  heave at the pipe  = {1e3 * (final.d - final.w):.2f} mm")
print(f"  effective depth d  = {1e3 * final.d:.2f} mm")
print(f"  volume error       = {a['error']:.1e} m² per half (pipe {a['pipe']:.3e} m²)")
