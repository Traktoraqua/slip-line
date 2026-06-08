"""
main.py — entry point for the Prandtl slip line analysis.

Usage:
  python main.py

Two foundation geometries available:
  1. Flat strip footing
  2. Partly-penetrated circular disk
"""

import numpy as np
from slip_line.boundary import FoundationBoundary, DiskFoundationBoundary
from slip_line.solver import PrandtlSolver, DiskPrandtlSolver
from slip_line.plot import plot_slip_lines, plot_disk_slip_lines


def get_float(prompt, default=None, lo=-np.inf, hi=np.inf):
    while True:
        raw = input(prompt).strip()
        if raw == '' and default is not None:
            return default
        try:
            val = float(raw)
            if lo <= val <= hi:
                return val
            print(f"  Value must be in [{lo}, {hi}].  Try again.")
        except ValueError:
            print("  Invalid input — please enter a number.")


def get_int(prompt, default, lo=1):
    while True:
        raw = input(prompt).strip()
        if raw == '':
            return default
        try:
            val = int(raw)
            if val >= lo:
                return val
            print(f"  Value must be ≥ {lo}.  Try again.")
        except ValueError:
            print("  Invalid input — please enter an integer.")


def run_flat():
    print()
    cu    = get_float("  Undrained shear strength  cu  [kPa, default 1.0]: ", 1.0, 1e-6)
    r     = get_float("  Roughness ratio           r   [0–1, default 1.0]: ", 1.0, 0.0, 1.0)
    B     = get_float("  Foundation half-width     B   [m,   default 1.0]: ", 1.0, 1e-6)
    n_fan = get_int  ("  Fan increments            n   [int, default  20]: ", 20)

    print("\n  Building slip line field …", end='', flush=True)
    fdn = FoundationBoundary(cu=cu, r=r, B=B)
    solver = PrandtlSolver(foundation=fdn, n_fan=n_fan)
    q = solver.solve()
    print(" done.\n")

    print(f"  ┌───────────────────────────────────────────┐")
    print(f"  │  Collapse load  q   = {q:>10.4f} kPa       │")
    print(f"  │  q / cu             = {q/cu:>10.4f}           │")
    print(f"  │  Prandtl exact      = {np.pi + 2:>10.4f}  (π+2)    │")
    print(f"  │  Fan angle          = {np.degrees(solver.delta_psi):>9.2f}°           │")
    print(f"  └───────────────────────────────────────────┘\n")

    plot_slip_lines(solver, q)


def run_disk():
    print()
    cu     = get_float("  Undrained shear strength  cu  [kPa, default 1.0]: ", 1.0, 1e-6)
    D      = get_float("  Disk diameter             D   [m,   default 1.0]: ", 1.0, 1e-6)
    p      = get_float("  Penetration ratio  p = z/D  [0 – 0.5, default 0.25]: ", 0.25, 1e-6, 0.5)
    r      = get_float("  Roughness ratio           r   [0–1, default 0.0]: ", 0.0, 0.0, 1.0)
    n_fan  = get_int  ("  Fan increments            n   [int, default  20]: ", 20)
    n_disk = get_int  ("  Disk surface nodes        m   [int, default  30]: ", 30)

    print("\n  Building slip line field …", end='', flush=True)
    disk = DiskFoundationBoundary(cu=cu, D=D, p=p, r=r)
    solver = DiskPrandtlSolver(disk=disk, n_fan=n_fan, n_disk=n_disk)
    result = solver.solve()
    print(" done.\n")

    print(f"  ┌────────────────────────────────────────────────┐")
    print(f"  │  p (penetration ratio)   = {p:>8.4f}             │")
    print(f"  │  B_s (contact half-width) = {result['B_s']:>8.4f} m           │")
    print(f"  │  θ_max                   = {result['theta_max_deg']:>8.2f}°            │")
    print(f"  │  Fan angle               = {result['fan_angle_deg']:>8.2f}°            │")
    print(f"  │                                                │")
    print(f"  │  σ_n / cu  at edge       = {result['sigma_n_at_edge']/cu:>8.4f}             │")
    print(f"  │  σ_n / cu  at bottom     = {result['sigma_n_at_bottom']/cu:>8.4f}             │")
    print(f"  │                                                │")
    print(f"  │  Q (per unit length)     = {result['Q_per_unit_length']:>8.4f} kN/m         │")
    print(f"  │  Q / D                   = {result['Q_over_D']:>8.4f} kPa           │")
    print(f"  │  Q / (D·cu)              = {result['Q_over_D_over_cu']:>8.4f}             │")
    print(f"  └────────────────────────────────────────────────┘\n")

    plot_disk_slip_lines(solver)


def main():
    print("=" * 64)
    print("  Slip Line Analysis — Tresca, weightless soil")
    print("=" * 64)
    print("    1) Flat strip footing")
    print("    2) Partly-penetrated circular disk")
    print()
    choice = input("  Select geometry [1/2, default 1]: ").strip() or '1'

    if choice == '2':
        run_disk()
    else:
        run_flat()


if __name__ == '__main__':
    main()
