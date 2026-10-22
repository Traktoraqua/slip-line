# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Slip-line (method of characteristics) analysis of bearing capacity for foundations on weightless/gravity soil. Python + numpy + matplotlib. `manual.tex` is the LaTeX user manual (theory, API, figures) — consult it for the derivations behind the solvers.

Managed with `uv` (Python 3.13; numpy and matplotlib declared in `pyproject.toml`). There is no linter config or pytest suite yet.

## Commands

- Interactive CLI (strip footing or partly-penetrated disk, Tresca only): `uv run main.py`
- The `test_*.py` files in the repo root are **not pytest tests** — they are standalone scripts that print tables and/or save PNGs (e.g. `test_mc.py` writes `plot_mc_phi*.png`). Run one with `uv run test_mc_disk_parallel.py`.
- `plot_*.py` scripts (and some `test_*.py`) regenerate the PNG figures in `plots/` (they use the `Agg` backend, so they save rather than show).
- `dbg.py` dumps Region I mesh nodes of a disk solve to `dbg.out`.
- `uv run run_penetration.py` simulates a smooth pipe sinking under its own weight with heave and prints the steps.
- Format/lint: `ruff format`, `ruff check`.

## Architecture

Package `slip_line/`:

- `yield_function.py` — `TrescaYield`, `MohrCoulombYield`: boundary stress states and Hencky/Mohr-circle relations.
- `boundary.py` — `FoundationBoundary` (strip, roughness `r`, half-width `B`) and `DiskFoundationBoundary` (diameter `D`, penetration ratio `p = z/D`, roughness `r`).
- `slip_line_field.py` — `Node(x, y, sigma, psi)` and the characteristic-net primitives: `alpha_slope`/`beta_slope`, `riemann_step` (Tresca), `*_mc` variants, and `riemann_step_mc_gravity`.
- `solver.py` (~2000 lines) — the core. Six solver classes, each building the field region by region (`_build_region_III` → `_build_region_II` (centred fan) → `_build_region_I`) and exposing `all_alpha_lines()` / `all_beta_lines()` for plotting:
  - `PrandtlSolver` — Tresca strip; `solve()` returns collapse load `q` (exact check: `cu·(π+2)` for smooth).
  - `MohrCoulombPrandtlSolver` — c-φ strip, weightless; checked against `analytical_Nc()`.
  - `MohrCoulombGravityPrandtlSolver` — c-φ strip with self-weight (γ>0), `solve(fit_foundation=...)`, `contact_pressure_curve()`.
  - `DiskPrandtlSolver` / `MohrCoulombDiskPrandtlSolver` — partly-penetrated disk (Tresca / MC). `solve()` returns a **dict** (`B_s`, `theta_max_deg`, `fan_angle_deg`, `Q_per_unit_length`, `Q_over_D`, `sigma_n_at_*`, …), and they provide `forces_on_half_disk()`. Both have `_resize_fan_and_wedge_to_match_region_I()`: Region I is built from the disk surface first, then the fan and Region III are resized to match it. The MC disk solver has a "parallel-β zone" mechanism for deep penetration/high roughness (`parallel_beta_zone` in the result).
- `velocity.py` — Geiringer velocity field on a Tresca net: `geiringer_step`, `PrandtlVelocityField` (strip) and `DiskVelocityField` (disk), both pushed down at `v0` (the disk also horizontally toward +x with `horizontal=True`, right half only, U = 0; for r > 0 `interface_slip()` gives the slip work against the interface shear, which is in `energy_balance()`); the outer α-line is the velocity discontinuity to rigid soil. `plastic_work_rate` / `DiskVelocityField.energy_balance()` check the sign of the plastic work and the virtual-work balance. `plot_velocity_strip.py` / `plot_velocity_disk.py` plot them with a hodograph.
- `penetration.py` — `PipePenetration`: quasi-static self-weight sinking of a smooth pipe in small steps `dw`. Each step uses a flat-surface `DiskPrandtlSolver` at the effective embedment (invert below the contact level), raises the surface by the Geiringer upward surface velocity × `dw` (vertical only), scales that increment so the soil volume is conserved, and stops when the collapse load reaches the weight `W`.
- `plot.py` — matplotlib plotting (`plot_slip_lines`, `plot_disk_slip_lines`, `plot_mc_slip_lines`, …), one per solver type.

Conventions that matter across files:

- **y points downward (depth)**; angles are measured CCW in (x, y_down) space, which looks clockwise on screen. `plot.py` flips y for display.
- **Compression positive.** `ψ` is the angle of the major principal stress σ₁; `σ` is mean stress. Along α-lines `σ − 2·cu·ψ = const`, along β-lines `σ + 2·cu·ψ = const` (Tresca).
- Only the right half (x ≥ 0) is computed; plotting mirrors it.

## Repo notes

- Generated PNGs are tracked in `plots/`. `Pipestab.numbers` is an unrelated, untracked file.
