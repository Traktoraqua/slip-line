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
- `uv run run_penetration.py` simulates a pipe sinking under its own weight with heave (smooth, and r = 0.75) and prints the steps.
- `uv run run_lateral.py` sinks a rough pipe (r = 0.8, W = 2.8 kN/m) and then pushes it 0.2 m sideways, free to move vertically (~4 min).
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
  - `DiskPrandtlSolver` / `MohrCoulombDiskPrandtlSolver` — partly-penetrated disk (Tresca / MC). `solve()` returns a **dict** (`B_s`, `theta_max_deg`, `fan_angle_deg`, `Q_per_unit_length`, `Q_over_D`, `sigma_n_at_*`, …), and they provide `forces_on_half_disk()`. Both have `_resize_fan_and_wedge_to_match_region_I()`: Region I is built from the disk surface first, then the fan and Region III are resized to match it. The MC disk solver has a "parallel-β zone" mechanism for deep penetration/high roughness (`parallel_beta_zone` in the result). The Tresca disk's interface shear has the direction of a horizontal push by default (ψ_disk = π/2 − φ − ω, so `collapse_load_Q()` ≠ `2·forces_on_half_disk()['F_v']` for r > 0); `vertical_shear=True` gives that of vertical penetration (ω changes sign, friction on the disk acts upward, and the two agree).
- `velocity.py` — Geiringer velocity field on a Tresca net: `geiringer_step`, `PrandtlVelocityField` (strip) and `DiskVelocityField` (disk), both pushed down at `v0` (the disk also horizontally toward +x with `horizontal=True`, right half only, U = 0; for r > 0 `interface_slip()` gives the slip work against the interface shear, which is in `energy_balance()`); the outer α-line is the velocity discontinuity to rigid soil. `plastic_work_rate` / `DiskVelocityField.energy_balance()` check the sign of the plastic work and the virtual-work balance (keys `F_stress`, `F_virtual_work`, `F_upper_bound`, `negative_area_fraction`, `negative_interface_fraction`). For the vertical disk d < 0 in ~30 % of the plastic area even for r = 0, so the field is not strictly kinematically admissible; `F_virtual_work` is ~1 % below `F_stress` at `n_fan=18, n_disk=30` and converges to it with a finer mesh. `plot_velocity_strip.py` / `plot_velocity_disk.py` plot them with a hodograph.
- `penetration.py` — `PipePenetration`: quasi-static self-weight sinking of a pipe (roughness `r`) in small steps `dw`. Each step uses a flat-surface `DiskPrandtlSolver(vertical_shear=True)` at the effective embedment (invert below the contact level), raises the surface by the Geiringer upward surface velocity × `dw` (vertical only), scales that increment so the soil volume is conserved, and stops when the collapse load reaches the weight `W`.
- `lateral.py` — `PipeLateral`: lateral push of a pipe that is free to move vertically (V = W constant), starting from a `PipePenetration` end state (`after_penetration`). Kinematic (upper bound): each step picks the direction β (below horizontal) and net angle ω that minimise H = (D − W·sin β)/cos β, where D = Σ|d|·A + cu·U·L + r·cu·∫|slip| on the Geiringer field of a `OneSidedDiskSolver` net. That net's contact chain runs from the trailing separation φ = −β (or the left contact angle on the original surface) to the right edge, all soil flows forward, and the free surface in front slopes at the fitted berm angle α. The trench behind stays open. The berm rises by the normal surface velocity, and the soil the pipe moves into is removed (volume-scaled). H starts slightly negative because this mechanism is a few % weaker than step 1's flat-surface stress solution.
- `plot.py` — matplotlib plotting (`plot_slip_lines`, `plot_disk_slip_lines`, `plot_mc_slip_lines`, …), one per solver type.

Conventions that matter across files:

- **y points downward (depth)**; angles are measured CCW in (x, y_down) space, which looks clockwise on screen. `plot.py` flips y for display.
- **Compression positive.** `ψ` is the angle of the major principal stress σ₁; `σ` is mean stress. Along α-lines `σ − 2·cu·ψ = const`, along β-lines `σ + 2·cu·ψ = const` (Tresca).
- Only the right half (x ≥ 0) is computed; plotting mirrors it.
- **Signed roughness in requests:** the user writes a negative `r` (e.g. r = −0.75) for the interface shear direction of vertical penetration. `DiskFoundationBoundary` only accepts 0 ≤ r ≤ 1, so this means `r=abs(r)` with `DiskPrandtlSolver(..., vertical_shear=True)`.
- **Contact stresses are shown as stresses from the pipe on the soil** (the user is a geotechnical engineer): σ_n as compression into the soil, and for penetration τ on the soil pointing toward the invert. `stress_distribution_curve()` returns σ_n and |τ|; `forces_on_half_disk()` gives forces on the disk; `interface_slip()["tau"]` is the signed shear on the soil along t = (cos φ, −sin φ).

## Repo notes

- Generated PNGs are tracked in `plots/`. `Pipestab.numbers` is an unrelated, untracked file.
- The local changes to `plots/mc_gravity_contact_pressure.png` and `plots/mc_gravity_phi30.png` are unrelated: never commit them.
