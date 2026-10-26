# slip-line

Slip-line (method of characteristics) analysis of the bearing capacity of foundations, written in Python with numpy and matplotlib.

## Features

- Strip footing on weightless Tresca soil (Prandtl solution, q = cu·(π+2) for a smooth base)
- Strip footing on Mohr-Coulomb soil, with and without self-weight
- Partly-penetrated circular disk, for both Tresca and Mohr-Coulomb soil; the rough Tresca disk takes the interface shear of either a horizontal push or vertical penetration (`vertical_shear=True`)
- Velocity fields (Geiringer) on the Tresca nets, with checks of the plastic work and the virtual-work balance
- A pipe sinking into the soil under its own weight, with heave of the surface beside it
- A lateral push of the pipe while it is free to move vertically (upper bound), followed by a push back
- Plots of the slip-line fields

The theory and API are described in `manual.tex` (LaTeX).

## Setup

Requires [uv](https://docs.astral.sh/uv/). Dependencies are installed automatically on first run:

```
uv run main.py
```

## Usage

`main.py` is an interactive program (Tresca soil): choose a strip footing or a disk, enter the parameters and get the collapse load plus a plot.

The pipe simulations print their steps as tables:

```
uv run run_penetration.py   # pipe sinking under its own weight, smooth and r = 0.75
uv run run_lateral.py       # sink, push 0.2 m sideways, then 40 mm back
```

`run_lateral.py` has a `REVERSING` switch at the top: `True` (default) makes the net carry the interface friction, reversing where the slip along the pipe changes sign; `False` lets the net angle be free.

The `test_*.py` and `plot_*.py` files in the root are standalone scripts that print result tables or save PNG figures in `plots/`, for example:

```
uv run test_mc.py
uv run test_mc_disk_parallel.py
uv run plot_velocity_disk.py
```

The solvers can also be used directly:

```python
from slip_line.boundary import FoundationBoundary
from slip_line.solver import PrandtlSolver

solver = PrandtlSolver(FoundationBoundary(cu=1.0, r=0.0, B=1.0), n_fan=20)
print(solver.solve())  # ≈ π + 2
```

## Code layout

- `slip_line/yield_function.py` – Tresca and Mohr-Coulomb stress relations
- `slip_line/boundary.py` – foundation boundary conditions
- `slip_line/slip_line_field.py` – nodes and characteristic-net steps
- `slip_line/solver.py` – the solvers
- `slip_line/velocity.py` – velocity fields, plastic work and energy balance
- `slip_line/penetration.py` – pipe sinking under its own weight (`PipePenetration`)
- `slip_line/lateral.py` – lateral push of the pipe (`PipeLateral`)
- `slip_line/plot.py` – plotting

Conventions: y points downward (depth), compression is positive.
