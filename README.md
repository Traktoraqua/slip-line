# slip-line

Slip-line (method of characteristics) analysis of the bearing capacity of foundations, written in Python with numpy and matplotlib.

## Features

- Strip footing on weightless Tresca soil (Prandtl solution, q = cu·(π+2) for a smooth base)
- Strip footing on Mohr-Coulomb soil, with and without self-weight
- Partly-penetrated circular disk, for both Tresca and Mohr-Coulomb soil
- Plots of the slip-line fields

The theory and API are described in `manual.tex` (LaTeX).

## Setup

Requires [uv](https://docs.astral.sh/uv/). Dependencies are installed automatically on first run:

```
uv run main.py
```

## Usage

`main.py` is an interactive program (Tresca soil): choose a strip footing or a disk, enter the parameters and get the collapse load plus a plot.

The `test_*.py` and `plot_*.py` files in the root are standalone scripts that print result tables or save PNG figures, for example:

```
uv run test_mc.py
uv run test_mc_disk_parallel.py
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
- `slip_line/plot.py` – plotting

Conventions: y points downward (depth), compression is positive.
