# 2D2V Wigner-Poisson Adaptive-Rank Solver

This repository contains the 2D2V MATLAB code used with the working manuscript
`Wigner_Poisson_Adaptive_Rank_highD`. The solver advances the
Wigner-Poisson equation in phase space `(x, vx, y, vy)` using Strang splitting,
semi-Lagrangian adaptive-rank transport, HTACA compression, FFT-based Poisson
solves, and conservation corrections.

## Entry Points

- `WP_str_SLphi3.m` runs the default strong-Landau test with the leading-mode adaptive-weight
  conservation correction.
- `QiuWP_str_SLphi1.m` runs the default two-stream benchmark with the spatially
  uniform conservation correction.

Each script sets the problem, mesh size, final time, CFL number, HT tree, rank
bounds, and output filename near the top of the file.

## Repository Layout

- `SetUp Initialization/`: grid construction and initial conditions.
- `SG/`: semi-Lagrangian spatial transport, WENO interpolation, and clamped
  Wigner force updates.
- `Methed Related/`: Poisson solves, moment diagnostics, and conservation
  corrections. The folder name is kept as-is for path compatibility with the
  MATLAB drivers.
- `ht_toolbox/`: Hierarchical Tucker tensor utilities used by the solver.

## Running

Open MATLAB from the repository root and run one of the driver scripts:

```matlab
WP_str_SLphi3
QiuWP_str_SLphi1
```

The scripts add the required helper folders to the MATLAB path and save run data
to `.mat` files named from the selected problem, mesh size, CFL, final time, and
`H` parameter.

## Problem Selectors

- `prob = 1`: weak Landau damping
- `prob = 2`: strong Landau damping
- `prob = 4`: two-stream benchmark
- `prob = 5, 6, 7`: bump-on-tail variants

The two main drivers use only the selector ranges documented in their headers.
