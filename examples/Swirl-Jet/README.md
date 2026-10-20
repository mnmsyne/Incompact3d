# Coaxial Swirling Jet

This case simulates an axial core jet surrounded by a swirling annular jet on a Cartesian grid. The default `input.i3d` defines a full-size, cold-start research configuration.

## Default configuration

| Parameter | Setting |
| --- | --- |
| Case selection | `itype=733` (`swirljet`) |
| Grid | `641 × 400 × 400` points (102,560,000 total) |
| Domain | `32 × 12 × 12`, with the jet axis along x |
| MPI decomposition | `p_row=0`, `p_col=0`; selected automatically by the decomposition library |
| Radii | `rai=0.5`, `rao=1.0`; both must be finite and satisfy `0 < rai < rao` |
| Flow parameters | `re=1000`, `u1=0.5`, `u2=0.3` |
| Time integration | Third-order Adams–Bashforth, `dt=0.002`, steps 1–10000 |
| Initialization | `irestart=0`; no restart data required |
| Perturbations | `init_noise=0.02`, `inflow_noise=0`; `iin=1` uses a clock-based random seed |
| Boundary conditions | Prescribed inflow and convective outflow in x; periodic in y and z |
| Passive scalars | Two scalars, each with `Sc=1` and `cp=1` |
| Output | Snapshots every 2500 steps, checkpoints every 5000 steps, log diagnostics every 10 steps |
| Visualization sampling | `nvisu=2`, spatial coarsening of visualization output |
| Spatial discretization | Sixth-order compact first derivatives and hyperviscous second derivatives; `nu0nu=4`, `cnu=0.44` |
| LES model | `ilesmod=0` |

The nominal simulation duration is 20 nondimensional time units. `nvisu` controls spatial sampling, not the snapshot interval. The case-specific postprocessing hook is currently empty, so `ipost=1` does not provide integral swirl-number or momentum-budget diagnostics.

## Inlet model

The radial coordinate is measured from the center of the transverse domain:

```text
r = sqrt((y - Ly/2)^2 + (z - Lz/2)^2)
```

The mean velocity profile has three regions:

- **Core, `r <= rai`:** a parabolic axial jet, `ux = 2*(1 - r^2/rai^2)`, with zero mean transverse velocity. Its peak axial velocity is fixed at 2 and does not depend on `u1`.
- **Annulus, `rai < r <= rao`:** a laminar annular axial profile scaled by `u1`, plus an axial coflow of `u1/100`. The tangential velocity is `u_theta = u2*r/rao`.
- **Exterior, `r > rao`:** an axial coflow of `u1/100`, with zero mean transverse velocity.

The default ratio `u2/u1=0.6` is a ratio of input velocity scales, rather than a specified integral swirl number. Changing `u1` also changes the velocity ratio between the annulus and the fixed core jet. An integral swirl number requires a defined reference radius and momentum-flux convention, followed by integration of the inlet profile.

Odd-numbered scalars mark the core and even-numbered scalars mark the annulus. Their radial profiles use error-function transitions with a width parameter of 0.1. Initial velocity noise is concentrated near the jet axis; inlet noise is disabled in the default configuration.

The transverse boundaries are periodic. Domain-size and grid-sensitivity checks are needed when interpreting the solution as an isolated jet.

## Running the case

Build the solver with compatible MPI and 2DECOMP&FFT installations. From this case directory, specify the executable path and input file explicitly:

```sh
mpirun -np <num_ranks> /path/to/build/bin/xcompact3d input.i3d > run.log 2>&1
```

The default grid is intended for a suitably provisioned parallel computing environment. Use a separate reduced configuration for startup checks before running the full case.
