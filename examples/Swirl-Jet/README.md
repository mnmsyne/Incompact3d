# Coaxial Swirling Jet

This case simulates an axial core jet surrounded by a swirling annular jet on a Cartesian grid. The default `input.i3d` defines a full-size, cold-start analytical configuration. `input_LES.i3d` provides a separate file-driven LES setup. `input_smoke.i3d` provides a reduced eight-rank integration test.

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

## Selecting the velocity inlet

For this case, `iin=0`, `1`, or `2` retains the analytical velocity profile. Values 1 and 2 retain clock-seeded and fixed-seed initial perturbations, respectively. Use **`iin=4`** for a steady Cartesian table containing all three nondimensional velocity components. `iin=3` remains reserved for the solver's precursor workflow and is not supported by this case.

The file mode reads `swirl_inlet.dat` from the existing `inflowpath` directory in `&InOutParam` (an empty path means the working directory). Rank 0 reads and validates the file, broadcasts it, and each rank bilinearly interpolates it onto its local X-pencil inlet coordinates, including stretched y coordinates. The local plane is cached; restart runs load it on their first boundary update. Coordinates outside the table cause a nonzero abort, rather than extrapolation. Table radii must match `rai` and `rao`.

The solver reads all three supplied Cartesian components without reflection, rescaling, or rotation reversal. Preprocessing can retain selected cylindrical modes, as described below. File mode requires `inflow_noise=0`. Optional `init_noise` perturbs only the initial volume; the steady boundary itself remains the supplied plane. A cold start extrudes that plane along x before the solver's pressure correction. Passive scalars retain their geometric analytical profiles because no scalar data are supplied.

The analytical `u1` and `u2` parameters do not rescale file velocities. In file mode, velocity and scalar outflow convection use the maximum axial speed on the penultimate plane, independent of those parameters.

The table defines the case's `bxx1/bxy1/bxz1` boundary targets. The common pressure-correction algorithm adds the previous transverse pressure-gradient contribution before projection and subtracts the current contribution afterwards. Consequently, snapshot `uy/uz` at the inlet can differ slightly from the table even when its reading, broadcast and interpolation are correct; compare the three targets before this correction, and inspect projected divergence separately. Do not overwrite the projected transverse velocities simply to force a table match.

## Preparing a Cartesian inlet

Install NumPy and SciPy in the preprocessing Python environment, then run from this directory:

```sh
python prepare_inlet.py --csv /path/to/plane.csv --input input_LES.i3d --output swirl_inlet.dat
```

The CSV must contain a `[Data]` section with x/y/z in metres and Cartesian velocity u/v/w in metres per second. Defaults are inner radius 0.0075 m, outer radius 0.015 m, viscosity 0.889e-6 m2/s, and a transverse center at (0, 0). These quantities can be changed through the script's command-line options.

The script computes the center-disk bulk velocity by area integration at the sampled plane and uses the inner diameter as the length scale. For the analyzed plane at x=2 mm, the scales are:

- `Di = 0.015 m`, `Ub = 0.441409596 m/s`.
- `Re = Ub*Di/nu = 7447.856`; `rai*=0.5`, `rao*=1`.
- Computational inlet x=0 represents the supplied physical plane at x=2 mm, rather than the original nozzle face.

Here `D=2*Ri=Ro=0.015 m`; the outer-radius and inner-diameter reference lengths coincide. The small difference from previous rounded estimates of bulk speed/Re comes from converged area integration of the actual scattered data.

Scattered data first use linear triangulation. By default, `--reconstruction fourier` samples cylindrical `ux,ur,utheta` on 720 angles and radial spacing at most one quarter of the transverse mesh width. It retains **m=0,12,24**, applies a Gaussian to their complex radial coefficients, and uses cubic splines for radial evaluation. The default `--smooth-sigma 0.25` is 0.25 of the minimum transverse cell width (0.078125 mm here). It does not force zero velocities at either nominal nozzle radius. Nonaxisymmetric spatial lobes are not a time-dependent turbulence model.

The convention is `theta=atan2(y,z)`, `ur=uy*sin(theta)+uz*cos(theta)`, `utheta=uy*cos(theta)-uz*sin(theta)`. Reconstruction uses the inverse of these formulas; the measured swirl sign is retained. Near the axis, a Cartesian Gaussian patch joins to the polar reconstruction over `r/D=0.1–0.2` with a quintic weight. This preserves a finite, continuous Cartesian center without imposing artificial cylindrical values at `r=0`.

A second quintic weight joins the distant flow to a common axial background over `r/D=3–3.75`. The constant exterior makes both periodic pairs continuous with flat neighborhoods. The axial background is derived by matching the unprocessed Cartesian rectangle's axial flux; the measured mean in the fringe is also recorded. For this CSV these are about `0.01717 Ub` and `0.01999 Ub`, respectively. Transverse background velocities are zero. Core and annular velocities are never rescaled to compensate for smoothing.

To compare with all-mode Cartesian smoothing, use `--reconstruction cartesian`. To generate the original unsmoothed interpolation as a reference:

```sh
python prepare_inlet.py --csv /path/to/plane.csv --input input_LES.i3d \
  --output raw_inlet.dat --reconstruction cartesian --smooth-sigma 0 --no-periodic-fringe
```

The table extent follows the input configuration. Both endpoints are included for periodic directions, giving a **385×385** table for the 384×384 mesh; the solver omits the duplicate endpoint. No extrapolation or nearest-neighbor filling is permitted. The source CSV is unchanged. The polar option also requires coverage of its sampling disk. The reader's bilinear interpolation is continuous in value; it does not guarantee globally continuous derivatives between every table cell.

A JSON sidecar records source/configuration/script/table hashes, scales, quadrature convergence, conditioning settings, component extrema, core flux error, periodic endpoint differences, and cylindrical mode energies/phase errors before and after processing. Regenerate after changing the configuration, data, geometry, or preprocessing options. Generated tables and metadata remain ignored by Git.

For this CSV, pure m=0/12/24 projection gives annular cylindrical R² of approximately 0.99927, 0.99905, 0.99921, compared with 0.85505, 0.10678, 0.82912 for pure m=0. These are area-weighted annular comparisons, not claims about every local point. The final grid retains 97.0–98.9% of m=12 amplitude and 93.9–95.6% of m=24 amplitude, with phase changes below 0.22 degrees (modal phase, not lobe orientation). Its core flux changes by -0.199%; inside `r/D<=1`, axial momentum and signed angular momentum change by approximately -0.391% and -0.585%. No-pressure swirl number changes from -0.39983 to -0.39905 using the above convention and reference radius 1.

The thin inner mixing layer remains resolution-sensitive: the azimuthal axial minimum rises from about 0.272 to 0.309 Ub; raw grid interpolation alone gives 0.300 Ub. The outer velocity drop also broadens on this mesh. Core transverse velocities are small, so relative transverse errors and relative core angular-momentum errors can be large even when absolute errors are small. Inspect local errors as well as integral agreement; use the proposed finer mesh for sensitivity tests.

## LES configuration

`input_LES.i3d` uses `iin=4`, `re=7450`, **iSVV** (`ilesmod=4`, `jles=4`), and the hyperviscous compact operator (`isecondder=5`, `nu0nu=4`, `cnu=0.44`). In this version `ilesmod` is a nonzero enable flag; `jles` selects the model. The coefficients are reference settings, not validated jet-specific optima. `jles=4` adds no explicit SGS stress. Hyperviscosity is controlled by the derivative operator even when the LES flag is zero; that combination should not be described as strict DNS.

| Parameter | File-driven LES setting |
| --- | --- |
| Domain | `24 × 8 × 8`, or `360 × 120 × 120 mm` |
| Mesh | `513 × 384 × 384` (75,644,928 points) |
| Spacing | `dx=0.046875 D`, `dy=dz=0.0208333 D`, or `0.703125 × 0.3125 × 0.3125 mm` |
| Reference | `Ub=0.441409596 m/s`, measured `Re=7447.856`; configured Re differs by 0.029% |
| Time integration | AB3, `dt=0.001`, approximately 0.03398 ms |
| Duration | 10000 steps = 10 reference time units, less than one `Lx/Ub` flow-through |

The initial duration is for startup inspection; extend it deliberately for downstream development and statistics. Output intervals remain specified in steps. The original analytical `input.i3d` is unchanged.

```sh
mpirun -np <num_ranks> /path/to/build/bin/xcompact3d input_LES.i3d > run.log 2>&1
```

This is a steady inlet plane. A complete solver run has not validated the model or grid. Check pressure projection, divergence, flow balance, runtime CFL, dissipation, scalar extrema/budgets, and statistical convergence. The inlet convective CFL estimate is about 0.0715; this is not a runtime stability proof. Suggested comparisons are `769×512×512` at the same extent and a wider transverse domain. Actual CFL and boundary stability determine the usable time step. ADIOS-specific output registration is outside this change.

## Eight-rank smoke test

`input_smoke.i3d` uses `65×72×72` points, an `8×8×8` domain, `Re=1000`, `dt=0.001`, 100 steps, iSVV, two passive scalars, and a fixed `2×4` MPI decomposition. Initial noise is disabled for reproducibility. The transverse extent matches the LES configuration, so reuse its generated **385×385** inlet table instead of smoothing the data again on the coarse smoke mesh. The 72-point transverse mesh also exercises interpolation between table nodes.

```sh
python prepare_inlet.py --csv /path/to/plane.csv --input input_LES.i3d --output swirl_inlet.dat
mpirun -np 8 /path/to/build/bin/xcompact3d input_smoke.i3d > smoke.log 2>&1
```

Use a scheduler allocation on an HPC system. Keep each run in a separate directory with its input and a valid `inflowpath`. This test writes snapshots and checkpoints every 50 steps. Check normal completion, finite fields, all three boundary targets against the table, scalar extrema, projection divergence and restart behavior. Apply the pressure-correction distinction above when inspecting snapshots. A low-Re coarse-grid smoke is not an LES or mixing-layer accuracy test. Deployment scripts, job IDs, raw logs and inspection results belong in the local `.hpc_log/` directory.

## Table format

`SWIRL_INLET_V1` is an ASCII format:

```text
SWIRL_INLET_V1
Ny_table Nz_table
rai_star rao_star
y_1 ... y_Ny
z_1 ... z_Nz
ux uy uz
... one triple per point, y varying fastest and z outermost ...
```

Coordinates are centered and nondimensional; the solver applies its existing yly/2 and zlz/2 centering convention. Axes must be finite and strictly increasing; all velocities must be finite. Each numeric record must contain exactly the required number of whitespace-separated values. Null fields, slash termination, repeat counts, extra values, invalid headers/dimensions/radii, truncated data, and trailing records are rejected. Values are compatible with both single- and double-precision builds.

## Passive scalar settings

Both configurations transport two geometric labels: scalar 1 marks the core, and scalar 2 marks the annulus. Their inlet and initial profiles retain the analytical error-function transitions with nondimensional width parameter 0.1 (1.5 mm with the current length scale). With both `cp=1`, their sum is a smoothed outer-disk indicator, not a constraint that the two scalars sum to one everywhere. Additional scalar pairs repeat these odd/even profiles; `numscalar=0` disables scalar transport.

The CSV provides velocity only. The file-driven case therefore assigns geometric labels at the supplied plane; these are not measured mixture fractions carried from the upstream nozzle. Use scalar data if physical mixing history at x=2 mm is required.

In the LES setup, `Sc(1)=Sc(2)=1`, `Ri(1)=Ri(2)=0` and `uset(1)=uset(2)=0`: no buoyancy coupling or settling is applied. The scalar equation uses `nu/Sc` with the same hyperviscous second-derivative choice as momentum. With `jles=4` no explicit `nut/Sc` SGS scalar term is added; this is not a WALE scalar model. There is no independently specified turbulent Schmidt number here.

Scalar x boundaries use prescribed inflow and the case's convective outflow, while y/z remain periodic. Scalar convection retains its existing non-skew setting. No hard [0,1] clipping is enabled: monitor minima/maxima and scalar budgets, since high-order transport can overshoot. Clipping would alter scalar conservation and is not introduced by this inlet change.
