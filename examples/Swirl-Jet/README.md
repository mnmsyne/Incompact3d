# Coaxial Swirling Jet

This case models an axial core jet surrounded by a swirling annular jet on a Cartesian grid. The jet axis is x and the inlet is centered in the transverse plane. It uses prescribed inflow and convective outflow in x, with periodic boundaries in y and z. The case is selected by `itype=733`.

## Flow and inlet model

The radii `rai` and `rao` define the core and outer annular edge, with `0 < rai < rao`. Two velocity inlet models are available:

- **Analytical inlet:** a parabolic axial core, an annular axial profile with tangential velocity proportional to radius, and a weak axial coflow outside the jet. `u1` scales the annular axial flow and coflow; `u2` scales rotation. The core peak is fixed by the profile rather than scaled by `u1`, and `u2/u1` is not an integral swirl number.
- **File-driven inlet:** a steady plane of all three Cartesian velocity components reconstructed from CSV data. It retains the supplied rotation sign and spatial nonuniformity; the analytical velocity scales do not rescale it.

For the CSV workflow, length is normalized by the inner diameter and velocity by the area-averaged core axial velocity. The preprocessing metadata records the dimensional scales and corresponding Reynolds number. The computational inlet represents the supplied measurement plane, which may be downstream of the nozzle.

Periodic transverse boundaries approximate an isolated jet only when the domain is sufficiently wide. Domain size, mixing-layer resolution and downstream development require separate assessment.

## Inputs and execution

- [input.i3d](input.i3d): analytical inlet configuration.
- [input_LES.i3d](input_LES.i3d): file-driven LES configuration.
- [input_smoke.i3d](input_smoke.i3d): reduced integration check.

The input files specify the mesh, flow parameters, time stepping and output settings. Check them before each run; numerical settings and historical test results are not duplicated here.

Build the solver following the repository instructions, then run from a directory containing the selected input and its inlet data:

```sh
mpirun -np <num_ranks> /path/to/build/bin/xcompact3d input_LES.i3d > run.log 2>&1
```

Use the site's scheduler on HPC systems. Keep separate run directories with their inputs and inlet metadata.

## Numerical method

The supplied LES setup uses sixth-order compact first derivatives, hyperviscous compact second derivatives for iSVV, pressure projection and third-order Adams–Bashforth time integration. iSVV provides dissipation through the derivative operator; it does not add an explicit SGS eddy-viscosity stress or scalar flux. The scalar equation scales diffusion with the molecular Schmidt number and uses the same second-derivative choice.

High-order accuracy does not ensure scalar boundedness or adequate resolution of thin mixing layers. A reduced integration run checks execution and data handling; it does not establish LES accuracy or grid convergence.

## CSV velocity inlet

With NumPy and SciPy installed, prepare the inlet from a CSV containing a `[Data]` section with Cartesian coordinates in metres and velocity components in metres per second:

```sh
python prepare_inlet.py --csv /path/to/plane.csv --input input_LES.i3d --output swirl_inlet.dat
python prepare_inlet.py --help
```

The script nondimensionalizes the data, reconstructs the inlet and writes a `SWIRL_INLET_V1` table plus JSON metadata. The default reconstruction retains selected cylindrical Fourier modes and smooths their radial coefficients; Cartesian reconstruction is also available. A regular Cartesian patch handles the axis, and a far-field join makes the periodic edges compatible. These spatial modes do not supply time-dependent turbulent inflow.

Check geometry, data coverage and preprocessing options for each dataset. Metadata records source hashes, scales and reconstruction diagnostics. Regenerate when changing the source data, physical domain, radii or intended preprocessing. Generated files are ignored by Git.

For the current grid-sensitivity study, prepare the table once from the 601x400x400 baseline in `input_LES.i3d`, then copy that exact table and its metadata to every run directory. The transverse table has 401x401 points, including both periodic endpoints. Reuse it for 801x400x400 and 601x480x480: the solver interpolates the same processed velocity field onto each mesh. Do not rerun preprocessing for each mesh, because the smoothing width depends on the transverse grid spacing. Compare table SHA-256 values before running.

The current baseline is a fresh 5000-step run with `dt=0.001` and checkpoint interval 5000. All three meshes use the same physical and numerical settings and inlet table. For `iin=4`, initial random perturbations use a fixed seed base of 42 while retaining rank-dependent seeds. `input_LES.i3d` fixes `p_row=6` and `p_col=8` for 48 MPI ranks. Repeating the initial random field requires the same mesh, decomposition and executable/runtime; changing the mesh or decomposition changes its spatial mapping. Bitwise reproducibility after startup projection and time integration still requires a repeated-run comparison.

Use `iin=4` to read the table from `inflowpath` (or the working directory when empty), with inlet noise disabled. Rank zero reads the table, broadcasts it, and each rank bilinearly interpolates and caches its inlet plane. Table coordinates must cover the inlet and its radii must match the case; extrapolation is rejected. Analytical inlet modes are `iin=0`, `1` and `2`; the precursor mode is not supported by this case.

The table contains centered, nondimensional y/z axes followed by velocity triples, with y varying fastest:

```text
SWIRL_INLET_V1
Ny_table Nz_table
rai_star rao_star
y_1 ... y_Ny
z_1 ... z_Nz
ux uy uz
...
```

## Initialization and scalars

File-driven cold starts project the initial velocity using the existing pressure-correction flow, retain the projected velocity, and clear the artificial pressure and gradient histories. Pressure correction can leave the inlet transverse velocities slightly different from the table targets.

The CSV supplies velocity only. Scalar 1 marks the core and scalar 2 marks the annulus, using smooth analytical error-function transitions at the radial interfaces. They are geometric labels rather than measured mixture fractions and need not sum to one everywhere. Scalar inflow is prescribed, outflow is convective, and transverse boundaries are periodic. Scalar convection retains the existing non-skew form; no scalar clipping is applied.

## Diagnostics

Enable scalar postprocessing with `ipost`; `iprocessing` sets its sampling interval. Rank zero writes:

- `swirl_scalar_bounds.csv`: raw extrema, volume fractions below -0.001 or above 1.001, and integrals of the raw excess below 0 or above 1. The integrals are not thresholded or divided by domain volume.
- `swirl_scalar_divergence.csv`: full-domain and interior scalar-operator divergence statistics every tenth bounds sample. The optional `interior` rows omit four x planes at each end; y/z remain fully periodic. Use the `all` rows for the current grid comparisons.

Scalar-operator divergence differs from the divergence constrained by pressure projection. These diagnostics do not change the fields. Start fresh diagnostic files when their column schema changes. Primary fields and checkpoints follow the input settings; `nvisu` controls spatial sampling of visualization output, not its time interval. Use full-resolution restart data for small-scale analysis.
