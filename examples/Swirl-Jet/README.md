# Coaxial Swirling Jet

This case models an axial core jet surrounded by a swirling annular jet on a Cartesian grid. The jet axis is x and the inlet is centered in the transverse plane. It uses prescribed inflow and convective outflow in x, with periodic boundaries in y and z.

## Flow and inlet model

The radii `rai` and `rao` define the core and outer annular edge, with `0 < rai < rao`. Two velocity inlet models are available:

- **Analytical inlet:** a parabolic axial core, an annular axial profile with tangential velocity proportional to radius, and a weak axial coflow outside the jet. `u1` scales the annular axial flow and coflow; `u2` scales rotation. The core peak is fixed by the profile rather than scaled by `u1`, and `u2/u1` is not an integral swirl number.
- **File-driven inlet:** a steady plane of all three Cartesian velocity components reconstructed from CSV data. It retains the supplied rotation sign and spatial nonuniformity; the analytical velocity scales do not rescale it.

For the CSV workflow, length is normalized by the inner diameter and velocity by the area-averaged core axial velocity. The preprocessing metadata records the dimensional scales and corresponding Reynolds number. The computational inlet represents the supplied measurement plane, which may be downstream of the nozzle.

Periodic transverse boundaries approximate an isolated jet only when the domain is sufficiently wide. Domain size, mixing-layer resolution and downstream development require separate assessment.

## Inputs and execution

- [input.i3d](input.i3d) uses the analytical inlet.
- [input_LES.i3d](input_LES.i3d) uses the file-driven LES inlet.

Read the selected input for mesh, domain, time stepping, MPI decomposition and output settings. For a reduced check, generate a small-grid input and a matching inlet.

Build the solver following the repository instructions. Run from a directory containing the selected input and any inlet data it requires. Set `input` to the input filename and choose `ranks` to match its decomposition.

```sh
mpirun -np "$ranks" /path/to/build/bin/xcompact3d "$input" > run.log 2>&1
```

Use the site's scheduler on HPC systems. Keep separate run directories with their inputs and inlet metadata.

## Numerical method

iSVV adds dissipation through the compact second-derivative operator. The scalar equation scales diffusion with the molecular Schmidt number and uses the selected second-derivative operator. Consult the case input for the discretization, time integrator and LES model settings.

High-order accuracy does not ensure scalar boundedness or adequate resolution of thin mixing layers. A reduced integration run checks execution and data handling; it does not establish LES accuracy or grid convergence.

## CSV velocity inlet

With NumPy and SciPy installed, prepare the inlet from a CSV containing a `[Data]` section with Cartesian coordinates in metres and velocity components in metres per second:

```sh
python prepare_inlet.py --csv /path/to/plane.csv --input input_LES.i3d --output swirl_inlet.dat
python prepare_inlet.py --help
```

The script nondimensionalizes the data, reconstructs the inlet and writes a `SWIRL_INLET_V1` table plus JSON metadata. The default reconstruction retains selected cylindrical Fourier modes and smooths their radial coefficients; Cartesian reconstruction is also available. A regular Cartesian patch handles the axis, and a far-field join makes the periodic edges compatible. These spatial modes do not supply time-dependent turbulent inflow.

Check geometry, data coverage and preprocessing options for each dataset. Metadata records source hashes, scales and reconstruction diagnostics. Regenerate when changing the source data, physical domain, radii or intended preprocessing. Generated files are ignored by Git.

For a matched grid comparison, prepare the table once for the chosen baseline and use identical table bytes and metadata on every mesh with the same domain and radii. The solver interpolates the table onto each mesh. Preprocessing again changes the physical smoothing width when transverse grid spacing changes. Verify the table SHA-256 when staging runs.

`prepare_inlet.py` derives the transverse table size from the input mesh and includes both end coordinates in each periodic direction. Use the generated metadata to identify the table and preprocessing settings for a run.

File-driven cold starts use fixed, rank-dependent seeds for initial random perturbations. Repeating the initial random field requires the same mesh, decomposition and executable/runtime. Check bitwise reproducibility after startup projection and time integration with a repeated run.

The file-driven mode reads the table from `inflowpath` or the working directory when that path is empty. It requires inlet noise to be disabled. Rank zero reads and broadcasts the table. Each rank interpolates and caches its inlet plane. Table coordinates must cover the inlet and its radii must match the case. See the example inputs for the inlet selection.

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

The CSV supplies velocity only. Scalar 1 marks the core and scalar 2 marks the annulus, using smooth analytical error-function transitions at the radial interfaces. Each is a geometric label with its own initial and inlet profile. Compare the same scalar index across numerical schemes.

Set `scalar_scheme` per tracer in `ScalarParam`. The startup summary shows the resolved scheme for each tracer and, when conservative transport is active, the FCT interval and the selected mapping backend.

| Scheme | Scalar update |
| --- | --- |
| `conv` | Original advective form without scalar clipping. |
| `div` | Unlimited conservative flux form for comparison with FCT. |
| `fct` | Conservative flux form that limits internal and outlet antidiffusive fluxes on the complete RK3 step. |

DIV and FCT require explicit native RK3, a uniform grid, sixth-order compact D1 and iSVV D2, and open x with periodic y and z for velocity and scalars. The prescribed inlet node is external to the active scalar control volumes. For a strict unit-interval test, explicitly set `scalar_fct_lower` and `scalar_fct_upper` to the endpoints of `[0,1]` and supply admissible initial and inlet data.

### Auxiliary face velocities

CG and Fourier solve the same auxiliary face-velocity mapping problem for DIV and FCT. The correction acts on the internal scalar faces and preserves the physical mass flux at each x boundary. Both backends require compatible net physical flux and abort on physical inlet or outlet backflow beyond `scalar_mapping_tol`. The `scalar_mapping_solver` option selects the backend. Fourier requires an Intel Fortran build linked against the CMake target `MKL::MKL`, with `X3D_SCALAR_MAPPING_MKL` enabled in the compile command.

Set `scalar_mapping_reference_velocity` to a positive velocity scale for the relative correction diagnostics. Set `scalar_mapping_max_correction` to the justified maximum absolute component correction at any internal face, in solver velocity units. Exceeding this budget aborts the run. The reference scale and diagnostic velocity floor do not relax it.

### Physical outlet flux

DIV integrates the physical outlet advective flux using each RK3 stage's outward mass velocity and outlet scalar. FCT uses the same RK integral as its high-order candidate. It limits the correction from the monotone baseline using the advected scalar interval and the last control volume's capacity. The outlet and adjacent internal faces share that capacity. Both schemes use the accepted face integral in the scalar update and physical mass balance. The inlet retains its prescribed scalar and step-initial diffusion closure.

To assess outlet time accuracy, advect the same smooth bounded scalar packet through the outlet on a fixed mesh. Compare the outlet flux history, cumulative discharged mass, and near-outlet field at common times as the time step is refined. Record limiter activity alongside the errors.

## Diagnostics

Enable scalar postprocessing with `ipost`; `iprocessing` sets its sampling interval. Rank zero writes:

- `swirl_scalar_bounds.csv` records raw extrema, volume fractions outside the diagnostic thresholds, and unnormalized integrals of the raw excess outside the unit interval.
- `swirl_scalar_divergence.csv` records divergence from the nodal scalar derivative operators at a coarser cadence. The optional `interior` rows exclude planes near the x boundaries. Use the same row selection across mesh comparisons.

The nodal scalar-operator divergence, pressure-projection divergence and scalar face divergence use different discrete operators. For mapping quality, read `SCALAR_MAP` and compare `raw_div` and `mapped_div` (maximum absolute face divergence before and after mapping) along with `max_du` (maximum absolute component correction).

For DIV and FCT, `SCALAR_STAGE` gives the unlimited stage candidate extrema and `SCALAR_BASE` gives the complete low-order baseline extrema and joint CFL. Assess complete-step bounds using `SCALAR_ACCEPT` and inspect the active-volume mass balance with `SCALAR_MASS`. For a limiter-activation check, inspect `limited_faces`, `removed_L1` and `removed_max` together in `SCALAR_INTERNAL` (the last log field is labeled `max`). A counted face must exceed the reduction-factor roundoff threshold and have nonzero removed flux. `theta_min` is the minimum factor over counted faces, or 1 when the count is zero. The removal norms are the sum and maximum of the absolute area- and time-integrated antidiffusive flux reductions, including reductions below the count threshold.

`SCALAR_OUTLET` reports area- and time-integrated physical fluxes, with positive values in the increasing-x direction. It gives the accepted inlet flux and the accepted, unlimited RK, and low-order outlet fluxes. For FCT, `SCALAR_BOUNDARY` gives the minimum fraction of the RK outlet correction accepted and the signed sum, absolute sum, and maximum magnitude of the applied correction relative to the low-order outlet baseline. `SCALAR_OUTLET_LIMIT` reports the sum and maximum of the absolute differences from the unlimited RK outlet integral, after weighting each face by area and before summing across faces.

Record the input, build and results with each verification run. Distinguish configuration, compilation and linking from successful execution and from numerical acceptance.

Use a fresh directory for each run so diagnostic files start with the expected schema. Primary fields and checkpoints follow the input settings. `nvisu` selects spatial sampling for visualization output. Use full-resolution restart data for small-scale analysis.
