# Repository Guidelines

## Project Structure & Module Organization

- `src/` contains the Fortran solver, numerical operators, and `Case-*.f90` modules. Register new sources in `src/CMakeLists.txt` and connect case routines in `src/case.f90`.
- `src/module_param.f90` defines case identifiers; `src/parameters.f90` reads namelists.
- `examples/` holds simulation inputs. `examples/Swirl-Jet/` documents the coaxial swirling jet, selected by `itype=733`.
- `tests/` contains CTest configurations, reference data, and Python verification scripts. `cmake/` contains dependency and compiler configuration; `docs/` contains documentation and figures.

## Build, Test, and Development Commands

Use CMake 3.20+, a modern Fortran compiler, matching MPI wrappers, and Python 3 with NumPy. The dependency configuration defaults to 2DECOMP&FFT v2.0.4.

```sh
export FC=mpif90
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build -j 4
ctest --test-dir build --output-on-failure
```

These commands configure, compile, and run the default Taylor–Green vortex verification. Enable additional cases at configuration time:

```sh
cmake -S . -B build -DBUILD_TESTING_FULL=ON
```

For an installed dependency, supply `-Ddecomp2d_DIR=/path/to/lib/decomp2d`. Run a case from its input directory:

```sh
mpirun -np 4 /absolute/path/to/build/bin/xcompact3d input.i3d
```

The default swirl grid is large; use a reduced input for startup checks.

## Coding Style & Naming Conventions

Use free-form Fortran, spaces, and the surrounding indentation style. Keep case modules named `Case-<name>.f90`, routine names descriptive, and constants precision-aware, such as `0.5_mytype`. Preserve copyright headers. Avoid unrelated formatting changes. `indent.sh` uses Emacs batch indentation but also references legacy directories; review its scope before use. Run `git diff --check` before review.

## Testing Guidelines

CTest runs MPI cases; TGV compares results against reference data using Python. Add case tests under `tests/<Case>/` with `CMakeLists.txt` and an `input_test.i3d`-style input. No percentage coverage threshold is specified. For numerical changes, report build configuration, process count, stability diagnostics, and relevant reference comparisons. Distinguish syntax checks from compiled, parallel, and physical validation.

## Commit & Pull Request Guidelines

History uses short imperative subjects such as “Fix…” and “Add…”, often followed by PR numbers; Conventional Commits prefixes are not mandatory. Keep changes focused. PR descriptions should explain changed behavior, affected cases, validation commands/results, and unresolved limitations; link related issues when applicable. Include figures for visualization changes. Exclude generated outputs and preserve case ignore rules. Commit or push only when explicitly authorized.
