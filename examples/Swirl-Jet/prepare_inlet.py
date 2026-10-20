#!/usr/bin/env python3
"""Convert a dimensional CSV plane to a nondimensional Cartesian swirl inlet.

Requires NumPy and SciPy. Reconstructs cylindrical modes 0/12/24 with mild radial
smoothing and a distant periodic fringe. Cartesian smoothing is also available.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicSpline, LinearNDInterpolator, RegularGridInterpolator
from scipy.ndimage import gaussian_filter, gaussian_filter1d
from scipy.spatial import Delaunay


def read_csv(path):
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    try:
        start = lines.index("[Data]") + 1
    except ValueError as exc:
        raise ValueError("CSV must contain a [Data] section") from exc
    if start >= len(lines):
        raise ValueError("CSV [Data] section is missing its field header")
    fields = [field.strip() for field in lines[start].split(",")]
    expected = ["x [ m ]", "y [ m ]", "z [ m ]", "Velocity u [ m s^-1 ]",
                "Velocity v [ m s^-1 ]", "Velocity w [ m s^-1 ]"]
    if fields != expected:
        raise ValueError("Expected coordinates in m and u/v/w in m s^-1")
    data = np.loadtxt(lines[start + 1:], delimiter=",")
    if data.ndim != 2 or data.shape[1] != 6 or not np.isfinite(data).all():
        raise ValueError("CSV must contain six finite numeric columns")
    if np.ptp(data[:, 0]) > 1e-9:
        raise ValueError("CSV is not a plane normal to the x direction")
    if len(np.unique(data[:, 1:3], axis=0)) != len(data):
        raise ValueError("Duplicate transverse points require explicit reconciliation")
    return data


def read_settings(path):
    text = "\n".join(line.split("!")[0] for line in path.read_text(encoding="utf-8-sig").splitlines())
    settings = {}
    for key in ("ny", "nz", "yly", "zlz", "ncly1", "nclyn", "nclz1", "nclzn", "rai", "rao", "re"):
        match = re.search(r"\b" + key + r"\s*=\s*([-+\d.eEdD]+)", text, re.I)
        if not match:
            raise ValueError(f"Missing {key} in input configuration")
        settings[key] = float(match.group(1).replace("D", "e").replace("d", "e"))
    if not all(np.isfinite(value) for value in settings.values()):
        raise ValueError("Non-finite grid or physical setting")
    for key in ("ny", "nz", "ncly1", "nclyn", "nclz1", "nclzn"):
        if not settings[key].is_integer():
            raise ValueError(f"{key} must be an integer")
        settings[key] = int(settings[key])
    if min(settings["ny"], settings["nz"]) < 2 or min(settings["yly"], settings["zlz"], settings["re"]) <= 0:
        raise ValueError("Grid sizes, domain lengths and Reynolds number must be positive")
    for axis in ("y", "z"):
        if (settings[f"ncl{axis}1"] == 0) != (settings[f"ncl{axis}n"] == 0):
            raise ValueError("Inconsistent periodic boundary flags")
    return settings


def bulk_velocity(interpolator, radius, nr, nt):
    s, weights = leggauss(nr)
    s, weights = (s + 1) / 2, weights / 2
    angles = np.arange(nt) * 2 * np.pi / nt
    radial = radius * np.sqrt(s)
    points = np.stack([radial[:, None] * np.cos(angles),
                       radial[:, None] * np.sin(angles)], axis=-1)
    velocity = interpolator(points)
    if not np.isfinite(velocity).all():
        raise ValueError("CSV/grid does not cover the complete reference disk")
    return float(weights @ velocity[:, :, 0].mean(axis=1))


def smooth_join(t):
    t = np.clip(t, 0, 1)
    return 10 * t**3 - 15 * t**4 + 6 * t**5


def fourier_plane(interpolator, yy, zz, sigma, radial_step, max_radius, core_field,
                  outer_radius):
    """Cylindrical reconstruction with a regular Cartesian center patch.

    The center patch avoids imposing false cylindrical conditions on a measured
    Cartesian vector at r=0. It joins with zero first/second blend derivatives.
    """
    radial = np.linspace(0, max_radius, int(np.ceil(max_radius / radial_step)) + 1)
    angles = np.arange(720) * (2 * np.pi / 720)
    sn, cs = np.sin(angles), np.cos(angles)
    points = np.stack([radial[:, None] * sn, radial[:, None] * cs], axis=-1)
    cart = interpolator(points)
    if not np.isfinite(cart).all():
        raise ValueError("CSV convex hull does not cover the Fourier sampling disk; extrapolation is disabled")
    cyl = np.stack([cart[..., 0], cart[..., 1] * sn + cart[..., 2] * cs,
                    cart[..., 1] * cs - cart[..., 2] * sn], axis=-1)
    modes = np.array([0, 12, 24])
    coefficients = np.fft.rfft(cyl, axis=1)[:, modes] / len(angles)
    if sigma:
        coefficients = gaussian_filter1d(coefficients, sigma / (radial[1] - radial[0]),
                                         axis=0, mode="nearest", truncate=4.0)
    radius = np.hypot(yy, zz)
    theta = np.arctan2(yy, zz)
    # Samples beyond the reconstruction disk lie in the constant fringe.
    at_radius = CubicSpline(radial, coefficients, axis=0)(np.minimum(radius, max_radius))
    factors = np.exp(1j * theta[..., None] * modes)
    factors[..., 1:] *= 2
    cyl = np.real(np.sum(at_radius * factors[..., None], axis=-2))
    st, ct = np.sin(theta), np.cos(theta)
    field = np.stack([cyl[..., 0], cyl[..., 1] * st + cyl[..., 2] * ct,
                      cyl[..., 1] * ct - cyl[..., 2] * st], axis=-1)
    join_start, join_end = 0.1 * outer_radius, 0.2 * outer_radius
    weight = smooth_join((radius - join_start) / (join_end - join_start))
    field = core_field + weight[..., None] * (field - core_field)
    return field, {"retained_cylindrical_modes": modes.tolist(),
                   "radial_sampling_step_over_Di": float(radial[1] - radial[0]),
                   "theta_convention": "atan2(y,z); utheta=(uy*z-uz*y)/r",
                   "Cartesian_center_join_radii_over_Di": [join_start, join_end]}


def condition_plane(interpolator, y, z, sigma_cells=0.25, periodic_fringe=True,
                    outer_radius=1.0, reconstruction="fourier"):
    """Condition three-component vectors without a polar-axis singularity.

    sigma is the same physical distance in both directions. Oversampling avoids
    an almost-identity Gaussian on the solver mesh at subcell widths. The fringe
    has zero first/second derivatives at both joins and is constant near edges.
    """
    if not np.isfinite(sigma_cells) or sigma_cells < 0:
        raise ValueError("Smoothing sigma must be finite and nonnegative")
    if reconstruction not in ("fourier", "cartesian"):
        raise ValueError("Reconstruction must be fourier or cartesian")
    step = min(y[1] - y[0], z[1] - z[0])
    sigma = sigma_cells * step
    fine_y = np.linspace(y[0], y[-1], 2 * len(y) - 1) if sigma else y
    fine_z = np.linspace(z[0], z[-1], 2 * len(z) - 1) if sigma else z
    yy, zz = np.meshgrid(fine_y, fine_z, indexing="ij")
    field = interpolator(np.stack([yy, zz], axis=-1))
    if not np.isfinite(field).all():
        raise ValueError("CSV convex hull does not cover the Cartesian inlet; extrapolation is disabled")
    reference_field = field[::2, ::2] if sigma else field
    if sigma:
        field = gaussian_filter(field, (sigma / (fine_y[1] - fine_y[0]),
                                        sigma / (fine_z[1] - fine_z[0]), 0),
                                mode="nearest", truncate=4.0)
        field = field[::2, ::2]
    details = {"method": "isotropic Cartesian Gaussian; no angular averaging",
               "sigma_solver_cells": sigma_cells, "sigma_over_Di": float(sigma),
               "fine_grid_factor": 2 if sigma else 1, "periodic_fringe": periodic_fringe}
    yy, zz = np.meshgrid(y, z, indexing="ij")
    half_width = min(y[-1], z[-1], -y[0], -z[0])
    start, end = 0.75 * half_width, 0.9375 * half_width
    if periodic_fringe and (start <= outer_radius + 4 * sigma or end - start < 4 * step):
        raise ValueError("Periodic fringe must lie beyond the jet and span at least four cells")
    if reconstruction == "fourier":
        max_radius = end + 4 * sigma if periodic_fringe else float(np.hypot(yy, zz).max())
        field, polar = fourier_plane(interpolator, yy, zz, sigma, step / 4,
                                     max_radius, field, outer_radius)
        details.update(polar)
        details["method"] = "cylindrical Fourier 0/12/24; Gaussian radial coefficients; Cartesian center"
    if periodic_fringe:
        # Record the measured distant coflow, then determine a common axial
        # background by preserving the raw Cartesian rectangle's axial flux.
        # Only the fringe changes; core/annular velocities are never rescaled.
        angles = np.arange(720) * (2 * np.pi / 720)
        radii = np.linspace(start, end, 65)
        ring = np.stack([radii[:, None] * np.cos(angles),
                         radii[:, None] * np.sin(angles)], axis=-1)
        measured_background = float(interpolator(ring)[..., 0].mean())
        weight = 1 - smooth_join((np.hypot(yy, zz) - start) / (end - start))
        wy, wz = np.ones(len(y)), np.ones(len(z))
        wy[[0, -1]], wz[[0, -1]] = 0.5, 0.5
        area_weight = wy[:, None] * wz
        target = np.sum(area_weight * reference_field[..., 0])
        background_ux = ((target - np.sum(area_weight * weight * field[..., 0])) /
                         np.sum(area_weight * (1 - weight)))
        background = np.array([background_ux, 0., 0.])
        field = background + weight[..., None] * (field - background)
        details.update(fringe_radii_over_Di=[float(start), float(end)],
                       background_velocity_Ub=background.tolist(),
                       measured_fringe_mean_ux_Ub=measured_background,
                       background_rule="preserve raw Cartesian rectangle axial flux by trapezoidal quadrature")
    return field, details


def azimuthal_diagnostics(before, after, outer_radius):
    """Area-weighted cylindrical Fourier energy and complex-mode error on disk.

    theta=atan2(y,z), matching Case-Swirl-Jet. No r=0 sample is used. The
    complex error captures amplitude and phase changes; tiny modes are null.
    """
    s, weights = leggauss(192)
    radial = outer_radius * np.sqrt((s + 1) / 2)
    weights = weights / 2
    angles = np.arange(720) * (2 * np.pi / 720)
    sn, cs = np.sin(angles), np.cos(angles)
    points = np.stack([radial[:, None] * sn, radial[:, None] * cs], axis=-1)
    coefficients = []
    for sample in (before, after):
        cart = sample(points)
        cyl = np.stack([cart[..., 0], cart[..., 1] * sn + cart[..., 2] * cs,
                        cart[..., 1] * cs - cart[..., 2] * sn], axis=-1)
        coefficients.append(np.fft.rfft(cyl, axis=1)[:, :33] / len(angles))
    raw, filtered = coefficients
    energies = [np.einsum("r,rmc->mc", weights, abs(v)**2) for v in coefficients]
    energies = [v * np.r_[1, np.full(32, 2)][:, None] for v in energies]
    difference = np.einsum("r,rmc->mc", weights, abs(filtered - raw)**2)
    reference = np.einsum("r,rmc->mc", weights, abs(raw)**2)
    result = {"disk_radius_over_Di": outer_radius, "theta_convention": "atan2(y,z)"}
    for component, name in enumerate(("ux", "ur", "utheta")):
        threshold = 1e-8 * energies[0][:, component].sum()
        cross = np.einsum("r,rm,rm->m", weights, filtered[..., component],
                          raw[..., component].conj())
        result[name] = {
            "area_mean_before_Ub": float(weights @ raw[:, 0, component].real),
            "area_mean_after_Ub": float(weights @ filtered[:, 0, component].real),
            "energy_before_m0_to_m32": energies[0][:, component].tolist(),
            "energy_after_m0_to_m32": energies[1][:, component].tolist(),
            "phase_change_degrees_m0_to_m32": [float(np.angle(cross[m], deg=True))
                if energies[0][m, component] > threshold else None for m in range(33)],
            "complex_relative_error_m0_to_m32": [
                float(np.sqrt(difference[m, component] / reference[m, component]))
                if energies[0][m, component] > threshold else None for m in range(33)]}
    return result


def prepare(csv_path, input_path, output, ri=0.0075, ro=0.015, nu=0.889e-6,
            center_y=0.0, center_z=0.0, sigma_cells=0.25, periodic_fringe=True,
            reconstruction="fourier"):
    if not all(np.isfinite(v) for v in (ri, ro, nu, center_y, center_z)):
        raise ValueError("Physical scales and center must be finite")
    if not 0 < ri < ro or nu <= 0:
        raise ValueError("Require 0 < inner radius < outer radius and nu > 0")
    outputs = {output.resolve(), output.with_suffix(output.suffix + ".json").resolve()}
    if outputs.intersection({csv_path.resolve(), input_path.resolve()}):
        raise ValueError("Output must not overwrite a source file")
    data = read_csv(csv_path)
    cfg = read_settings(input_path)
    length = 2 * ri
    if not np.allclose([cfg["rai"], cfg["rao"]], [ri / length, ro / length], rtol=0, atol=1e-10):
        raise ValueError("Configuration radii do not match the physical radii/inner-diameter scale")
    points = data[:, 1:3] - [center_y, center_z]
    interpolation = LinearNDInterpolator(Delaunay(points), data[:, 3:6])
    coarse = bulk_velocity(interpolation, ri, 256, 720)
    speed = bulk_velocity(interpolation, ri, 512, 1440)
    if speed <= 0 or not np.isfinite(speed):
        raise ValueError("Center-disk average axial velocity must be positive")
    if abs(speed - coarse) / speed > 1e-4:
        raise ValueError("Reference disk integration has not converged")

    # Include both bounding coordinates even for periodic solver meshes. The
    # solver omits the duplicate endpoint; the table still covers the full box.
    axes = []
    for axis, length_key in (("y", "yly"), ("z", "zlz")):
        count = cfg[f"n{axis}"] + int(cfg[f"ncl{axis}1"] == 0)
        axes.append(np.linspace(-cfg[length_key] / 2, cfg[length_key] / 2, count))
    y, z = axes
    yy, zz = np.meshgrid(y, z, indexing="ij")
    def nondimensional(points):
        return interpolation(points * length) / speed

    raw_velocity = nondimensional(np.stack([yy, zz], axis=-1))
    if not np.isfinite(raw_velocity).all():
        raise ValueError("CSV convex hull does not cover the Cartesian inlet; extrapolation is disabled")
    periodic = cfg["ncly1"] == 0 or cfg["nclz1"] == 0
    velocity, conditioning = condition_plane(nondimensional, y, z, sigma_cells,
                                              periodic_fringe and periodic, cfg["rao"], reconstruction)
    if not np.isfinite(velocity).all():
        raise ValueError("Conditioned inlet contains a non-finite velocity")
    gridded = RegularGridInterpolator((y, z), velocity)
    grid_bulk = bulk_velocity(gridded, cfg["rai"], 256, 720)
    raw_gridded = RegularGridInterpolator((y, z), raw_velocity)
    raw_bulk = bulk_velocity(raw_gridded, cfg["rai"], 256, 720)
    modes = azimuthal_diagnostics(nondimensional, gridded, cfg["rao"])
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="ascii", newline="\n") as stream:
        stream.write(f"SWIRL_INLET_V1\n{len(y)} {len(z)}\n")
        stream.write(f"{cfg['rai']:.17e} {cfg['rao']:.17e}\n")
        stream.write(" ".join(f"{value:.17e}" for value in y) + "\n")
        stream.write(" ".join(f"{value:.17e}" for value in z) + "\n")
        # j (y) varies fastest, k (z) outermost; triples are ux, uy, uz.
        np.savetxt(stream, velocity.transpose(1, 0, 2).reshape(-1, 3), fmt="%.17e")
    metadata = {
        "format": "SWIRL_INLET_V1", "csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
        "preprocessor_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "table_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "source_rows": len(data), "source_x_m": float(data[0, 0]),
        "source_x_over_Di": float(data[0, 0] / length), "center_yz_m": [center_y, center_z],
        "inner_radius_m": ri, "outer_radius_m": ro, "Di_m": length, "Ubulk_m_s": speed,
        "nu_m2_s": nu, "Re_measured": speed * length / nu, "Re_configured": cfg["re"],
        "Q_core_L_min": speed * np.pi * ri**2 * 60000,
        "reference_quadrature_relative_change": abs(speed - coarse) / speed,
        "table_dimensions_yz": [len(y), len(z)],
        "table_step_yz": [float(y[1] - y[0]), float(z[1] - z[0])],
        "cartesian_core_bulk_Ub": grid_bulk, "cartesian_core_flux_relative_error": grid_bulk - 1,
        "unconditioned_cartesian_core_bulk_Ub": raw_bulk,
        "conditioning": conditioning, "azimuthal_modes": modes,
        "conditioning_only_modes": azimuthal_diagnostics(raw_gridded, gridded, cfg["rao"]),
        "component_min_Ub": velocity.min(axis=(0, 1)).tolist(),
        "component_max_Ub": velocity.max(axis=(0, 1)).tolist(),
        "unconditioned_component_max_Ub": raw_velocity.max(axis=(0, 1)).tolist(),
        "periodic_y_endpoint_jump_Ub": np.max(np.abs(velocity[-1] - velocity[0]), axis=0).tolist(),
        "periodic_z_endpoint_jump_Ub": np.max(np.abs(velocity[:, -1] - velocity[:, 0]), axis=0).tolist(),
        "notes": ["All Cartesian components divided by the same center-disk bulk speed.",
                  "No imposed rotation reversal, full-plane axisymmetrization or interior velocity rescaling.",
                  "The sampled physical x plane becomes computational inlet x=0.",
                  "Core bulk is defined at the sampled plane, not inferred inside the upstream pipe.",
                  "Smoothing changes small scales; inspect cylindrical mode energies and complex errors.",
                  "Periodic fringe modifies only the distant exterior; no solenoidal projection is performed."]}
    output.with_suffix(output.suffix + ".json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=Path("input_LES.i3d"))
    parser.add_argument("--output", type=Path, default=Path("swirl_inlet.dat"))
    parser.add_argument("--inner-radius", type=float, default=0.0075, help="m")
    parser.add_argument("--outer-radius", type=float, default=0.015, help="m")
    parser.add_argument("--nu", type=float, default=0.889e-6, help="m2/s")
    parser.add_argument("--center-y", type=float, default=0.0, help="m")
    parser.add_argument("--center-z", type=float, default=0.0, help="m")
    parser.add_argument("--smooth-sigma", type=float, default=0.25,
                        help="Gaussian sigma in minimum transverse solver-cell widths; 0 disables smoothing")
    parser.add_argument("--reconstruction", choices=("fourier", "cartesian"), default="fourier",
                        help="cylindrical 0/12/24 modes or all-mode Cartesian interpolation")
    parser.add_argument("--no-periodic-fringe", action="store_true",
                        help="disable the constant far-field join (periodic endpoints may then differ)")
    args = parser.parse_args()
    try:
        result = prepare(args.csv, args.input, args.output, args.inner_radius, args.outer_radius,
                         args.nu, args.center_y, args.center_z, args.smooth_sigma,
                         not args.no_periodic_fringe, args.reconstruction)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
