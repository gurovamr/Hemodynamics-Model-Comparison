#!/usr/bin/env python3
"""Reproduce the CFD-versus-ROM workflow audit without changing any notebook."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "workflow_audit"
OUT.mkdir(parents=True, exist_ok=True)

RHO = 1060.0
MU = 3.5e-3
PA_TO_MMHG = 1.0 / 133.322

CENTERLINE = ROOT / "data/mr_limited/geometry/segment_test_centerline_candidate.vtu"
CFD_SURFACE = ROOT / "openfoam/segment_test_V2/constant/triSurface/segment_test_closed.stl"
FOAM = ROOT / "openfoam/segment_test_V2/segment_test.foam"

PLANES = {
    "inlet": {
        "origin_m": np.array(
            [0.005926580083397241, -0.021677941320206097, -0.006866637498613615]
        ),
        "normal": np.array(
            [0.5103037804570406, -0.8596550644734353, -0.024149985019174598]
        ),
    },
    "midslice": {
        "origin_m": np.array(
            [0.010303987000313618, -0.02773153503893052, -0.007040573269073331]
        ),
        "normal": np.array(
            [0.6324414356546973, -0.7702129056188429, 0.08240091313331853]
        ),
    },
    "outlet": {
        "origin_m": np.array(
            [0.014722232293941245, -0.03462298076987406, -0.005146773400980403]
        ),
        "normal": np.array(
            [0.4278871251313104, -0.9017025116588478, 0.06200958486385685]
        ),
    },
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_centerline():
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(CENTERLINE))
    reader.Update()
    data = reader.GetOutput()
    points = vtk_to_numpy(data.GetPoints().GetData()).astype(float)
    original_point_ids = vtk_to_numpy(
        data.GetPointData().GetArray("vtkOriginalPointIds")
    )
    ce_radius = vtk_to_numpy(data.GetCellData().GetArray("ce_radius")).astype(float)
    labels = vtk_to_numpy(data.GetCellData().GetArray("labels")).astype(int)
    original_cell_ids = vtk_to_numpy(
        data.GetCellData().GetArray("vtkOriginalCellIds")
    )

    adjacency = {i: [] for i in range(data.GetNumberOfPoints())}
    edge_cell = {}
    for cell_index in range(data.GetNumberOfCells()):
        cell = data.GetCell(cell_index)
        p0, p1 = cell.GetPointId(0), cell.GetPointId(1)
        adjacency[p0].append(p1)
        adjacency[p1].append(p0)
        edge_cell[tuple(sorted((p0, p1)))] = cell_index

    terminals = [node for node, neighbors in adjacency.items() if len(neighbors) == 1]
    pair = max(
        (
            (a, b)
            for index, a in enumerate(terminals)
            for b in terminals[index + 1 :]
        ),
        key=lambda pair_: np.linalg.norm(points[pair_[0]] - points[pair_[1]]),
    )
    queue = [pair[0]]
    parent = {pair[0]: None}
    for node in queue:
        if node == pair[1]:
            break
        for neighbor in adjacency[node]:
            if neighbor not in parent:
                parent[neighbor] = node
                queue.append(neighbor)
    path = []
    node = pair[1]
    while node is not None:
        path.append(node)
        node = parent[node]
    path = path[::-1]

    path_points = points[path]
    path_s_mm = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(path_points, axis=0), axis=1))]
    path_cells = np.array(
        [edge_cell[tuple(sorted(edge))] for edge in zip(path[:-1], path[1:])],
        dtype=int,
    )
    path_radius = ce_radius[path_cells]
    node_radius = np.r_[
        path_radius[0], 0.5 * (path_radius[:-1] + path_radius[1:]), path_radius[-1]
    ]
    return {
        "data": data,
        "points": points,
        "path": np.asarray(path),
        "path_points": path_points,
        "path_s_mm": path_s_mm,
        "path_cells": path_cells,
        "path_radius_mm": path_radius,
        "node_radius_mm": node_radius,
        "labels": labels,
        "original_point_ids": original_point_ids,
        "original_cell_ids": original_cell_ids,
        "terminals": terminals,
    }


def intersect_plane(points_mm, s_mm, origin_m, normal):
    points_m = points_mm * 1e-3
    normal = normal / np.linalg.norm(normal)
    signed = (points_m - origin_m) @ normal
    hits = []
    for index, (d0, d1) in enumerate(zip(signed[:-1], signed[1:])):
        if d0 == 0.0 or d0 * d1 <= 0.0:
            fraction = d0 / (d0 - d1) if d0 != d1 else 0.0
            point_m = points_m[index] + fraction * (points_m[index + 1] - points_m[index])
            tangent = points_m[index + 1] - points_m[index]
            tangent /= np.linalg.norm(tangent)
            angle_deg = np.degrees(
                np.arccos(np.clip(abs(np.dot(tangent, normal)), -1.0, 1.0))
            )
            position_mm = s_mm[index] + fraction * (s_mm[index + 1] - s_mm[index])
            hits.append(
                {
                    "edge_index": index,
                    "fraction": fraction,
                    "global_s_mm": position_mm,
                    "point_m": point_m,
                    "radial_offset_mm": np.linalg.norm(point_m - origin_m) * 1e3,
                    "normal_tangent_angle_deg": angle_deg,
                    "extrapolated": False,
                }
            )
    if hits:
        return min(hits, key=lambda row: row["radial_offset_mm"])

    index = len(points_m) - 2
    d0, d1 = signed[-2], signed[-1]
    fraction = d0 / (d0 - d1)
    point_m = points_m[index] + fraction * (points_m[index + 1] - points_m[index])
    tangent = points_m[-1] - points_m[-2]
    tangent /= np.linalg.norm(tangent)
    return {
        "edge_index": index,
        "fraction": fraction,
        "global_s_mm": s_mm[index] + fraction * (s_mm[-1] - s_mm[index]),
        "point_m": point_m,
        "radial_offset_mm": np.linalg.norm(point_m - origin_m) * 1e3,
        "normal_tangent_angle_deg": np.degrees(
            np.arccos(np.clip(abs(np.dot(tangent, normal)), -1.0, 1.0))
        ),
        "extrapolated": True,
    }


def load_cfd_grid():
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(str(FOAM))
    reader.UpdateInformation()
    info = reader.GetExecutive().GetOutputInformation(0)
    info.Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(), 384.0)
    reader.GetExecutive().Update()
    block = reader.GetOutput().GetBlock(0)
    cell_to_point = vtk.vtkCellDataToPointData()
    cell_to_point.SetInputData(block)
    cell_to_point.Update()
    return cell_to_point.GetOutput()


def probe_pressure(grid, points_m):
    polydata = vtk.vtkPolyData()
    vtk_points = vtk.vtkPoints()
    vtk_points.SetData(numpy_to_vtk(np.asarray(points_m), deep=True))
    polydata.SetPoints(vtk_points)
    probe = vtk.vtkProbeFilter()
    probe.SetInputData(polydata)
    probe.SetSourceData(grid)
    probe.Update()
    result = probe.GetOutput()
    valid = vtk_to_numpy(
        result.GetPointData().GetArray("vtkValidPointMask")
    ).astype(bool)
    pressure = vtk_to_numpy(result.GetPointData().GetArray("p")) * RHO
    return pressure, valid


def equivalent_radius_cuts(surface, sample_s_mm, sample_points_mm, full_s_mm, full_points_mm):
    radius_mm = []
    offsets_mm = []
    for position_mm, point_mm in zip(sample_s_mm, sample_points_mm):
        index = int(np.searchsorted(full_s_mm, position_mm))
        index = int(np.clip(index, 1, len(full_s_mm) - 1))
        if (
            index < len(full_s_mm) - 1
            and np.isclose(position_mm, full_s_mm[index], atol=1e-8)
        ):
            tangent = full_points_mm[index + 1] - full_points_mm[index - 1]
        else:
            tangent = full_points_mm[index] - full_points_mm[index - 1]
        tangent = tangent / np.linalg.norm(tangent)

        plane = vtk.vtkPlane()
        plane.SetOrigin(*point_mm)
        plane.SetNormal(*tangent)
        cutter = vtk.vtkCutter()
        cutter.SetCutFunction(plane)
        cutter.SetInputData(surface)
        cutter.Update()
        stripper = vtk.vtkStripper()
        stripper.SetInputData(cutter.GetOutput())
        stripper.JoinContiguousSegmentsOn()
        stripper.Update()
        contours = stripper.GetOutput()

        candidates = []
        for cell_index in range(contours.GetNumberOfCells()):
            cell = contours.GetCell(cell_index)
            polygon = np.array(
                [
                    contours.GetPoint(cell.GetPointId(point_index))
                    for point_index in range(cell.GetNumberOfPoints())
                ]
            )
            if len(polygon) < 3:
                continue
            if np.linalg.norm(polygon[0] - polygon[-1]) < 1e-5:
                polygon = polygon[:-1]
            area_vector = 0.5 * np.sum(
                np.cross(polygon, np.roll(polygon, -1, axis=0)), axis=0
            )
            area_mm2 = abs(np.dot(area_vector, tangent))
            candidates.append(
                (
                    np.linalg.norm(polygon.mean(axis=0) - point_mm),
                    np.sqrt(area_mm2 / np.pi),
                )
            )
        if not candidates:
            raise RuntimeError(f"No closed surface contour at s={position_mm:.6f} mm")
        offset, radius = min(candidates)
        offsets_mm.append(offset)
        radius_mm.append(radius)
    return np.asarray(radius_mm), np.asarray(offsets_mm)


def pressure_method_audit(slice_hits, grid):
    centerline_points = [slice_hits[name]["point_m"] for name in PLANES]
    centerline_pressure, valid = probe_pressure(grid, centerline_points)
    rows = []
    for index, name in enumerate(PLANES):
        integrated = read_rows(
            ROOT / f"openfoam/segment_test_V2/data/{name}_integration_variable.csv"
        )[0]
        area_average = float(integrated["p"]) / float(integrated["Area"]) * RHO

        raw = read_rows(ROOT / f"openfoam/segment_test_V2/data/{name}_integration.csv")
        raw_mean = np.mean([float(row["p"]) for row in raw]) * RHO

        point_name = "midsclice_point_u.csv" if name == "midslice" else f"{name}_point_u.csv"
        point_rows = read_rows(ROOT / "data_parafoam" / point_name)
        point_mean = np.mean([float(row["p"]) for row in point_rows]) * RHO
        rows.append(
            {
                "slice": name,
                "area_weighted_pressure_Pa": area_average,
                "raw_slice_unweighted_mean_Pa": raw_mean,
                "point_export_unweighted_mean_Pa": point_mean,
                "centreline_probe_Pa": centerline_pressure[index] if valid[index] else np.nan,
            }
        )
    return rows


def make_dependency_graph():
    fig, ax = plt.subplots(figsize=(12, 5.2))
    ax.axis("off")
    boxes = [
        (0.02, 0.65, "OpenFOAM\nsimpleFoam t=384"),
        (0.20, 0.78, "ParaView Slice +\nIntegrate Variables"),
        (0.20, 0.45, "VTK cell-to-point +\nProbeFilter"),
        (0.39, 0.78, "Area-averaged\np, Q, A CSVs"),
        (0.39, 0.45, "Centreline point\np(s), U(s), valid mask"),
        (0.02, 0.15, "Graph/candidate VTU\ncoordinates + ce_radius"),
        (0.39, 0.15, "Path ordering,\n0D/1D construction"),
        (0.62, 0.62, "Notebook 04\nmodel summary"),
        (0.82, 0.62, "Report resistance /\npressure-drop figures"),
        (0.62, 0.25, "Pressure-profile\nalignment / interpolation"),
        (0.82, 0.25, "Final pressure-profile\nfigure"),
    ]
    for x, y, text in boxes:
        ax.text(
            x,
            y,
            text,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.4", "facecolor": "white", "edgecolor": "#444"},
        )
    arrows = [
        ((0.15, 0.67), (0.20, 0.78)),
        ((0.15, 0.64), (0.20, 0.47)),
        ((0.33, 0.78), (0.39, 0.78)),
        ((0.33, 0.47), (0.39, 0.47)),
        ((0.20, 0.18), (0.39, 0.18)),
        ((0.53, 0.78), (0.62, 0.65)),
        ((0.53, 0.18), (0.62, 0.58)),
        ((0.74, 0.62), (0.82, 0.62)),
        ((0.53, 0.47), (0.62, 0.28)),
        ((0.53, 0.18), (0.62, 0.28)),
        ((0.74, 0.28), (0.82, 0.28)),
    ]
    for start, end in arrows:
        ax.annotate(
            "",
            xy=end,
            xytext=start,
            xycoords=ax.transAxes,
            textcoords=ax.transAxes,
            arrowprops={"arrowstyle": "->", "lw": 1.3, "color": "#555"},
        )
    ax.set_title("CFD-to-ROM Pressure and Resistance Dependency Graph", fontsize=14)
    fig.savefig(OUT / "workflow_dependency_graph.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    vtk.vtkLogger.SetStderrVerbosity(vtk.vtkLogger.VERBOSITY_OFF)
    centerline = load_centerline()
    path_points = centerline["path_points"]
    path_s = centerline["path_s_mm"]
    slice_hits = {
        name: intersect_plane(
            path_points,
            path_s,
            plane["origin_m"],
            plane["normal"],
        )
        for name, plane in PLANES.items()
    }

    geometry_rows = []
    for name, hit in slice_hits.items():
        rom_endpoint = path_points[0] if name == "inlet" else path_points[-1]
        geometry_rows.append(
            {
                "slice": name,
                "CFD_plane_origin_x_m": PLANES[name]["origin_m"][0],
                "CFD_plane_origin_y_m": PLANES[name]["origin_m"][1],
                "CFD_plane_origin_z_m": PLANES[name]["origin_m"][2],
                "centreline_plane_x_m": hit["point_m"][0],
                "centreline_plane_y_m": hit["point_m"][1],
                "centreline_plane_z_m": hit["point_m"][2],
                "global_s_mm": hit["global_s_mm"],
                "plane_origin_to_centreline_mm": hit["radial_offset_mm"],
                "normal_tangent_angle_deg": hit["normal_tangent_angle_deg"],
                "original_ROM_endpoint_distance_mm": np.linalg.norm(
                    rom_endpoint * 1e-3 - PLANES[name]["origin_m"]
                )
                * 1e3,
                "plane_intersection_extrapolated": hit["extrapolated"],
            }
        )
    write_rows(
        OUT / "geometry_alignment.csv", list(geometry_rows[0]), geometry_rows
    )

    grid = load_cfd_grid()
    pressure_rows = pressure_method_audit(slice_hits, grid)
    write_rows(
        OUT / "pressure_method_comparison.csv",
        list(pressure_rows[0]),
        pressure_rows,
    )

    s_in = slice_hits["inlet"]["global_s_mm"]
    s_out = min(slice_hits["outlet"]["global_s_mm"], path_s[-1])
    sample_s = np.r_[
        s_in,
        path_s[(path_s > s_in) & (path_s < s_out)],
        s_out,
    ]
    sample_points = np.column_stack(
        [np.interp(sample_s, path_s, path_points[:, axis]) for axis in range(3)]
    )
    rom_radius = np.interp(sample_s, path_s, centerline["node_radius_mm"])

    surface_reader = vtk.vtkSTLReader()
    surface_reader.SetFileName(str(CFD_SURFACE))
    surface_reader.Update()
    cfd_radius, centerline_offset = equivalent_radius_cuts(
        surface_reader.GetOutput(), sample_s, sample_points, path_s, path_points
    )
    radius_rows = []
    for values in zip(sample_s, cfd_radius, rom_radius, centerline_offset):
        s_value, cfd_r, rom_r, offset = values
        radius_rows.append(
            {
                "global_s_mm": s_value,
                "radius_CFD_equivalent_mm": cfd_r,
                "radius_ROM_ce_mm": rom_r,
                "difference_CFD_minus_ROM_mm": cfd_r - rom_r,
                "relative_difference_pct": (cfd_r - rom_r) / rom_r * 100.0,
                "centreline_to_section_centroid_mm": offset,
            }
        )
    write_rows(OUT / "radius_consistency.csv", list(radius_rows[0]), radius_rows)

    flow_rows = read_rows(ROOT / "output/06_clear_results/flow_rate_summary.csv")
    pressure_summary = {
        row["Slice"]: row
        for row in read_rows(ROOT / "output/06_clear_results/pressure_slice_summary.csv")
    }
    flow_audit = []
    for row in flow_rows:
        name = row["Slice"]
        area_m2 = float(row["Area [mm²]"]) * 1e-6
        flow_m3_s = float(row["Flow rate [m³/s]"])
        bulk_velocity = flow_m3_s / area_m2
        pressure_pa = float(pressure_summary[name]["Average pressure [Pa]"])
        flow_audit.append(
            {
                "slice": name,
                "area_mm2": area_m2 * 1e6,
                "flow_mL_s": flow_m3_s * 1e6,
                "bulk_normal_velocity_m_s": bulk_velocity,
                "static_pressure_Pa": pressure_pa,
                "bulk_dynamic_pressure_Pa": 0.5 * RHO * bulk_velocity**2,
                "bulk_total_pressure_Pa": pressure_pa
                + 0.5 * RHO * bulk_velocity**2,
            }
        )
    write_rows(OUT / "flow_audit.csv", list(flow_audit[0]), flow_audit)

    q_mean = np.mean([row["flow_mL_s"] for row in flow_audit]) * 1e-6
    cfd_dp = (
        pressure_rows[0]["area_weighted_pressure_Pa"]
        - pressure_rows[-1]["area_weighted_pressure_Pa"]
    )
    r_cfd = cfd_dp / q_mean
    r_surface = (
        8.0
        * MU
        / np.pi
        * np.trapz((cfd_radius * 1e-3) ** -4, sample_s * 1e-3)
    )
    r_rom_integral = (
        8.0
        * MU
        / np.pi
        * np.trapz((rom_radius * 1e-3) ** -4, sample_s * 1e-3)
    )
    common_rows = read_rows(
        ROOT / "output/reduced_order_models/model_comparison_common_domain.csv"
    )
    r_rom_native = float(
        next(row for row in common_rows if row["model"] == "Native 0D")["R_Pa_s_m3"]
    )
    resistance_rows = [
        {
            "method": "A: CFD area-averaged DeltaP / mean Q",
            "R_Pa_s_m3": r_cfd,
            "R_mmHg_s_mL": r_cfd * PA_TO_MMHG * 1e-6,
            "relative_to_CFD_pct": 0.0,
        },
        {
            "method": "B1: Poiseuille integral using CFD surface-cut equivalent radius",
            "R_Pa_s_m3": r_surface,
            "R_mmHg_s_mL": r_surface * PA_TO_MMHG * 1e-6,
            "relative_to_CFD_pct": (r_surface - r_cfd) / r_cfd * 100.0,
        },
        {
            "method": "B2: Poiseuille integral using ROM ce_radius",
            "R_Pa_s_m3": r_rom_integral,
            "R_mmHg_s_mL": r_rom_integral * PA_TO_MMHG * 1e-6,
            "relative_to_CFD_pct": (r_rom_integral - r_cfd) / r_cfd * 100.0,
        },
        {
            "method": "C: clipped native 0D resistor sum",
            "R_Pa_s_m3": r_rom_native,
            "R_mmHg_s_mL": r_rom_native * PA_TO_MMHG * 1e-6,
            "relative_to_CFD_pct": (r_rom_native - r_cfd) / r_cfd * 100.0,
        },
    ]
    write_rows(
        OUT / "resistance_audit.csv", list(resistance_rows[0]), resistance_rows
    )

    transforms = [
        {
            "stage": "OpenFOAM pressure",
            "transformation": "kinematic p multiplied by rho=1060 kg/m3",
            "risk": "low; unit conversion is correct",
        },
        {
            "stage": "Slice pressure",
            "transformation": "ParaView p integral divided by slice area",
            "risk": "low; area-weighted cross-sectional average",
        },
        {
            "stage": "Slice flow",
            "transformation": "integrated velocity vector dotted with supplied normal",
            "risk": "low; 0.591% spread across slices",
        },
        {
            "stage": "Centreline pressure",
            "transformation": "cell-to-point conversion then vtkProbeFilter linear interpolation",
            "risk": "medium; centreline pressure is not cross-sectional average",
        },
        {
            "stage": "Invalid CFD probes",
            "transformation": "invalid mask removes 62/124 main-path points",
            "risk": "high if arc length is recomputed after filtering",
        },
        {
            "stage": "Path selection",
            "transformation": "farthest Euclidean terminal pair plus unweighted graph shortest path",
            "risk": "medium; candidate contains a side branch",
        },
        {
            "stage": "Network CFD labels",
            "transformation": "bbox overlap selects labels 1,2,8, then full labels are summed",
            "risk": "critical; includes geometry outside CFD and side branch label 8",
        },
        {
            "stage": "0D coarsening",
            "transformation": "harmonic averaging of r^4 within bins",
            "risk": "low for total resistance; nodal positions need explicit provenance",
        },
        {
            "stage": "Pressure reference",
            "transformation": "pressure differences remove OpenFOAM gauge offset",
            "risk": "low when the same two sections are used",
        },
    ]
    write_rows(
        OUT / "interpolation_transformations.csv",
        list(transforms[0]),
        transforms,
    )

    # Geometry overlay.
    fig = plt.figure(figsize=(9.2, 7.2))
    ax = fig.add_subplot(111, projection="3d")
    surface_points = vtk_to_numpy(
        surface_reader.GetOutput().GetPoints().GetData()
    )
    ax.scatter(
        surface_points[::20, 0],
        surface_points[::20, 1],
        surface_points[::20, 2],
        s=1,
        alpha=0.12,
        color="#777777",
        label="CFD surface",
    )
    ax.plot(*path_points.T, color="#777777", lw=1.5, label="full ROM main path")
    common = (path_s >= s_in) & (path_s <= s_out)
    ax.plot(
        *path_points[common].T,
        color="#0072B2",
        lw=3.0,
        label="slice-to-slice common path",
    )
    side_ids = np.setdiff1d(np.arange(len(centerline["points"])), centerline["path"])
    if len(side_ids):
        ax.scatter(
            *centerline["points"][side_ids].T,
            color="#E69F00",
            s=25,
            label="candidate side branch",
        )
    for name, plane in PLANES.items():
        origin_mm = plane["origin_m"] * 1e3
        normal = plane["normal"] / np.linalg.norm(plane["normal"])
        ax.scatter(*origin_mm, s=34, label=f"{name} plane origin")
        ax.quiver(*origin_mm, *(normal * 2.0), length=1.0, normalize=False)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    ax.set_title("CFD Surface, Candidate Centreline, and Slice Definitions")
    ax.legend(fontsize=7, loc="best")
    fig.savefig(OUT / "geometry_centreline_plane_overlay.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.2, 6.5), sharex=True)
    ax1.plot(sample_s, cfd_radius, label="CFD surface-cut equivalent radius", lw=2.1)
    ax1.plot(sample_s, rom_radius, label="ROM ce_radius", lw=1.9)
    ax1.set_ylabel("Radius [mm]")
    ax1.legend()
    ax1.grid(alpha=0.25)
    difference = cfd_radius - rom_radius
    ax2.plot(sample_s, difference, color="#CC6677", lw=2)
    ax2.axhline(0.0, color="black", lw=1)
    ax2.set_xlabel("Global centreline arc length [mm]")
    ax2.set_ylabel("CFD - ROM radius [mm]")
    ax2.grid(alpha=0.25)
    fig.suptitle("Radius Consistency on the Common CFD Domain")
    fig.savefig(OUT / "radius_consistency.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    x = np.arange(len(pressure_rows))
    methods = [
        ("area_weighted_pressure_Pa", "area weighted"),
        ("raw_slice_unweighted_mean_Pa", "raw slice mean"),
        ("point_export_unweighted_mean_Pa", "point-export mean"),
        ("centreline_probe_Pa", "centreline probe"),
    ]
    for key, label in methods:
        ax.plot(x, [row[key] for row in pressure_rows], marker="o", label=label)
    ax.set_xticks(x, [row["slice"] for row in pressure_rows])
    ax.set_ylabel("Pressure [Pa]")
    ax.set_title("Pressure Extraction Methods at Identical Sections")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(OUT / "pressure_extraction_methods.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    labels = [row["method"].split(":")[0] for row in resistance_rows]
    values = [row["R_Pa_s_m3"] * 1e-8 for row in resistance_rows]
    bars = ax.bar(labels, values, color=["#0072B2", "#CC6677", "#E69F00", "#009E73"])
    for bar, row in zip(bars, resistance_rows):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.03,
            f"{row['relative_to_CFD_pct']:+.1f}%",
            ha="center",
        )
    ax.set_ylabel("Resistance [10^8 Pa s/m^3]")
    ax.set_title("Independent Hydraulic Resistance Checks")
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(OUT / "resistance_methods.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    make_dependency_graph()

    file_rows = []
    for path in [
        ROOT / "data/graph_mr_025.vtp",
        ROOT / "data/mr_limited/graph_mr_025.vtp",
        CENTERLINE,
        ROOT / "data/mr_limited/geometry/segment_test.stl",
        CFD_SURFACE,
    ]:
        file_rows.append(
            {
                "file": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    write_rows(OUT / "geometry_file_hashes.csv", list(file_rows[0]), file_rows)

    print(f"Audit outputs written to {OUT}")
    print(f"Common domain: {s_in:.6f}-{s_out:.6f} mm")
    print(f"R_CFD={r_cfd:.6e}, R_surface={r_surface:.6e}, R_ROM={r_rom_native:.6e}")


if __name__ == "__main__":
    main()
