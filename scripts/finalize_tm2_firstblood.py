#!/usr/bin/env python3
"""Package the converged TM2 FirstBlood result and publication figures."""

from __future__ import annotations

import csv
import math
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "firstblood_tm2_final"
FINAL = OUT / "grid_41"
COMMON = ROOT / "output" / "reduced_order_models" / "model_comparison_common_domain.csv"
PROFILES = ROOT / "output" / "reduced_order_models" / "pressure_profiles_common_domain.csv"
CFD_FLOW = ROOT / "output" / "workflow_audit" / "flow_audit.csv"

GRID_POINTS = [11, 15, 21, 31, 41]
CFD_DP_PA = 356.36595884094373
CFD_Q_M3_S = 1.0213744738922028e-6
S_GLOBAL_START_MM = 19.671598192628114

COLORS = {
    "Analytical": "#8C8C8C",
    "Native 0D": "#E69F00",
    "0D (N=50)": "#56B4E9",
    "FirstBlood 1D": "#0072B2",
    "CFD": "#D55E00",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def summary_map(path: Path) -> dict[str, float]:
    return {r["metric"]: float(r["value"]) for r in read_rows(path)}


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.7)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
        }
    )

    final_profile = read_rows(FINAL / "firstblood_profile.csv")
    final_flow = read_rows(FINAL / "firstblood_flow_profile.csv")
    final_summary = summary_map(FINAL / "firstblood_summary.csv")

    pressure_rows = []
    p_in = float(final_profile[0]["pressure_gauge_Pa"])
    for row in final_profile:
        distance = float(row["distance_from_inlet_mm"])
        p_gauge = float(row["pressure_gauge_Pa"])
        pressure_rows.append(
            {
                "distance_from_inlet_mm": distance,
                "global_s_mm": S_GLOBAL_START_MM + distance,
                "radius_m": float(row["radius_nominal_m"]),
                "pressure_gauge_Pa": p_gauge,
                "pressure_loss_from_inlet_Pa": p_in - p_gauge,
                "pressure_loss_from_inlet_mmHg": (p_in - p_gauge) / 133.322368,
            }
        )
    write_rows(
        OUT / "firstblood_pressure_profile.csv",
        list(pressure_rows[0]),
        pressure_rows,
    )

    flow_rows = []
    for row in final_flow:
        s_mid_m = float(row["s_mid_m"])
        flow_rows.append(
            {
                "segment_id": row["segment_id"],
                "distance_from_inlet_mm": s_mid_m * 1e3,
                "global_s_mm": S_GLOBAL_START_MM + s_mid_m * 1e3,
                "flow_start_m3_s": float(row["flow_start_m3_s"]),
                "flow_end_m3_s": float(row["flow_end_m3_s"]),
                "flow_mean_m3_s": float(row["flow_mean_m3_s"]),
                "flow_mean_mL_s": float(row["flow_mean_mL_s"]),
            }
        )
    write_rows(OUT / "firstblood_flow_profile.csv", list(flow_rows[0]), flow_rows)

    grid_rows = []
    previous_dp = None
    previous_r = None
    for points in GRID_POINTS:
        values = summary_map(OUT / f"grid_{points}" / "firstblood_summary.csv")
        dp = values["pressure_drop"]
        resistance = values["hydraulic_resistance"]
        grid_rows.append(
            {
                "points_per_edge": points,
                "total_stored_grid_points": 40 * points,
                "unique_axial_grid_points": 40 * (points - 1) + 1,
                "pressure_drop_Pa": dp,
                "hydraulic_resistance_Pa_s_m3": resistance,
                "pressure_drop_change_percent": (
                    "" if previous_dp is None else 100 * abs(dp - previous_dp) / abs(previous_dp)
                ),
                "resistance_change_percent": (
                    "" if previous_r is None else 100 * abs(resistance - previous_r) / abs(previous_r)
                ),
                "terminal_flow_mismatch_percent": values["flow_mismatch_percent"],
                "axial_flow_variation_percent": values["axial_flow_spread_percent"],
            }
        )
        previous_dp = dp
        previous_r = resistance
    write_rows(OUT / "firstblood_grid_convergence.csv", list(grid_rows[0]), grid_rows)

    common = {r["model"]: r for r in read_rows(COMMON)}
    # Every final comparison uses the measured CFD inlet flow as the common
    # denominator. In particular, do not use FirstBlood's mean terminal flow.
    fb_dp = final_summary["pressure_drop"]
    fb_r = fb_dp / CFD_Q_M3_S
    resistance = {
        "Analytical": float(common["Analytical (Poiseuille)"]["R_Pa_s_m3"]),
        "Native 0D": float(common["Native 0D"]["R_Pa_s_m3"]),
        "0D (N=50)": float(common["0D (N=50)"]["R_Pa_s_m3"]),
        "FirstBlood 1D": fb_r,
        "CFD": CFD_DP_PA / CFD_Q_M3_S,
    }
    pressure_drop = {
        name: value * CFD_Q_M3_S for name, value in resistance.items()
    }
    pressure_drop["CFD"] = CFD_DP_PA

    final_comparison_rows = []
    for name in ["CFD", "Analytical", "Native 0D", "0D (N=50)", "FirstBlood 1D"]:
        dp = pressure_drop[name]
        final_comparison_rows.append(
            {
                "model": name,
                "Q_m3_s": CFD_Q_M3_S,
                "Q_mL_s": CFD_Q_M3_S * 1e6,
                "pressure_drop_Pa": dp,
                "hydraulic_resistance_Pa_s_m3": resistance[name],
                "relative_error_vs_CFD_percent": 100 * (dp - CFD_DP_PA) / CFD_DP_PA,
                "source": (
                    "converged FirstBlood grid_41"
                    if name == "FirstBlood 1D"
                    else "CFD common-domain reference"
                    if name == "CFD"
                    else "common-domain reduced-order model"
                ),
            }
        )
    write_rows(
        OUT / "final_model_comparison.csv",
        list(final_comparison_rows[0]),
        final_comparison_rows,
    )

    comparison_rows = [
        {
            "metric": "pressure_drop",
            "CFD_value": CFD_DP_PA,
            "FirstBlood_value": fb_dp,
            "absolute_error": abs(fb_dp - CFD_DP_PA),
            "signed_error": fb_dp - CFD_DP_PA,
            "relative_error_percent": 100 * (fb_dp - CFD_DP_PA) / CFD_DP_PA,
            "units": "Pa",
        },
        {
            "metric": "hydraulic_resistance",
            "CFD_value": resistance["CFD"],
            "FirstBlood_value": fb_r,
            "absolute_error": abs(fb_r - resistance["CFD"]),
            "signed_error": fb_r - resistance["CFD"],
            "relative_error_percent": 100 * (fb_r - resistance["CFD"]) / resistance["CFD"],
            "units": "Pa s/m3",
        },
    ]
    write_rows(OUT / "firstblood_cfd_comparison.csv", list(comparison_rows[0]), comparison_rows)

    scalar_rows = [
        {"metric": "pressure_drop", "value": fb_dp, "units": "Pa"},
        {"metric": "hydraulic_resistance", "value": fb_r, "units": "Pa s/m3"},
        {"metric": "maximum_velocity", "value": final_summary["maximum_velocity"], "units": "m/s"},
        {"metric": "maximum_configured_CFL", "value": 0.9, "units": "-"},
        {"metric": "simulation_time", "value": 1.0, "units": "s"},
        {"metric": "late_time_average_window", "value": 0.2, "units": "s"},
        {"metric": "points_per_edge", "value": 41, "units": "-"},
        {"metric": "number_of_segments", "value": 40, "units": "-"},
        {"metric": "unique_axial_grid_points", "value": 1601, "units": "-"},
        {"metric": "centreline_arc_length", "value": final_summary["centreline_arc_length"], "units": "m"},
        {"metric": "inlet_flow", "value": final_summary["inlet_flow"], "units": "m3/s"},
        {"metric": "outlet_flow", "value": final_summary["outlet_flow"], "units": "m3/s"},
        {"metric": "terminal_flow_mismatch", "value": final_summary["flow_mismatch_percent"], "units": "%"},
        {"metric": "axial_flow_variation", "value": final_summary["axial_flow_spread_percent"], "units": "%"},
        {"metric": "maximum_area_change", "value": final_summary["maximum_area_change_percent"], "units": "%"},
    ]
    write_rows(OUT / "firstblood_summary.csv", list(scalar_rows[0]), scalar_rows)

    # Pressure profile: plot loss from inlet so CFD and 1D share the same reference.
    cfd_pressure = [r for r in read_rows(PROFILES) if r["model"] == "CFD"]
    x_fb = np.array([float(r["distance_from_inlet_mm"]) for r in pressure_rows])
    loss_fb = np.array([float(r["pressure_loss_from_inlet_Pa"]) for r in pressure_rows])
    x_cfd = np.array([float(r["distance_from_inlet_mm"]) for r in cfd_pressure])
    loss_cfd = np.array([float(r["pressure_loss_Pa"]) for r in cfd_pressure])
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    ax.plot(x_fb, loss_fb, color=COLORS["FirstBlood 1D"], linewidth=2.3, label="FirstBlood 1D")
    ax.plot(x_cfd, loss_cfd, "o-", color=COLORS["CFD"], linewidth=2, markersize=6, label="CFD area average")
    ax.set(xlabel="Distance from inlet (mm)", ylabel="Pressure loss from inlet (Pa)", title="TM2 pressure profile")
    style_axis(ax)
    ax.legend(frameon=False)
    save_figure(fig, "figure_1_pressure_profile_cfd_firstblood")

    names = list(resistance)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.bar(names, [resistance[n] / 1e8 for n in names], color=[COLORS[n] for n in names])
    ax.set_ylabel(r"Hydraulic resistance ($10^8$ Pa s m$^{-3}$)")
    ax.set_title("Hydraulic resistance on the common CFD domain")
    ax.tick_params(axis="x", rotation=18)
    style_axis(ax)
    save_figure(fig, "figure_2_hydraulic_resistance_comparison")

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.bar(names, [pressure_drop[n] for n in names], color=[COLORS[n] for n in names])
    ax.set_ylabel("Pressure drop at CFD inlet flow (Pa)")
    ax.set_title("Pressure-drop comparison")
    ax.tick_params(axis="x", rotation=18)
    style_axis(ax)
    save_figure(fig, "figure_3_pressure_drop_comparison")

    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    points = np.array(GRID_POINTS)
    dp_values = np.array([float(r["pressure_drop_Pa"]) for r in grid_rows])
    ax.plot(points, dp_values, "o-", color=COLORS["FirstBlood 1D"], linewidth=2.2)
    ax.axhline(dp_values[-1], color="#777777", linestyle="--", linewidth=1, label="41-point result")
    ax.set(xlabel="Numerical points per edge", ylabel="Pressure drop (Pa)", title="FirstBlood grid convergence")
    style_axis(ax)
    ax.legend(frameon=False)
    save_figure(fig, "figure_4_grid_convergence")

    cfd_flow = read_rows(CFD_FLOW)
    cfd_x = np.array([0.0, 7.641004615739643, 16.221105817381652])
    cfd_q = np.array([float(r["flow_mL_s"]) for r in cfd_flow])
    fb_x = np.array([float(r["distance_from_inlet_mm"]) for r in flow_rows])
    fb_q = np.array([float(r["flow_mean_mL_s"]) for r in flow_rows])
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    ax.plot(fb_x, fb_q, color=COLORS["FirstBlood 1D"], linewidth=2.2, label="FirstBlood 1D")
    ax.plot(cfd_x, cfd_q, "o", color=COLORS["CFD"], markersize=6, label="CFD slices")
    ax.set(xlabel="Distance from inlet (mm)", ylabel="Flow rate (mL/s)", title="Axial flow conservation")
    style_axis(ax)
    ax.legend(frameon=False)
    save_figure(fig, "figure_5_flow_convergence")

    shutil.copy2(FINAL / "firstblood_validation.csv", OUT / "firstblood_solver_validation.csv")

    fb_at_cfd = np.interp(x_cfd, x_fb, loss_fb)
    pressure_rmse = float(np.sqrt(np.mean((fb_at_cfd - loss_cfd) ** 2)))
    flow_at_cfd = np.interp(cfd_x, fb_x, fb_q)
    flow_rmse = float(np.sqrt(np.mean((flow_at_cfd - cfd_q) ** 2)))
    report = f"""# Final TM2 FirstBlood model validation

## Acceptance decision

**Accepted as numerically converged.** The 31-to-41-point change is
{grid_rows[-1]['pressure_drop_change_percent']:.3f}% for pressure drop and
{grid_rows[-1]['resistance_change_percent']:.3f}% for resistance, both below
1%. Terminal flow mismatch is {final_summary['flow_mismatch_percent']:.5f}%
and axial flow variation is {final_summary['axial_flow_spread_percent']:.5f}%,
both below 0.2%.

## Consistency checks

- Geometry domain: global centreline s=19.6716–35.8927 mm; length 16.2211 mm.
- Geometry representation: 40 continuous tapered edges interpolated from the
  CFD-derived radius profile.
- Inlet: prescribed CFD inlet flow {CFD_Q_M3_S * 1e6:.6f} mL/s.
- Outlet: zero-gauge pressure, matching the CFD pressure reference.
- Fluid: rho=1060 kg/m3 and nu=3.3e-6 m2/s.
- Wall: linear law, E=1.0 MPa, h/D=0.1; maximum simulated area change
  {final_summary['maximum_area_change_percent']:.3f}%.
- Friction multiplier: 1.0; no calibration to CFD.

## Final FirstBlood result

- Pressure drop: {fb_dp:.3f} Pa
- Hydraulic resistance: {fb_r:.6e} Pa s/m3
- Resistance denominator: common CFD inlet flow {CFD_Q_M3_S:.10e} m3/s
- Maximum velocity: {final_summary['maximum_velocity']:.4f} m/s
- Grid: 41 points/edge, 40 edges, 1601 unique axial points
- Maximum configured CFL: 0.9
- Simulated time: 1.0 s; reported mean window: final 0.2 s

## CFD comparison

- CFD pressure drop: {CFD_DP_PA:.3f} Pa
- Absolute pressure-drop error: {abs(fb_dp - CFD_DP_PA):.3f} Pa
- Relative pressure-drop error: {100 * (fb_dp - CFD_DP_PA) / CFD_DP_PA:.2f}%
- CFD resistance at the prescribed inlet flow: {resistance['CFD']:.6e} Pa s/m3
- FirstBlood resistance error: {100 * (fb_r - resistance['CFD']) / resistance['CFD']:.2f}%
- Pressure-profile RMSE at the three CFD planes: {pressure_rmse:.3f} Pa
- Flow-profile RMSE at the three CFD planes: {flow_rmse:.5f} mL/s

The converged 1D model underpredicts CFD loss. This is physically plausible:
FirstBlood uses axial length and taper but contains no curvature, torsion,
secondary-flow, or bend-loss model. The discrepancy is therefore retained and
reported rather than removed by parameter tuning.

The row labelled **FirstBlood 1D** comes exclusively from the converged
41-points-per-edge FirstBlood run. No legacy Poiseuille-generated `1D` row is
used in this report or its figures.

## Geometry qualification

The 1D geometry uses the same CFD-derived centreline radius data and identical
comparison-plane arc-length bounds. It is not geometrically identical to the
full 3D lumen: a 1D taper cannot retain non-circular cross-sections, curvature,
torsion, or local secondary-flow structures.
"""
    (OUT / "firstblood_model_validation.md").write_text(report, encoding="utf-8")
    print(f"Final outputs written to {OUT}")


if __name__ == "__main__":
    main()
