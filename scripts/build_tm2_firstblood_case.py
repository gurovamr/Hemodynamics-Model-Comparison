#!/usr/bin/env python3
"""Generate the standalone FirstBlood input files for the TM2 vessel.

The geometry source is the existing CFD-derived FirstBlood parameter table.
No Abel_ref2 circulation data are used.
"""

from __future__ import annotations

import csv
import bisect
import argparse
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "output" / "reduced_order_models" / "firstblood_segment_parameters.csv"
COMMON_DOMAIN = ROOT / "output" / "reduced_order_models" / "model_comparison_common_domain.csv"

RHO_KG_M3 = 1060.0
NU_M2_S = 3.3e-6
Q_INLET_ML_S = 1.0213744738922028
OUTLET_GAUGE_PA = 0.0

# A cerebral-artery wall is commonly O(0.1D). This value is not measured in
# TM2, so a constant h/D ratio is used and explicitly treated as an assumption.
WALL_THICKNESS_TO_DIAMETER = 0.10

# This is a stiff physiological-scale modulus. It limits expected
# pressure-driven area changes to a negligible level while avoiding the tiny
# time steps produced by an effectively infinite stiffness.
YOUNGS_MODULUS_PA = 1.0e6

# FirstBlood's viscous source already has the circular Poiseuille coefficient.
FRICTION_MULTIPLIER = 1.0
MODEL_SEGMENTS = 40
DIVISION_POINTS_PER_SEGMENT = 5
SIMULATION_TIME_S = 1.0


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-name", default="TM2_tortuous_standalone")
    parser.add_argument("--density", type=float, default=RHO_KG_M3)
    parser.add_argument("--nu", type=float, default=NU_M2_S)
    parser.add_argument("--wall-ratio", type=float, default=WALL_THICKNESS_TO_DIAMETER)
    parser.add_argument("--youngs-modulus", type=float, default=YOUNGS_MODULUS_PA)
    parser.add_argument("--friction-multiplier", type=float, default=FRICTION_MULTIPLIER)
    parser.add_argument("--division-points", type=int, default=DIVISION_POINTS_PER_SEGMENT)
    parser.add_argument("--simulation-time", type=float, default=SIMULATION_TIME_S)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    case = ROOT / "first_blood" / "models" / args.case_name
    rows = read_rows(SOURCE)
    if not rows:
        raise RuntimeError(f"No geometry rows found in {SOURCE}")
    domain_rows = read_rows(COMMON_DOMAIN)
    cfd_domain = next(row for row in domain_rows if row["model"] == "CFD")
    domain_start = float(cfd_domain["s_in_global_mm"]) * 1.0e-3
    domain_end = float(cfd_domain["s_out_global_mm"]) * 1.0e-3

    case.mkdir(parents=True, exist_ok=True)

    main_csv = (
        "run,forward\n"
        f"time,{args.simulation_time:.9g}\n"
        "material,linear\n"
        "solver,maccormack\n"
        "\n"
        "type,name,main node,model node,...\n"
        "moc,tm2_vessel\n"
    )
    (case / "main.csv").write_text(main_csv, encoding="utf-8")

    source_s = [0.0]
    source_d = [float(rows[0]["diameter_prox_m"])]
    for row in rows:
        source_s.append(float(row["s_end_m"]))
        source_d.append(float(row["diameter_dist_m"]))
    if domain_start < source_s[0] or domain_end > source_s[-1]:
        raise RuntimeError("CFD comparison domain is outside the geometry profile.")
    total_length = domain_end - domain_start

    def diameter_at(s_query: float) -> float:
        if s_query <= source_s[0]:
            return source_d[0]
        if s_query >= source_s[-1]:
            return source_d[-1]
        j = bisect.bisect_right(source_s, s_query)
        s0, s1 = source_s[j - 1], source_s[j]
        d0, d1 = source_d[j - 1], source_d[j]
        weight = (s_query - s0) / (s1 - s0)
        return d0 + weight * (d1 - d0)

    model_bounds = [
        domain_start + total_length * i / MODEL_SEGMENTS
        for i in range(MODEL_SEGMENTS + 1)
    ]

    header = (
        "type,ID,name,start_node,end_node,start_diameter[SI],"
        "end_diameter[SI],start_thickness[SI],end_thickness[SI],"
        "length[SI],division_points,elasticity[SI],res_start[SI],"
        "res_end[SI],visc_fact[1],k1[SI],k2[SI],k3[SI]"
    )
    lines = [header]
    for i in range(MODEL_SEGMENTS):
        d0 = diameter_at(model_bounds[i])
        d1 = diameter_at(model_bounds[i + 1])
        length = model_bounds[i + 1] - model_bounds[i]
        lines.append(
            ",".join(
                [
                    "vis_f",
                    f"V{i:03d}",
                    f"TM2 segment {i:03d}",
                    f"n{i:03d}",
                    f"n{i + 1:03d}",
                    f"{d0:.16g}",
                    f"{d1:.16g}",
                    f"{args.wall_ratio * d0:.16g}",
                    f"{args.wall_ratio * d1:.16g}",
                    f"{length:.16g}",
                    str(args.division_points),
                    f"{args.youngs_modulus:.16g}",
                    "0",
                    "0",
                    f"{args.friction_multiplier:.16g}",
                    "2e6",
                    "-2253",
                    "8.65e4",
                ]
            )
        )

    lines.extend(["", "type,ID,name,value,parameter,file name"])
    lines.append("heart,n000,0,Q,inlet_flow")
    for i in range(1, MODEL_SEGMENTS):
        lines.append(f"node,n{i:03d},0,0")
    lines.append(f"perif,n{MODEL_SEGMENTS:03d},0,0")
    lines.append("")
    (case / "tm2_vessel.csv").write_text("\n".join(lines), encoding="utf-8")

    # FirstBlood expects Q in mL/s and repeats this table periodically.
    waveform = (
        f"0.0,{Q_INLET_ML_S:.16g}\n"
        f"1.0,{Q_INLET_ML_S:.16g}\n"
    )
    (case / "inlet_flow.csv").write_text(waveform, encoding="utf-8")

    assumptions = (
        "parameter,value,units,provenance\n"
        f"density,{args.density},kg/m3,TM2 CFD or sensitivity override\n"
        f"kinematic_viscosity,{args.nu},m2/s,TM2 CFD or sensitivity override\n"
        f"inlet_flow,{Q_INLET_ML_S},mL/s,CFD inlet integration\n"
        f"outlet_gauge_pressure,{OUTLET_GAUGE_PA},Pa,TM2 CFD boundary\n"
        f"wall_thickness_to_diameter,{args.wall_ratio},-,assumed\n"
        f"youngs_modulus,{args.youngs_modulus},Pa,rigid-wall approximation\n"
        f"friction_multiplier,{args.friction_multiplier},-,unmodified circular friction\n"
        f"division_points_per_segment,{args.division_points},-,numerical assumption\n"
        f"simulation_time,{args.simulation_time},s,numerical assumption\n"
        f"source_geometry_intervals,{len(rows)},-,CFD-derived profile\n"
        f"model_segments,{MODEL_SEGMENTS},-,continuous linear resampling\n"
        f"global_s_start,{domain_start:.16g},m,CFD comparison inlet\n"
        f"global_s_end,{domain_end:.16g},m,CFD comparison outlet\n"
        f"centreline_arc_length,{total_length:.16g},m,CFD-derived profile\n"
    )
    (case / "assumptions.csv").write_text(assumptions, encoding="utf-8")

    print(f"Generated {case}")
    print(f"Source intervals: {len(rows)}")
    print(f"Model segments: {MODEL_SEGMENTS}")


if __name__ == "__main__":
    main()
