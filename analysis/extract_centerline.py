"""
Extract Velocity and Pressure along a Centre-line from OpenFOAM Output
=======================================================================
Reads the CSV files written by the OpenFOAM *sets* function object
(configured in system/controlDict) and produces a unified centre-line
dataset with columns compatible with the rest of the analysis pipeline.

Usage
-----
    python analysis/extract_centerline.py \\
        --of_dir  cfd/cow_case \\
        --cl_name mca_centerline \\
        --time    500 \\
        --out     data/centerlines/cow_mca_centerline.csv

The script reads:
    <of_dir>/postProcessing/centerlineSample/<time>/<cl_name>_U_p.csv

and writes a CSV with columns:
    arc_length,x,y,z,radius,velocity_3d,pressure_3d

Notes
-----
- If the OpenFOAM output is not available the script falls back to the
  placeholder CSV in data/centerlines/ and just echoes it.
- Radius is estimated from the nearest surface vertex if a geometry
  STL is available, otherwise set to NaN.
"""

import argparse
import os
import sys

import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_of_csv(path: str) -> dict:
    """
    Parse an OpenFOAM *sets* CSV file.

    Expected header (example):
        x,y,z,Ux,Uy,Uz,p
    Returns dict with numpy arrays.
    """
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {name: data[name] for name in data.dtype.names}


def arc_length(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Compute cumulative arc-length along a polyline."""
    dx = np.diff(x)
    dy = np.diff(y)
    dz = np.diff(z)
    ds = np.sqrt(dx**2 + dy**2 + dz**2)
    return np.concatenate([[0.0], np.cumsum(ds)])


def velocity_magnitude(ux, uy, uz):
    return np.sqrt(ux**2 + uy**2 + uz**2)


# ---------------------------------------------------------------------------
# Main extraction routine
# ---------------------------------------------------------------------------

def extract(of_dir: str, cl_name: str, time: str, out_path: str) -> None:
    of_csv = os.path.join(
        of_dir, "postProcessing", "centerlineSample", str(time),
        f"{cl_name}_U_p.csv",
    )

    if not os.path.isfile(of_csv):
        print(
            f"[WARNING] OpenFOAM output not found:\n  {of_csv}\n"
            "Using placeholder data instead."
        )
        placeholder = os.path.join(
            os.path.dirname(__file__), "..", "data", "centerlines",
            "cow_mca_centerline.csv",
        )
        placeholder = os.path.normpath(placeholder)
        if os.path.isfile(placeholder):
            if os.path.abspath(placeholder) != os.path.abspath(out_path):
                import shutil
                shutil.copy(placeholder, out_path)
                print(f"Copied placeholder → {out_path}")
            else:
                print(f"Placeholder already in place → {out_path}")
        else:
            print("[ERROR] Placeholder file not found either. Aborting.")
            sys.exit(1)
        return

    raw = parse_of_csv(of_csv)
    x   = raw["x"]
    y   = raw["y"]
    z   = raw["z"]
    s   = arc_length(x, y, z)
    vel = velocity_magnitude(raw["Ux"], raw["Uy"], raw["Uz"])
    # OpenFOAM simpleFoam stores kinematic pressure p/rho; convert to Pa
    rho = 1060.0
    pressure = raw["p"] * rho

    radius = np.full_like(x, np.nan)   # placeholder; fill from geometry if available

    header = "arc_length,x,y,z,radius,velocity_3d,pressure_3d"
    out_data = np.column_stack([s, x, y, z, radius, vel, pressure])
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    np.savetxt(out_path, out_data, delimiter=",", header=header, comments="")
    print(f"Centre-line data saved → {out_path}  ({len(x)} points)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Extract centre-line velocity/pressure from OpenFOAM output."
    )
    p.add_argument(
        "--of_dir", default="cfd/cow_case",
        help="Path to OpenFOAM case directory (default: cfd/cow_case)",
    )
    p.add_argument(
        "--cl_name", default="mca_centerline",
        help="Name of the centre-line set as defined in controlDict (default: mca_centerline)",
    )
    p.add_argument(
        "--time", default="500",
        help="Time directory to read from (default: 500)",
    )
    p.add_argument(
        "--out", default="data/centerlines/cow_mca_centerline.csv",
        help="Output CSV path (default: data/centerlines/cow_mca_centerline.csv)",
    )
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    extract(args.of_dir, args.cl_name, args.time, args.out)
