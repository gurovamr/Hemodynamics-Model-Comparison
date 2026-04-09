"""
Plot Velocity and Pressure Profiles along a Centre-line
=======================================================
Reads a centre-line CSV (produced by extract_centerline.py or the ROMs)
and generates publication-quality figures.

Usage
-----
    python analysis/plot_centerline.py \\
        --csv  data/centerlines/cow_mca_centerline.csv \\
        --out  data/centerlines/cow_centerline_profiles.png

Expected CSV columns (any subset is accepted):
    arc_length, velocity_3d, pressure_3d,
    velocity_1d, pressure_1d, pressure_0d
"""

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_csv(path: str) -> dict:
    """Load a comma-separated file with a header row. Returns dict of arrays."""
    if not os.path.isfile(path):
        print(f"[ERROR] File not found: {path}")
        sys.exit(1)
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {name: data[name] for name in data.dtype.names}


def column_present(data: dict, col: str) -> bool:
    return col in data and not np.all(np.isnan(data[col]))


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_profiles(data: dict, out_path: str) -> None:
    s = data["arc_length"] * 1e3   # m → mm

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # ---- Velocity ----
    ax = axes[0]
    if column_present(data, "velocity_3d"):
        ax.plot(s, data["velocity_3d"], label="3-D CFD", color="tab:blue",
                linewidth=2)
    if column_present(data, "velocity_1d"):
        ax.plot(s, data["velocity_1d"], label="1-D model", color="tab:orange",
                linestyle="--", linewidth=1.5)
    ax.set_xlabel("Arc length  [mm]")
    ax.set_ylabel("Velocity  [m/s]")
    ax.set_title("Centre-line velocity – MCA")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)

    # ---- Pressure ----
    ax = axes[1]
    mmhg = 133.322
    if column_present(data, "pressure_3d"):
        ax.plot(s, data["pressure_3d"] / mmhg, label="3-D CFD", color="tab:blue",
                linewidth=2)
    if column_present(data, "pressure_1d"):
        ax.plot(s, data["pressure_1d"] / mmhg, label="1-D model",
                color="tab:orange", linestyle="--", linewidth=1.5)
    if column_present(data, "pressure_0d"):
        ax.plot(
            [s[0], s[-1]],
            [data["pressure_0d"][0] / mmhg, data["pressure_0d"][-1] / mmhg],
            label="0-D Windkessel", color="tab:green",
            linestyle=":", linewidth=1.5,
            marker="o", markersize=5,
        )
    ax.set_xlabel("Arc length  [mm]")
    ax.set_ylabel("Pressure  [mmHg]")
    ax.set_title("Centre-line pressure – MCA")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)

    plt.suptitle("Circle of Willis – Multiscale model comparison", fontsize=13,
                 fontweight="bold")
    plt.tight_layout()

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    plt.savefig(out_path, dpi=150)
    print(f"Figure saved → {out_path}")
    plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Plot centre-line velocity and pressure profiles."
    )
    p.add_argument(
        "--csv", default="data/centerlines/cow_mca_centerline.csv",
        help="Path to centre-line CSV file",
    )
    p.add_argument(
        "--out", default="data/centerlines/cow_centerline_profiles.png",
        help="Output figure path (PNG)",
    )
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    data = load_csv(args.csv)
    plot_profiles(data, args.out)
