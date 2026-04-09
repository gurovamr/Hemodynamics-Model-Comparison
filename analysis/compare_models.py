"""
Comparison of 3-D CFD, 1-D and 0-D Blood-Flow Models
=====================================================
Loads centre-line results from all three modelling levels and performs:

  1. Qualitative overlay plots (velocity, pressure along arc-length)
  2. Quantitative error metrics (RMSE, max error, R²) vs. the 3-D CFD baseline

Usage
-----
    python analysis/compare_models.py \\
        --cfd3d  data/centerlines/cow_mca_centerline.csv \\
        --rom1d  data/centerlines/cow_1d_results.csv \\
        --rom0d  data/centerlines/cow_0d_results.csv \\
        --out    data/centerlines/comparison_summary.png

The script can also be run without any arguments; it will auto-detect
available files and produce whatever comparison is possible.
"""

import argparse
import os
import sys
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_csv(path: str) -> Optional[dict]:
    """Load CSV with header. Returns None if file does not exist."""
    if not os.path.isfile(path):
        print(f"[INFO] File not found (skipping): {path}")
        return None
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {name: data[name] for name in data.dtype.names}


def safe_col(data: Optional[dict], col: str) -> Optional[np.ndarray]:
    """Return column array or None if not available."""
    if data is None or col not in data:
        return None
    arr = data[col]
    if np.all(np.isnan(arr)):
        return None
    return arr


# ---------------------------------------------------------------------------
# Error metrics
# ---------------------------------------------------------------------------

def rmse(ref: np.ndarray, pred: np.ndarray) -> float:
    return float(np.sqrt(np.mean((ref - pred) ** 2)))


def max_error(ref: np.ndarray, pred: np.ndarray) -> float:
    return float(np.max(np.abs(ref - pred)))


def r_squared(ref: np.ndarray, pred: np.ndarray) -> float:
    ss_res = np.sum((ref - pred) ** 2)
    ss_tot = np.sum((ref - np.mean(ref)) ** 2)
    return float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")


def interpolate_to(x_target, x_src, y_src):
    """Linearly interpolate y_src(x_src) onto x_target."""
    return np.interp(x_target, x_src, y_src)


# ---------------------------------------------------------------------------
# Comparison routine
# ---------------------------------------------------------------------------

def compare(cfd3d_path: str, rom1d_path: str, rom0d_path: str, out_path: str) -> None:
    d3  = load_csv(cfd3d_path)
    d1  = load_csv(rom1d_path)
    d0  = load_csv(rom0d_path)

    if d3 is None:
        print("[WARNING] No 3-D CFD data available. Comparison requires a baseline.")

    # Arc-length from 3-D data (or 1-D if 3-D is missing)
    s_ref = safe_col(d3, "arc_length") if d3 else safe_col(d1, "arc_length_m")
    if s_ref is None:
        print("[ERROR] No arc-length data found in any input file.")
        sys.exit(1)

    s_mm = s_ref * 1e3  # m → mm

    # Collect velocity and pressure from each model
    vel3d  = safe_col(d3, "velocity_3d")
    pres3d = safe_col(d3, "pressure_3d")

    # 1-D data may use different column names
    vel1d_raw  = safe_col(d1, "velocity_1d_ms") if safe_col(d1, "velocity_1d_ms") is not None else safe_col(d1, "velocity_1d")
    pres1d_raw = safe_col(d1, "pressure_1d_Pa") if safe_col(d1, "pressure_1d_Pa") is not None else safe_col(d1, "pressure_1d")
    s1d        = safe_col(d1, "arc_length_m")   if safe_col(d1, "arc_length_m")   is not None else safe_col(d1, "arc_length")

    # Interpolate 1-D onto 3-D arc-length if needed
    vel1d  = interpolate_to(s_ref, s1d, vel1d_raw)  if (vel1d_raw  is not None and s1d is not None) else None
    pres1d = interpolate_to(s_ref, s1d, pres1d_raw) if (pres1d_raw is not None and s1d is not None) else None

    # 0-D gives a single time-series; extract mean pressure as a scalar
    pres0d_scalar = None
    if d0 is not None:
        p0_arr = safe_col(d0, "inlet_pressure_Pa") if safe_col(d0, "inlet_pressure_Pa") is not None else safe_col(d0, "compliance_pressure_Pa")
        if p0_arr is not None:
            pres0d_scalar = float(np.mean(p0_arr))

    # ---- print metrics ----
    mmhg = 133.322
    print("\n=== Model Comparison ===")
    print(f"{'Metric':<30} {'1-D model':>14} {'0-D model':>14}")
    print("-" * 60)

    if vel3d is not None and vel1d is not None:
        print(f"{'Velocity RMSE [m/s]':<30} {rmse(vel3d, vel1d):>14.4f} {'N/A':>14}")
        print(f"{'Velocity max error [m/s]':<30} {max_error(vel3d, vel1d):>14.4f} {'N/A':>14}")
        print(f"{'Velocity R²':<30} {r_squared(vel3d, vel1d):>14.4f} {'N/A':>14}")

    if pres3d is not None and pres1d is not None:
        p0_str = f"{abs(np.mean(pres3d) - pres0d_scalar)/mmhg:.2f} mmHg" if pres0d_scalar else "N/A"
        print(f"{'Pressure RMSE [Pa]':<30} {rmse(pres3d, pres1d):>14.2f} {'N/A':>14}")
        print(f"{'Mean pressure error [mmHg]':<30} {'N/A':>14} {p0_str:>14}")
        print(f"{'Pressure R²':<30} {r_squared(pres3d, pres1d):>14.4f} {'N/A':>14}")

    # ---- plots ----
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # -- velocity --
    ax = axes[0]
    if vel3d is not None:
        ax.plot(s_mm, vel3d, label="3-D CFD", color="tab:blue", linewidth=2)
    if vel1d is not None:
        ax.plot(s_mm, vel1d, label="1-D model", color="tab:orange",
                linestyle="--", linewidth=1.5)
    ax.set_xlabel("Arc length  [mm]")
    ax.set_ylabel("Velocity  [m/s]")
    ax.set_title("Velocity comparison")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)

    # -- pressure --
    ax = axes[1]
    if pres3d is not None:
        ax.plot(s_mm, pres3d / mmhg, label="3-D CFD", color="tab:blue", linewidth=2)
    if pres1d is not None:
        ax.plot(s_mm, pres1d / mmhg, label="1-D model", color="tab:orange",
                linestyle="--", linewidth=1.5)
    if pres0d_scalar is not None:
        ax.axhline(pres0d_scalar / mmhg, label="0-D mean pressure",
                   color="tab:green", linestyle=":", linewidth=1.5)
    ax.set_xlabel("Arc length  [mm]")
    ax.set_ylabel("Pressure  [mmHg]")
    ax.set_title("Pressure comparison")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)

    plt.suptitle("Circle of Willis – 3D / 1D / 0D comparison", fontsize=13,
                 fontweight="bold")
    plt.tight_layout()

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    plt.savefig(out_path, dpi=150)
    print(f"\nFigure saved → {out_path}")
    plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compare 3-D CFD, 1-D and 0-D blood-flow model results."
    )
    p.add_argument(
        "--cfd3d", default="data/centerlines/cow_mca_centerline.csv",
        help="3-D CFD centre-line CSV (default: data/centerlines/cow_mca_centerline.csv)",
    )
    p.add_argument(
        "--rom1d", default="data/centerlines/cow_1d_results.csv",
        help="1-D ROM results CSV (default: data/centerlines/cow_1d_results.csv)",
    )
    p.add_argument(
        "--rom0d", default="data/centerlines/cow_0d_results.csv",
        help="0-D ROM results CSV (default: data/centerlines/cow_0d_results.csv)",
    )
    p.add_argument(
        "--out", default="data/centerlines/comparison_summary.png",
        help="Output comparison figure (default: data/centerlines/comparison_summary.png)",
    )
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    compare(args.cfd3d, args.rom1d, args.rom0d, args.out)
