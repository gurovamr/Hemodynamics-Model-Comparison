#!/usr/bin/env python3
"""Prototype steady 1D rigid-wall solver for single-vessel CFD-ROM workflow.

This script reuses exported centerline geometry and slice pressure/flow summaries.
It solves the steady 1D rigid-wall momentum balance on a prescribed area profile:

    dQ/ds = 0
    dp/ds = alpha * rho * Q^2 / A^3 * dA/ds - 8 * pi * mu * Q / A^2

Unknown Q is obtained from inlet/outlet pressure boundary conditions by bisection.
Then p(s) is reconstructed and compared against:
- existing analytical Poiseuille profile from pressure_profiles_common_domain.csv
- CFD area-averaged slice pressures on the common domain
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "workflow_audit"
OUT.mkdir(parents=True, exist_ok=True)

GEOM_PROFILE = ROOT / "output" / "reduced_order_models" / "segment_1D_geometry_profile.csv"
COMMON_DOMAIN = ROOT / "output" / "reduced_order_models" / "model_comparison_common_domain.csv"
PRESSURE_PROFILES = ROOT / "output" / "reduced_order_models" / "pressure_profiles_common_domain.csv"
PRESSURE_SLICES = ROOT / "output" / "06_clear_results" / "pressure_slice_summary.csv"
FLOW_SLICES = ROOT / "output" / "06_clear_results" / "flow_rate_summary.csv"

MMHG_PER_PA = 1.0 / 133.322368


@dataclass
class GeometryDomain:
    s: np.ndarray
    r: np.ndarray
    a: np.ndarray


@dataclass
class CfdSlices:
    s_global_mm: np.ndarray
    p_pa: np.ndarray
    q_m3_s: float


def read_rows(path: Path) -> List[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, fieldnames: Iterable[str], rows: Iterable[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        writer.writerows(rows)


def load_common_domain() -> Tuple[float, float]:
    rows = read_rows(COMMON_DOMAIN)
    cfd = next(row for row in rows if row["model"] == "CFD")
    return float(cfd["s_in_global_mm"]) * 1e-3, float(cfd["s_out_global_mm"]) * 1e-3


def load_geometry_domain(s_in: float, s_out: float) -> GeometryDomain:
    rows = read_rows(GEOM_PROFILE)
    s = np.array([float(r["s_m"]) for r in rows], dtype=float)
    r = np.array([float(r["r_m"]) for r in rows], dtype=float)

    # Some exported comparison tables extend the outlet by a tiny interpolation
    # offset; clamp if the mismatch is within sub-mm tolerance.
    tol = 5e-4  # m
    if s_in < s[0] - tol or s_out > s[-1] + tol:
        raise ValueError("Common-domain bounds are outside geometry profile range.")
    s_in = max(s_in, s[0])
    s_out = min(s_out, s[-1])

    mask = (s >= s_in) & (s <= s_out)
    s_mid = s[mask]
    r_mid = r[mask]

    # Ensure exact endpoints are represented by interpolation.
    s_use = np.r_[s_in, s_mid[(s_mid > s_in) & (s_mid < s_out)], s_out]
    r_use = np.interp(s_use, s, r)
    a_use = math.pi * r_use * r_use
    return GeometryDomain(s=s_use, r=r_use, a=a_use)


def load_cfd_slices(s_in: float) -> CfdSlices:
    rows = read_rows(PRESSURE_PROFILES)
    cfd_rows = [r for r in rows if r["model"] == "CFD"]
    cfd_rows.sort(key=lambda r: float(r["global_s_mm"]))

    s_global_mm = np.array([float(r["global_s_mm"]) for r in cfd_rows], dtype=float)
    p_loss_pa = np.array([float(r["pressure_loss_Pa"]) for r in cfd_rows], dtype=float)

    p_slice = read_rows(PRESSURE_SLICES)
    p_map = {r["Slice"]: float(r["Average pressure [Pa]"]) for r in p_slice}
    p_in = p_map["inlet"]
    p_out = p_map["outlet"]
    # pressure_profiles_common_domain uses inlet as reference p_loss=0.
    # Recover absolute pressures in the common domain from inlet pressure and loss.
    p_pa = p_in - p_loss_pa
    # Force last point to outlet pressure for consistency with summary table.
    p_pa[-1] = p_out

    q_rows = read_rows(FLOW_SLICES)
    q_m3_s = float(np.mean([float(r["Flow rate [m³/s]"]) for r in q_rows]))

    return CfdSlices(s_global_mm=s_global_mm, p_pa=p_pa, q_m3_s=q_m3_s)


def solve_pressure_profile_for_q(
    geom: GeometryDomain,
    q_m3_s: float,
    p_in_pa: float,
    rho: float,
    mu: float,
    alpha: float,
) -> np.ndarray:
    n = len(geom.s)
    p = np.zeros(n, dtype=float)
    p[0] = p_in_pa

    for i in range(n - 1):
        ds = geom.s[i + 1] - geom.s[i]
        ai = geom.a[i]
        aj = geom.a[i + 1]
        ui = q_m3_s / ai
        uj = q_m3_s / aj
        r_seg = 0.5 * (geom.r[i] + geom.r[i + 1])

        # Convective (kinetic-energy) correction from area change.
        dp_conv = 0.5 * alpha * rho * (ui * ui - uj * uj)
        # Viscous loss on each segment (Poiseuille law, rigid tube).
        dp_visc = 8.0 * mu * q_m3_s * ds / (math.pi * r_seg**4)

        p[i + 1] = p[i] + dp_conv - dp_visc

    return p


def predicted_delta_p(
    geom: GeometryDomain,
    q_m3_s: float,
    rho: float,
    mu: float,
    alpha: float,
) -> float:
    p = solve_pressure_profile_for_q(
        geom=geom,
        q_m3_s=q_m3_s,
        p_in_pa=0.0,
        rho=rho,
        mu=mu,
        alpha=alpha,
    )
    return p[0] - p[-1]


def solve_q_from_delta_p(
    geom: GeometryDomain,
    dp_target_pa: float,
    rho: float,
    mu: float,
    alpha: float,
) -> float:
    def f(q: float) -> float:
        return predicted_delta_p(geom, q, rho, mu, alpha) - dp_target_pa

    # Adaptive bracket search on a log grid. The model can be non-monotonic
    # over very wide flow ranges because the convective term scales with Q^2.
    q_grid = np.logspace(-10, -3, 200)
    f_grid = [f(float(q)) for q in q_grid]
    q_lo = None
    q_hi = None
    flo = None
    fhi = None
    for i in range(len(q_grid) - 1):
        if f_grid[i] == 0.0:
            return float(q_grid[i])
        if f_grid[i] * f_grid[i + 1] <= 0.0:
            q_lo = float(q_grid[i])
            q_hi = float(q_grid[i + 1])
            flo = float(f_grid[i])
            fhi = float(f_grid[i + 1])
            break

    if q_lo is None or q_hi is None:
        raise RuntimeError("Could not bracket flow solution for steady 1D equation.")

    for _ in range(120):
        q_mid = 0.5 * (q_lo + q_hi)
        fm = f(q_mid)
        if abs(fm) < 1e-9:
            return q_mid
        if flo * fm <= 0:
            q_hi = q_mid
            fhi = fm
        else:
            q_lo = q_mid
            flo = fm
    return 0.5 * (q_lo + q_hi)


def interpolate_at_positions(x: np.ndarray, y: np.ndarray, xq: np.ndarray) -> np.ndarray:
    xq_clip = np.clip(xq, x[0], x[-1])
    return np.interp(xq_clip, x, y)


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a - b) ** 2)))


def main() -> None:
    # Fluid assumptions for prototype steady 1D solver.
    rho = 1060.0
    mu = 3.5e-3
    alpha = 4.0 / 3.0  # laminar profile correction

    s_in, s_out = load_common_domain()
    geom = load_geometry_domain(s_in=s_in, s_out=s_out)
    cfd = load_cfd_slices(s_in=s_in)

    p_in = float(cfd.p_pa[0])
    p_out = float(cfd.p_pa[-1])
    dp_cfd = p_in - p_out
    r_cfd = dp_cfd / cfd.q_m3_s

    # Solve unknown Q from pressure BCs.
    q_1d = solve_q_from_delta_p(geom, dp_target_pa=dp_cfd, rho=rho, mu=mu, alpha=alpha)
    p_1d = solve_pressure_profile_for_q(geom, q_m3_s=q_1d, p_in_pa=p_in, rho=rho, mu=mu, alpha=alpha)
    dp_1d = p_1d[0] - p_1d[-1]
    r_1d = dp_1d / q_1d

    # Also evaluate the same 1D equation at the CFD-measured flow rate.
    p_1d_qcfd = solve_pressure_profile_for_q(
        geom, q_m3_s=cfd.q_m3_s, p_in_pa=p_in, rho=rho, mu=mu, alpha=alpha
    )
    dp_1d_qcfd = p_1d_qcfd[0] - p_1d_qcfd[-1]
    r_1d_qcfd = dp_1d_qcfd / cfd.q_m3_s

    # Current Poiseuille baseline from existing exported comparison table.
    common_rows = read_rows(COMMON_DOMAIN)
    pois_row = next(r for r in common_rows if r["model"] == "Analytical (Poiseuille)")
    r_pois = float(pois_row["R_Pa_s_m3"])
    dp_pois_at_qcfd = r_pois * cfd.q_m3_s

    # Profile agreement against CFD (3 measured area-averaged slices on common domain).
    s_cfd_m = cfd.s_global_mm * 1e-3
    p_1d_at_cfd = interpolate_at_positions(geom.s, p_1d, s_cfd_m)
    p_1d_qcfd_at_cfd = interpolate_at_positions(geom.s, p_1d_qcfd, s_cfd_m)

    # Reconstruct Poiseuille profile at same positions from exported file.
    pp_rows = [r for r in read_rows(PRESSURE_PROFILES) if r["model"] == "Analytical (Poiseuille)"]
    s_pp = np.array([float(r["global_s_mm"]) * 1e-3 for r in pp_rows], dtype=float)
    p_loss_pp = np.array([float(r["pressure_loss_Pa"]) for r in pp_rows], dtype=float)
    p_pp = p_in - p_loss_pp
    p_pois_at_cfd = interpolate_at_positions(s_pp, p_pp, s_cfd_m)

    rmse_1d = rmse(p_1d_at_cfd, cfd.p_pa)
    rmse_1d_qcfd = rmse(p_1d_qcfd_at_cfd, cfd.p_pa)
    rmse_pois = rmse(p_pois_at_cfd, cfd.p_pa)

    # Export profile and scalar summary.
    profile_rows = []
    for si, ri, ai, pi in zip(geom.s, geom.r, geom.a, p_1d):
        profile_rows.append(
            {
                "s_m": si,
                "distance_from_inlet_mm": (si - geom.s[0]) * 1e3,
                "radius_m": ri,
                "area_m2": ai,
                "pressure_1d_pa": pi,
                "pressure_1d_mmhg": pi * MMHG_PER_PA,
            }
        )
    write_rows(
        OUT / "steady_1d_rigid_profile.csv",
        fieldnames=list(profile_rows[0].keys()),
        rows=profile_rows,
    )

    cfd_pts_rows = []
    for s_mm, p_cfd, p1, pp in zip(cfd.s_global_mm, cfd.p_pa, p_1d_at_cfd, p_pois_at_cfd):
        cfd_pts_rows.append(
            {
                "global_s_mm": s_mm,
                "p_cfd_pa": p_cfd,
                "p_1d_pa": p1,
                "p_poiseuille_pa": pp,
                "err_1d_pa": p1 - p_cfd,
                "err_poiseuille_pa": pp - p_cfd,
            }
        )
    write_rows(
        OUT / "steady_1d_vs_cfd_points.csv",
        fieldnames=list(cfd_pts_rows[0].keys()),
        rows=cfd_pts_rows,
    )

    summary = [
        {
            "metric": "Q_CFD_m3_s",
            "value": cfd.q_m3_s,
            "units": "m3/s",
        },
        {
            "metric": "Q_1D_solved_m3_s",
            "value": q_1d,
            "units": "m3/s",
        },
        {
            "metric": "Q_error_pct_vs_CFD",
            "value": (q_1d - cfd.q_m3_s) / cfd.q_m3_s * 100.0,
            "units": "%",
        },
        {
            "metric": "DeltaP_CFD_pa",
            "value": dp_cfd,
            "units": "Pa",
        },
        {
            "metric": "DeltaP_1D_pa",
            "value": dp_1d,
            "units": "Pa",
        },
        {
            "metric": "DeltaP_1D_at_Qcfd_pa",
            "value": dp_1d_qcfd,
            "units": "Pa",
        },
        {
            "metric": "DeltaP_poiseuille_at_Qcfd_pa",
            "value": dp_pois_at_qcfd,
            "units": "Pa",
        },
        {
            "metric": "R_CFD_pa_s_m3",
            "value": r_cfd,
            "units": "Pa s/m3",
        },
        {
            "metric": "R_1D_pa_s_m3",
            "value": r_1d,
            "units": "Pa s/m3",
        },
        {
            "metric": "R_1D_at_Qcfd_pa_s_m3",
            "value": r_1d_qcfd,
            "units": "Pa s/m3",
        },
        {
            "metric": "R_poiseuille_pa_s_m3",
            "value": r_pois,
            "units": "Pa s/m3",
        },
        {
            "metric": "R_error_1D_pct_vs_CFD",
            "value": (r_1d - r_cfd) / r_cfd * 100.0,
            "units": "%",
        },
        {
            "metric": "R_error_1D_at_Qcfd_pct_vs_CFD",
            "value": (r_1d_qcfd - r_cfd) / r_cfd * 100.0,
            "units": "%",
        },
        {
            "metric": "R_error_poiseuille_pct_vs_CFD",
            "value": (r_pois - r_cfd) / r_cfd * 100.0,
            "units": "%",
        },
        {
            "metric": "pressure_rmse_1D_vs_CFD_pa",
            "value": rmse_1d,
            "units": "Pa",
        },
        {
            "metric": "pressure_rmse_1D_at_Qcfd_vs_CFD_pa",
            "value": rmse_1d_qcfd,
            "units": "Pa",
        },
        {
            "metric": "pressure_rmse_poiseuille_vs_CFD_pa",
            "value": rmse_pois,
            "units": "Pa",
        },
        {
            "metric": "alpha",
            "value": alpha,
            "units": "-",
        },
    ]
    write_rows(
        OUT / "steady_1d_rigid_summary.csv",
        fieldnames=["metric", "value", "units"],
        rows=summary,
    )

    print("steady_1d_rigid_profile.csv written")
    print("steady_1d_vs_cfd_points.csv written")
    print("steady_1d_rigid_summary.csv written")


if __name__ == "__main__":
    main()
