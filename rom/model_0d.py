"""
0-D Windkessel Model for the Circle of Willis
==============================================
A lumped-parameter (Windkessel) representation of the cerebral circulation.

The standard three-element Windkessel (RCR) model:

    Q_in(t) --> R_p --> +--C--+ --> R_d --> P_out
                        |     |
                        +--p--+

Governing ODE:
    C * dp/dt = Q_in(t) - p / R_d        (diastolic decay)

with the total pressure drop:
    P_in(t) = p(t) + R_p * Q_in(t)

Parameters are derived from physiological data for the MCA segment.

Usage
-----
    python rom/model_0d.py

Output
------
    data/centerlines/cow_0d_results.csv   -- time series of p(t) and Q(t)
    Plots displayed interactively
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ---------------------------------------------------------------------------
# Physiological parameters (MCA, approximate values)
# ---------------------------------------------------------------------------
RHO = 1060.0          # blood density [kg/m^3]
HEART_RATE = 75       # beats per minute
T_CYCLE = 60.0 / HEART_RATE  # cardiac cycle duration [s]

# Three-element Windkessel parameters for one MCA segment
R_PROXIMAL = 1.0e7    # proximal resistance  [Pa·s/m^3]
R_DISTAL   = 8.0e7    # distal resistance    [Pa·s/m^3]
C_COMPLIANCE = 2.0e-10  # arterial compliance [m^3/Pa]

# Mean inlet flow rate
Q_MEAN = 2.5e-6       # [m^3/s]  (~150 mL/min, both MCAs)
Q_AMPLITUDE = 1.5e-6  # pulsatile amplitude [m^3/s]


# ---------------------------------------------------------------------------
# Inlet flow waveform  Q(t)  – simplified sinusoidal pulsatile profile
# ---------------------------------------------------------------------------
def inlet_flow(t: float) -> float:
    """Return pulsatile inlet flow rate at time t [m^3/s]."""
    phase = 2 * np.pi * t / T_CYCLE
    return Q_MEAN + Q_AMPLITUDE * (
        0.8 * np.sin(phase)
        + 0.2 * np.sin(2 * phase - 0.5)
    )


# ---------------------------------------------------------------------------
# ODE integration using forward Euler
# ---------------------------------------------------------------------------
def run_windkessel(
    n_cycles: int = 5,
    dt: float = 1e-4,
) -> dict:
    """
    Integrate the three-element Windkessel ODE over *n_cycles* cardiac cycles.

    Returns a dict with keys: t, Q, p_c, p_in
    """
    t_end = n_cycles * T_CYCLE
    t = np.arange(0.0, t_end + dt, dt)

    p_c = np.zeros_like(t)    # compliance pressure
    p_in = np.zeros_like(t)   # inlet pressure
    Q = np.zeros_like(t)

    # Initial condition: steady state pressure
    p_c[0] = Q_MEAN * R_DISTAL

    for i in range(1, len(t)):
        q = inlet_flow(t[i - 1])
        Q[i - 1] = q
        # dp/dt = (Q - p_c / R_d) / C
        dp_dt = (q - p_c[i - 1] / R_DISTAL) / C_COMPLIANCE
        p_c[i] = p_c[i - 1] + dp_dt * dt
        p_in[i - 1] = p_c[i - 1] + R_PROXIMAL * q

    # Last point
    Q[-1] = inlet_flow(t[-1])
    p_in[-1] = p_c[-1] + R_PROXIMAL * Q[-1]

    return {"t": t, "Q": Q, "p_c": p_c, "p_in": p_in}


# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------
def save_results(results: dict, output_dir: str = "data/centerlines") -> None:
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "cow_0d_results.csv")

    header = "time_s,flow_rate_m3s,compliance_pressure_Pa,inlet_pressure_Pa"
    data = np.column_stack(
        [results["t"], results["Q"], results["p_c"], results["p_in"]]
    )
    np.savetxt(out_path, data, delimiter=",", header=header, comments="")
    print(f"Saved 0-D results → {out_path}")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_results(results: dict) -> None:
    t = results["t"]
    # Show only the last two cycles (settled)
    t_start = t[-1] - 2 * T_CYCLE
    mask = t >= t_start
    t_plot = t[mask] - t_start  # shift to start from 0

    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    axes[0].plot(t_plot, results["Q"][mask] * 1e6, color="tab:blue", linewidth=1.5)
    axes[0].set_ylabel("Flow rate  Q  [mL/s]")
    axes[0].set_title("0-D Windkessel Model – MCA segment")
    axes[0].grid(True, linestyle="--", alpha=0.5)

    axes[1].plot(
        t_plot, results["p_in"][mask] / 133.322, color="tab:red", label="Inlet"
    )
    axes[1].plot(
        t_plot, results["p_c"][mask] / 133.322, color="tab:orange",
        linestyle="--", label="Compliance node"
    )
    axes[1].set_ylabel("Pressure  [mmHg]")
    axes[1].set_xlabel("Time  [s]")
    axes[1].legend()
    axes[1].grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig("data/centerlines/cow_0d_pressure_flow.png", dpi=150)
    plt.show()
    print("Plot saved → data/centerlines/cow_0d_pressure_flow.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Running 0-D Windkessel model …")
    results = run_windkessel(n_cycles=6, dt=5e-5)
    save_results(results)
    plot_results(results)
    print("Done.")
