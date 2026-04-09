"""
1-D Pulse-Wave Propagation Model for the Circle of Willis
==========================================================
Implements the linearised 1-D blood-flow equations (Euler equations with
viscous damping) in a single straight vessel segment.

Governing equations (linearised, elastic wall):
    ∂A/∂t + ∂(AU)/∂x = 0                     (mass)
    ∂U/∂t + U ∂U/∂x + 1/ρ ∂p/∂x = -f*U      (momentum, friction term f)

Tube law (elastic wall):
    p = p_ext + β/A₀ * (√A - √A₀)            with β = √(π) * E * h / (1-ν²)

This script uses a simple MacCormack finite-difference scheme.

Usage
-----
    python rom/model_1d.py

Output
------
    data/centerlines/cow_1d_results.csv   -- spatial profiles at each saved time
    Plots displayed interactively
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ---------------------------------------------------------------------------
# Physical and geometric parameters (MCA, approximate)
# ---------------------------------------------------------------------------
RHO   = 1060.0       # blood density [kg/m^3]
MU    = 3.5e-3       # dynamic viscosity [Pa·s]
N_VIS = 9            # Womersley friction profile order

# Vessel segment (MCA)
L         = 0.06     # vessel length [m]
R0        = 0.002    # reference radius [m]
A0        = np.pi * R0**2  # reference cross-section [m^2]
E_WALL    = 5.0e5    # Young's modulus of vessel wall [Pa]
H_WALL    = 3.0e-4   # wall thickness [m]
NU_POISSON = 0.5     # Poisson's ratio

BETA = np.sqrt(np.pi) * E_WALL * H_WALL / ((1 - NU_POISSON**2) * A0)

# Friction coefficient (Poiseuille):  f = 8*mu / (rho * R0^2)
# Equivalently: 8*pi*mu / (rho * A0)  since A0 = pi * R0^2
FRICTION_COEFF = 8 * np.pi * MU / (RHO * A0)  # [1/s]

# Pulse-wave speed (Moens-Korteweg)
C0 = np.sqrt(BETA * np.sqrt(A0) / (2 * RHO))

# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------
HEART_RATE = 75
T_CYCLE    = 60.0 / HEART_RATE
N_CYCLES   = 4
T_END      = N_CYCLES * T_CYCLE

NX    = 100          # spatial grid points
DX    = L / (NX - 1)
CFL   = 0.8
DT    = CFL * DX / (C0 + 0.5)  # Courant-limited time step

# ---------------------------------------------------------------------------
# Tube law and wave speed
# ---------------------------------------------------------------------------

def tube_law(A: np.ndarray) -> np.ndarray:
    """Pressure from cross-sectional area via elastic tube law [Pa]."""
    return BETA * (np.sqrt(A) - np.sqrt(A0))


def wave_speed(A: np.ndarray) -> np.ndarray:
    """Local pulse-wave speed c = sqrt(A/rho * dp/dA) [m/s]."""
    return np.sqrt(BETA * np.sqrt(A) / (2 * RHO))


# ---------------------------------------------------------------------------
# Inlet flow waveform
# ---------------------------------------------------------------------------
Q_MEAN = 2.5e-6       # mean flow [m^3/s]
Q_AMP  = 1.5e-6       # pulsatile amplitude [m^3/s]


def inlet_flow(t: float) -> float:
    """Pulsatile inlet flow rate [m^3/s]."""
    phase = 2 * np.pi * t / T_CYCLE
    return Q_MEAN + Q_AMP * (0.8 * np.sin(phase) + 0.2 * np.sin(2 * phase - 0.5))


# ---------------------------------------------------------------------------
# MacCormack scheme helper
# ---------------------------------------------------------------------------

def fluxes(A: np.ndarray, U: np.ndarray):
    """Return flux vectors F1 = A*U, F2 = A*U^2 + A*p/rho."""
    p = tube_law(A)
    F1 = A * U
    F2 = A * U**2 + A * p / RHO
    return F1, F2


def source(A: np.ndarray, U: np.ndarray):
    """Return source term for momentum: -f*A*U (friction)."""
    return -FRICTION_COEFF * A * U


# ---------------------------------------------------------------------------
# Run the simulation
# ---------------------------------------------------------------------------

def run_1d() -> dict:
    """
    Integrate the 1-D blood-flow equations using the MacCormack scheme.

    Returns a dict with keys: x, t_saved, A_saved, U_saved, p_saved
    """
    x = np.linspace(0.0, L, NX)

    # Initial conditions: rest
    A = np.full(NX, A0)
    U = np.zeros(NX)

    t = 0.0
    save_every = T_CYCLE / 20   # save 20 snapshots per cycle
    next_save  = 0.0

    t_saved: list = []
    A_saved: list = []
    U_saved: list = []
    p_saved: list = []

    while t < T_END:
        # --- predictor step ---
        F1, F2 = fluxes(A, U)
        S2     = source(A, U)

        A_star = A.copy()
        U_star = U.copy()

        A_star[:-1] = A[:-1] - DT / DX * (F1[1:] - F1[:-1])
        U_star[:-1] = (
            U[:-1]
            - DT / DX * (F2[1:] - F2[:-1]) / A[:-1]
            + DT * S2[:-1] / A[:-1]
        )

        # --- corrector step ---
        F1s, F2s = fluxes(A_star, U_star)
        S2s      = source(A_star, U_star)

        A_new = 0.5 * (A + A_star)
        U_new = 0.5 * (U + U_star)

        A_new[1:] -= 0.5 * DT / DX * (F1s[1:] - F1s[:-1])
        U_new[1:] -= (
            0.5 * DT / DX * (F2s[1:] - F2s[:-1]) / A_star[1:]
            - 0.5 * DT * S2s[1:] / A_star[1:]
        )

        # --- boundary conditions ---
        # Inlet: prescribed flow rate → velocity = Q / A
        q_in    = inlet_flow(t + DT)
        A_new[0] = A0                  # fixed cross-section at inlet
        U_new[0] = q_in / A_new[0]

        # Outlet: non-reflecting (zero-gradient)
        A_new[-1] = A_new[-2]
        U_new[-1] = U_new[-2]

        A[:] = A_new
        U[:] = U_new
        t   += DT

        # --- save snapshot ---
        if t >= next_save:
            t_saved.append(t)
            A_saved.append(A.copy())
            U_saved.append(U.copy())
            p_saved.append(tube_law(A))
            next_save += save_every

    return {
        "x": x,
        "t": np.array(t_saved),
        "A": np.array(A_saved),
        "U": np.array(U_saved),
        "p": np.array(p_saved),
    }


# ---------------------------------------------------------------------------
# Save results (last cycle centre-line profile)
# ---------------------------------------------------------------------------

def save_results(results: dict, output_dir: str = "data/centerlines") -> None:
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "cow_1d_results.csv")

    x = results["x"]
    # Use the snapshot closest to peak systole of the last cycle
    U_peak_idx = np.argmax(results["U"][:, 0])  # peak at inlet
    U_prof = results["U"][U_peak_idx]
    p_prof = results["p"][U_peak_idx]
    A_prof = results["A"][U_peak_idx]
    r_prof = np.sqrt(A_prof / np.pi)

    header = "arc_length_m,radius_m,velocity_1d_ms,pressure_1d_Pa"
    data   = np.column_stack([x, r_prof, U_prof, p_prof])
    np.savetxt(out_path, data, delimiter=",", header=header, comments="")
    print(f"Saved 1-D results → {out_path}")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_results(results: dict) -> None:
    x = results["x"] * 1e3   # m → mm

    # Spatial profiles at peak systole and end diastole
    U_last = results["U"][-20:]  # last 20 snapshots
    p_last = results["p"][-20:]

    i_peak  = np.argmax(U_last[:, 0])
    i_diast = np.argmin(U_last[:, 0])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].plot(x, U_last[i_peak],  label="Peak systole",   color="tab:red")
    axes[0].plot(x, U_last[i_diast], label="End diastole",   color="tab:blue", linestyle="--")
    axes[0].set_xlabel("Arc length  [mm]")
    axes[0].set_ylabel("Cross-sectional mean velocity  [m/s]")
    axes[0].set_title("1-D Model – Velocity profile along MCA")
    axes[0].legend()
    axes[0].grid(True, linestyle="--", alpha=0.5)

    axes[1].plot(x, p_last[i_peak]  / 133.322, label="Peak systole",  color="tab:red")
    axes[1].plot(x, p_last[i_diast] / 133.322, label="End diastole",  color="tab:blue", linestyle="--")
    axes[1].set_xlabel("Arc length  [mm]")
    axes[1].set_ylabel("Pressure  [mmHg]")
    axes[1].set_title("1-D Model – Pressure profile along MCA")
    axes[1].legend()
    axes[1].grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig("data/centerlines/cow_1d_profiles.png", dpi=150)
    plt.show()
    print("Plot saved → data/centerlines/cow_1d_profiles.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print(f"Running 1-D model  (Δt = {DT*1e6:.1f} µs, NX = {NX}) …")
    results = run_1d()
    save_results(results)
    plot_results(results)
    print("Done.")
