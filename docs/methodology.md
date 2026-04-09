# Methodology

## Overview

This project compares three levels of fidelity for modelling blood flow in the
**Circle of Willis** (CoW):

| Level | Name | Equations | Spatial dims | Time |
|-------|------|-----------|--------------|------|
| 3-D | Full Navier-Stokes CFD | Navier-Stokes | 3 | steady / transient |
| 1-D | Pulse-wave propagation | Euler + tube law | 1 (axial) | transient |
| 0-D | Windkessel | ODE | 0 (lumped) | transient |

---

## 3-D CFD (OpenFOAM)

### Governing equations

Incompressible Navier-Stokes:

```
∇ · u = 0
∂u/∂t + (u · ∇)u = -1/ρ ∇p + ν ∇²u
```

### Assumptions

- Blood treated as a Newtonian fluid (valid for large cerebral arteries, D > 1 mm)
- Rigid vessel walls (compliant-wall extension possible via fluid-structure interaction)
- Laminar flow (Re ~ 200–500 in MCA; turbulence model optional)

### Solver

`simpleFoam` (steady-state) or `pimpleFoam` (transient, pulsatile inlet)

### Physical properties

| Property | Symbol | Value |
|----------|--------|-------|
| Density | ρ | 1060 kg/m³ |
| Dynamic viscosity | μ | 3.5 × 10⁻³ Pa·s |
| Kinematic viscosity | ν = μ/ρ | 3.3 × 10⁻⁶ m²/s |

---

## 1-D Model (Pulse-wave propagation)

### Governing equations

The 1-D Euler equations for an axisymmetric compliant vessel:

```
∂A/∂t + ∂(AU)/∂x = 0                     (continuity)
∂U/∂t + U ∂U/∂x + 1/ρ ∂p/∂x = -f U     (momentum)
```

Tube law (elastic wall):

```
p = β/A₀ · (√A − √A₀)
β = √π · E · h / (1 − ν²)
```

Pulse-wave speed (Moens-Korteweg):

```
c = √(A/ρ · dp/dA) = √(β √A / (2ρ))
```

### Numerical method

MacCormack predictor-corrector finite-difference scheme (second-order in space
and time for smooth solutions).

### Boundary conditions

- **Inlet**: prescribed pulsatile flow rate Q(t)
- **Outlet**: non-reflecting (zero-gradient A and U)

---

## 0-D Model (Three-element Windkessel)

### Circuit analogy

```
Q_in(t) → R_p → ┬─ C ─┬ → R_d → P_out = 0
                 │      │
                 └──────┘
                   p_c
```

### Governing ODE

```
C dp_c/dt = Q_in(t) − p_c / R_d
P_in(t)   = p_c(t) + R_p · Q_in(t)
```

### Parameters (MCA, approximate)

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Proximal resistance | R_p | 1×10⁷ Pa·s/m³ |
| Distal resistance | R_d | 8×10⁷ Pa·s/m³ |
| Compliance | C | 2×10⁻¹⁰ m³/Pa |

---

## Centre-line comparison

Results are compared along vessel centre-lines extracted from the 3-D mesh.
The OpenFOAM `sets` function object samples U and p along a user-defined line.
The 1-D model provides cross-sectional mean velocity and pressure along the same
arc-length.  The 0-D model provides a single pressure value (mean or
pulse-averaged) representing the vessel segment.

### Error metrics

- **RMSE**: root-mean-square error
- **Max error**: maximum absolute pointwise error
- **R²**: coefficient of determination (Pearson)

---

## Limitations and future work

- Non-Newtonian rheology (Carreau-Yasuda model) for small vessels
- Compliant walls (FSI coupling with OpenFOAM solid mechanics)
- Patient-specific geometry from medical imaging (MRI / CTA)
- Bifurcation handling in the 1-D network model
- Uncertainty quantification (Monte Carlo over vascular geometry)
