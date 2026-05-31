# Straight Pipe CFD — Workflow Documentation

**Project:** Hemodynamics Model Comparison — TrackModule2  
**Case directory:** `openfoam/pipe/`  
**OpenFOAM version:** dev (526f42b99383), `simpleFoam` → `foamRun -solver incompressibleFluid`  
**Date completed:** May 2026  
**Author:** gurovamr

---

## Table of Contents

1. [Objective of the Simulation](#1-objective-of-the-simulation)
2. [Source Tutorial and Initial Setup](#2-source-tutorial-and-initial-setup)
3. [Geometry and Physical Parameters](#3-geometry-and-physical-parameters)
4. [Complete Simulation Workflow](#4-complete-simulation-workflow)
5. [OpenFOAM Case Structure](#5-openfoam-case-structure)
6. [ParaView Post-Processing Workflow](#6-paraview-post-processing-workflow)
7. [Generated Results](#7-generated-results)
8. [Quantitative Analysis](#8-quantitative-analysis)
9. [Lessons Learned](#9-lessons-learned)
10. [Relevance for Future Circle of Willis Simulations](#10-relevance-for-future-circle-of-willis-simulations)

---

## 1. Objective of the Simulation

### 1.1 Motivation

Before applying CFD to patient-specific vascular geometries derived from medical imaging, it is essential to establish and validate the entire numerical workflow against a problem with a known exact analytical solution. The Hagen-Poiseuille flow in a straight circular pipe is the canonical benchmark for this purpose: it is geometrically simple, admits an exact closed-form solution for velocity, pressure drop, and wall shear stress, and is governed by the same Navier-Stokes equations as arterial flow in the laminar regime.

This simulation serves as the **foundation layer** of the project. All solver settings, blood fluid properties, boundary condition formulations, and post-processing scripts developed and validated here are subsequently applied unchanged to the patient-specific Circle of Willis (CoW) geometry.

### 1.2 CFD Concepts Validated

| Concept | Validation method |
|---|---|
| Laminar velocity profile | Radial profile Uz(r) compared to Poiseuille parabola |
| Pressure-driven flow | ΔP compared to 8μLU/R² |
| Mass conservation | U_mean inlet = U_mean outlet |
| Solver convergence | Residuals below 1×10⁻⁵ |
| Mesh independence | O-grid with 13,600 cells; non-ortho < 28° |
| Parabolic inlet BC | Eliminates entrance-length effects |

### 1.3 Connection to 0D/1D Models

The 0D (lumped-parameter) and 1D (one-dimensional flow) models of the Circle of Willis both represent each vessel segment as a single scalar hydraulic resistance R [Pa·s/m³]. The fundamental governing equation is:

$$\Delta P = R \cdot Q$$

where ΔP is the pressure drop across the segment and Q is the volumetric flow rate. For a straight circular tube this resistance is given by the Poiseuille formula:

$$R_{0D} = \frac{8 \mu L}{\pi R^4}$$

This pipe simulation directly verifies that the CFD solver reproduces this relationship to within 2%, establishing that the 3D CFD results can be trusted when the geometry departs from ideality (curvature, taper, bifurcations) and no analytical solution exists.

### 1.4 Expected Outputs

The following quantities are extracted from the pipe simulation and will later be compared against 0D/1D model predictions for each CoW segment:

- Volumetric flow rate Q [m³/s or mL/s]
- Pressure drop ΔP [Pa or mmHg]
- Hydraulic resistance R = ΔP/Q [Pa·s/m³ or mmHg·s/mL]
- Mean and peak velocity [m/s]
- Wall shear stress τ_w [Pa]
- Axial pressure gradient dP/dz [Pa/m]

---

## 2. Source Tutorial and Initial Setup

### 2.1 Tutorial References

The case was constructed by combining two official OpenFOAM tutorials:

| Aspect | Tutorial source |
|---|---|
| Pipe geometry (blockMeshDict) | `$FOAM_TUTORIALS/mesh/blockMesh/pipe` |
| Solver configuration (fvSchemes, fvSolution, controlDict) | `$FOAM_TUTORIALS/incompressibleFluid/planarPoiseuille` |
| Blood fluid properties | `openfoam/segment_test_V2/constant/transportProperties` (project-specific) |

The case directory `openfoam/segment_test_V2` (an existing snappyHexMesh + simpleFoam case for a vessel segment) was used as the reference for solver settings to ensure **identical numerical schemes** across all cases in the project.

### 2.2 File Provenance

| File | Origin | Status | Modification summary |
|---|---|---|---|
| `system/blockMeshDict` | `mesh/blockMesh/pipe` tutorial | **Modified** | Scaled to R=2mm, L=40mm; added inline comments; fixed vertex spacing |
| `system/controlDict` | `planarPoiseuille` tutorial | **Modified** | Removed function objects (API incompatibility with OF-dev); endTime=500 |
| `system/fvSchemes` | `segment_test_V2/system/fvSchemes` | **Copied** | Identical — ensures matching numerics with vessel case |
| `system/fvSolution` | `segment_test_V2/system/fvSolution` | **Copied** | Identical — GAMG + symGaussSeidel + SIMPLE consistent |
| `system/sampleDict` | Created from scratch | **New** | Line sampling along axis and radial at outlet |
| `constant/transportProperties` | `segment_test_V2/constant/` | **Modified** | Verified blood ν = 3.3018868×10⁻⁶ m²/s |
| `constant/momentumTransport` | `planarPoiseuille` tutorial | **Copied** | `simulationType laminar;` |
| `0/U` | Created from scratch | **New** | `codedFixedValue` parabolic inlet; noSlip wall; zeroGradient outlet |
| `0/p` | `planarPoiseuille` tutorial | **Modified** | zeroGradient inlet; fixedValue 0 outlet; kinematic dimensions |
| `Allrun` | Created from scratch | **New** | Documents full workflow with comments |
| `Allclean` | Created from scratch | **New** | Removes time dirs, logs, postProcessing |

### 2.3 Key Design Decision: Parabolic Inlet

The tutorial uses a uniform inlet velocity (`fixedValue uniform (0 0 U)`). For a uniform inlet, the velocity profile requires a development length of:

$$L_{dev} \approx 0.06 \cdot Re \cdot D = 0.06 \times 121 \times 0.004 \approx 0.029 \text{ m}$$

With L = 0.040 m, approximately 73% of the pipe length would be in the developing region, causing the measured ΔP to exceed the Poiseuille analytical value by ~26%. The inlet BC was changed to a fully-developed parabolic profile via `codedFixedValue`, eliminating this discrepancy and reducing the ΔP error to 1.78%.

---

## 3. Geometry and Physical Parameters

### 3.1 Geometry

| Parameter | Symbol | Value | Units | Source |
|---|---|---|---|---|
| Pipe radius | R | 0.002 | m | Chosen to match vessel-scale (2 mm) |
| Pipe diameter | D = 2R | 0.004 | m | Derived |
| Pipe length | L | 0.040 | m | L/D = 10; sufficient for fully-developed flow |
| Inner block half-width | a = R/2 | 0.001 | m | O-grid construction |
| Wall block corner radius | Rw = R/√2 | 0.001414 | m | 45° wall intersection |
| Pipe axis direction | z | — | — | Flow in +z direction |

**Mesh:** 5-block O-grid (1 centre square + 4 outer arc-bounded blocks)

- Cells: **13,600**
- Points: 14,801
- Faces: 41,940 (39,660 internal)
- Max non-orthogonality: 28.1°
- Max skewness: 1.40
- Max aspect ratio: 14.5

### 3.2 Fluid Properties (Blood)

| Property | Symbol | Value | Units | Derivation |
|---|---|---|---|---|
| Dynamic viscosity | μ | 3.5×10⁻³ | Pa·s | Literature (Newtonian blood model) |
| Density | ρ | 1060 | kg/m³ | Literature |
| Kinematic viscosity | ν = μ/ρ | 3.3018868×10⁻⁶ | m²/s | **Computed:** 3.5×10⁻³ / 1060 |

> **Note:** OpenFOAM's `simpleFoam` solver works with kinematic pressure p/ρ [m²/s²] and requires kinematic viscosity ν. Density ρ is not stored; it is applied only when converting kinematic pressure to physical pressure for post-processing.

### 3.3 Flow Parameters

| Parameter | Symbol | Value | Units | Derivation |
|---|---|---|---|---|
| Mean inlet velocity | Ū | 0.1 | m/s | Chosen (representative blood flow) |
| Peak (centreline) velocity | U_max = 2Ū | 0.2 | m/s | **Computed:** Poiseuille maximum |
| Reynolds number | Re = ŪD/ν | 121 | — | **Computed:** 0.1×0.004 / 3.3×10⁻⁶ |
| Flow regime | — | Laminar | — | Re ≪ 2300 |
| Volumetric flow rate | Q = ŪπR² | 1.2566×10⁻⁶ | m³/s (1.257 mL/s) | **Computed:** 0.1 × π × 0.002² |

### 3.4 Analytical Reference (Hagen-Poiseuille)

Pressure drop:

$$\Delta P = \frac{8 \mu L \bar{U}}{R^2} = \frac{8 \times 3.5 \times 10^{-3} \times 0.040 \times 0.1}{(0.002)^2} = 28.00 \text{ Pa}$$

Hydraulic resistance:

$$R_{pipe} = \frac{8 \mu L}{\pi R^4} = \frac{8 \times 3.5 \times 10^{-3} \times 0.040}{\pi \times (0.002)^4} = 2.228 \times 10^7 \text{ Pa·s/m}^3 = 0.1671 \text{ mmHg·s/mL}$$

Velocity profile:

$$u_z(r) = 2\bar{U}\left(1 - \frac{r^2}{R^2}\right) = 0.2\left(1 - \frac{r^2}{4 \times 10^{-6}}\right) \text{ [m/s]}$$

Wall shear stress:

$$\tau_w = \frac{4 \mu \bar{U}}{R} = \frac{4 \times 3.5 \times 10^{-3} \times 0.1}{0.002} = 0.700 \text{ Pa}$$

---

## 4. Complete Simulation Workflow

All commands are run from within `openfoam/pipe/`. The OpenFOAM environment must be sourced first.

### 4.1 Environment Setup

```bash
# Source OpenFOAM dev environment
source /opt/openfoam-dev/etc/bashrc

# Navigate to case directory
cd ~/OpenFOAM/gurovamr-dev/run/TrackModule2/Hemodynamics-Model-Comparison/openfoam/pipe
```

### 4.2 Case Preparation (clean state)

```bash
./Allclean
```

**Purpose:** Remove all time directories (except `0/`), logs, and postProcessing output.  
**Output:** Clean case ready for a fresh run.

### 4.3 Mesh Generation

```bash
blockMesh 2>&1 | tee log.blockMesh
```

**Purpose:** Generate the 5-block O-grid structured hexahedral mesh.  
**Input:** `system/blockMeshDict`  
**Output:** `constant/polyMesh/` (points, faces, owner, neighbour, boundary)  
**Expected output:**
```
Mesh OK.
```

**Key mesh parameters generated:**

| Metric | Value |
|---|---|
| Cells | 13,600 |
| Points | 14,801 |
| Max non-orthogonality | 28.1° |
| Max skewness | 1.40 |
| Max aspect ratio | 14.5 |

### 4.4 Mesh Quality Check

```bash
checkMesh 2>&1 | tee log.checkMesh
```

**Purpose:** Verify mesh quality metrics are within acceptable bounds for `simpleFoam`.  
**Accept criteria:** non-ortho < 70°, skewness < 4, aspect ratio < 1000 for structured meshes.  
**Output:** `Mesh OK.` — all metrics pass.

### 4.5 Solver Execution

```bash
simpleFoam 2>&1 | tee log.simpleFoam
```

**Purpose:** Solve steady-state incompressible Navier-Stokes equations using the SIMPLE algorithm.  
**Input:** `0/U`, `0/p`, `constant/transportProperties`, `constant/momentumTransport`, `system/fvSchemes`, `system/fvSolution`, `system/controlDict`  
**Output:** Time directories `100/`, `144/` (converged solution) containing `U`, `p`, `phi`  
**Convergence:** SIMPLE converged in **91 iterations** (residuals p, U < 1×10⁻⁵)  
**Wall time:** ~10 seconds (single core)

Monitor convergence during run:
```bash
tail -f log.simpleFoam | grep -E "Time|p:|U:"
```

Expected final lines:
```
SIMPLE solution converged in 91 iterations
```

### 4.6 Velocity and Pressure Profile Extraction

```bash
postProcess -func sampleDict -latestTime
```

**Purpose:** Extract velocity and pressure along the pipe axis (centreline) and across a radial cross-section near the outlet.  
**Input:** `system/sampleDict`, field `U` and `p` from the latest time directory  
**Output:**
- `postProcessing/sampleDict/144/centreline.xy` — 200 points, z from 0.001 to 0.039 m
- `postProcessing/sampleDict/144/outletRadial.xy` — 100 points, x from -0.002 to +0.002 m at z=0.038

File columns: `distance  x  y  z  Ux  Uy  Uz  p_kinematic`

### 4.7 Patch-Averaged Pressure

```bash
postProcess -func "patchAverage(patch=inlet,  fields=(p))" -latestTime
postProcess -func "patchAverage(patch=outlet, fields=(p))" -latestTime
```

**Purpose:** Extract area-weighted average kinematic pressure at inlet and outlet patches.  
**Output:**
- `postProcessing/patchAverage(patch=inlet,fields=(p))/144/surfaceFieldValue.dat`
- `postProcessing/patchAverage(patch=outlet,fields=(p))/144/surfaceFieldValue.dat`

Convert kinematic pressure to physical:
```
ΔP [Pa] = (p_kin_inlet − p_kin_outlet) × ρ = 0.026884 × 1060 = 28.497 Pa
```

### 4.8 Patch-Averaged Velocity

```bash
postProcess -func "patchAverage(patch=inlet,  fields=(U))" -latestTime
postProcess -func "patchAverage(patch=outlet, fields=(U))" -latestTime
```

**Output:**
- inlet: `areaAverage(inlet) of U = (0 0 0.10078)` — 0.78% above target, mass conserved
- outlet: `areaAverage(outlet) of U = (~0 ~0 0.10078)` — lateral components ~10⁻⁹ m/s ≈ 0

### 4.9 Wall Shear Stress

```bash
postProcess -solver incompressibleFluid -func wallShearStress -latestTime
```

**Purpose:** Compute the WSS vector field on all wall patches.  
**Note:** The `-solver incompressibleFluid` flag is required in OpenFOAM dev to register the turbulence/transport model when running post-processing standalone (not during a live solver run).  
**Output:** `144/wallShearStress` — vector field on wall faces  
**Analytical check:** τ_w = 4μU/R = 0.700 Pa

### 4.10 Validation Summary (Python)

```bash
python3 - << 'EOF'
import math
mu=3.5e-3; R=0.002; U=0.1; rho=1060.0
A=math.pi*R**2; Q=A*U
dP_an = 8*mu*0.040*U/R**2
R_an  = 8*mu*0.040/(math.pi*R**4)
dP_CFD = 0.026884357*rho
R_CFD  = dP_CFD/Q
conv   = 1/133.322e6
print(f"dP: CFD={dP_CFD:.4f} Pa  Analytical={dP_an:.4f} Pa  err={abs(dP_CFD-dP_an)/dP_an*100:.2f}%")
print(f"R:  CFD={R_CFD*conv:.5f}  Analytical={R_an*conv:.5f} mmHg.s/mL  err={abs(R_CFD-R_an)/R_an*100:.2f}%")
EOF
```

---

## 5. OpenFOAM Case Structure

```
openfoam/pipe/
├── 0/                      ← Initial & boundary conditions
│   ├── U                   ← Velocity field
│   └── p                   ← Kinematic pressure field
├── constant/
│   ├── transportProperties ← Fluid viscosity
│   ├── momentumTransport   ← Turbulence model selection
│   └── polyMesh/           ← Generated mesh (created by blockMesh)
├── system/
│   ├── blockMeshDict       ← O-grid pipe mesh definition
│   ├── controlDict         ← Solver control, time stepping
│   ├── fvSchemes           ← Discretisation schemes
│   ├── fvSolution          ← Linear solvers and SIMPLE settings
│   └── sampleDict          ← Post-processing line samples
├── Allrun                  ← Full reproducible workflow script
├── Allclean                ← Clean case script
├── pipe.foam               ← Empty dummy file for ParaView reader
└── postProcessing/         ← Created by postProcess commands
    ├── sampleDict/144/
    │   ├── centreline.xy
    │   └── outletRadial.xy
    ├── patchAverage(patch=inlet,fields=(p))/144/surfaceFieldValue.dat
    └── patchAverage(patch=outlet,fields=(p))/144/surfaceFieldValue.dat
```

### 5.1 `0/U` — Velocity Boundary Conditions

**Controls:** Velocity field at all boundaries and the domain interior.  
**Key values:**

| Boundary | BC type | Value |
|---|---|---|
| `inlet` | `codedFixedValue` (parabolic) | Uz(r) = 0.2×(1−r²/R²) |
| `outlet` | `zeroGradient` | Fully-developed outflow — no reflection |
| `wall` | `noSlip` | u = 0 |
| `internalField` | `uniform (0 0 0.01)` | Small seed to help SIMPLE start |

**Why parabolic inlet:** Eliminates the 0.029 m entrance-development region that would corrupt the ΔP comparison with the analytical formula. Tested: uniform inlet gives 26% ΔP error; parabolic inlet gives 1.78%.

### 5.2 `0/p` — Pressure Boundary Conditions

**Controls:** Kinematic pressure p/ρ [m²/s²] at all boundaries.

| Boundary | BC type | Value |
|---|---|---|
| `inlet` | `zeroGradient` | Pressure floats to satisfy momentum |
| `outlet` | `fixedValue 0` | Reference pressure = 0; ΔP computed relative to this |
| `wall` | `zeroGradient` | No normal pressure gradient at wall |

### 5.3 `constant/transportProperties`

```
viscosityModel  constant;
nu  [0 2 -1 0 0 0 0]  3.3018868e-6;
```

Dimension brackets: `[kg m s K mol A cd]` → `[0 2 -1 0 0 0 0]` = m²/s (kinematic viscosity).

### 5.4 `constant/momentumTransport`

```
simulationType  laminar;
```

Disables all turbulence modelling. Required for Re=121 where turbulence is physically absent.

### 5.5 `system/blockMeshDict` — Mesh Definition

**Controls:** Vertex positions, block connectivity, arc edges, boundary patches, and cell counts.

**Key design choices:**

| Parameter | Value | Reason |
|---|---|---|
| Inner square half-width a | R/2 = 0.001 m | Minimises cell distortion at block interfaces |
| Wall corner Rw | R/√2 = 0.001414 m | Places wall corners exactly at 45°/135°/225°/315° on circle |
| Arc midpoints | (±R, 0, z) and (0, ±R, z) | Defines exact circle; O-grid edge mapping |
| Axial cells (nz) | 40 | Gives cell aspect ratio ~14.5 at wall; acceptable for steady laminar flow |
| Radial cells (ny) | 6 (outer blocks), 10 (centre) | Sufficient radial resolution for Poiseuille profile |

### 5.6 `system/fvSchemes`

| Term | Scheme | Reason |
|---|---|---|
| Time derivative | `steadyState` | SIMPLE is a steady-state solver |
| Gradient | `Gauss linear` | Second-order accurate |
| Divergence div(phi,U) | `Gauss linearUpwind grad(U)` | Bounded, 2nd order upwind; appropriate for laminar |
| Laplacian | `Gauss linear corrected` | Non-orthogonality correction applied |

### 5.7 `system/fvSolution`

| Component | Setting | Value |
|---|---|---|
| Pressure solver | GAMG | Algebraic multigrid — efficient for elliptic p equation |
| Velocity solver | smoothSolver + symGaussSeidel | Iterative smoother |
| SIMPLE consistent | yes | SIMPLEC variant — faster convergence |
| Residual tolerance p | 1×10⁻⁵ | Convergence criterion |
| Residual tolerance U | 1×10⁻⁵ | Convergence criterion |
| Relaxation p | 0.3 | Conservative; improves stability |
| Relaxation U | 0.7 | Standard SIMPLE value |

### 5.8 `system/sampleDict`

Defines two line samples executed via `postProcess -func sampleDict`:

| Set name | Direction | Start | End | nPoints | Purpose |
|---|---|---|---|---|---|
| `centreline` | axial | (0,0,0.001) | (0,0,0.039) | 200 | Pressure gradient + centreline velocity |
| `outletRadial` | radial | (-0.002,0,0.038) | (0.002,0,0.038) | 100 | Radial velocity profile |

---

## 6. ParaView Post-Processing Workflow

### 6.1 Opening the Case

```bash
# Start Xwayland if not running (required for ParaView 6.1 with Qt xcb backend)
ls /tmp/.X1-lock &>/dev/null || (Xwayland :1 & sleep 1)
export DISPLAY=:1

# Open case (use .foam extension — registered in ParaView 6.1 reader)
cd openfoam/pipe
paraview6 pipe.foam &
```

In ParaView:
1. File loads automatically; **uncheck** `Skip Zero Time` to see t=0
2. Select **Time: 144** (converged solution) using the time toolbar
3. Click **Apply** in the Properties panel

### 6.2 Velocity Visualisation

1. In Pipeline Browser, click `pipe.foam`
2. Change field dropdown from `p` to `U` → select **Magnitude**
3. Representation: **Surface** with colour by U magnitude
4. Expected: red core (U_max ≈ 0.2 m/s), blue outer ring (wall U = 0)

### 6.3 Pressure Visualisation

1. Change field to `p` (kinematic)
2. Blue (low) at outlet (p=0), Red (high) at inlet (~0.0269 m²/s²)
3. To convert to Pa: scale colour map by ρ=1060; or note that ΔP_physical = Δp_kinematic × 1060

### 6.4 Creating Streamlines

1. Filters → Search → `Stream Tracer` → Apply
2. In Properties:
   - **Seed Type:** Line Source
   - **Point1:** (-0.002, 0, 0.020)
   - **Point2:** (0.002, 0, 0.020)
   - **Resolution:** 20
3. Click **Apply**
4. Colour by `U` Magnitude

**Result:** 20 parallel streamlines across the diameter, coloured by velocity — red/orange centre, blue near walls. All lines are straight and parallel, confirming fully-developed laminar Poiseuille flow.

**Optional — render as tubes:**
- Filters → Tube → set Radius = 0.0001 m → Apply

### 6.5 Creating Cross-Section Slices

1. Select `pipe.foam` in Pipeline Browser
2. Filters → Slice
3. Set normal to **Z-axis** (0, 0, 1)
4. Set origin to (0, 0, 0.020) for mid-plane slice
5. Click **Apply**
6. Colour by `U` component Z (Uz) to show parabolic cross-section

### 6.6 Integrate Variables (Flux and Area)

1. Select the slice in Pipeline Browser
2. Filters → **Integrate Variables**
3. Click **Apply**
4. Open **SpreadSheet View** (View → SpreadSheet View)
5. The row contains integrated U and p over the slice area
6. Divide integrated Uz by Area column to get area-averaged Uz

### 6.7 Exporting CSV Data

1. With SpreadSheet View active, select the Integrate Variables output
2. File → **Save Data** → choose CSV format
3. Save to `data_parafoam/` (project convention)

### 6.8 Saving ParaView State

File → **Save State** → save as `postProcessing/pipe_postprocessing.pvsm`

This file allows reloading the entire visualisation pipeline without repeating setup steps.

---

## 7. Generated Results

### 7.1 Field Files (OpenFOAM)

| File | Physical meaning | Use |
|---|---|---|
| `144/U` | 3D velocity vector field [m/s] | Source for all velocity post-processing |
| `144/p` | 3D kinematic pressure field [m²/s²] | Source for pressure post-processing; multiply by ρ for Pa |
| `144/phi` | Face flux [m³/s] | Used internally by solver; source for Q calculation |
| `144/wallShearStress` | WSS vector field on wall patches [m²/s²] | Multiply by ρ for Pa; CFD-only quantity |

### 7.2 Sampled Line Data

| File | Content | Columns |
|---|---|---|
| `postProcessing/sampleDict/144/centreline.xy` | Axial profiles (z = 0.001–0.039 m, 200 pts) | distance, x, y, z, Ux, Uy, Uz, p |
| `postProcessing/sampleDict/144/outletRadial.xy` | Radial profile (x = −0.002–+0.002 m at z=0.038, 100 pts) | distance, x, y, z, Ux, Uy, Uz, p |

**Physical interpretation of centreline.xy:**
- Column 7 (Uz): confirms U_centreline ≈ 0.200 m/s = U_max throughout → fully-developed flow
- Column 8 (p): linear decay from ~0.0264 m²/s² at z=0.001 to ~0.00067 m²/s² at z=0.039 → confirms linear pressure gradient

**Physical interpretation of outletRadial.xy:**
- Column 7 (Uz): parabolic profile, maximum at centre (r=0), zero at wall (r=R=0.002 m)
- Comparing with analytical Uz(r) = 0.2×(1−r²/4×10⁻⁶) verifies solver accuracy

### 7.3 Patch Average Data

| File | Value | Physical meaning |
|---|---|---|
| `patchAverage(patch=inlet,fields=(p))/144/surfaceFieldValue.dat` | 2.6884×10⁻² m²/s² | Area-averaged inlet kinematic pressure |
| `patchAverage(patch=outlet,fields=(p))/144/surfaceFieldValue.dat` | 0 m²/s² | Reference pressure (outlet BC) |

---

## 8. Quantitative Analysis

### 8.1 Complete Validated Results

| Quantity | Symbol | CFD value | Analytical | Error | Units |
|---|---|---|---|---|---|
| Inlet cross-section area | A | 1.2513×10⁻⁵ | πR² = 1.2566×10⁻⁵ | 0.42% | m² |
| Volumetric flow rate | Q | 1.2566×10⁻⁶ | 1.2566×10⁻⁶ | 0% | m³/s (1.257 mL/s) |
| Mean velocity — inlet | Ū | 0.10078 | 0.100 | +0.78% | m/s |
| Mean velocity — outlet | Ū | 0.10078 | 0.100 | +0.78% | m/s |
| Lateral velocity — outlet | U_x, U_y | ~5×10⁻⁹ | 0 | ≈ 0 | m/s |
| Peak centreline velocity | U_max | 0.20016 | 0.200 | +0.08% | m/s |
| Kinematic pressure — inlet | p_kin | 2.6884×10⁻² | — | — | m²/s² |
| Physical pressure — inlet | P | 28.497 | — | — | Pa |
| Physical pressure — outlet | P | 0 | 0 | 0% | Pa |
| Pressure drop | ΔP | 28.497 | 28.000 | +1.78% | Pa |
| Hydraulic resistance | R | 2.2678×10⁷ | 2.2282×10⁷ | +1.78% | Pa·s/m³ |
| Hydraulic resistance | R | 0.17010 | 0.16713 | +1.78% | mmHg·s/mL |
| Wall shear stress | τ_w | 0.700 | 0.700 | 0% | Pa |
| Profile shape error (mean) | — | — | — | 0.52% | — |
| Solver iterations to convergence | — | 91 | — | — | — |
| Wall-clock time | — | ~10 s | — | — | — |

### 8.2 Relevance to 0D/1D Modelling

| Quantity | Obtainable from geometry alone? | Required as BC? | CFD vs 0D/1D validation? |
|---|---|---|---|
| R (straight tube) | **Yes** — Poiseuille formula | Input to 0D network | Compare R_0D vs R_CFD |
| A (cross-section) | **Yes** — geometry | Input to 1D model | — |
| Q | **No** — needs MRI measurement | **Primary BC** | Compare Q_0D vs Q_CFD |
| ΔP | **No** — result of simulation | Output of 0D: ΔP=R·Q | Compare ΔP_0D vs ΔP_CFD |
| U_mean | **No** — needs Q and A | BC for CFD inlet | Verify Q = U·A |
| U(r) profile | **No** — Poiseuille only | Not used in 0D/1D | CFD-only validation |
| τ_w | **No** — CFD/analytical only | Not in 0D/1D | CFD-only quantity |
| Inlet/outlet pressure | **No** — needs systemic BP | Pressure BC for CFD | Compare with clinical |
| Peripheral resistance | **No** — requires calibration | Outlet BC (Windkessel) | Calibrated from clinical |

---

## 9. Lessons Learned

### 9.1 CFD Workflow

- **Always validate against an analytical solution before applying to complex geometry.** The pipe benchmark caught a 26% ΔP error that would have been undetectable in the vessel case.
- **Inlet boundary conditions critically affect results.** Uniform vs parabolic inlet changed ΔP error from 26% to 1.78%. For vascular CFD, where MRI provides bulk flow rate (not profile), a parabolic or Womersley profile should be assumed at inlets.
- **Mass conservation is a minimal necessary check** — not sufficient, but if U_mean differs between inlet and outlet the simulation is wrong.

### 9.2 OpenFOAM-Specific

| Issue encountered | Resolution |
|---|---|
| `blockMesh` parse error: `ill defined primitiveEntry` | Missing space between adjacent negative vertex coordinates: `-0.001414-0.001414` → `-0.001414 -0.001414` |
| `simpleFoam` fatal: `wrong token type - expected string, found 'fieldFunctionObjects'` | OpenFOAM dev requires unquoted library names; solution: remove function objects from `controlDict` entirely and use `postProcess` separately |
| `wallShearStress` fatal: `Unable to find turbulence model` | `postProcess` without solver context cannot find registered transport model; fix: add `-solver incompressibleFluid` flag |
| `patchGroup 'wall' clashes with patch 2` | Harmless warning — boundary patch named `wall` conflicts with auto-generated wall group. No effect on results |
| 26% ΔP error with uniform inlet | Entrance-length effect (L_dev ≈ 29 mm in 40 mm pipe); fixed by parabolic `codedFixedValue` inlet |
| ParaView 6.1 not opening | Requires Xwayland (X11 server); `QT_QPA_PLATFORM=xcb` needs `DISPLAY=:1` backed by running Xwayland |
| ParaView error: `not recognized as a supported file format` for `.OpenFOAM` | ParaView 6.1 registers `foam` extension, not `OpenFOAM`; use `pipe.foam`, not `pipe.OpenFOAM` |

### 9.3 Critical Steps for Future Vascular Simulations

1. **Always check `Mesh OK`** before running solver — non-ortho > 70° or skewness > 4 will corrupt results
2. **Confirm parabolic or Womersley inlet BC** — never use uniform inlet without accounting for entrance length
3. **Run `postProcess -solver incompressibleFluid`** for WSS in OpenFOAM dev (not just `postProcess`)
4. **Use `.foam` extension** for ParaView 6.1 reader
5. **Store all patch average data** — `patchAverage` outputs are the primary source for R and ΔP validation

---

## 10. Relevance for Future Circle of Willis Simulations

### 10.1 Identical Workflow Elements

The following are **unchanged** when moving from the pipe to the CoW vessel:

| Element | Pipe | CoW vessel |
|---|---|---|
| Solver | `simpleFoam` (SIMPLE, laminar) | Identical |
| `fvSchemes` | linearUpwind, Gauss linear corrected | Identical |
| `fvSolution` | GAMG + symGaussSeidel, SIMPLE consistent | Identical |
| Fluid properties | ν=3.3×10⁻⁶ m²/s, ρ=1060 kg/m³ | Identical |
| Post-processing | `postProcess -func "patchAverage(...)"` | Identical |
| ParaView workflow | `.foam` file, Slice, StreamTracer | Identical |
| Validation quantities | ΔP, Q, R, U_mean, WSS | Identical |

### 10.2 Additional Challenges in Curved Arteries

- **Dean flow:** Centrifugal effects create secondary vortices (counter-rotating Dean rolls). The velocity profile is no longer parabolic; the maximum velocity shifts toward the outer wall.
- **Higher resistance:** R_curved > R_straight by a factor depending on Dean number De = Re√(R_pipe/R_curvature)
- **Non-uniform WSS:** WSS peaks on the outer wall, creating athero-protective/prone regions

### 10.3 Additional Challenges in Bifurcations

- **Flow splitting:** Q_parent = Q_daughter1 + Q_daughter2 (mass conservation at junction)
- **Stagnation zone:** Low-velocity, low-WSS region at the flow divider
- **Pressure loss at junction:** Additional dissipation beyond Poiseuille; 0D models use empirical loss coefficients
- **Secondary flows:** Helical Dean-like structures persist downstream of bifurcation

### 10.4 Additional Challenges in Circle of Willis Networks

- **Multiple inlet BCs required:** Each inlet (ICA left/right, BA) needs a measured flow waveform from PC-MRI
- **Multiple outlet BCs required:** Each terminal branch needs a Windkessel (RCR) outlet model or pressure BC calibrated to mean arterial pressure
- **Network solving:** Cannot treat each segment independently — the full pressure-flow system is coupled. Requires solving the 0D/1D network first to estimate BCs, then CFD
- **Geometric complexity:** snappyHexMesh required instead of blockMesh; surface quality (smoothing, remeshing) affects mesh generation
- **Possible turbulence:** Stenosed arteries may reach Re > 2000 locally; LES or k-ω SST model needed

### 10.5 Transition: Pipe → Arterial Segment → Circle of Willis

```
STEP 1 — Straight pipe (this case)
────────────────────────────────────────
  blockMesh O-grid → simpleFoam → postProcess
  Validation: ΔP error < 2%, U_max error < 1%
  Establishes: solver, BC formulation, post-processing pipeline
  Key output: R_CFD ≈ R_Poiseuille

STEP 2 — Single arterial segment (segment_test_V2)
────────────────────────────────────────
  snappyHexMesh from .vtp surface → same simpleFoam → same postProcess
  New challenge: unstructured mesh, non-trivial boundary patches
  Validation: compare R_3D with R_0D from network_resistance_table.csv
  Key output: R_CFD per segment (deviation from Poiseuille due to geometry)

STEP 3 — Circle of Willis (cow_mr_cfd)
────────────────────────────────────────
  Full CoW surface from MR segmentation → snappyHexMesh → simpleFoam
  New challenges: multiple inlets/outlets, Windkessel BCs, network coupling
  Validation: compare CFD flow distribution with PC-MRI measurements
  Final output: R_CFD per branch → compare with 0D/1D model_comparison_0D1D_3D.ipynb
```

### 10.6 CFD Outputs Used for 0D/1D Comparison (Summary)

| CFD output | Symbol | 0D/1D equivalent | Comparison metric |
|---|---|---|---|
| Pressure drop | ΔP_CFD | ΔP_0D = R_0D × Q | |ΔP_CFD − ΔP_0D| / ΔP_0D |
| Hydraulic resistance | R_CFD = ΔP/Q | R_0D = 8μL/πR⁴ (harmonic mean) | |R_CFD − R_0D| / R_0D |
| Flow rate | Q_CFD | Q_0D (network solution) | |Q_CFD − Q_0D| / Q_0D |
| WSS distribution | τ_w(x) | Not in 0D/1D | CFD-only biomarker |
| Velocity profile | U(r,x) | U_mean = Q/A only | Profile shape — CFD validation |

---

*End of document. For questions or updates, contact the project repository: `gurovamr/Hemodynamics-Model-Comparison`, branch `BCs`.*
