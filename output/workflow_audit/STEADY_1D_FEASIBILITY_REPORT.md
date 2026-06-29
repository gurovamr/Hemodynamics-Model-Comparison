# Steady 1D Rigid-Wall Upgrade Assessment

Date: 2026-06-29
Repository: Hemodynamics-Model-Comparison

## 1. Inventory (Current Workflow Data Availability)

### Already available

- Ordered centerline path and segment connectivity:
  - scripts/03_reduced_order_models_0D_1D.ipynb (ordered path reconstruction and ordered cells)
  - scripts/model_comparison_0D1D_3D.ipynb (graph loading and edge endpoint extraction)
- Arc-length coordinate:
  - scripts/03_reduced_order_models_0D_1D.ipynb (s_pts_m, s_cell_m)
  - output/reduced_order_models/segment_1D_geometry_profile.csv
- Local radius profile:
  - scripts/03_reduced_order_models_0D_1D.ipynb (ce_radius, mis_radius, voreen_radius)
  - output/reduced_order_models/segment_1D_geometry_profile.csv
- Cross-sectional area:
  - scripts/03_reduced_order_models_0D_1D.ipynb (A_i = pi r^2)
  - output/reduced_order_models/segment_1D_geometry_profile.csv (area_m2)
- Flow rate (slice-integrated CFD):
  - output/06_clear_results/flow_rate_summary.csv
- Inlet/outlet pressure (slice area-averaged CFD):
  - output/06_clear_results/pressure_slice_summary.csv
- Fluid properties:
  - scripts/03_reduced_order_models_0D_1D.ipynb (MU)
  - scripts/audit_cfd_rom_workflow.py and scripts/figures_for_report.ipynb (RHO, MU)

### Missing for a fully validated hemodynamic 1D simulation

- Continuous area-averaged pressure profile along the vessel (not only 3 slices).
- Friction closure calibration beyond Poiseuille laminar rigid-wall assumptions.
- High-fidelity boundary-condition consistency (only 3 pressure points for this domain in current common-domain comparison).
- Branch/junction treatment (for network-level genuine 1D, outside this single-vessel scope).

## 2. Feasibility (Single-Vessel, Steady, Rigid)

A genuine **steady 1D rigid-wall** solver is feasible with current data for a single vessel.

Minimal governing system:

- Continuity:
  - dQ/ds = 0
- Momentum (cross-section averaged, rigid wall):
  - dp/ds = alpha * rho * Q^2 / A^3 * dA/ds - 8*pi*mu*Q/A^2

Equivalent segment form used in prototype:

- dp_conv = 0.5 * alpha * rho * (u_i^2 - u_{i+1}^2)
- dp_visc = 8 * mu * Q * ds / (pi * r^4)
- p_{i+1} = p_i + dp_conv - dp_visc

This is a real 1D momentum solve on variable area (not just R=sum 8muL/(pi r^4)), but still a simplified steady rigid model.

## 3. Minimum Implementation Scope

Implemented with reuse of existing exports (no project rewrite):

- New file:
  - scripts/steady_1d_rigid_solver_prototype.py
- New equations:
  - 2 governing equations (steady continuity + momentum)
  - 1 nonlinear scalar solve for Q from pressure BCs (bisection)
- New functions:
  - 11 helper/solver functions (CSV I/O, domain loading, pressure march, nonlinear solve, interpolation, metrics)
- Code size:
  - ~370 lines total file; core solver logic ~140 lines.

Complexity assessment:

- Low-to-moderate coding complexity.
- High scientific-risk complexity (closure/boundary consistency dominate error).

## 4. Prototype Outputs

Generated files:

- output/workflow_audit/steady_1d_rigid_profile.csv
- output/workflow_audit/steady_1d_vs_cfd_points.csv
- output/workflow_audit/steady_1d_rigid_summary.csv

Reference baseline used:

- output/reduced_order_models/model_comparison_common_domain.csv
- output/reduced_order_models/pressure_profiles_common_domain.csv

## 5. Validation Results

From output/workflow_audit/steady_1d_rigid_summary.csv:

- CFD reference:
  - Q_CFD = 1.0251e-06 m3/s
  - DeltaP_CFD = 356.37 Pa
  - R_CFD = 3.4764e+08 Pa s/m3

- Existing Poiseuille baseline (at Q_CFD):
  - DeltaP = 145.53 Pa
  - R = 1.4197e+08 Pa s/m3
  - Resistance error vs CFD = -59.16%
  - Pressure RMSE vs CFD (3 points) = 125.89 Pa

- New steady 1D rigid prototype, evaluated at Q_CFD:
  - DeltaP = 166.08 Pa
  - R = 1.6202e+08 Pa s/m3
  - Resistance error vs CFD = -53.39%
  - Pressure RMSE vs CFD (3 points) = 109.90 Pa

- New steady 1D rigid prototype, solving Q from inlet/outlet pressure:
  - Q_solved = 2.2978e-06 m3/s
  - Flow error vs CFD = +124.15%
  - Pressure RMSE vs CFD (3 points) = 210.06 Pa

Interpretation:

- The added convective-area term gives only a modest improvement in profile RMSE at fixed Q.
- Hydraulic resistance remains far from CFD (still >50% low).
- Solving Q from the current reduced model with present closures gives nonphysical mismatch.

## 6. Scientific Decision

A small in-notebook steady 1D rigid solver is technically feasible and has been prototyped.

However, with current available validation data and closure assumptions, improvement is limited and not enough to claim robust hemodynamic upgrade over the existing Poiseuille-based ROM.

Therefore, the current project should remain primarily Poiseuille-based for defensible reporting unless further development is added.

## 7. What is required for a defensible genuine 1D hemodynamic model

Minimum additional development:

1. Better pressure-profile truth data along centerline (not just 3 area-averaged slices).
2. Consistent 1D friction/energy closure calibration for curved, tapering vessel segments.
3. Controlled inlet/outlet BC strategy for 1D model (avoid over-constrained mismatch).
4. If moving beyond single-vessel steady case:
   - junction coupling,
   - outlet model closures,
   - and preferably transient validation.

External software path (if desired):

- Use a dedicated 1D solver stack (for example SimVascular svOneDSolver or a validated hemodynamics 1D code) and import this workflow's centerline/radius as geometry input.
- Keep current notebooks as geometry extraction + CFD comparison harness.
