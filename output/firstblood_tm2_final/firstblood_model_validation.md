# Final TM2 FirstBlood model validation

## Acceptance decision

**Accepted as numerically converged.** The 31-to-41-point change is
0.283% for pressure drop and
0.283% for resistance, both below
1%. Terminal flow mismatch is 0.00139%
and axial flow variation is 0.00240%,
both below 0.2%.

## Consistency checks

- Geometry domain: global centreline s=19.6716–35.8927 mm; length 16.2211 mm.
- Geometry representation: 40 continuous tapered edges interpolated from the
  CFD-derived radius profile.
- Inlet: prescribed CFD inlet flow 1.021374 mL/s.
- Outlet: zero-gauge pressure, matching the CFD pressure reference.
- Fluid: rho=1060 kg/m3 and nu=3.3e-6 m2/s.
- Wall: linear law, E=1.0 MPa, h/D=0.1; maximum simulated area change
  0.125%.
- Friction multiplier: 1.0; no calibration to CFD.

## Final FirstBlood result

- Pressure drop: 166.122 Pa
- Hydraulic resistance: 1.626453e+08 Pa s/m3
- Resistance denominator: common CFD inlet flow 1.0213744739e-06 m3/s
- Maximum velocity: 0.6179 m/s
- Grid: 41 points/edge, 40 edges, 1601 unique axial points
- Maximum configured CFL: 0.9
- Simulated time: 1.0 s; reported mean window: final 0.2 s

## CFD comparison

- CFD pressure drop: 356.366 Pa
- Absolute pressure-drop error: 190.244 Pa
- Relative pressure-drop error: -53.38%
- CFD resistance at the prescribed inlet flow: 3.489082e+08 Pa s/m3
- FirstBlood resistance error: -53.38%
- Pressure-profile RMSE at the three CFD planes: 110.234 Pa
- Flow-profile RMSE at the three CFD planes: 0.00459 mL/s

The converged 1D model underpredicts CFD loss. This is physically plausible:
FirstBlood uses axial length and taper but contains no curvature, torsion,
secondary-flow, or bend-loss model. The discrepancy is therefore retained and
reported rather than removed by parameter tuning.

The row labelled **FirstBlood 1D** comes exclusively from the converged
41-points-per-edge FirstBlood run. No legacy Poiseuille-generated `1D` row is
used in this report or its figures.

## Geometry qualification

The 1D geometry uses the same CFD-derived centreline radius data and identical
comparison-plane arc-length bounds. It is not geometrically identical to the
full 3D lumen: a 1D taper cannot retain non-circular cross-sections, curvature,
torsion, or local secondary-flow structures.
