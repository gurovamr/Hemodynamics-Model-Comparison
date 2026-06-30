# TM2 standalone FirstBlood parameter sensitivity

> **Superseded-domain warning:** this sensitivity study used the full
> 35.893 mm source centreline. The audited CFD comparison domain is only
> global s=19.672–35.893 mm (16.221 mm). These results diagnose the old model
> but must not be used as the final CFD comparison.

## Baseline

The completed TM2 baseline used 40 tapered edges, five numerical points per
edge, linear wall mechanics, `h/D=0.1`, `E=1.0 MPa`, friction multiplier 1.0,
`rho=1060 kg/m3`, and `nu=3.3e-6 m2/s`.

Its late-time result was:

- pressure drop: 659.386 Pa
- hydraulic resistance: 6.41243e8 Pa s/m3
- terminal flow mismatch: 1.346%
- axial flow spread: approximately 2.3%

The CFD pressure drop is 356.366 Pa, giving a baseline FirstBlood excess of
303.020 Pa.

## Abel_ref2 comparison

Abel_ref2 uses the same linear wall law, MacCormack solver, `h/D=0.1`,
Poisson ratio 0.5, CFL 0.9, and `k1=2e6`, `k2=-2253`, `k3=8.65e4`.
The Olufsen constants are parsed but inactive under the linear wall law.

Across its 103 arteries, Abel_ref2 has:

- Young's modulus: 0.292–2.59 MPa; arithmetic mean 1.01005 MPa
- friction multiplier: 2.75
- division points: 5–53; arithmetic mean 11.466
- density: 1055 kg/m3
- kinematic viscosity: 3.0e-6 m2/s

The TM2 and Abel edge lengths differ greatly, so the Abel arithmetic mean of
11 points is a reference substitution, not a claim that equal point counts
give equal physical grid spacing.

## OAT results

| Rank | One changed parameter | Pressure drop (Pa) | Change from baseline | Resistance (Pa s/m3) |
|---:|---|---:|---:|---:|
| 1 | Points per edge: 5 -> 11 | 383.422 | -275.964 (-41.85%) | 3.74881e8 |
| 2 | Friction multiplier: 1 -> 2.75 | 1188.860 | +529.475 (+80.30%) | 1.15603e9 |
| 3 | Viscosity: 3.3e-6 -> 3.0e-6 m2/s | 631.746 | -27.640 (-4.19%) | 6.14368e8 |
| 4 | Density: 1060 -> 1055 kg/m3 | 657.259 | -2.127 (-0.32%) | 6.39165e8 |
| 5 | E: 1.0 -> 1.01005 MPa | 661.470 | +2.084 (+0.32%) | 6.43247e8 |

Identical Abel substitutions for `h/D`, `k1-k3`, wall law, and solver produce
zero input change and therefore zero OAT response.

## Attribution of the 659 Pa result

The five-point mesh is the primary cause. Changing only to 11 points reduces
the FirstBlood–CFD difference from 303.020 Pa to 27.056 Pa, removing 91.1% of
the original gap. It also improves terminal flow mismatch to 0.276% and axial
flow spread to 0.498%, both within the chosen 0.5% acceptance threshold.

Neither wall thickness nor Young's modulus explains the gap: `h/D` is
identical to Abel and the TM2 modulus is already essentially Abel's mean.
Density is negligible and viscosity accounts for only 27.6 Pa on the
under-resolved baseline.

Abel's friction multiplier is highly influential but is not responsible for
the current excess because the TM2 model uses 1.0. Replacing it with 2.75
would increase the discrepancy substantially.

The remaining 27 Pa refined-grid difference can contain discretisation error,
the compliant-wall approximation, radius resampling, and differences between
the axial 1D equations and the curved 3D CFD flow. A formal grid-convergence
study above 11 points per edge is required before assigning that residual to
physics.
