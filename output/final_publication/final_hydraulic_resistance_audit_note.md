# Final Hydraulic Resistance Audit Note

## Why the old and new plots conflicted

The old and new figures came from different pipelines and are not directly comparable.

1. Unit system difference:
- Old figure often plotted mmHg s/mL.
- New figure plotted 10^8 Pa s/m^3.
- These are convertible, but mixing plots without explicit conversion creates apparent disagreement.

2. Model identity difference:
- Old pipeline included an older pseudo-1D/network-1D row (e.g., "1D" in model_comparison_common_domain.csv), not the final converged FirstBlood standalone run.
- Final report now uses only FirstBlood grid_41 for the 1D entry.

3. Domain-definition difference:
- Historical scripts also used full-domain and label-based selections in some outputs.
- Final table here enforces the common CFD domain only: s=19.671598 to 35.892704 mm.

4. CFD extraction method and denominator difference:
- Final CFD uses area-averaged inlet/outlet slice pressures and Q_mean from slice-integrated flow.
- Some prior summaries used alternate denominator conventions (e.g., inlet-only or other table flow values).
- Final table uses one denominator policy for comparison rows: Q_used = CFD mean slice flow.

5. FirstBlood resistance definition difference:
- Raw FirstBlood resistance uses its own mean terminal flow (R_raw = DeltaP/Q_FirstBloodMean).
- Final comparison resistance uses common CFD denominator (R_cmp = DeltaP/Q_CFD_mean) for apples-to-apples across all models.

## Final policy applied in this notebook
- One equation: R = DeltaP / Q
- DeltaP in Pa, Q in m^3/s, R in Pa s/m^3
- Common CFD domain for all models
- One authoritative CSV: final_hydraulic_resistance_audit.csv
- One resistance figure: final_hydraulic_resistance_common_domain.png
- One pressure-drop figure: final_pressure_drop_common_domain.png
