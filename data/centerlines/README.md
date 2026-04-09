# Centerlines

CSV files with centre-line data extracted from the vessel geometry or CFD results.

## Expected column format

```
arc_length,x,y,z,radius,velocity_3d,pressure_3d,velocity_1d,pressure_1d,pressure_0d
```

| Column | Unit | Description |
|--------|------|-------------|
| `arc_length` | m | Cumulative arc-length from inlet |
| `x`, `y`, `z` | m | Spatial coordinates |
| `radius` | m | Local vessel radius |
| `velocity_3d` | m/s | Centre-line velocity from 3-D CFD |
| `pressure_3d` | Pa | Centre-line pressure from 3-D CFD |
| `velocity_1d` | m/s | Cross-sectional mean velocity from 1-D model |
| `pressure_1d` | Pa | Pressure from 1-D model |
| `pressure_0d` | Pa | Pressure from 0-D Windkessel model |

## Extraction

Run `analysis/extract_centerline.py` to populate these files from CFD output.
The placeholder file `cow_mca_centerline.csv` is provided as a minimal example.
