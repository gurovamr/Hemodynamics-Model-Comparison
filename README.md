# Hemodynamics Model Comparison – Circle of Willis

A scientific computing project for multiscale blood flow modelling in the **Circle of Willis**.
Three levels of fidelity are compared along vessel centrelines:

| Model | Fidelity | Tool |
|-------|----------|------|
| 3-D CFD | Full Navier-Stokes | OpenFOAM |
| 1-D | Cross-sectional averaged | Python (numpy) |
| 0-D | Lumped-parameter (Windkessel) | Python (numpy) |

---

## Project goal

Understand how well reduced-order models (1D, 0D) reproduce the velocity and pressure fields
predicted by a full 3-D CFD simulation in the cerebral vasculature.
Results are compared along vessel centre-lines extracted from the 3-D geometry.

---

## Repository layout

```
Hemodynamics-Model-Comparison/
├── data/
│   ├── geometry/          # STL / VTK surface and volume meshes
│   └── centerlines/       # Centre-line CSV files (x, y, z, radius, …)
├── cfd/
│   └── cow_case/          # OpenFOAM case template for the Circle of Willis
│       ├── 0/             # Initial & boundary conditions
│       ├── constant/      # Physical properties and mesh
│       └── system/        # Solver settings (controlDict, fvSolution, …)
├── rom/
│   ├── model_1d.py        # 1-D reduced-order model (pulse-wave propagation)
│   └── model_0d.py        # 0-D lumped-parameter Windkessel model
├── analysis/
│   ├── extract_centerline.py   # Extract velocity/pressure along a centre-line
│   ├── plot_centerline.py      # Plot centre-line profiles
│   └── compare_models.py       # Side-by-side comparison of 3D / 1D / 0D
├── docs/
│   ├── methodology.md     # Governing equations and modelling assumptions
│   └── references.md      # Key literature
└── README.md
```

---

## Workflow

```
1. Geometry preparation
   └─ Place STL/VTK files in data/geometry/

2. Mesh generation
   └─ Run blockMesh / snappyHexMesh inside cfd/cow_case/

3. 3-D CFD simulation
   └─ cd cfd/cow_case && foamRun   (or: simpleFoam / pimpleFoam)

4. Centre-line extraction
   └─ python analysis/extract_centerline.py

5. Reduced-order models
   └─ python rom/model_1d.py
   └─ python rom/model_0d.py

6. Comparison & plotting
   └─ python analysis/compare_models.py
   └─ python analysis/plot_centerline.py
```

---

## Quick start

### Requirements

* **OpenFOAM** v9+ (for 3-D CFD)
* **Python 3.9+** with `numpy` and `matplotlib`

```bash
pip install numpy matplotlib
```

### Run the 0-D model (no CFD required)

```bash
python rom/model_0d.py
```

### Run the 1-D model (no CFD required)

```bash
python rom/model_1d.py
```

### Run the OpenFOAM case

```bash
cd cfd/cow_case
blockMesh
foamRun          # or: simpleFoam / pimpleFoam depending on your OF version
```

### Compare all models

```bash
python analysis/compare_models.py
```

---

## Data

Centre-line CSV files should have the columns:

```
arc_length, x, y, z, radius, velocity_3d, pressure_3d
```

Example placeholder files are provided in `data/centerlines/`.

---

## Contributing

1. Fork the repository
2. Create a feature branch
3. Submit a pull request with a clear description of your changes

---

## License

[MIT](LICENSE)
