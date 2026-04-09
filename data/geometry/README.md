# Geometry

Place surface and volume mesh files here.

Supported formats:
- `*.stl`  – triangulated surface (STL ASCII or binary)
- `*.vtk`  – VTK legacy or XML format
- `*.vtp`  – VTK PolyData (surface)

## Recommended naming convention

| File | Contents |
|------|----------|
| `cow_surface.stl` | Full Circle-of-Willis surface mesh |
| `cow_clipped.stl` | Surface with inlet/outlet patches clipped |
| `cow_volume.vtk` | Volumetric mesh (post blockMesh/snappyHexMesh) |

## Obtaining geometry

A suitable patient-averaged Circle-of-Willis geometry can be downloaded from:
- [Vascular Model Repository](https://www.vascularmodel.com)
- [DICOM to STL](https://github.com/KitwareMedical/SlicerVMTK) via 3D Slicer + VMTK
