# OpenFOAM Vascular CFD — Boundary Patch Creation

**Case:** `openfoam/segment_test_V2`  
**Date:** 2026-05-26  
**Stage:** Post-snappyHexMesh → CFD boundary condition preparation

---

## Final Mesh State

```
constant/polyMesh/boundary:
  artery   wall   24974 faces  ← vessel wall (no-slip)
  inlet    patch   1310 faces  ← inlet opening (~4 × 3 mm)
  outlet   patch    350 faces  ← outlet opening (~1.7 × 1.8 mm)
```

Total boundary faces: 26634 (= 24974 + 1310 + 350 ✓)

---

## What Happened and Why — Full Conceptual Explanation

### Boundary patches in OpenFOAM

Every face in an OpenFOAM mesh has a global index. Faces are sorted:

```
index 0 … 226294    → internal faces (shared between two cells)
index 226295 … 252928 → boundary faces (one cell only; in a named patch)
```

`createPatch` can only reassign faces whose index ≥ `nInternalFaces`. If you
hand it an internal face index, it throws:

```
Face XXXX specified in set inletFaces is not an external face of the mesh.
```

### Why the previous attempt failed

Four independent bugs stacked up:

| # | Bug | Effect |
|---|-----|--------|
| 1 | Stray `*/` after the `FoamFile {}` block in `topoSetDict` | OpenFOAM's dictionary parser read it as a `*` token, reported "found punctuation token `*`", never parsed `actions` → `topoSet` exited 1 |
| 2 | Old `sourceInfo {}` wrapper syntax (removed in Foundation dev) | Even without bug 1 this would have silently produced an empty set |
| 3 | Box coordinates from ParaView pointed at vessel *wall* faces (x ≈ 4.8–15.7), not at the open ends (x ≈ 3.76 and x ≈ 16.3) | Both sets had size 0 even when the syntax was right |
| 4 | Stale on-disk faceSets from a prior pure-`boxToFace` run (no `patchToFace` guard) | Those sets contained internal faces → `createPatch` rejected them |

### Why `patchToFace` + `subset boxToFace` is the correct pattern

```
patchToFace(artery)   →  {face 226295, 226296, …, 252928}   (all external ✓)
boxToFace(inlet_box)  →  {face 0, 1, …, 252928}             (internal + external)
SUBSET (intersection) →  {face ≥ 226295 AND in box}         (external only ✓)
```

Pure `boxToFace` without the `patchToFace` guard returns every face whose
centre is in the box — including internal faces — which caused the original
`createPatch` error.

---

## Correct Workflow

### 1. `system/topoSetDict`

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  dev
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    object      topoSetDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

actions
(
    // ---- Inlet: all artery boundary faces ∩ inlet box ----
    {
        name    inletFaces;
        type    faceSet;
        action  new;
        source  patchToFace;
        patch   artery;
    }
    {
        name    inletFaces;
        type    faceSet;
        action  subset;
        source  boxToFace;
        box     (3.74 -21.5 -7.8)(3.78 -17.3 -4.5);
    }

    // ---- Outlet: all artery boundary faces ∩ outlet box ----
    {
        name    outletFaces;
        type    faceSet;
        action  new;
        source  patchToFace;
        patch   artery;
    }
    {
        name    outletFaces;
        type    faceSet;
        action  subset;
        source  boxToFace;
        box     (16.25 -36.1 -6.8)(16.52 -34.3 -4.9);
    }
);

// ************************************************************************* //
```

Key syntax rules for OpenFOAM Foundation dev:
- **No** `sourceInfo {}` sub-dictionary — parameters go directly in the action block
- **No** stray `*/` after `FoamFile {}`
- `patch artery;` not `sourceInfo { patch artery; }`
- `box (min)(max);` not `sourceInfo { box ...; }`

### 2. `system/createPatchDict`

```cpp
FoamFile
{
    format      ascii;
    class       dictionary;
    object      createPatchDict;
}

pointSync false;

patches
(
    {
        name        inlet;
        patchInfo   { type patch; }
        constructFrom set;
        set         inletFaces;
    }
    {
        name        outlet;
        patchInfo   { type patch; }
        constructFrom set;
        set         outletFaces;
    }
);
```

The `artery` patch retains all faces not assigned to `inlet` or `outlet` and
automatically becomes the vessel wall.

### 3. Run sequence

```bash
# Remove any stale faceSets from previous runs
rm -f constant/polyMesh/sets/inletFaces \
      constant/polyMesh/sets/outletFaces \
      constant/polyMesh/sets/arteryPatchFaces

# Create faceSets (boundary-only intersection)
topoSet

# Split artery patch into inlet / outlet / wall
createPatch -overwrite

# Verify
checkMesh | grep -E "Boundary|patches|OK"
```

---

## How the Opening Coordinates Were Found

The original ParaView estimates were pointing at vessel wall faces. The actual
opening locations were determined by reading the polyMesh directly:

```
Mesh bounding box:  (3.75 -36.99 -8.52)  →  (16.51 -17.22 -3.49)

Inlet  (x-min end):  1269 faces at x ≈ 3.755–3.770
                     y: -21.46 .. -17.37   z: -7.72 .. -4.57
                     → use box (3.74 -21.5 -7.8)(3.78 -17.3 -4.5)

Outlet (x-max end):  326  faces at x ≈ 16.27–16.51
                     y: -36.08 .. -34.37   z: -6.73 .. -4.96
                     → use box (16.25 -36.1 -6.8)(16.52 -34.3 -4.9)
```

Diagnostic Python snippet (pure stdlib, no numpy):

```python
import re

def parse_foam_boundary(path):
    with open(path) as f: txt = f.read()
    txt = re.sub(r'/\*.*?\*/', '', txt, flags=re.DOTALL)
    txt = re.sub(r'//[^\n]*', '', txt)
    m = re.search(r'artery\s*\{([^}]*)\}', txt, re.DOTALL)
    b = m.group(1)
    return (int(re.search(r'nFaces\s+(\d+)', b).group(1)),
            int(re.search(r'startFace\s+(\d+)', b).group(1)))

def read_foam_list(path):
    with open(path) as f: txt = f.read()
    txt = re.sub(r'/\*.*?\*/', '', txt, flags=re.DOTALL)
    txt = re.sub(r'//[^\n]*', '', txt)
    m = re.search(r'(\d+)\s*\((.+)\)\s*$', txt, re.DOTALL)
    return m.group(2).strip(), int(m.group(1))

nFaces_b, startFace = parse_foam_boundary('constant/polyMesh/boundary')
raw_p, _ = read_foam_list('constant/polyMesh/points')
coords = re.findall(r'\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)', raw_p)
pts = [(float(a), float(b), float(c)) for a, b, c in coords]
raw_f, _ = read_foam_list('constant/polyMesh/faces')
face_groups = re.findall(r'\d+\(([^)]+)\)', raw_f)

def centre(g):
    ids = list(map(int, g.split()))
    return tuple(sum(pts[i][k] for i in ids)/len(ids) for k in range(3))

artery = face_groups[startFace : startFace + nFaces_b]
centres = [centre(g) for g in artery]

# Print x-extremes to locate openings
by_x = sorted(centres, key=lambda c: c[0])
for label, subset in [("Inlet (low-x)", by_x[:300]), ("Outlet (high-x)", by_x[-300:])]:
    print(f"{label}:")
    print(f"  x: {min(c[0] for c in subset):.4f} .. {max(c[0] for c in subset):.4f}")
    print(f"  y: {min(c[1] for c in subset):.4f} .. {max(c[1] for c in subset):.4f}")
    print(f"  z: {min(c[2] for c in subset):.4f} .. {max(c[2] for c in subset):.4f}")
```

---

## Boundary Conditions (`0/`)

### `0/U`

```cpp
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0.01 0 0);

boundaryField
{
    // Vessel wall — no-slip
    artery
    {
        type    noSlip;
    }

    // Inlet — steady plug-flow; adjust direction and speed as needed
    inlet
    {
        type    fixedValue;
        value   uniform (0.1 0 0);   // 0.1 m/s in +x; refine with actual normal
    }

    // Outlet — zero-gradient (outflow)
    outlet
    {
        type    zeroGradient;
    }
}
```

### `0/p`

```cpp
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;

boundaryField
{
    artery
    {
        type    zeroGradient;
    }

    // Inlet — zero-gradient pressure (velocity-driven)
    inlet
    {
        type    zeroGradient;
    }

    // Outlet — fixed reference pressure
    outlet
    {
        type    fixedValue;
        value   uniform 0;
    }
}
```

---

## Next Steps

1. **Inlet velocity direction** — the value `(0.1 0 0)` assumes the vessel axis
   is aligned with +x at the inlet. Apply *Surface Normals* on the inlet patch
   in ParaView and use the inward normal as the velocity direction.

2. **`transportProperties`** — set kinematic viscosity to blood:
   ```
   nu  [0 2 -1 0 0 0 0]  3.5e-6;   // m²/s (Newtonian approximation)
   ```

3. **Minor mesh quality** — 2 face-pyramid errors and 3 skew faces are
   pre-existing snappyHexMesh artefacts near the cut planes. If the solver
   diverges near boundaries, add to `fvSolution`:
   ```cpp
   relaxed { maxNonOrtho 75; }
   ```

4. **`topoSet` deprecation** — in OpenFOAM-dev, `topoSet` is superseded by
   `createZones`. The `topoSet` + `createPatch` workflow remains fully
   functional but will print a deprecation warning.
