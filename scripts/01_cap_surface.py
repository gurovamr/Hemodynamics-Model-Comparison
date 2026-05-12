#!/usr/bin/env python3
"""
01_cap_surface.py
─────────────────
Prepares the CoW MR surface mesh for OpenFOAM CFD by:
  1. Identifying "dome" triangles at each vessel inlet/outlet (sphere + half-plane
     selection relative to the centerline endpoint and local vessel direction)
  2. Removing dome triangles from the wall mesh
  3. Filling each removed dome with a flat cap (fan triangulation from centroid)
  4. Exporting per-patch STL files + combined multi-solid STL for snappyHexMesh

Coordinate convention: VTP files are in mm; STL output is in metres (SI).

Run with the tm2 conda environment (has vtk 9.2.6):
    /home/gurovamr/miniforge3/envs/tm2/bin/python3 scripts/01_cap_surface.py
"""

import os, json, collections
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

REPO     = "/home/gurovamr/OpenFOAM/gurovamr-dev/run/TrackModule2/Hemodynamics-Model-Comparison"
DATA_DIR = os.path.join(REPO, "data")
OF_CASE  = os.path.join(REPO, "openfoam", "cow_mr_cfd")
STL_DIR  = os.path.join(OF_CASE, "constant", "triSurface")

MESH_FILE    = os.path.join(DATA_DIR, "mesh_mr_025.vtp")
GRAPH_FILE   = os.path.join(DATA_DIR, "graph_mr_025.vtp")
NODE_FILE    = os.path.join(DATA_DIR, "node_mr_025.json")
FEATURE_FILE = os.path.join(DATA_DIR, "features_mr_025.json")

MM_TO_M = 1e-3

LABEL_TO_PATCH = {
    ("1",  "BA start"):  "inlet_BA",
    ("2",  "PCA end"):   "outlet_R_PCA",
    ("3",  "PCA end"):   "outlet_L_PCA",
    ("4",  "ICA start"): "inlet_R_ICA",
    ("6",  "ICA start"): "inlet_L_ICA",
    ("7",  "MCA end"):   "outlet_L_MCA",
    ("11", "ACA end"):   "outlet_R_ACA",
    ("12", "ACA end"):   "outlet_L_ACA",
}

SPHERE_MULT   = 3.0
INWARD_THRESH = 1.0


def load_vtp(path):
    r = vtk.vtkXMLPolyDataReader()
    r.SetFileName(path)
    r.Update()
    return r.GetOutput()


def load_endpoints(node_file, graph_file, feature_file):
    with open(node_file) as f:
        nodes = json.load(f)
    with open(feature_file) as f:
        feat = json.load(f)

    g = load_vtp(graph_file)
    gpts = vtk_to_numpy(g.GetPoints().GetData())
    graph_adj = collections.defaultdict(list)
    for c in range(g.GetNumberOfCells()):
        cell = g.GetCell(c)
        p0, p1 = cell.GetPointId(0), cell.GetPointId(1)
        graph_adj[p0].append(gpts[p1])
        graph_adj[p1].append(gpts[p0])

    endpoints = {}
    for seg_id, seg_nodes in nodes.items():
        for node_type, node_list in seg_nodes.items():
            key = (seg_id, node_type)
            if key not in LABEL_TO_PATCH:
                continue
            patch = LABEL_TO_PATCH[key]
            n = node_list[0]
            nid = n["id"]
            coord = np.array(n["coords"], dtype=float)

            neighbours = graph_adj.get(nid, [])
            if neighbours:
                inward = np.array(neighbours[0]) - coord
                inward /= np.linalg.norm(inward)
            else:
                inward = np.array([0., 0., 1.])

            radius = 2.0
            seg_feat = feat.get(seg_id, {})
            for _, seg_list in seg_feat.items():
                if isinstance(seg_list, list) and seg_list and "radius" in seg_list[0]:
                    radius = float(seg_list[0]["radius"]["mean"])
                    break

            endpoints[patch] = {
                "coord":      coord,
                "inward_dir": inward,
                "radius":     radius,
            }
    return endpoints


def select_dome_triangles(all_pts, all_tris, ep):
    coord   = ep["coord"]
    inward  = ep["inward_dir"]
    R       = SPHERE_MULT * ep["radius"]
    max_dep = INWARD_THRESH * ep["radius"]

    centroids = (all_pts[all_tris[:, 0]] +
                 all_pts[all_tris[:, 1]] +
                 all_pts[all_tris[:, 2]]) / 3.0

    dists = np.linalg.norm(centroids - coord, axis=1)
    depth = np.einsum("ij,j->i", centroids - coord, inward)
    mask  = (dists < R) & (depth < max_dep)
    return np.where(mask)[0]


def build_boundary_loops(all_pts, all_tris, removed_indices):
    removed_set = set(int(i) for i in removed_indices)

    edge_tris = collections.defaultdict(list)
    for ti, tri in enumerate(all_tris):
        for j in range(3):
            e = tuple(sorted([int(tri[j]), int(tri[(j+1) % 3])]))
            edge_tris[e].append(ti)

    boundary_edges = []
    for e, tlist in edge_tris.items():
        in_removed = [t in removed_set for t in tlist]
        if len(tlist) == 2 and any(in_removed) and not all(in_removed):
            boundary_edges.append(e)
        elif len(tlist) == 1 and tlist[0] in removed_set:
            boundary_edges.append(e)

    if not boundary_edges:
        return []

    adj = collections.defaultdict(list)
    for e in boundary_edges:
        adj[e[0]].append(e[1])
        adj[e[1]].append(e[0])

    visited_pts = set()
    loops = []
    for start in sorted(adj.keys()):
        if start in visited_pts:
            continue
        loop = [start]
        visited_pts.add(start)
        prev, curr = -1, start
        while True:
            nbs = [nb for nb in adj[curr] if nb != prev]
            if not nbs:
                break
            nxt = nbs[0]
            if nxt == start:
                break
            if nxt in visited_pts:
                break
            visited_pts.add(nxt)
            loop.append(nxt)
            prev, curr = curr, nxt
        if len(loop) >= 3:
            loops.append(np.array([all_pts[pid] for pid in loop]))
    return loops


def make_flat_cap(loop_coords, inward_dir):
    centroid = loop_coords.mean(axis=0)
    n = len(loop_coords)
    cap_pts = np.vstack([loop_coords, centroid])
    c_idx = n

    tris = np.array([[i, (i + 1) % n, c_idx] for i in range(n)])

    v0 = cap_pts[tris[:, 1]] - cap_pts[tris[:, 0]]
    v1 = cap_pts[tris[:, 2]] - cap_pts[tris[:, 0]]
    normals = np.cross(v0, v1)
    cap_normal = normals.sum(axis=0)
    norm = np.linalg.norm(cap_normal)
    if norm > 0:
        cap_normal /= norm

    if np.dot(cap_normal, inward_dir) > 0:
        tris = tris[:, [0, 2, 1]]

    return cap_pts, tris


def solid_to_stl_text(pts_mm, tris, solid_name):
    pts_m = pts_mm * MM_TO_M
    lines = [f"solid {solid_name}"]
    for tri in tris:
        v = pts_m[tri]
        e1 = v[1] - v[0]; e2 = v[2] - v[0]
        n = np.cross(e1, e2)
        nm = np.linalg.norm(n)
        n = n / nm if nm > 0 else np.array([0., 0., 1.])
        lines.append(f"  facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}")
        lines.append("    outer loop")
        for vid in tri:
            vv = pts_m[vid]
            lines.append(f"      vertex {vv[0]:.8e} {vv[1]:.8e} {vv[2]:.8e}")
        lines.append("    endloop")
        lines.append("  endfacet")
    lines.append(f"endsolid {solid_name}")
    return "\n".join(lines)


def write_stl_one(pts_mm, tris, solid_name, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(solid_to_stl_text(pts_mm, tris, solid_name) + "\n")
    print(f"  {os.path.basename(path):38s}  {len(tris):6d} triangles")


def write_stl_multi(solids, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for pts_mm, tris, solid_name in solids:
            f.write(solid_to_stl_text(pts_mm, tris, solid_name) + "\n")
    total = sum(len(t) for _, t, _ in solids)
    print(f"  {os.path.basename(path):38s}  {len(solids)} solids, {total} tris total")


def main():
    os.makedirs(STL_DIR, exist_ok=True)

    print("Loading surface mesh ...")
    pd = load_vtp(MESH_FILE)
    all_pts = vtk_to_numpy(pd.GetPoints().GetData()).astype(float)
    all_tris_raw = []
    for c in range(pd.GetNumberOfCells()):
        cell = pd.GetCell(c)
        all_tris_raw.append([cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())])
    all_tris = np.array(all_tris_raw, dtype=int)
    print(f"  {len(all_pts)} pts, {len(all_tris)} tris  bounds x["
          f"{all_pts[:,0].min():.1f},{all_pts[:,0].max():.1f}] "
          f"y[{all_pts[:,1].min():.1f},{all_pts[:,1].max():.1f}] "
          f"z[{all_pts[:,2].min():.1f},{all_pts[:,2].max():.1f}] mm")

    print("\nLoading endpoint data ...")
    endpoints = load_endpoints(NODE_FILE, GRAPH_FILE, FEATURE_FILE)
    for patch, ep in endpoints.items():
        print(f"  {patch:20s}: r={ep['radius']:.2f} mm  "
              f"search_R={SPHERE_MULT*ep['radius']:.2f} mm  "
              f"coord={ep['coord'].round(1)}")

    print("\nSelecting dome triangles ...")
    removed_global = set()
    patch_removed  = {}
    for patch, ep in endpoints.items():
        sel = select_dome_triangles(all_pts, all_tris, ep)
        sel_clean = [i for i in sel if i not in removed_global]
        patch_removed[patch] = set(int(i) for i in sel_clean)
        removed_global.update(patch_removed[patch])
        print(f"  {patch:20s}: {len(sel_clean):5d} dome triangles")

    wall_indices = [i for i in range(len(all_tris)) if i not in removed_global]
    wall_tris    = all_tris[wall_indices]
    print(f"\nWall triangles remaining: {len(wall_tris)} / {len(all_tris)}")

    print("\nBuilding flat caps ...")
    cap_solids = []
    for patch, ep in endpoints.items():
        removed_for_this = patch_removed[patch]
        if not removed_for_this:
            print(f"  WARNING {patch}: 0 dome tris — skipping")
            continue

        loops = build_boundary_loops(all_pts, all_tris, removed_for_this)
        if not loops:
            print(f"  WARNING {patch}: no boundary loop found")
            continue

        loops.sort(key=lambda l: -len(l))
        loop_pts = loops[0]
        if len(loops) > 1:
            print(f"  {patch}: {len(loops)} loops, using largest ({len(loop_pts)} pts)")
        cap_pts, cap_tris = make_flat_cap(loop_pts, ep["inward_dir"])
        cap_solids.append((cap_pts, cap_tris, patch))
        print(f"  {patch:20s}: cap {len(cap_tris)} tris  loop {len(loop_pts)} pts")

    print("\nWriting STL files ...")
    all_solids = [(all_pts, wall_tris, "wall")]
    write_stl_one(all_pts, wall_tris, "wall", os.path.join(STL_DIR, "wall.stl"))
    for cap_pts, cap_tris, patch_name in cap_solids:
        write_stl_one(cap_pts, cap_tris, patch_name,
                      os.path.join(STL_DIR, f"{patch_name}.stl"))
        all_solids.append((cap_pts, cap_tris, patch_name))
    write_stl_multi(all_solids, os.path.join(STL_DIR, "cow_mr_all.stl"))

    bc_summary = {}
    for cap_pts, cap_tris, patch_name in cap_solids:
        centroid = cap_pts[:-1].mean(axis=0)
        ep = endpoints[patch_name]
        bc_summary[patch_name] = {
            "type":        "inlet" if "inlet" in patch_name else "outlet",
            "centroid_mm": [float(round(c, 4)) for c in centroid],
            "radius_mm":   float(round(ep["radius"], 4)),
            "inward_dir":  [float(round(d, 6)) for d in ep["inward_dir"]],
            "n_cap_tris":  int(len(cap_tris)),
        }
    os.makedirs(OF_CASE, exist_ok=True)
    with open(os.path.join(OF_CASE, "patch_summary.json"), "w") as f:
        json.dump(bc_summary, f, indent=2)

    print(f"\nDone. STL files in: {STL_DIR}")
    print("Run the OpenFOAM pipeline with: cd openfoam/cow_mr_cfd && ./Allrun")


if __name__ == "__main__":
    main()
