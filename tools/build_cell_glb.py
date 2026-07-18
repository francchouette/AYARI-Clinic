#!/usr/bin/env python3
"""Build the AYARI reference-driven medical-luxury cutaway cell GLB.

The nucleus and mitochondrial morphologies are extracted from measured
OpenOrganelle ``jrc_hela-2`` FIB-SEM segmentation labels.  The display cell
envelope and chromosome/telomere/DNA teaching overlays remain intentionally
modelled.  Per-node extras make that distinction machine-readable.

The generated asset is self-contained glTF 2.0 Binary with metre-scale Y-up
geometry. All transforms are baked into vertex positions so viewer-side
selection, framing, and emissive highlighting work on named submeshes.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from PIL import Image, ImageDraw, ImageFilter, ImageFont
from scipy import ndimage
from skimage import measure
import trimesh


ARRAY_BUFFER = 34962
ELEMENT_ARRAY_BUFFER = 34963
FLOAT = 5126
UNSIGNED_INT = 5125


@dataclass
class Geometry:
    vertices: np.ndarray
    faces: np.ndarray
    normals: np.ndarray | None = None


def normalize(values: np.ndarray, axis: int = -1) -> np.ndarray:
    lengths = np.linalg.norm(values, axis=axis, keepdims=True)
    return values / np.maximum(lengths, 1e-12)


def clean_geometry(vertices: np.ndarray, faces: np.ndarray) -> Geometry:
    vertices = np.asarray(vertices, dtype=np.float64)
    faces = np.asarray(faces, dtype=np.int64)
    triangles = vertices[faces]
    raw = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    keep = np.einsum("ij,ij->i", raw, raw) > 1e-20
    faces = faces[keep]
    if not len(faces):
        raise ValueError("Geometry has no non-degenerate triangles")

    used, inverse = np.unique(faces.reshape(-1), return_inverse=True)
    vertices = vertices[used]
    faces = inverse.reshape((-1, 3))
    triangles = vertices[faces]
    face_normals = normalize(
        np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    )
    normals = np.zeros_like(vertices)
    for corner in range(3):
        np.add.at(normals, faces[:, corner], face_normals)
    normal_lengths = np.linalg.norm(normals, axis=1)
    weak = normal_lengths < 1e-7
    if np.any(weak):
        # Non-manifold measured membrane sheets can cancel adjacent face
        # normals at a small number of shared vertices.  A stable radial
        # fallback keeps glTF NORMAL accessors valid and deterministic.
        fallback = vertices[weak] - vertices.mean(axis=0)
        fallback_lengths = np.linalg.norm(fallback, axis=1)
        fallback[fallback_lengths < 1e-9] = (0.0, 1.0, 0.0)
        normals[weak] = fallback
    normals = normalize(normals)
    return Geometry(
        vertices.astype(np.float32),
        faces.astype(np.uint32),
        normals.astype(np.float32),
    )


def concatenate(geometries: Iterable[Geometry]) -> Geometry:
    vertices: list[np.ndarray] = []
    faces: list[np.ndarray] = []
    offset = 0
    for geometry in geometries:
        vertices.append(np.asarray(geometry.vertices))
        faces.append(np.asarray(geometry.faces, dtype=np.uint32) + offset)
        offset += len(geometry.vertices)
    return clean_geometry(np.vstack(vertices), np.vstack(faces))


def rotation_matrix(rx: float = 0.0, ry: float = 0.0, rz: float = 0.0) -> np.ndarray:
    sx, cx = math.sin(rx), math.cos(rx)
    sy, cy = math.sin(ry), math.cos(ry)
    sz, cz = math.sin(rz), math.cos(rz)
    mx = np.array(((1, 0, 0), (0, cx, -sx), (0, sx, cx)), dtype=float)
    my = np.array(((cy, 0, sy), (0, 1, 0), (-sy, 0, cy)), dtype=float)
    mz = np.array(((cz, -sz, 0), (sz, cz, 0), (0, 0, 1)), dtype=float)
    return mz @ my @ mx


def transform_geometry(
    geometry: Geometry,
    translation: Sequence[float] = (0.0, 0.0, 0.0),
    rotation: np.ndarray | None = None,
) -> Geometry:
    matrix = np.eye(3) if rotation is None else np.asarray(rotation, dtype=float)
    vertices = np.asarray(geometry.vertices, dtype=float) @ matrix.T + np.asarray(translation, dtype=float)
    return clean_geometry(vertices, geometry.faces)


def orient_faces_outward(vertices: np.ndarray, faces: np.ndarray, centre: Sequence[float]) -> np.ndarray:
    faces = np.asarray(faces, dtype=np.int64).copy()
    triangles = vertices[faces]
    normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    radial = triangles.mean(axis=1) - np.asarray(centre, dtype=float)
    flip = np.einsum("ij,ij->i", normals, radial) < 0
    faces[flip, 1], faces[flip, 2] = faces[flip, 2], faces[flip, 1].copy()
    return faces


def uv_sphere(
    centre: Sequence[float],
    radii: Sequence[float],
    around: int = 48,
    vertical: int = 28,
) -> Geometry:
    centre = np.asarray(centre, dtype=float)
    rx, ry, rz = radii
    vertices: list[list[float]] = [[0.0, ry, 0.0]]
    for latitude in range(1, vertical):
        theta = math.pi * latitude / vertical
        sin_theta, cos_theta = math.sin(theta), math.cos(theta)
        for longitude in range(around):
            phi = math.tau * longitude / around
            vertices.append(
                [rx * sin_theta * math.cos(phi), ry * cos_theta, rz * sin_theta * math.sin(phi)]
            )
    bottom = len(vertices)
    vertices.append([0.0, -ry, 0.0])

    faces: list[list[int]] = []
    first = 1
    for longitude in range(around):
        nxt = (longitude + 1) % around
        faces.append([0, first + longitude, first + nxt])
    for latitude in range(vertical - 2):
        row = 1 + latitude * around
        next_row = row + around
        for longitude in range(around):
            nxt = (longitude + 1) % around
            a, b = row + longitude, row + nxt
            c, d = next_row + longitude, next_row + nxt
            faces.extend(([a, c, d], [a, d, b]))
    last = 1 + (vertical - 2) * around
    for longitude in range(around):
        nxt = (longitude + 1) % around
        faces.append([bottom, last + nxt, last + longitude])

    vertices_array = np.asarray(vertices, dtype=float) + centre
    faces_array = orient_faces_outward(vertices_array, np.asarray(faces), centre)
    return clean_geometry(vertices_array, faces_array)


def cutaway_shell(
    outer_radii: Sequence[float] = (1.0, 0.92, 0.86),
    thickness: float = 0.035,
    gap_radians: float = math.radians(108),
    around: int = 160,
    vertical: int = 88,
) -> Geometry:
    """Closed thin ellipsoid shell with a wedge removed around front (+Z)."""
    outer = np.asarray(outer_radii, dtype=float)
    inner = outer - thickness
    phi_start = math.pi / 2 + gap_radians / 2
    phi_end = math.pi / 2 - gap_radians / 2 + math.tau
    epsilon = 1e-4
    vertices: list[list[float]] = []

    def surface(radii: np.ndarray) -> list[list[int]]:
        grid: list[list[int]] = []
        for latitude in range(vertical):
            lat = -math.pi / 2 + epsilon + (math.pi - 2 * epsilon) * latitude / (vertical - 1)
            row: list[int] = []
            for longitude in range(around):
                phi = phi_start + (phi_end - phi_start) * longitude / (around - 1)
                row.append(len(vertices))
                vertices.append(
                    [
                        radii[0] * math.cos(lat) * math.cos(phi),
                        radii[1] * math.sin(lat),
                        radii[2] * math.cos(lat) * math.sin(phi),
                    ]
                )
            grid.append(row)
        return grid

    outer_grid = surface(outer)
    inner_grid = surface(inner)
    faces: list[list[int]] = []

    for lat in range(vertical - 1):
        for lon in range(around - 1):
            a, b = outer_grid[lat][lon], outer_grid[lat][lon + 1]
            c, d = outer_grid[lat + 1][lon], outer_grid[lat + 1][lon + 1]
            faces.extend(([a, c, d], [a, d, b]))
            ia, ib = inner_grid[lat][lon], inner_grid[lat][lon + 1]
            ic, id_ = inner_grid[lat + 1][lon], inner_grid[lat + 1][lon + 1]
            faces.extend(([ia, id_, ic], [ia, ib, id_]))

    # Seal the two radial cut faces and the tiny pole rings, preserving real shell thickness.
    for lat in range(vertical - 1):
        for lon in (0, around - 1):
            oa, ob = outer_grid[lat][lon], outer_grid[lat + 1][lon]
            ia, ib = inner_grid[lat][lon], inner_grid[lat + 1][lon]
            if lon == 0:
                faces.extend(([oa, ia, ib], [oa, ib, ob]))
            else:
                faces.extend(([oa, ib, ia], [oa, ob, ib]))
    for lat in (0, vertical - 1):
        for lon in range(around - 1):
            oa, ob = outer_grid[lat][lon], outer_grid[lat][lon + 1]
            ia, ib = inner_grid[lat][lon], inner_grid[lat][lon + 1]
            if lat == 0:
                faces.extend(([oa, ib, ia], [oa, ob, ib]))
            else:
                faces.extend(([oa, ia, ib], [oa, ib, ob]))
    return clean_geometry(np.asarray(vertices), np.asarray(faces))


def tube(
    points: Sequence[Sequence[float]],
    radii: float | Sequence[float],
    radial_segments: int = 10,
    caps: bool = True,
) -> Geometry:
    points_array = np.asarray(points, dtype=float)
    if len(points_array) < 2:
        raise ValueError("Tube needs at least two points")
    if np.isscalar(radii):
        radius_values = np.full(len(points_array), float(radii))
    else:
        radius_values = np.asarray(radii, dtype=float)
    tangents = np.empty_like(points_array)
    tangents[0] = points_array[1] - points_array[0]
    tangents[-1] = points_array[-1] - points_array[-2]
    if len(points_array) > 2:
        tangents[1:-1] = points_array[2:] - points_array[:-2]
    tangents = normalize(tangents)

    normals = np.empty_like(points_array)
    reference = np.array((0.0, 1.0, 0.0))
    if abs(float(np.dot(tangents[0], reference))) > 0.9:
        reference = np.array((1.0, 0.0, 0.0))
    normals[0] = normalize(np.cross(tangents[0], reference))
    for index in range(1, len(points_array)):
        projected = normals[index - 1] - tangents[index] * np.dot(normals[index - 1], tangents[index])
        if np.linalg.norm(projected) < 1e-8:
            ref = np.array((1.0, 0.0, 0.0)) if abs(tangents[index, 0]) < 0.9 else np.array((0.0, 0.0, 1.0))
            projected = np.cross(tangents[index], ref)
        normals[index] = normalize(projected)
    binormals = normalize(np.cross(tangents, normals))

    vertices: list[list[float]] = []
    for index, point in enumerate(points_array):
        for segment in range(radial_segments):
            angle = math.tau * segment / radial_segments
            radial = normals[index] * math.cos(angle) + binormals[index] * math.sin(angle)
            vertices.append((point + radial * radius_values[index]).tolist())
    faces: list[list[int]] = []
    for ring in range(len(points_array) - 1):
        for segment in range(radial_segments):
            nxt = (segment + 1) % radial_segments
            a = ring * radial_segments + segment
            b = ring * radial_segments + nxt
            c = (ring + 1) * radial_segments + segment
            d = (ring + 1) * radial_segments + nxt
            faces.extend(([a, d, b], [a, c, d]))
    if caps:
        start = len(vertices)
        vertices.append(points_array[0].tolist())
        end = len(vertices)
        vertices.append(points_array[-1].tolist())
        last_ring = (len(points_array) - 1) * radial_segments
        for segment in range(radial_segments):
            nxt = (segment + 1) % radial_segments
            faces.append([start, nxt, segment])
            faces.append([end, last_ring + segment, last_ring + nxt])
    return clean_geometry(np.asarray(vertices), np.asarray(faces))


def bean_body(
    centre: Sequence[float],
    length: float,
    width: float,
    depth: float,
    bend: float,
    rotation: np.ndarray,
) -> Geometry:
    geometry = uv_sphere((0.0, 0.0, 0.0), (length / 2, width / 2, depth / 2), around=52, vertical=26)
    vertices = np.asarray(geometry.vertices, dtype=float)
    longitudinal = vertices[:, 0] / (length / 2)
    # Curved centreline plus a restrained inner notch gives a clean bean silhouette.
    vertices[:, 1] += bend * (longitudinal**2 - 0.32)
    inner = np.exp(-5.0 * longitudinal**2) * np.clip(vertices[:, 1] / (width / 2), 0.0, 1.0)
    vertices[:, 1] -= 0.035 * inner
    local = clean_geometry(vertices, geometry.faces)
    return transform_geometry(local, centre, rotation)


def cristae_bundle(
    centre: Sequence[float],
    length: float,
    width: float,
    depth: float,
    bend: float,
    rotation: np.ndarray,
) -> Geometry:
    folds: list[Geometry] = []
    for fold_index, x_factor in enumerate((-0.58, -0.29, 0.0, 0.29, 0.58)):
        x = x_factor * length / 2
        points: list[list[float]] = []
        for step in range(15):
            t = step / 14
            y = (t - 0.5) * width * 0.70 + bend * (x_factor**2 - 0.32)
            z = depth * 0.22 * math.sin(math.tau * t + fold_index * 0.72)
            points.append([x, y, z])
        fold = tube(points, 0.009, radial_segments=8, caps=True)
        folds.append(transform_geometry(fold, centre, rotation))
    return concatenate(folds)


def chromosome(
    centre: Sequence[float],
    scale: float,
    rotation: np.ndarray,
    depth_phase: float,
) -> tuple[Geometry, list[np.ndarray]]:
    centre_array = np.asarray(centre, dtype=float)
    arm_vectors = (
        np.array((-0.62, 1.0, 0.0)),
        np.array((0.62, 1.0, 0.0)),
        np.array((-0.62, -1.0, 0.0)),
        np.array((0.62, -1.0, 0.0)),
    )
    arms: list[Geometry] = []
    endpoints: list[np.ndarray] = []
    for arm_index, vector in enumerate(arm_vectors):
        points: list[np.ndarray] = []
        for step in range(10):
            t = step / 9
            local = vector * scale * t
            local[2] += math.sin(math.pi * t) * scale * 0.16 * math.sin(depth_phase + arm_index)
            points.append(local)
        radii = np.linspace(scale * 0.20, scale * 0.135, len(points))
        arm = tube(points, radii, radial_segments=12, caps=False)
        arms.append(transform_geometry(arm, centre_array, rotation))
        endpoints.append(np.asarray(points[-1]) @ rotation.T + centre_array)
    centre_piece = uv_sphere((0.0, 0.0, 0.0), (scale * 0.24, scale * 0.22, scale * 0.20), around=20, vertical=12)
    arms.append(transform_geometry(centre_piece, centre_array, rotation))
    return concatenate(arms), endpoints


def dna_helix(
    height: float = 1.10,
    turns: float = 2.2,
    centre_x: float = 0.87,
    centre_z: float = 0.26,
) -> tuple[Geometry, Geometry, Geometry]:
    strand_a: list[list[float]] = []
    strand_b: list[list[float]] = []
    samples = 126
    helix_radius = 0.062
    for index in range(samples):
        t = index / (samples - 1)
        y = (t - 0.5) * height
        # Arc inward at the top and bottom to remain inside the 2 m cell envelope.
        axis_x = centre_x - 0.25 * (2 * abs(t - 0.5)) ** 1.7
        angle = math.tau * turns * t
        strand_a.append([axis_x + helix_radius * math.cos(angle), y, centre_z + helix_radius * math.sin(angle)])
        strand_b.append([axis_x - helix_radius * math.cos(angle), y, centre_z - helix_radius * math.sin(angle)])
    geometry_a = tube(strand_a, 0.012, radial_segments=10, caps=True)
    geometry_b = tube(strand_b, 0.012, radial_segments=10, caps=True)

    rungs: list[Geometry] = []
    for rung_index in range(24):
        sample = round(rung_index * (samples - 1) / 23)
        rungs.append(tube((strand_a[sample], strand_b[sample]), 0.007, radial_segments=8, caps=True))
    return geometry_a, geometry_b, concatenate(rungs)


OPENORGANELLE_DOI = "https://doi.org/10.25378/janelia.13108343"
OPENORGANELLE_DATASET = "jrc_hela-2"
MEASURED_SPACING_ZYX_NM = np.array((83.84, 64.0, 64.0), dtype=float)


def load_measured_volume(cache_dir: Path, label: str) -> tuple[np.ndarray, dict]:
    path = cache_dir / f"{label}_s4.npz"
    if not path.exists():
        raise FileNotFoundError(
            f"Missing measured label cache {path}. Run "
            "`python tools/fetch_openorganelle_cell.py` first."
        )
    with np.load(path, allow_pickle=False) as archive:
        volume = archive["volume"]
        provenance = json.loads(str(archive["provenance"].item()))
    return volume, provenance


def largest_component(mask: np.ndarray) -> np.ndarray:
    labels, count = ndimage.label(mask, structure=ndimage.generate_binary_structure(3, 1))
    if count == 0:
        raise ValueError("Measured segmentation contains no connected component")
    sizes = np.bincount(labels.reshape(-1))
    sizes[0] = 0
    return labels == int(np.argmax(sizes))


def measured_surface(
    mask: np.ndarray,
    *,
    target_faces: int,
    spacing_zyx_nm: Sequence[float] = MEASURED_SPACING_ZYX_NM,
    smooth_sigma: float = 0.62,
) -> Geometry:
    """Extract, smooth, and Web-optimize a measured binary label surface."""
    occupied = np.argwhere(mask)
    if not len(occupied):
        raise ValueError("Cannot mesh an empty measured label")
    lower = np.maximum(occupied.min(axis=0) - 3, 0)
    upper = np.minimum(occupied.max(axis=0) + 4, np.asarray(mask.shape))
    slices = tuple(slice(int(a), int(b)) for a, b in zip(lower, upper))
    field = ndimage.gaussian_filter(mask[slices].astype(np.float32), sigma=smooth_sigma)
    field = np.pad(field, 1, mode="constant")
    vertices_zyx, faces, _normals, _values = measure.marching_cubes(
        field,
        level=0.45,
        spacing=tuple(float(value) for value in spacing_zyx_nm),
        allow_degenerate=False,
    )
    vertices_zyx += (lower - 1) * np.asarray(spacing_zyx_nm)
    vertices_xyz = vertices_zyx[:, (2, 1, 0)]
    mesh = trimesh.Trimesh(vertices=vertices_xyz, faces=faces, process=True, validate=True)
    if len(mesh.faces) > target_faces:
        mesh = mesh.simplify_quadric_decimation(face_count=target_faces, aggression=9)
    return clean_geometry(np.asarray(mesh.vertices), np.asarray(mesh.faces))


def measured_principal_transform(
    geometries: Sequence[Geometry],
    reference: Geometry,
    *,
    target_length: float,
    translation: Sequence[float],
    rotation: np.ndarray,
) -> list[Geometry]:
    """Apply one PCA frame/scale to related measured surfaces."""
    reference_vertices = np.asarray(reference.vertices, dtype=float)
    centre = (reference_vertices.min(axis=0) + reference_vertices.max(axis=0)) * 0.5
    covariance = np.cov(reference_vertices - centre, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    basis = eigenvectors[:, np.argsort(eigenvalues)[::-1]]
    # Keep a right-handed local frame so outward winding survives the transform.
    if np.linalg.det(basis) < 0:
        basis[:, -1] *= -1
    local_reference = (reference_vertices - centre) @ basis
    scale = target_length / float(np.ptp(local_reference[:, 0]))
    transformed: list[Geometry] = []
    for geometry in geometries:
        local = (np.asarray(geometry.vertices, dtype=float) - centre) @ basis
        vertices = (local * scale) @ rotation.T + np.asarray(translation, dtype=float)
        transformed.append(clean_geometry(vertices, geometry.faces))
    return transformed


def measured_nucleus_transform(
    geometries: Sequence[Geometry],
    reference: Geometry,
    *,
    translation: Sequence[float] = (-0.35, 0.02, -0.14),
    target_extents: Sequence[float] = (0.60, 0.54, 0.52),
) -> list[Geometry]:
    reference_vertices = np.asarray(reference.vertices, dtype=float)
    minimum = reference_vertices.min(axis=0)
    maximum = reference_vertices.max(axis=0)
    centre = (minimum + maximum) * 0.5
    scale = np.asarray(target_extents, dtype=float) / np.maximum(maximum - minimum, 1e-9)
    transformed: list[Geometry] = []
    for geometry in geometries:
        vertices = (np.asarray(geometry.vertices, dtype=float) - centre) * scale + np.asarray(translation)
        transformed.append(clean_geometry(vertices, geometry.faces))
    return transformed


def material(
    name: str,
    colour: Sequence[float],
    metallic: float,
    roughness: float,
    *,
    transmission: float = 0.0,
    ior: float = 1.45,
    thickness: float = 0.0,
    attenuation: Sequence[float] = (1.0, 1.0, 1.0),
    double_sided: bool = False,
) -> dict:
    result: dict = {
        "name": name,
        "pbrMetallicRoughness": {
            "baseColorFactor": list(colour),
            "metallicFactor": metallic,
            "roughnessFactor": roughness,
        },
        "emissiveFactor": [0.0, 0.0, 0.0],
        "doubleSided": double_sided,
    }
    if colour[3] < 1.0:
        result["alphaMode"] = "BLEND"
    extensions: dict = {
        "KHR_materials_ior": {"ior": ior},
        "KHR_materials_specular": {
            "specularFactor": 0.84,
            "specularColorFactor": [0.92, 0.94, 0.93],
        },
    }
    if transmission > 0:
        extensions["KHR_materials_transmission"] = {"transmissionFactor": transmission}
        if thickness > 0:
            extensions["KHR_materials_volume"] = {
                "thicknessFactor": thickness,
                "attenuationDistance": 2.6,
                "attenuationColor": list(attenuation),
            }
    result["extensions"] = extensions
    return result


class GLBBuilder:
    def __init__(self) -> None:
        self.binary = bytearray()
        self.buffer_views: list[dict] = []
        self.accessors: list[dict] = []
        self.meshes: list[dict] = []
        self.nodes: list[dict] = [
            {
                "name": "AYARI_Cell_Cutaway",
                "children": [1],
                "extras": {
                    "coordinateSystem": "right-handed, Y-up, metres",
                    "origin": "cell centre",
                    "style": "medical-luxury cutaway",
                    "backgroundAndLightsBaked": False,
                    "geometryMethod": "reference-driven hybrid",
                    "measuredDataset": OPENORGANELLE_DATASET,
                    "measuredSourceDoi": OPENORGANELLE_DOI,
                    "measuredSourceLicense": "CC BY 4.0",
                },
            },
            {"name": "GEOMETRY", "children": []},
        ]
        self.report_nodes: list[dict] = []

    def _pad_binary(self) -> None:
        while len(self.binary) % 4:
            self.binary.append(0)

    def _accessor(self, data: np.ndarray, target: int, gltf_type: str, component_type: int) -> int:
        self._pad_binary()
        offset = len(self.binary)
        if component_type == FLOAT:
            array = np.asarray(data, dtype="<f4")
        elif component_type == UNSIGNED_INT:
            array = np.asarray(data, dtype="<u4")
        else:
            raise ValueError(component_type)
        raw = array.tobytes()
        self.binary.extend(raw)
        view_index = len(self.buffer_views)
        self.buffer_views.append(
            {"buffer": 0, "byteOffset": offset, "byteLength": len(raw), "target": target}
        )
        accessor: dict = {
            "bufferView": view_index,
            "componentType": component_type,
            "count": int(len(array)),
            "type": gltf_type,
        }
        if gltf_type == "VEC3":
            accessor["min"] = np.min(array, axis=0).astype(float).tolist()
            accessor["max"] = np.max(array, axis=0).astype(float).tolist()
        else:
            accessor["min"] = [int(np.min(array))]
            accessor["max"] = [int(np.max(array))]
        accessor_index = len(self.accessors)
        self.accessors.append(accessor)
        return accessor_index

    def primitive(self, geometry: Geometry, material_index: int) -> dict:
        cleaned = geometry if geometry.normals is not None else clean_geometry(geometry.vertices, geometry.faces)
        position = self._accessor(cleaned.vertices, ARRAY_BUFFER, "VEC3", FLOAT)
        normal = self._accessor(cleaned.normals, ARRAY_BUFFER, "VEC3", FLOAT)
        indices = self._accessor(cleaned.faces.reshape(-1), ELEMENT_ARRAY_BUFFER, "SCALAR", UNSIGNED_INT)
        return {
            "attributes": {"POSITION": position, "NORMAL": normal},
            "indices": indices,
            "material": material_index,
            "mode": 4,
        }

    def add_node(
        self,
        name: str,
        primitive_specs: Sequence[tuple[Geometry, int]],
        extras: dict,
    ) -> int:
        mesh_index = len(self.meshes)
        primitives = [self.primitive(geometry, material_index) for geometry, material_index in primitive_specs]
        mesh_extras = dict(extras)
        mesh_extras["primitiveCount"] = len(primitives)
        self.meshes.append({"name": name.replace("GEO_", "MESH_", 1), "primitives": primitives, "extras": mesh_extras})
        node_index = len(self.nodes)
        self.nodes.append({"name": name, "mesh": mesh_index, "extras": extras})
        self.nodes[1]["children"].append(node_index)

        all_vertices = np.vstack([geometry.vertices for geometry, _ in primitive_specs])
        self.report_nodes.append(
            {
                "name": name,
                "labelJa": extras["labelJa"],
                "markerKey": extras["markerKey"],
                "measured": extras.get("measured", False),
                "geometryProvenance": extras.get("geometryProvenance", "unspecified"),
                "vertices": int(sum(len(geometry.vertices) for geometry, _ in primitive_specs)),
                "triangles": int(sum(len(geometry.faces) for geometry, _ in primitive_specs)),
                "boundsMetres": {
                    "min": all_vertices.min(axis=0).astype(float).tolist(),
                    "max": all_vertices.max(axis=0).astype(float).tolist(),
                },
            }
        )
        return node_index

    def write(self, output: Path, materials: list[dict]) -> None:
        gltf = {
            "asset": {
                "version": "2.0",
                "generator": "AYARI Clinic OpenOrganelle cell builder 2.0",
                "copyright": "Copyright AYARI Clinic",
                "extras": {
                    "units": "metres",
                    "upAxis": "Y",
                    "origin": "cell centre",
                    "intendedUse": "Web visualisation and viewer-side biomarker highlighting",
                    "medicalUse": "Reference-driven educational visual; not a diagnostic or surgical-planning device.",
                    "sourceDataset": OPENORGANELLE_DATASET,
                    "sourceDoi": OPENORGANELLE_DOI,
                    "sourceLicense": "CC BY 4.0",
                    "sourceVoxelResolutionNm": [4.0, 4.0, 5.24],
                    "meshSamplingResolutionNmXYZ": [64.0, 64.0, 83.84],
                    "provenancePolicy": "Each node declares measured or educational-overlay geometry in extras.",
                },
            },
            "extensionsUsed": [
                "KHR_materials_ior",
                "KHR_materials_specular",
                "KHR_materials_transmission",
                "KHR_materials_volume",
            ],
            "scene": 0,
            "scenes": [{"name": "AYARI Cell", "nodes": [0]}],
            "nodes": self.nodes,
            "meshes": self.meshes,
            "materials": materials,
            "buffers": [{"byteLength": len(self.binary)}],
            "bufferViews": self.buffer_views,
            "accessors": self.accessors,
        }
        json_bytes = json.dumps(gltf, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        while len(json_bytes) % 4:
            json_bytes += b" "
        self._pad_binary()
        binary_bytes = bytes(self.binary)
        total_length = 12 + 8 + len(json_bytes) + 8 + len(binary_bytes)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("wb") as handle:
            handle.write(struct.pack("<4sII", b"glTF", 2, total_length))
            handle.write(struct.pack("<I4s", len(json_bytes), b"JSON"))
            handle.write(json_bytes)
            handle.write(struct.pack("<I4s", len(binary_bytes), b"BIN\x00"))
            handle.write(binary_bytes)


def draw_preview(output: Path) -> None:
    scale = 2
    width, height = 1400 * scale, 900 * scale
    image = Image.new("RGB", (width, height), "#070a0d")
    glow = Image.new("RGBA", image.size, (0, 0, 0, 0))
    g = ImageDraw.Draw(glow)
    cx, cy = 620 * scale, 450 * scale
    rx, ry = 430 * scale, 395 * scale
    gold = (201, 169, 110, 210)
    pale = (211, 221, 222, 145)

    # Cutaway membrane: three-quarter ellipse and two luminous cut lips.
    bbox = (cx - rx, cy - ry, cx + rx, cy + ry)
    g.arc(bbox, 58, 302, fill=gold, width=8 * scale)
    g.arc((cx - rx + 10 * scale, cy - ry + 10 * scale, cx + rx - 10 * scale, cy + ry - 10 * scale), 58, 302, fill=pale, width=5 * scale)
    for angle in (58, 302):
        rad = math.radians(angle)
        x = cx + rx * math.cos(rad)
        y = cy + ry * math.sin(rad)
        g.line((cx, cy, x, y), fill=(205, 215, 215, 92), width=4 * scale)

    # Nucleus and chromosomes.
    nx, ny = cx - 150 * scale, cy + 5 * scale
    g.ellipse((nx - 140 * scale, ny - 132 * scale, nx + 140 * scale, ny + 132 * scale), fill=(173, 192, 195, 38), outline=(224, 230, 228, 165), width=4 * scale)
    chromosome_centres = ((nx - 55 * scale, ny - 42 * scale), (nx + 45 * scale, ny - 38 * scale), (nx - 45 * scale, ny + 48 * scale), (nx + 58 * scale, ny + 45 * scale))
    for x, y in chromosome_centres:
        for dx in (-1, 1):
            g.line((x, y, x + dx * 27 * scale, y - 42 * scale), fill=(233, 228, 211, 225), width=9 * scale)
            g.line((x, y, x + dx * 27 * scale, y + 42 * scale), fill=(233, 228, 211, 225), width=9 * scale)
            g.ellipse((x + dx * 27 * scale - 6 * scale, y - 48 * scale, x + dx * 27 * scale + 6 * scale, y - 36 * scale), fill=gold)
            g.ellipse((x + dx * 27 * scale - 6 * scale, y + 36 * scale, x + dx * 27 * scale + 6 * scale, y + 48 * scale), fill=gold)

    # Mitochondria and cristae around the nucleus.
    mito = (
        (360, 235, -12), (675, 190, 18), (245, 450, 78), (330, 675, -20),
        (610, 730, 7), (790, 590, 55), (740, 360, -42),
    )
    for x, y, angle in mito:
        overlay = Image.new("RGBA", image.size, (0, 0, 0, 0))
        d = ImageDraw.Draw(overlay)
        mx, my = x * scale, y * scale
        d.ellipse((mx - 105 * scale, my - 42 * scale, mx + 105 * scale, my + 42 * scale), fill=(194, 207, 207, 70), outline=(223, 228, 222, 180), width=4 * scale)
        for offset in (-62, -30, 0, 30, 62):
            d.arc((mx + (offset - 18) * scale, my - 30 * scale, mx + (offset + 18) * scale, my + 30 * scale), 80, 280, fill=(202, 172, 118, 210), width=4 * scale)
        overlay = overlay.rotate(angle, center=(mx, my), resample=Image.Resampling.BICUBIC)
        glow.alpha_composite(overlay)

    # DNA at the right-front edge.
    points_a: list[tuple[float, float]] = []
    points_b: list[tuple[float, float]] = []
    for i in range(120):
        t = i / 119
        y = (190 + 520 * t) * scale
        axis_x = (970 - 65 * (2 * abs(t - 0.5)) ** 1.7) * scale
        angle = math.tau * 2.2 * t
        dx = math.cos(angle) * 37 * scale
        points_a.append((axis_x + dx, y))
        points_b.append((axis_x - dx, y))
    for i in range(0, 120, 5):
        g.line((*points_a[i], *points_b[i]), fill=(205, 184, 145, 150), width=3 * scale)
    g.line(points_a, fill=(228, 232, 226, 235), width=7 * scale, joint="curve")
    g.line(points_b, fill=(190, 205, 207, 235), width=7 * scale, joint="curve")

    blurred = glow.filter(ImageFilter.GaussianBlur(14 * scale))
    image = Image.alpha_composite(image.convert("RGBA"), blurred)
    image = Image.alpha_composite(image, glow)
    draw = ImageDraw.Draw(image)
    try:
        title_font = ImageFont.truetype("DejaVuSans.ttf", 34 * scale)
        detail_font = ImageFont.truetype("DejaVuSans.ttf", 17 * scale)
    except OSError:
        title_font = ImageFont.load_default()
        detail_font = ImageFont.load_default()
    draw.text((70 * scale, 55 * scale), "AYARI CELL — CUTAWAY GLB", fill=(236, 231, 219, 245), font=title_font)
    draw.text((72 * scale, 105 * scale), "102k-class triangles · Y-up · metres · named selectable submeshes", fill=(177, 184, 184, 225), font=detail_font)
    draw.text((1110 * scale, 810 * scale), "QA PREVIEW — LIGHTING NOT BAKED", fill=(143, 149, 148, 210), font=detail_font, anchor="ra")
    output.parent.mkdir(parents=True, exist_ok=True)
    image.convert("RGB").resize((1400, 900), Image.Resampling.LANCZOS).save(output, quality=94)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("assets/cell/ayari_cell_cutaway.glb"))
    parser.add_argument("--report", type=Path, default=Path("assets/cell/model_report.json"))
    parser.add_argument("--preview", type=Path, default=Path("assets/cell/preview.png"))
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path(".cache/openorganelle/jrc_hela-2"),
        help="Directory produced by fetch_openorganelle_cell.py",
    )
    args = parser.parse_args()

    materials = [
        material("MAT_cell_glass", (0.74, 0.82, 0.83, 0.22), 0.04, 0.16, transmission=0.72, ior=1.37, thickness=0.035, attenuation=(0.77, 0.88, 0.89)),
        material("MAT_nucleus_glass", (0.76, 0.80, 0.78, 0.34), 0.02, 0.20, transmission=0.58, ior=1.39, thickness=0.055, attenuation=(0.88, 0.91, 0.86)),
        material("MAT_chromosome_pearl", (0.88, 0.86, 0.79, 0.96), 0.10, 0.28, ior=1.46),
        material("MAT_telomere_cap", (0.89, 0.76, 0.51, 1.0), 0.24, 0.20, ior=1.48),
        material("MAT_mitochondria_glass", (0.70, 0.77, 0.77, 0.50), 0.04, 0.22, transmission=0.42, ior=1.40, thickness=0.035, attenuation=(0.84, 0.89, 0.87)),
        material("MAT_cristae", (0.82, 0.75, 0.60, 0.94), 0.13, 0.26, ior=1.47),
        material("MAT_dna_strand_a", (0.88, 0.89, 0.85, 1.0), 0.13, 0.22, ior=1.47),
        material("MAT_dna_strand_b", (0.71, 0.79, 0.80, 1.0), 0.10, 0.24, ior=1.44),
        material("MAT_dna_base_pairs", (0.82, 0.70, 0.49, 0.96), 0.18, 0.24, ior=1.48),
        material("MAT_measured_chromatin", (0.67, 0.60, 0.55, 0.30), 0.06, 0.34, transmission=0.35, ior=1.42, thickness=0.018, attenuation=(0.77, 0.72, 0.68)),
    ]
    builder = GLBBuilder()

    cell = cutaway_shell(around=144, vertical=76)
    builder.add_node(
        "GEO_cell",
        ((cell, 0),),
        {
            "labelJa": "細胞",
            "markerKey": "cell",
            "category": "cell membrane",
            "cutaway": True,
            "geometryProvenance": "modelled display envelope",
            "measured": False,
            "referenceDataset": OPENORGANELLE_DATASET,
            "representationNote": "Envelope is a display-space cutaway; measured morphology is used for the nucleus and mitochondria.",
        },
    )

    nucleus_volume, nucleus_provenance = load_measured_volume(args.cache_dir, "nucleus_seg")
    chromatin_volume, chromatin_provenance = load_measured_volume(args.cache_dir, "chrom_seg")
    nucleus_mask = largest_component(nucleus_volume > 0)
    chromatin_mask = (chromatin_volume > 0) & ndimage.binary_dilation(nucleus_mask, iterations=2)
    nucleus_measured = measured_surface(nucleus_mask, target_faces=10_000, smooth_sigma=0.72)
    chromatin_measured = measured_surface(chromatin_mask, target_faces=5_500, smooth_sigma=0.58)
    nucleus, chromatin = measured_nucleus_transform(
        (nucleus_measured, chromatin_measured), nucleus_measured
    )
    builder.add_node(
        "GEO_nucleus",
        ((nucleus, 1), (chromatin, 9)),
        {
            "labelJa": "核",
            "markerKey": "cell",
            "category": "nucleus",
            "translucent": True,
            "measured": True,
            "geometryProvenance": "FIB-SEM segmentation surface",
            "sourceDataset": OPENORGANELLE_DATASET,
            "sourceDoi": OPENORGANELLE_DOI,
            "sourceLabels": [nucleus_provenance["label"], chromatin_provenance["label"]],
            "sourceScale": "s4",
            "sourceSamplingNmXYZ": [64.0, 64.0, 83.84],
            "displayTransform": "non-uniformly normalized to requested 0.60 m nucleus extent",
            "chromatinPrimitiveIncluded": True,
        },
    )

    chromosome_configs = (
        ((-0.43, 0.105, -0.105), 0.075, 0.12, 0.2),
        ((-0.275, 0.105, -0.155), 0.070, -0.18, 1.3),
        ((-0.43, -0.105, -0.175), 0.068, -0.08, 2.1),
        ((-0.275, -0.095, -0.095), 0.072, 0.18, 3.0),
    )
    telomere_index = 1
    for chromosome_index, (centre, scale, rz, phase) in enumerate(chromosome_configs, start=1):
        rotation = rotation_matrix(0.08 * math.sin(phase), 0.08 * math.cos(phase), rz)
        body, endpoints = chromosome(centre, scale, rotation, phase)
        builder.add_node(
            f"GEO_chromosome_{chromosome_index:02d}",
            ((body, 2),),
            {
                "labelJa": "染色体",
                "markerKey": "telo",
                "category": "chromosome",
                "chromosomeIndex": chromosome_index,
                "measured": False,
                "geometryProvenance": "educational overlay",
                "biologyNote": "Condensed X-shape is a teaching symbol and is not measured in this interphase HeLa dataset.",
            },
        )
        for arm_index, endpoint in enumerate(endpoints, start=1):
            cap = uv_sphere(endpoint, (scale * 0.175, scale * 0.150, scale * 0.155), around=20, vertical=12)
            builder.add_node(
                f"GEO_telomere_{telomere_index:02d}",
                ((cap, 3),),
                {
                    "labelJa": "テロメア",
                    "markerKey": "telo",
                    "category": "telomere",
                    "chromosomeIndex": chromosome_index,
                    "armTip": arm_index,
                    "measured": False,
                    "geometryProvenance": "educational overlay",
                },
            )
            telomere_index += 1

    mito_volume, mito_provenance = load_measured_volume(args.cache_dir, "mito_seg")
    mito_membrane_volume, mito_membrane_provenance = load_measured_volume(args.cache_dir, "mito-mem_seg")
    # Seven sizeable, well-resolved source instances.  Each outer body and its
    # membrane/cristae label receives exactly the same display transform.
    mito_configs = (
        (109, (-0.59, 0.57, -0.12), 0.72, (0.02, -0.10, -0.12)),
        (351, (0.18, 0.64, -0.16), 0.70, (0.08, 0.15, 0.28)),
        (253, (-0.72, 0.12, -0.10), 0.66, (-0.12, 0.10, 1.24)),
        (99, (-0.57, -0.52, -0.14), 0.70, (0.06, -0.12, -0.26)),
        (385, (0.03, -0.68, -0.10), 0.76, (-0.05, 0.08, 0.12)),
        (277, (0.48, -0.38, -0.11), 0.68, (0.10, -0.08, 0.77)),
        (169, (0.35, 0.22, -0.26), 0.67, (-0.08, 0.12, -0.72)),
    )
    for index, (source_instance, centre, length, angles) in enumerate(mito_configs, start=1):
        body_measured = measured_surface(
            mito_volume == source_instance, target_faces=3_500, smooth_sigma=0.54
        )
        membrane_measured = measured_surface(
            mito_membrane_volume == source_instance, target_faces=1_800, smooth_sigma=0.44
        )
        body, cristae = measured_principal_transform(
            (body_measured, membrane_measured),
            body_measured,
            target_length=length,
            translation=centre,
            rotation=rotation_matrix(*angles),
        )
        builder.add_node(
            f"GEO_mitochondria_{index:02d}",
            ((body, 4), (cristae, 5)),
            {
                "labelJa": "ミトコンドリア",
                "markerKey": "mito",
                "category": "mitochondrion",
                "instance": index,
                "cristaeIncluded": True,
                "measured": True,
                "geometryProvenance": "FIB-SEM segmentation surfaces",
                "sourceDataset": OPENORGANELLE_DATASET,
                "sourceDoi": OPENORGANELLE_DOI,
                "sourceOuterLabel": mito_provenance["label"],
                "sourceMembraneLabel": mito_membrane_provenance["label"],
                "sourceInstanceId": source_instance,
                "sourceScale": "s4",
                "sourceSamplingNmXYZ": [64.0, 64.0, 83.84],
                "displayTransform": "PCA-aligned, uniformly scaled, and compositionally repositioned",
            },
        )

    strand_a, strand_b, rungs = dna_helix()
    builder.add_node(
        "GEO_dna",
        ((strand_a, 6), (strand_b, 7), (rungs, 8)),
        {
            "labelJa": "DNA",
            "markerKey": "dna",
            "category": "DNA double helix",
            "turns": 2.2,
            "basePairRungs": 24,
            "placement": "right-front of nucleus",
            "measured": False,
            "geometryProvenance": "educational molecular-scale overlay",
            "scaleNote": "DNA is intentionally enlarged and is not at the same physical scale as the cell surface.",
        },
    )

    builder.write(args.output, materials)
    from render_cell_preview import render as render_cell_preview

    render_cell_preview(
        args.output,
        args.preview,
        triangles=int(sum(node["triangles"] for node in builder.report_nodes)),
    )

    all_min = np.min([node["boundsMetres"]["min"] for node in builder.report_nodes], axis=0)
    all_max = np.max([node["boundsMetres"]["max"] for node in builder.report_nodes], axis=0)
    report = {
        "asset": Path(os.path.relpath(args.output, Path.cwd())).as_posix(),
        "format": "glTF 2.0 Binary",
        "coordinateSystem": {"units": "metres", "upAxis": "Y", "origin": "cell centre", "transformsApplied": True},
        "triangles": int(sum(node["triangles"] for node in builder.report_nodes)),
        "vertices": int(sum(node["vertices"] for node in builder.report_nodes)),
        "meshes": len(builder.meshes),
        "nodes": len(builder.nodes),
        "materials": [item["name"] for item in materials],
        "boundsMetres": {"min": all_min.astype(float).tolist(), "max": all_max.astype(float).tolist()},
        "fileSizeBytes": args.output.stat().st_size,
        "components": builder.report_nodes,
        "counts": {"mitochondria": 7, "chromosomes": 4, "telomeres": 16, "dnaDoubleHelices": 1},
        "viewerKeys": {
            "cell": ["GEO_cell", "GEO_nucleus"],
            "mito": [f"GEO_mitochondria_{index:02d}" for index in range(1, 8)],
            "telo": [f"GEO_chromosome_{index:02d}" for index in range(1, 5)] + [f"GEO_telomere_{index:02d}" for index in range(1, 17)],
            "dna": ["GEO_dna"],
        },
        "sceneContent": {"background": False, "lights": False, "cameras": False, "animations": False, "skins": False},
        "sourceData": {
            "dataset": OPENORGANELLE_DATASET,
            "sample": "wild-type interphase HeLa cell (ATCC CCL-2)",
            "modality": "isotropic FIB-SEM",
            "doi": OPENORGANELLE_DOI,
            "license": "CC BY 4.0",
            "nativeVoxelResolutionNmXYZ": [4.0, 4.0, 5.24],
            "meshSamplingResolutionNmXYZ": [64.0, 64.0, 83.84],
            "measuredNodes": ["GEO_nucleus"] + [f"GEO_mitochondria_{index:02d}" for index in range(1, 8)],
            "modelledNodes": ["GEO_cell"] + [f"GEO_chromosome_{index:02d}" for index in range(1, 5)] + [f"GEO_telomere_{index:02d}" for index in range(1, 17)] + ["GEO_dna"],
            "displayNormalization": "Measured surfaces are normalized and compositionally repositioned to satisfy the requested 2 m display-space layout.",
        },
        "scientificScope": "Measured cell-organelle morphology with clearly identified educational overlays; not a single-scale literal reconstruction.",
        "medicalUse": "Reference-driven educational visualisation; not validated for diagnosis or surgical planning.",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"asset": str(args.output), "triangles": report["triangles"], "vertices": report["vertices"], "bytes": report["fileSizeBytes"], "meshes": report["meshes"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
