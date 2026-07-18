#!/usr/bin/env python3
"""Build the AYARI mobile anatomy asset from BodyParts3D release 3.0.

The upstream model uses millimetres with Z as the superior axis. This builder
converts it to glTF's metre-scale Y-up convention and moves the origin to the
centre of the feet. Geometry is simplified to a configurable triangle budget.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import struct
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import trimesh
from pygltflib import (
    ARRAY_BUFFER,
    ELEMENT_ARRAY_BUFFER,
    FLOAT,
    GLTF2,
    UNSIGNED_INT,
    Accessor,
    Asset,
    Buffer,
    BufferView,
    Material,
    Mesh,
    Node,
    PbrMetallicRoughness,
    Primitive,
    Scene,
)


UPSTREAM_URL = "https://github.com/Kevin-Mattheus-Moerman/BodyParts3D.git"
UPSTREAM_COMMIT = "f0eeb6e843380cfe6b83797cf8c3e1af74de5e61"
UPSTREAM_ATTRIBUTION = (
    "BodyParts3D, (c) The Database Center for Life Science licensed under "
    "CC Attribution-Share Alike 2.1 Japan"
)

STL_ROOT = Path("assets/BodyParts3D_data/stl")
PARTS_LIST = Path("assets/BodyParts3D_data/parts_list_e.txt")
COMPOSITE_LIST = Path("assets/BodyParts3D_data/composite_parts.txt")
LICENSE_FILE = Path("assets/BodyParts3D_data/LICENSE_content")


@dataclass(frozen=True)
class ComponentSpec:
    key: str
    label_ja: str
    fma_ids: tuple[str, ...]
    target_triangles: int
    material: str
    marker_aliases: tuple[str, ...] = ()


FIXED_COMPONENTS = (
    ComponentSpec("skin", "皮膚", ("FMA7163",), 3_500, "skin", ("皮膚",)),
    ComponentSpec("heart", "心臓", ("FMA7274",), 1_400, "heart", ("心臓", "炎症")),
    ComponentSpec(
        "right_lung",
        "右肺",
        ("FMA7333", "FMA7383", "FMA7337"),
        1_250,
        "lung",
        ("右肺",),
    ),
    ComponentSpec(
        "left_lung",
        "左肺",
        ("FMA7370", "FMA7371"),
        1_250,
        "lung",
        ("左肺",),
    ),
    ComponentSpec("liver", "肝臓", ("FMA7197",), 1_500, "liver", ("肝臓",)),
    ComponentSpec("stomach", "胃", ("FMA7148",), 800, "digestive", ("胃",)),
    ComponentSpec("pancreas", "膵臓", ("FMA7198nsn",), 500, "pancreas", ("膵臓",)),
    ComponentSpec("spleen", "脾臓", ("FMA7196",), 500, "spleen", ("脾臓",)),
    ComponentSpec("right_kidney", "右腎臓", ("FMA7204",), 500, "kidney", ("右腎臓",)),
    ComponentSpec("left_kidney", "左腎臓", ("FMA7205",), 500, "kidney", ("左腎臓",)),
    ComponentSpec("right_adrenal", "右副腎", ("FMA15629",), 150, "adrenal", ("右副腎",)),
    ComponentSpec("left_adrenal", "左副腎", ("FMA15630",), 150, "adrenal", ("左副腎",)),
    ComponentSpec(
        "small_intestine",
        "小腸",
        ("FMA7206", "FMA7207", "FMA7208"),
        1_800,
        "intestine",
        ("小腸",),
    ),
    ComponentSpec(
        "large_intestine",
        "大腸",
        ("FMA14543nsn",),
        1_300,
        "colon",
        ("大腸",),
    ),
    ComponentSpec("bladder", "膀胱", ("FMA15900",), 500, "bladder", ("膀胱",)),
    ComponentSpec("prostate", "前立腺", ("FMA9600",), 300, "prostate", ("前立腺",)),
    ComponentSpec("trachea", "気管", ("FMA7394",), 500, "airway", ("気管",)),
    ComponentSpec("esophagus", "食道", ("FMA7131",), 500, "digestive", ("食道",)),
)


MATERIALS = {
    "skin": dict(base=(0.62, 0.69, 0.73, 0.16), metallic=0.0, roughness=0.48, alpha="BLEND", double=True),
    "bone": dict(base=(0.88, 0.84, 0.72, 1.0), metallic=0.05, roughness=0.72),
    "brain": dict(base=(0.77, 0.56, 0.60, 1.0), metallic=0.0, roughness=0.62),
    "heart": dict(base=(0.50, 0.08, 0.10, 1.0), metallic=0.02, roughness=0.48),
    "lung": dict(base=(0.66, 0.36, 0.43, 1.0), metallic=0.0, roughness=0.66),
    "liver": dict(base=(0.39, 0.17, 0.12, 1.0), metallic=0.0, roughness=0.65),
    "digestive": dict(base=(0.60, 0.33, 0.26, 1.0), metallic=0.0, roughness=0.68),
    "pancreas": dict(base=(0.77, 0.52, 0.28, 1.0), metallic=0.0, roughness=0.62),
    "spleen": dict(base=(0.42, 0.13, 0.21, 1.0), metallic=0.0, roughness=0.65),
    "kidney": dict(base=(0.45, 0.15, 0.13, 1.0), metallic=0.0, roughness=0.68),
    "adrenal": dict(base=(0.78, 0.60, 0.25, 1.0), metallic=0.0, roughness=0.60),
    "intestine": dict(base=(0.67, 0.42, 0.32, 1.0), metallic=0.0, roughness=0.72),
    "colon": dict(base=(0.58, 0.31, 0.25, 1.0), metallic=0.0, roughness=0.72),
    "bladder": dict(base=(0.72, 0.64, 0.42, 1.0), metallic=0.0, roughness=0.58),
    "prostate": dict(base=(0.66, 0.47, 0.34, 1.0), metallic=0.0, roughness=0.62),
    "airway": dict(base=(0.52, 0.38, 0.34, 1.0), metallic=0.0, roughness=0.66),
}


def run(*args: str, cwd: Path | None = None) -> str:
    result = subprocess.run(
        args,
        cwd=cwd,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return result.stdout


def ensure_checkout(cache: Path) -> Path:
    if not (cache / ".git").exists():
        cache.parent.mkdir(parents=True, exist_ok=True)
        run(
            "git",
            "clone",
            "--filter=blob:none",
            "--no-checkout",
            "--depth",
            "1",
            UPSTREAM_URL,
            str(cache),
        )

    current = run("git", "rev-parse", "HEAD", cwd=cache).strip()
    if current != UPSTREAM_COMMIT:
        raise RuntimeError(f"Unexpected BodyParts3D revision: {current}")

    run(
        "git",
        "checkout",
        "HEAD",
        "--",
        str(PARTS_LIST),
        str(COMPOSITE_LIST),
        str(LICENSE_FILE),
        cwd=cache,
    )
    return cache


def read_name_map(source: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for line in (source / PARTS_LIST).read_text(encoding="utf-8").splitlines():
        if "\t" not in line:
            continue
        fma_id, name = line.split("\t", 1)
        mapping[fma_id] = name.strip()
    return mapping


def available_fma_ids(source: Path) -> set[str]:
    paths = run("git", "ls-tree", "-r", "--name-only", "HEAD", str(STL_ROOT), cwd=source)
    return {Path(path).stem for path in paths.splitlines() if path.endswith(".stl")}


def select_skeleton_ids(names: dict[str, str], available: set[str]) -> list[str]:
    anatomy_terms = re.compile(
        r"(bone|vertebra|rib|sternum|sacrum|coccyx|clavicle|scapula|humerus|"
        r"radius|ulna|femur|patella|tibia|fibula|phalanx|metacarpal|metatarsal|"
        r"calcaneus|talus|navicular|cuneiform|cuboid|trapezium|trapezoid|"
        r"capitate|hamate|lunate|pisiform|scaphoid|triquetral|maxilla|mandible|"
        r"ethmoid|sphenoid|vomer|hyoid|temporal|parietal|frontal|occipital|"
        r"zygomatic|lacrimal)",
        re.IGNORECASE,
    )
    unpaired = {
        "FMA12519",  # atlas
        "FMA12520",  # axis
        "FMA16202",  # sacrum
        "FMA7485",   # sternum
        "FMA52736",  # sphenoid
        "FMA52740",  # ethmoid
        "FMA52748",  # mandible
        "FMA52749",  # hyoid
        "FMA9710",   # vomer
    }
    selected: set[str] = set()
    for fma_id, name in names.items():
        if fma_id not in available or not anatomy_terms.search(name):
            continue
        lower = name.lower()
        is_lateral = lower.startswith("right ") or lower.startswith("left ")
        is_specific_vertebra = bool(
            re.search(r"(first|second|third|fourth|fifth|sixth|seventh|eighth|ninth|tenth|eleventh|twelfth) .*vertebra$", lower)
        )
        if is_lateral or is_specific_vertebra or fma_id in unpaired:
            selected.add(fma_id)
    return sorted(selected)


def select_brain_ids(source: Path, names: dict[str, str], available: set[str]) -> list[str]:
    selected = {"FMA62004", "FMA67943", "FMA67944"} & available
    surface_terms = re.compile(r"(gyrus|lobule|cerebellum|pons|medulla oblongata)$", re.IGNORECASE)
    for line in (source / COMPOSITE_LIST).read_text(encoding="utf-8").splitlines():
        columns = line.split("\t")
        if len(columns) < 4 or columns[0] != "FMA50801":
            continue
        child_id, child_name = columns[2], columns[3]
        if child_id in available and surface_terms.search(child_name):
            selected.add(child_id)
    return sorted(selected)


def checkout_source_meshes(source: Path, fma_ids: Iterable[str]) -> None:
    paths = [str(STL_ROOT / f"{fma_id}.stl") for fma_id in sorted(set(fma_ids))]
    missing = [path for path in paths if not (source / path).exists()]
    if missing:
        run("git", "checkout", "HEAD", "--", *missing, cwd=source)


def load_group(source: Path, fma_ids: Iterable[str]) -> trimesh.Trimesh:
    meshes: list[trimesh.Trimesh] = []
    for fma_id in fma_ids:
        path = source / STL_ROOT / f"{fma_id}.stl"
        if not path.exists():
            raise FileNotFoundError(path)
        loaded = trimesh.load_mesh(path, file_type="stl", process=True, validate=False)
        if isinstance(loaded, trimesh.Scene):
            loaded = loaded.to_geometry()
        if not isinstance(loaded, trimesh.Trimesh):
            raise TypeError(f"Unsupported mesh type for {path}: {type(loaded)!r}")
        loaded.remove_unreferenced_vertices()
        meshes.append(loaded)
    merged = trimesh.util.concatenate(meshes)
    merged.merge_vertices()
    merged.remove_unreferenced_vertices()
    return merged


def simplify(mesh: trimesh.Trimesh, target: int) -> trimesh.Trimesh:
    if len(mesh.faces) <= target:
        return mesh
    simplified = mesh.simplify_quadric_decimation(face_count=target, aggression=7)
    simplified.remove_unreferenced_vertices()
    return simplified


def bodyparts_to_gltf(vertices_mm: np.ndarray, floor_mm: float) -> np.ndarray:
    # BodyParts3D: X=left, Y=posterior, Z=superior. glTF: Y=up, metres.
    result = np.empty_like(vertices_mm, dtype=np.float32)
    result[:, 0] = vertices_mm[:, 0] * 0.001
    result[:, 1] = (vertices_mm[:, 2] - floor_mm) * 0.001
    result[:, 2] = vertices_mm[:, 1] * 0.001
    return result


def pad4(data: bytearray) -> None:
    while len(data) % 4:
        data.append(0)


def clean_geometry(vertices: np.ndarray, faces: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Drop zero-area triangles and return compact vertices with unit normals."""
    triangles = vertices[faces]
    raw_face_normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    squared_lengths = np.einsum("ij,ij->i", raw_face_normals, raw_face_normals)
    faces = faces[squared_lengths > 1e-20]
    if not len(faces):
        raise RuntimeError("Mesh contains no non-degenerate triangles")

    used, inverse = np.unique(faces.reshape(-1), return_inverse=True)
    vertices = np.asarray(vertices[used], dtype=np.float32)
    faces = inverse.reshape((-1, 3)).astype(np.uint32)

    triangles = vertices[faces]
    face_normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]).astype(np.float64)
    face_lengths = np.linalg.norm(face_normals, axis=1)
    face_normals /= face_lengths[:, None]

    normals = np.zeros((len(vertices), 3), dtype=np.float64)
    for corner in range(3):
        np.add.at(normals, faces[:, corner], face_normals)
    lengths = np.linalg.norm(normals, axis=1)
    zero = lengths < 1e-12
    normals[~zero] /= lengths[~zero, None]
    if np.any(zero):
        fallback = vertices[zero].astype(np.float64) - vertices.mean(axis=0, dtype=np.float64)
        fallback_lengths = np.linalg.norm(fallback, axis=1)
        fallback[fallback_lengths > 1e-12] /= fallback_lengths[fallback_lengths > 1e-12, None]
        fallback[fallback_lengths <= 1e-12] = (0.0, 1.0, 0.0)
        normals[zero] = fallback

    return vertices, faces, normals.astype(np.float32)


def make_material(name: str, spec: dict) -> Material:
    material = Material(
        name=f"MAT_{name}",
        pbrMetallicRoughness=PbrMetallicRoughness(
            baseColorFactor=list(spec["base"]),
            metallicFactor=spec["metallic"],
            roughnessFactor=spec["roughness"],
        ),
        alphaMode=spec.get("alpha", "OPAQUE"),
        doubleSided=spec.get("double", False),
    )
    return material


def build_glb(
    components: list[tuple[ComponentSpec, trimesh.Trimesh]],
    floor_mm: float,
    output: Path,
) -> dict:
    gltf = GLTF2(
        asset=Asset(
            version="2.0",
            generator="AYARI Clinic anatomy builder 1.0",
            copyright=UPSTREAM_ATTRIBUTION,
            extras={
                "units": "metres",
                "upAxis": "Y",
                "origin": "centre of feet at floor",
                "source": "BodyParts3D release 3.0",
                "sourceCommit": UPSTREAM_COMMIT,
                "medicalUse": "Reference visualisation only; not a diagnostic or surgical-planning device.",
            },
        ),
        scene=0,
        scenes=[Scene(name="AYARI Anatomy", nodes=[0])],
        nodes=[],
        meshes=[],
        materials=[],
        buffers=[Buffer(byteLength=0)],
        bufferViews=[],
        accessors=[],
    )

    material_order = list(MATERIALS)
    gltf.materials = [make_material(name, MATERIALS[name]) for name in material_order]
    material_index = {name: idx for idx, name in enumerate(material_order)}

    # Root -> geometry and locators. All transforms are baked into metre-scale vertices.
    gltf.nodes.append(
        Node(
            name="AYARI_Human_Adult_Male",
            children=[1, 2],
            extras={
                "coordinateSystem": "right-handed, Y-up, metres",
                "origin": "foot-floor centre",
                "referenceSubject": "adult male atlas",
            },
        )
    )
    gltf.nodes.append(Node(name="GEOMETRY", children=[]))
    gltf.nodes.append(Node(name="LOCATORS", children=[]))

    binary = bytearray()
    report_components: list[dict] = []
    locator_positions: dict[str, list[float]] = {}

    for spec, mesh in components:
        vertices = bodyparts_to_gltf(np.asarray(mesh.vertices), floor_mm)
        faces = np.asarray(mesh.faces, dtype=np.uint32)
        vertices, faces, normals = clean_geometry(vertices, faces)

        position_offset = len(binary)
        position_bytes = vertices.astype("<f4", copy=False).tobytes()
        binary.extend(position_bytes)
        pad4(binary)
        position_view = len(gltf.bufferViews)
        gltf.bufferViews.append(
            BufferView(buffer=0, byteOffset=position_offset, byteLength=len(position_bytes), target=ARRAY_BUFFER)
        )
        position_accessor = len(gltf.accessors)
        gltf.accessors.append(
            Accessor(
                bufferView=position_view,
                byteOffset=0,
                componentType=FLOAT,
                count=len(vertices),
                type="VEC3",
                min=vertices.min(axis=0).astype(float).tolist(),
                max=vertices.max(axis=0).astype(float).tolist(),
            )
        )

        normal_offset = len(binary)
        normal_bytes = normals.astype("<f4", copy=False).tobytes()
        binary.extend(normal_bytes)
        pad4(binary)
        normal_view = len(gltf.bufferViews)
        gltf.bufferViews.append(
            BufferView(buffer=0, byteOffset=normal_offset, byteLength=len(normal_bytes), target=ARRAY_BUFFER)
        )
        normal_accessor = len(gltf.accessors)
        gltf.accessors.append(
            Accessor(
                bufferView=normal_view,
                byteOffset=0,
                componentType=FLOAT,
                count=len(normals),
                type="VEC3",
            )
        )

        index_data = faces.reshape(-1)
        index_offset = len(binary)
        index_bytes = index_data.astype("<u4", copy=False).tobytes()
        binary.extend(index_bytes)
        pad4(binary)
        index_view = len(gltf.bufferViews)
        gltf.bufferViews.append(
            BufferView(buffer=0, byteOffset=index_offset, byteLength=len(index_bytes), target=ELEMENT_ARRAY_BUFFER)
        )
        index_accessor = len(gltf.accessors)
        gltf.accessors.append(
            Accessor(
                bufferView=index_view,
                byteOffset=0,
                componentType=UNSIGNED_INT,
                count=len(index_data),
                type="SCALAR",
                min=[int(index_data.min())],
                max=[int(index_data.max())],
            )
        )

        mesh_index = len(gltf.meshes)
        gltf.meshes.append(
            Mesh(
                name=f"MESH_{spec.key}",
                primitives=[
                    Primitive(
                        attributes={"POSITION": position_accessor, "NORMAL": normal_accessor},
                        indices=index_accessor,
                        material=material_index[spec.material],
                    )
                ],
                extras={"fmaIds": list(spec.fma_ids), "labelJa": spec.label_ja},
            )
        )
        node_index = len(gltf.nodes)
        gltf.nodes.append(
            Node(
                name=f"GEO_{spec.key}",
                mesh=mesh_index,
                extras={
                    "category": "skin" if spec.key == "skin" else ("skeletal" if spec.key == "skeleton" else "organ"),
                    "fmaIds": list(spec.fma_ids),
                    "labelJa": spec.label_ja,
                },
            )
        )
        gltf.nodes[1].children.append(node_index)

        centre = ((vertices.min(axis=0) + vertices.max(axis=0)) * 0.5).astype(float).tolist()
        for alias in spec.marker_aliases:
            locator_positions[alias] = centre

        report_components.append(
            {
                "key": spec.key,
                "labelJa": spec.label_ja,
                "fmaIds": list(spec.fma_ids),
                "vertices": int(len(vertices)),
                "triangles": int(len(faces)),
                "boundsMetres": {
                    "min": vertices.min(axis=0).astype(float).tolist(),
                    "max": vertices.max(axis=0).astype(float).tolist(),
                },
            }
        )

    # Whole-lung and paired-organ aliases are derived from child centres.
    def midpoint(*names: str) -> list[float]:
        values = np.asarray([locator_positions[name] for name in names], dtype=float)
        return values.mean(axis=0).tolist()

    locator_positions["肺"] = midpoint("右肺", "左肺")
    locator_positions["腎臓"] = midpoint("右腎臓", "左腎臓")
    locator_positions["副腎"] = midpoint("右副腎", "左副腎")

    # '炎症' is a generic marker-category fallback, anchored at the heart. Apps
    # should prefer an organ-specific locator whenever the measurement supplies one.
    locator_positions["炎症"] = list(locator_positions["心臓"])

    for name in sorted(locator_positions):
        node_index = len(gltf.nodes)
        gltf.nodes.append(
            Node(
                name=name,
                translation=locator_positions[name],
                extras={
                    "nodeType": "locator",
                    "markerKey": name,
                    "emissiveColor": "#c9a96e",
                    "genericFallback": name == "炎症",
                },
            )
        )
        gltf.nodes[2].children.append(node_index)

    gltf.buffers[0].byteLength = len(binary)
    gltf.set_binary_blob(bytes(binary))
    output.parent.mkdir(parents=True, exist_ok=True)
    gltf.save_binary(output)

    all_min = np.min([c["boundsMetres"]["min"] for c in report_components], axis=0)
    all_max = np.max([c["boundsMetres"]["max"] for c in report_components], axis=0)
    return {
        "asset": Path(os.path.relpath(output, Path.cwd())).as_posix(),
        "source": {
            "name": "BodyParts3D release 3.0",
            "commit": UPSTREAM_COMMIT,
            "attribution": UPSTREAM_ATTRIBUTION,
            "license": "CC BY-SA 2.1 Japan",
        },
        "coordinateSystem": {
            "units": "metres",
            "upAxis": "Y",
            "origin": "centre of feet at floor",
            "sourceMapping": {"gltfX": "BodyParts3D X", "gltfY": "BodyParts3D Z", "gltfZ": "BodyParts3D Y"},
        },
        "triangles": sum(c["triangles"] for c in report_components),
        "vertices": sum(c["vertices"] for c in report_components),
        "boundsMetres": {"min": all_min.astype(float).tolist(), "max": all_max.astype(float).tolist()},
        "locators": locator_positions,
        "components": report_components,
        "medicalUse": "Reference visualisation only; not validated for diagnosis or surgical planning.",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-cache", type=Path, default=Path(".cache/BodyParts3D"))
    parser.add_argument("--output", type=Path, default=Path("assets/anatomy/ayari_human_anatomy.glb"))
    parser.add_argument("--report", type=Path, default=Path("assets/anatomy/model_report.json"))
    args = parser.parse_args()

    source = ensure_checkout(args.source_cache.resolve())
    names = read_name_map(source)
    available = available_fma_ids(source)

    skeleton_ids = select_skeleton_ids(names, available)
    brain_ids = select_brain_ids(source, names, available)
    dynamic_components = (
        ComponentSpec("skeleton", "骨格", tuple(skeleton_ids), 18_000, "bone", ("骨格", "骨")),
        ComponentSpec("brain", "脳", tuple(brain_ids), 3_500, "brain", ("脳",)),
    )
    specs = [FIXED_COMPONENTS[0], *dynamic_components, *FIXED_COMPONENTS[1:]]

    all_ids = [fma_id for spec in specs for fma_id in spec.fma_ids]
    missing_ids = sorted(set(all_ids) - available)
    if missing_ids:
        raise RuntimeError(f"BodyParts3D source is missing expected IDs: {missing_ids}")
    checkout_source_meshes(source, all_ids)

    # The skin is the only full-height component, so its inferior Z bound gives
    # a stable floor independent of later anatomy selection.
    skin_source = load_group(source, ("FMA7163",))
    floor_mm = float(skin_source.bounds[0, 2])

    built: list[tuple[ComponentSpec, trimesh.Trimesh]] = []
    for spec in specs:
        print(f"Building {spec.key}: {len(spec.fma_ids)} source mesh(es)", flush=True)
        mesh = skin_source.copy() if spec.key == "skin" else load_group(source, spec.fma_ids)
        mesh = simplify(mesh, spec.target_triangles)
        built.append((spec, mesh))

    report = build_glb(built, floor_mm, args.output.resolve())
    report["sourceSelection"] = {
        "skeletonFmaIds": skeleton_ids,
        "brainFmaIds": brain_ids,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    if report["triangles"] > 50_000:
        raise RuntimeError(f"Triangle budget exceeded: {report['triangles']}")
    print(json.dumps({"triangles": report["triangles"], "vertices": report["vertices"], "asset": report["asset"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
