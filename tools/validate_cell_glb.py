#!/usr/bin/env python3
"""Validate the AYARI cutaway cell GLB structure, budgets, and metadata."""

from __future__ import annotations

import argparse
import json
import struct
from pathlib import Path


REQUIRED_MARKERS = {
    "GEO_cell": ("細胞", "cell"),
    "GEO_nucleus": ("核", "cell"),
    "GEO_dna": ("DNA", "dna"),
}


def read_glb(path: Path) -> dict:
    raw = path.read_bytes()
    magic, version, total_length = struct.unpack_from("<4sII", raw, 0)
    assert magic == b"glTF", magic
    assert version == 2, version
    assert total_length == len(raw), (total_length, len(raw))
    offset = 12
    chunks: dict[bytes, bytes] = {}
    while offset < len(raw):
        length, kind = struct.unpack_from("<I4s", raw, offset)
        offset += 8
        chunks[kind] = raw[offset : offset + length]
        offset += length
    assert b"JSON" in chunks and b"BIN\x00" in chunks
    document = json.loads(chunks[b"JSON"].decode("utf-8"))
    assert document["buffers"][0]["byteLength"] <= len(chunks[b"BIN\x00"])
    return document


def triangle_count(document: dict) -> int:
    total = 0
    for mesh in document["meshes"]:
        for primitive in mesh["primitives"]:
            assert primitive.get("mode", 4) == 4
            accessor = document["accessors"][primitive["indices"]]
            assert accessor["count"] % 3 == 0
            total += accessor["count"] // 3
    return total


def bounds_for_node(document: dict, node: dict) -> tuple[list[float], list[float]]:
    minimum = [float("inf")] * 3
    maximum = [float("-inf")] * 3
    mesh = document["meshes"][node["mesh"]]
    for primitive in mesh["primitives"]:
        accessor = document["accessors"][primitive["attributes"]["POSITION"]]
        minimum = [min(a, b) for a, b in zip(minimum, accessor["min"])]
        maximum = [max(a, b) for a, b in zip(maximum, accessor["max"])]
    return minimum, maximum


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("asset", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    document = read_glb(args.asset)
    report = json.loads(args.report.read_text(encoding="utf-8"))
    assert document["asset"]["version"] == "2.0"
    assert document["asset"]["extras"]["units"] == "metres"
    assert document["asset"]["extras"]["upAxis"] == "Y"
    assert document["asset"]["extras"]["origin"] == "cell centre"
    assert document["asset"]["extras"]["sourceDataset"] == "jrc_hela-2"
    assert document["asset"]["extras"]["sourceLicense"] == "CC BY 4.0"
    assert document["scenes"][document["scene"]]["nodes"] == [0]
    assert document["nodes"][0]["name"] == "AYARI_Cell_Cutaway"
    assert document["nodes"][1]["name"] == "GEOMETRY"

    names = {node["name"]: node for node in document["nodes"]}
    for name, (label_ja, marker_key) in REQUIRED_MARKERS.items():
        assert name in names, name
        assert names[name]["extras"]["labelJa"] == label_ja
        assert names[name]["extras"]["markerKey"] == marker_key

    assert names["GEO_cell"]["extras"]["measured"] is False
    assert names["GEO_nucleus"]["extras"]["measured"] is True
    assert names["GEO_nucleus"]["extras"]["sourceDataset"] == "jrc_hela-2"
    assert len(document["meshes"][names["GEO_nucleus"]["mesh"]]["primitives"]) == 2

    mitochondria = [names[f"GEO_mitochondria_{index:02d}"] for index in range(1, 8)]
    chromosomes = [names[f"GEO_chromosome_{index:02d}"] for index in range(1, 5)]
    telomeres = [names[f"GEO_telomere_{index:02d}"] for index in range(1, 17)]
    for node in mitochondria:
        assert node["extras"]["labelJa"] == "ミトコンドリア"
        assert node["extras"]["markerKey"] == "mito"
        assert node["extras"]["cristaeIncluded"] is True
        assert node["extras"]["measured"] is True
        assert node["extras"]["sourceDataset"] == "jrc_hela-2"
        assert isinstance(node["extras"]["sourceInstanceId"], int)
        assert len(document["meshes"][node["mesh"]]["primitives"]) == 2
    for node in chromosomes:
        assert node["extras"]["labelJa"] == "染色体"
        assert node["extras"]["markerKey"] == "telo"
        assert node["extras"]["measured"] is False
        assert node["extras"]["geometryProvenance"] == "educational overlay"
    for node in telomeres:
        assert node["extras"]["labelJa"] == "テロメア"
        assert node["extras"]["markerKey"] == "telo"
        assert node["extras"]["measured"] is False

    # Four protective caps must map to each chromosome.
    for chromosome_index in range(1, 5):
        mapped = [node for node in telomeres if node["extras"]["chromosomeIndex"] == chromosome_index]
        assert len(mapped) == 4, (chromosome_index, len(mapped))
        assert {node["extras"]["armTip"] for node in mapped} == {1, 2, 3, 4}

    dna_mesh = document["meshes"][names["GEO_dna"]["mesh"]]
    assert len(dna_mesh["primitives"]) == 3  # strand A, strand B, base-pair rungs
    assert names["GEO_dna"]["extras"]["turns"] == 2.2
    assert names["GEO_dna"]["extras"]["basePairRungs"] == 24
    assert names["GEO_dna"]["extras"]["measured"] is False

    # All transforms are intentionally baked; geometry nodes carry no TRS or matrix.
    for node in document["nodes"]:
        assert not any(key in node for key in ("translation", "rotation", "scale", "matrix")), node["name"]

    triangles = triangle_count(document)
    assert triangles == report["triangles"]
    assert report["sourceData"]["dataset"] == "jrc_hela-2"
    assert report["sourceData"]["license"] == "CC BY 4.0"
    assert 50_000 <= triangles <= 150_000, triangles
    file_size = args.asset.stat().st_size
    assert file_size == report["fileSizeBytes"]
    assert 2_000_000 <= file_size <= 6_000_000, file_size

    cell_min, cell_max = bounds_for_node(document, names["GEO_cell"])
    assert abs(cell_min[0] + 1.0) < 0.01 and abs(cell_max[0] - 1.0) < 0.01
    assert abs(cell_min[1] + 0.92) < 0.01 and abs(cell_max[1] - 0.92) < 0.01
    nucleus_min, nucleus_max = bounds_for_node(document, names["GEO_nucleus"])
    nucleus_centre_x = (nucleus_min[0] + nucleus_max[0]) / 2
    nucleus_diameter_x = nucleus_max[0] - nucleus_min[0]
    assert abs(nucleus_centre_x + 0.35) < 0.005
    assert abs(nucleus_diameter_x - 0.60) < 0.005
    dna_min, dna_max = bounds_for_node(document, names["GEO_dna"])
    assert 1.08 <= dna_max[1] - dna_min[1] <= 1.14
    assert dna_max[0] >= 0.92  # near the requested x≈+0.95 right-front placement
    assert dna_min[2] > 0.15

    for material in document["materials"]:
        assert "pbrMetallicRoughness" in material
        assert "KHR_materials_ior" in material["extensions"]
        assert material.get("emissiveFactor") == [0.0, 0.0, 0.0]
    assert "KHR_materials_transmission" in document["materials"][0]["extensions"]
    assert "KHR_materials_transmission" in document["materials"][1]["extensions"]

    assert not document.get("cameras")
    assert not document.get("animations")
    assert not document.get("skins")
    assert "KHR_lights_punctual" not in document.get("extensionsUsed", [])

    result = {
        "valid": True,
        "gltfVersion": document["asset"]["version"],
        "triangles": triangles,
        "meshes": len(document["meshes"]),
        "fileSizeBytes": file_size,
        "counts": {"mitochondria": len(mitochondria), "chromosomes": len(chromosomes), "telomeres": len(telomeres), "dna": 1},
        "cellBoundsMetres": {"min": cell_min, "max": cell_max},
        "dnaBoundsMetres": {"min": dna_min, "max": dna_max},
        "lightsBaked": False,
        "transformsApplied": True,
    }
    print(json.dumps(result, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
