#!/usr/bin/env python3
"""Validate AYARI-specific structure and budgets in the generated GLB."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pygltflib import GLTF2


REQUIRED_LOCATORS = {
    "炎症",
    "心臓",
    "肺",
    "右肺",
    "左肺",
    "肝臓",
    "胃",
    "膵臓",
    "脾臓",
    "腎臓",
    "右腎臓",
    "左腎臓",
    "副腎",
    "右副腎",
    "左副腎",
    "小腸",
    "大腸",
    "膀胱",
    "前立腺",
    "気管",
    "食道",
    "脳",
    "骨格",
    "骨",
    "皮膚",
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("asset", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    gltf = GLTF2().load_binary(args.asset)
    report = json.loads(args.report.read_text(encoding="utf-8"))
    names = {node.name for node in gltf.nodes}
    missing = sorted(REQUIRED_LOCATORS - names)
    assert not missing, f"Missing locator nodes: {missing}"
    for node in gltf.nodes:
        if node.name in REQUIRED_LOCATORS:
            assert node.mesh is None, f"Locator {node.name} unexpectedly has geometry"
            assert node.extras and node.extras.get("nodeType") == "locator"

    assert gltf.asset.version == "2.0"
    assert report["triangles"] <= 50_000
    bounds = report["boundsMetres"]
    # Decimation may remove the original lowest skin vertex; allow 5 mm while
    # keeping the coordinate origin itself at the source model's floor.
    assert abs(bounds["min"][1]) < 0.005, bounds
    assert 1.5 < bounds["max"][1] < 2.1, bounds
    assert len(gltf.meshes) == len(report["components"])
    assert len(gltf.animations) == 0
    assert len(gltf.skins) == 0

    print(
        json.dumps(
            {
                "valid": True,
                "gltfVersion": gltf.asset.version,
                "triangles": report["triangles"],
                "meshes": len(gltf.meshes),
                "locators": len(REQUIRED_LOCATORS),
                "boundsMetres": bounds,
                "rigged": False,
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
