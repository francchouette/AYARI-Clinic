#!/usr/bin/env python3
"""Render a dark-studio QA preview from the actual cell GLB geometry."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import trimesh
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


MATERIAL_STYLE = {
    "MAT_cell_glass": ((0.62, 0.78, 0.80), 0.055),
    "MAT_nucleus_glass": ((0.74, 0.80, 0.77), 0.12),
    "MAT_measured_chromatin": ((0.66, 0.53, 0.43), 0.08),
    "MAT_chromosome_pearl": ((0.91, 0.88, 0.78), 0.92),
    "MAT_telomere_cap": ((0.91, 0.70, 0.35), 1.00),
    "MAT_mitochondria_glass": ((0.64, 0.77, 0.77), 0.30),
    "MAT_cristae": ((0.84, 0.66, 0.39), 0.92),
    "MAT_dna_strand_a": ((0.91, 0.91, 0.84), 0.96),
    "MAT_dna_strand_b": ((0.57, 0.75, 0.77), 0.96),
    "MAT_dna_base_pairs": ((0.88, 0.68, 0.35), 0.92),
}


def material_name(mesh: trimesh.Trimesh) -> str:
    material = getattr(mesh.visual, "material", None)
    return getattr(material, "name", "") or ""


def render(asset: Path, output: Path, triangles: int | None = None) -> None:
    scene = trimesh.load(asset, force="scene")
    geometry = list(scene.geometry.items())
    if triangles is None:
        triangles = sum(len(mesh.faces) for _name, mesh in geometry)

    def layer(item: tuple[str, trimesh.Trimesh]) -> int:
        mat = material_name(item[1])
        if mat == "MAT_cell_glass":
            return 0
        if mat in ("MAT_nucleus_glass", "MAT_mitochondria_glass"):
            return 1
        return 2

    fig = plt.figure(figsize=(14, 9), dpi=140, facecolor="#070a0d")
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("#070a0d")
    ax.set_proj_type("persp", focal_length=0.92)

    light = np.array((-0.35, 0.55, 0.76), dtype=float)
    light /= np.linalg.norm(light)
    for name, mesh in sorted(geometry, key=layer):
        mat = material_name(mesh)
        base, alpha = MATERIAL_STYLE.get(mat, ((0.76, 0.78, 0.74), 0.85))
        faces = np.asarray(mesh.faces)
        # The shell is translucent; a deterministic reduction keeps the preview
        # responsive without changing the delivered GLB.
        if mat == "MAT_cell_glass" and len(faces) > 22_000:
            faces = faces[np.linspace(0, len(faces) - 1, 22_000, dtype=int)]
        vertices = np.asarray(mesh.vertices)
        triangles_xyz = vertices[faces]
        raw_normals = np.cross(
            triangles_xyz[:, 1] - triangles_xyz[:, 0],
            triangles_xyz[:, 2] - triangles_xyz[:, 0],
        )
        raw_normals /= np.maximum(np.linalg.norm(raw_normals, axis=1, keepdims=True), 1e-12)
        lambert = np.clip(np.abs(raw_normals @ light), 0.0, 1.0)
        rim = np.power(1.0 - np.abs(raw_normals[:, 2]), 2.3)
        brightness = np.clip(0.38 + 0.48 * lambert + 0.32 * rim, 0.0, 1.18)
        rgb = np.clip(np.asarray(base)[None, :] * brightness[:, None], 0.0, 1.0)
        facecolors = np.column_stack((rgb, np.full(len(rgb), alpha)))
        # Matplotlib uses Z as up; the asset uses Y as up.
        plot_triangles = triangles_xyz[:, :, (0, 2, 1)]
        collection = Poly3DCollection(
            plot_triangles,
            facecolors=facecolors,
            edgecolors=(0.79, 0.66, 0.43, 0.018 if mat == "MAT_cell_glass" else 0.0),
            linewidths=0.05,
            zsort="average",
        )
        ax.add_collection3d(collection)

    ax.set_xlim(-1.08, 1.08)
    ax.set_ylim(-0.98, 0.92)
    ax.set_zlim(-1.0, 1.0)
    ax.set_box_aspect((2.16, 1.90, 2.0))
    ax.view_init(elev=7, azim=-83, roll=0)
    ax.set_axis_off()
    fig.suptitle(
        "AYARI CELL — OPENORGANELLE REFERENCE BUILD",
        color="#eee8dc",
        fontsize=17,
        x=0.055,
        y=0.955,
        ha="left",
    )
    fig.text(
        0.057,
        0.905,
        "Measured FIB-SEM nucleus, chromatin, mitochondria and membrane/cristae · jrc_hela-2",
        color="#aeb9b8",
        fontsize=9.5,
        ha="left",
    )
    fig.text(
        0.945,
        0.055,
        f"ACTUAL GLB QA RENDER · {triangles:,} TRIANGLES · LIGHTING NOT BAKED",
        color="#98886d",
        fontsize=8.5,
        ha="right",
    )
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("asset", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.asset, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
