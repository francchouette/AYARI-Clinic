#!/usr/bin/env python3
"""Render reproducible front and side QA previews of the anatomy GLB."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import trimesh
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


COLORS = {
    "skin": (0.55, 0.64, 0.68, 0.10),
    "skeleton": (0.85, 0.81, 0.68, 0.92),
    "brain": (0.78, 0.55, 0.60, 0.98),
    "heart": (0.55, 0.06, 0.09, 1.0),
    "lung": (0.69, 0.37, 0.45, 0.94),
    "liver": (0.39, 0.15, 0.10, 1.0),
    "kidney": (0.49, 0.13, 0.11, 1.0),
    "intestine": (0.69, 0.43, 0.31, 1.0),
    "default": (0.72, 0.45, 0.28, 1.0),
}


def color_for(name: str) -> tuple[float, float, float, float]:
    for key in ("skin", "skeleton", "brain", "heart", "lung", "liver", "kidney", "intestine"):
        if key in name:
            return COLORS[key]
    return COLORS["default"]


def add_scene(ax, scene: trimesh.Scene, azimuth: float) -> None:
    draw_order = sorted(scene.geometry, key=lambda name: name == "MESH_skin")
    for name in draw_order:
        mesh = scene.geometry[name]
        vertices = np.asarray(mesh.vertices)
        # Matplotlib's Z is vertical; glTF's Y is vertical.
        plot_vertices = vertices[:, [0, 2, 1]]
        collection = Poly3DCollection(
            plot_vertices[np.asarray(mesh.faces)],
            facecolor=color_for(name),
            edgecolor="none",
            linewidth=0,
            zsort="average",
        )
        ax.add_collection3d(collection)

    ax.set_xlim(-0.38, 0.38)
    ax.set_ylim(-0.30, 0.10)
    ax.set_zlim(0.0, 1.72)
    ax.set_box_aspect((0.76, 0.40, 1.72))
    ax.view_init(elev=0, azim=azimuth)
    ax.set_axis_off()
    ax.set_facecolor("white")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("asset", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    scene = trimesh.load(args.asset, force="scene")
    triangles = sum(len(mesh.faces) for mesh in scene.geometry.values())
    profile = "HIGH-RES"
    if args.report:
        report = json.loads(args.report.read_text(encoding="utf-8"))
        triangles = report["triangles"]
        profile = report.get("profile", {}).get("name", profile).upper()
    fig = plt.figure(figsize=(11.2, 7.2), dpi=160, facecolor="white")
    front = fig.add_subplot(1, 2, 1, projection="3d")
    side = fig.add_subplot(1, 2, 2, projection="3d")
    add_scene(front, scene, -90)
    add_scene(side, scene, 0)
    front.set_title("ANTERIOR", fontsize=10, color="#5d6670", pad=4)
    side.set_title("RIGHT LATERAL", fontsize=10, color="#5d6670", pad=4)
    fig.suptitle("AYARI HUMAN ANATOMY — MULTI-LOD GLB", fontsize=15, color="#26313a", y=0.97)
    fig.text(
        0.5,
        0.035,
        f"{profile} · {triangles:,} triangles · Y-up · metres · foot-centred origin",
        ha="center",
        fontsize=9,
        color="#7b858e",
    )
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.07, top=0.92, wspace=0.0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
