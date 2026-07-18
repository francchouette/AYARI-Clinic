#!/usr/bin/env python3
"""Fetch selected jrc_hela-2 N5 label volumes from Janelia OpenOrganelle.

The public source dataset is the automatic organelle segmentation of the
wild-type interphase HeLa cell ``jrc_hela-2``.  Only the coarse ``s4`` level
is used for the Web GLB build; it samples the measured FIB-SEM labels at
64 x 64 x 83.84 nm and keeps the generated asset tractable.

Source DOI: https://doi.org/10.25378/janelia.13108343
License: CC BY 4.0
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import struct
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np


BUCKET = "https://janelia-cosem-datasets.s3.us-east-1.amazonaws.com"
DATASET = "jrc_hela-2"
DEFAULT_LABELS = (
    "nucleus_seg",
    "chrom_seg",
    "mito_seg",
    "mito-mem_seg",
)
N5_DTYPES = {
    "uint8": np.dtype(">u1"),
    "uint16": np.dtype(">u2"),
    "uint32": np.dtype(">u4"),
    "uint64": np.dtype(">u8"),
    "int8": np.dtype(">i1"),
    "int16": np.dtype(">i2"),
    "int32": np.dtype(">i4"),
    "int64": np.dtype(">i8"),
    "float32": np.dtype(">f4"),
    "float64": np.dtype(">f8"),
}


def n5_prefix(label: str, scale: str) -> str:
    return f"{DATASET}/{DATASET}.n5/labels/{label}/{scale}/"


def fetch_bytes(url: str, timeout: int = 45) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "AYARI-Clinic-cell-builder/2.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def fetch_json(url: str) -> dict:
    return json.loads(fetch_bytes(url).decode("utf-8"))


def list_keys(prefix: str) -> list[str]:
    """List every public S3 object below a prefix."""
    keys: list[str] = []
    token: str | None = None
    namespace = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}
    while True:
        params = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if token:
            params["continuation-token"] = token
        xml = fetch_bytes(f"{BUCKET}/?{urllib.parse.urlencode(params)}")
        root = ET.fromstring(xml)
        keys.extend(item.text for item in root.findall("s3:Contents/s3:Key", namespace) if item.text)
        truncated = root.findtext("s3:IsTruncated", default="false", namespaces=namespace) == "true"
        if not truncated:
            return keys
        token = root.findtext("s3:NextContinuationToken", namespaces=namespace)
        if not token:
            raise RuntimeError(f"S3 listing for {prefix!r} was truncated without a continuation token")


def decode_n5_block(payload: bytes, dtype: np.dtype) -> np.ndarray:
    """Decode one gzip-compressed N5 primitive block to a Z/Y/X array."""
    if len(payload) < 4:
        raise ValueError("N5 block is shorter than its header")
    mode, ndim = struct.unpack(">HH", payload[:4])
    if mode not in (0, 1):
        raise ValueError(f"Unsupported N5 block mode {mode}")
    offset = 4
    dimensions = struct.unpack(">" + "I" * ndim, payload[offset : offset + 4 * ndim])
    offset += 4 * ndim
    if mode == 1:
        # Variable-length blocks carry their element count after the dimensions.
        element_count = struct.unpack(">I", payload[offset : offset + 4])[0]
        offset += 4
    else:
        element_count = math.prod(dimensions)
    raw = gzip.decompress(payload[offset:])
    values = np.frombuffer(raw, dtype=dtype, count=element_count)
    if values.size != math.prod(dimensions):
        raise ValueError(
            f"Variable-length N5 block has {values.size} values; "
            f"dense labels require {math.prod(dimensions)}"
        )
    # N5 stores X/Y/Z dimensions with X varying fastest.  A C-order Python
    # view is therefore Z/Y/X, which also matches the COSEM transform metadata.
    return values.reshape(tuple(reversed(dimensions)), order="C")


def fetch_volume(label: str, scale: str, workers: int) -> tuple[np.ndarray, dict]:
    prefix = n5_prefix(label, scale)
    attributes = fetch_json(f"{BUCKET}/{prefix}attributes.json")
    dtype_name = attributes["dataType"]
    if dtype_name not in N5_DTYPES:
        raise ValueError(f"Unsupported N5 dtype {dtype_name!r}")
    source_dtype = N5_DTYPES[dtype_name]
    native_dtype = source_dtype.newbyteorder("=")
    dimensions_xyz = tuple(int(value) for value in attributes["dimensions"])
    block_xyz = tuple(int(value) for value in attributes["blockSize"])
    volume = np.zeros(tuple(reversed(dimensions_xyz)), dtype=native_dtype)

    keys = [key for key in list_keys(prefix) if key.rsplit("/", 1)[-1] != "attributes.json"]

    def download(key: str) -> tuple[str, bytes]:
        return key, fetch_bytes(f"{BUCKET}/{key}")

    completed = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(download, key) for key in keys]
        for future in as_completed(futures):
            key, payload = future.result()
            chunk = decode_n5_block(payload, source_dtype).astype(native_dtype, copy=False)
            ix, iy, iz = (int(part) for part in key[len(prefix) :].split("/"))
            x0, y0, z0 = ix * block_xyz[0], iy * block_xyz[1], iz * block_xyz[2]
            z1 = min(z0 + chunk.shape[0], volume.shape[0])
            y1 = min(y0 + chunk.shape[1], volume.shape[1])
            x1 = min(x0 + chunk.shape[2], volume.shape[2])
            volume[z0:z1, y0:y1, x0:x1] = chunk[: z1 - z0, : y1 - y0, : x1 - x0]
            completed += 1
            if completed % 25 == 0 or completed == len(keys):
                print(f"{label}: {completed}/{len(keys)} N5 blocks", flush=True)

    provenance = {
        "dataset": DATASET,
        "label": label,
        "scale": scale,
        "source": f"s3://janelia-cosem-datasets/{DATASET}/{DATASET}.n5/labels/{label}/{scale}",
        "doi": "https://doi.org/10.25378/janelia.13108343",
        "license": "CC BY 4.0",
        "dimensionsXYZ": list(dimensions_xyz),
        "blockSizeXYZ": list(block_xyz),
        "pixelResolutionNmXYZ": attributes["pixelResolution"]["dimensions"],
        "dtype": dtype_name,
        "nonzeroVoxels": int(np.count_nonzero(volume)),
        "uniqueNonzeroLabels": int(max(0, np.unique(volume).size - 1)),
    }
    return volume, provenance


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path(".cache/openorganelle/jrc_hela-2"))
    parser.add_argument("--scale", default="s4")
    parser.add_argument("--labels", nargs="+", default=list(DEFAULT_LABELS))
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, dict] = {}
    for label in args.labels:
        file_label = label.replace("/", "-")
        output = args.output_dir / f"{file_label}_{args.scale}.npz"
        if output.exists() and not args.force:
            with np.load(output, allow_pickle=False) as archive:
                volume = archive["volume"]
                provenance = json.loads(str(archive["provenance"].item()))
            print(f"{label}: cached {list(volume.shape)} {volume.dtype}")
        else:
            volume, provenance = fetch_volume(label, args.scale, args.workers)
            np.savez_compressed(output, volume=volume, provenance=json.dumps(provenance, separators=(",", ":")))
            print(f"{label}: wrote {output} ({output.stat().st_size:,} bytes)")
        manifest[label] = provenance

    (args.output_dir / "manifest.json").write_text(
        json.dumps({"dataset": DATASET, "labels": manifest}, indent=2) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
