#!/usr/bin/env python3
"""Block-level fingerprint of a recession ray: A_b(r) per PSD sector.

Reads a --dump-blocks payload and a ray JSON (label -> value), assembles
each block's pairing matrix A_b(r) = sum_l r_l A_l^b, and reports trace,
Frobenius norm, rank and extremal eigenvalues per block.  This locates
the sectors through which an escape direction routes its PSD mass —
in particular whether the surviving theta-orthogonal escape lives in
the even two-root sectors (theta territory), the odd/orientation
sectors, the minors, or elsewhere.

Usage:
  .venv/bin/python theta_ray_blocks.py \
      --dump sdpa_runs/theta_ab_dump_nocuts.json \
      --ray sdpa_runs/theta_ab_ray_fingerprint_q1_2.json --variant both
"""

from __future__ import annotations

import argparse
import json

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dump", required=True)
    parser.add_argument("--ray", required=True)
    parser.add_argument(
        "--variant",
        help="results key when the ray file is a theta_ab output",
    )
    args = parser.parse_args()

    with open(args.dump) as handle:
        dump = json.load(handle)
    with open(args.ray) as handle:
        ray_data = json.load(handle)
    if args.variant:
        ray = ray_data["results"][args.variant]["ray"]
    else:
        ray = ray_data["ray"]

    print(f"ray support {len(ray)}; blocks {len(dump['blocks'])}")
    print(
        f"{'block':<28} {'size':>4} {'trace':>12} {'fro':>12} "
        f"{'min_eig':>12} {'max_eig':>12} {'rank@1e-6':>9}"
    )
    rows = []
    for name, label_matrices in dump["blocks"].items():
        size = len(next(iter(label_matrices.values())))
        matrix = np.zeros((size, size))
        for label, coefficients in label_matrices.items():
            value = ray.get(label)
            if value:
                matrix += value * np.array(coefficients)
        matrix = 0.5 * (matrix + matrix.T)
        eigenvalues = np.linalg.eigvalsh(matrix)
        rows.append(
            (
                name,
                size,
                float(np.trace(matrix)),
                float(np.linalg.norm(matrix)),
                float(eigenvalues[0]),
                float(eigenvalues[-1]),
                int(np.count_nonzero(eigenvalues > 1e-6)),
            )
        )
    rows.sort(key=lambda row: -row[3])
    for name, size, trace, fro, low, high, rank in rows:
        print(
            f"{name:<28} {size:>4} {trace:>12.4f} {fro:>12.4f} "
            f"{low:>12.2e} {high:>12.4f} {rank:>9}"
        )


if __name__ == "__main__":
    main()
