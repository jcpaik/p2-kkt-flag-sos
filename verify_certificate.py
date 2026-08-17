#!/usr/bin/env python3
"""Independent verification of a selector-certificate JSON.

Recomputes, from the exported PROBLEM.dat-s data and a certificate's
block matrices:

- max |F_i . Y - c_i|      (equality residual: target reproduction)
- F_0 . Y                  (certified raw level; T-bound = level + shift)
- min eigenvalue per block (PSD check)

Works for both Clarabel (frob*.json, key "blocks") and any future
selector output in the same layout.
"""

from __future__ import annotations

import argparse
import json

import numpy as np


def parse_dats(path):
    lines = [line.rstrip("\n") for line in open(path) if line.strip()]
    m = int(lines[0].split("=")[0])
    sizes = [
        int(token)
        for token in lines[2].split("=")[0].strip().strip("()").split(",")
    ]
    c = [float(token) for token in lines[3].strip()[1:-1].split(",")]
    entries = []
    for line in lines[4:]:
        parts = line.split()
        entries.append(
            (
                int(parts[0]),
                int(parts[1]),
                int(parts[2]),
                int(parts[3]),
                float(parts[4]),
            )
        )
    return m, sizes, c, entries


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("problem", help="PROBLEM.dat-s (the original, not the selector)")
    parser.add_argument("certificate", help="frob*.json with block matrices")
    # Objective shift of the exported problem (printed in the export
    # JSON): -4/3 for target E, 2/3 for target h2*E.  Legacy exports made
    # before the E-normalization used (3/16)E and shift -1/4 or 1/8.
    parser.add_argument("--shift", type=float, default=-4.0 / 3.0)
    args = parser.parse_args()

    m, sizes, c, entries = parse_dats(args.problem)
    document = json.load(open(args.certificate))
    blocks = [np.array([[float(v) for v in row] for row in mat])
              for mat in document["blocks"]]
    # The selector added a 1x1 slack block at the end; original problem
    # has len(sizes) blocks.
    if len(blocks) == len(sizes) + 1:
        slack = float(blocks[-1][0, 0])
        blocks = blocks[: len(sizes)]
    else:
        slack = None

    dot = np.zeros(m + 1)  # index 0 = F0
    for k, b, i, j, v in entries:
        Y = blocks[b - 1]
        contribution = v * Y[i - 1, j - 1]
        if i != j:
            contribution *= 2
        dot[k] += contribution

    residuals = np.abs(dot[1:] - np.array(c))
    eigenvalues = [float(np.linalg.eigvalsh(B)[0]) for B in blocks]
    level = dot[0]
    print(f"max equality residual : {residuals.max():.3e}")
    print(f"certified raw level   : {level:.12f}")
    print(f"certified T-bound     : {level + args.shift:.3e}")
    print(f"min block eigenvalue  : {min(eigenvalues):.3e}")
    if slack is not None:
        print(f"slack block value     : {slack:.3e}")


if __name__ == "__main__":
    main()
