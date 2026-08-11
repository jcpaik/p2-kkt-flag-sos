#!/usr/bin/env python3
"""Extract high-precision solutions from an SDPA-GMP result file.

The exporter (``sos_search.py --export-sdpa``) writes ``PROBLEM.dat-s`` and
a sidecar ``PROBLEM.dat-s.map.json``.  This script parses the SDPA-GMP
result file and reconstructs, at the printed precision (40 digits with the
repository ``param.p2``):

- the reduced solution vector ``z`` (SDPA's xVec);
- the moment vector ``y = y0 + sum_j z_j q_j / N_j`` on all labels;
- the certificate Gram matrices (SDPA's dual yMat), one per named block.

Output is a JSON document with decimal strings, suitable for degreewise
coefficient extrapolation and integer-relation detection.
"""

from __future__ import annotations

import argparse
import json
import re
from fractions import Fraction

from mpmath import mp, mpf

NUMBER = re.compile(r"[+-]\d\.\d+e[+-]\d+")


def parse_vector(text: str, marker: str) -> list[mpf]:
    # Anchor at line start to avoid collisions such as "xMatTime".
    start = text.index(f"\n{marker} = ")
    open_brace = text.index("{", start)
    depth = 0
    for position in range(open_brace, len(text)):
        if text[position] == "{":
            depth += 1
        elif text[position] == "}":
            depth -= 1
            if depth == 0:
                end = position
                break
    else:
        raise ValueError(f"Unbalanced braces after {marker}")
    return [mpf(token) for token in NUMBER.findall(text[open_brace:end])]


def parse_result(path: str) -> dict[str, object]:
    text = open(path).read()
    values: dict[str, object] = {}
    for line in text.splitlines():
        if line.startswith("objValPrimal"):
            values["objValPrimal"] = mpf(line.split("=")[1].strip())
        if line.startswith("objValDual"):
            values["objValDual"] = mpf(line.split("=")[1].strip())
        if line.startswith("phase.value"):
            values["phase"] = line.split("=")[1].strip()
    values["xVec"] = parse_vector(text, "xVec")
    values["yMat"] = parse_vector(text, "yMat")
    values["xMat"] = parse_vector(text, "xMat")
    return values


def split_blocks(flat: list[mpf], sizes: list[int]) -> list[list[list[mpf]]]:
    matrices: list[list[list[mpf]]] = []
    cursor = 0
    for size in sizes:
        block = [
            [flat[cursor + row * size + column] for column in range(size)]
            for row in range(size)
        ]
        cursor += size * size
        matrices.append(block)
    if cursor != len(flat):
        raise ValueError(
            f"Matrix payload mismatch: consumed {cursor} of {len(flat)}"
        )
    return matrices


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("problem", help="path to PROBLEM.dat-s")
    parser.add_argument("result", help="path to the SDPA-GMP result file")
    parser.add_argument("--output", help="write full JSON here")
    parser.add_argument("--digits", type=int, default=45)
    arguments = parser.parse_args()
    mp.dps = arguments.digits

    mapping = json.load(open(arguments.problem + ".map.json"))
    result = parse_result(arguments.result)

    sizes = [block["size"] for block in mapping["blocks"]]
    names = [block["name"] for block in mapping["blocks"]]

    z = result["xVec"]
    directions = mapping["directions"]
    if len(z) != len(directions):
        raise ValueError(
            f"xVec has {len(z)} entries, map has {len(directions)}"
        )

    moments: dict[str, mpf] = {
        label: mpf(Fraction(value).numerator) / mpf(Fraction(value).denominator)
        for label, value in mapping["base_point"].items()
    }
    for value, direction in zip(z, directions, strict=True):
        normalizer = Fraction(direction["normalizer"])
        scale = value / (
            mpf(normalizer.numerator) / mpf(normalizer.denominator)
        )
        for label, coefficient in direction["coefficients"].items():
            fraction = Fraction(coefficient)
            moments[label] = moments.get(label, mpf(0)) + scale * (
                mpf(fraction.numerator) / mpf(fraction.denominator)
            )

    shift = Fraction(mapping["objective_shift"])
    # objValPrimal is printed at 16 digits only; recompute c.z + shift at
    # full working precision from the 40-digit xVec.
    bound = mpf(shift.numerator) / mpf(shift.denominator)
    for value, coefficient in zip(z, mapping["objective"], strict=True):
        fraction = Fraction(coefficient)
        bound += value * (
            mpf(fraction.numerator) / mpf(fraction.denominator)
        )

    certificate = split_blocks(result["yMat"], sizes)
    moment_blocks = split_blocks(result["xMat"], sizes)

    summary = {
        "phase": result["phase"],
        "bound": str(bound),
        "bound_dual": str(
            result["objValDual"]
            + mpf(shift.numerator) / mpf(shift.denominator)
        ),
        "p2": str(moments.get("('pair', 2)", mpf(0))),
        "p4": str(moments.get("('pair', 4)", mpf(0))),
        "p6": str(moments.get("('pair', 6)", mpf(0))),
    }
    print(json.dumps(summary, indent=2))

    if arguments.output:
        document = {
            "summary": summary,
            "moments": {
                label: str(value) for label, value in moments.items()
            },
            "certificate_blocks": {
                name: [[str(value) for value in row] for row in matrix]
                for name, matrix in zip(names, certificate, strict=True)
            },
            "moment_blocks": {
                name: [[str(value) for value in row] for row in matrix]
                for name, matrix in zip(names, moment_blocks, strict=True)
            },
        }
        with open(arguments.output, "w") as handle:
            json.dump(document, handle)
        print(f"wrote {arguments.output}")


if __name__ == "__main__":
    main()
