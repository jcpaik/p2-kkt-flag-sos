#!/usr/bin/env python3
"""Build the trace-regularized canonical-certificate selector problem.

Given an exported problem ``PROBLEM.dat-s`` (SDPA-P form: min c'z with
X = sum z_i F_i - F0 >= 0, whose SDPA-D side max F0.Y, F_i.Y = c_i,
Y >= 0 carries the certificate), write ``SELECTOR.dat-s`` whose SDPA-D
side is

    max  -tr(Y)   subject to   F_i . Y = c_i,
                               F0 . Y - s = t0,   s >= 0,
                               Y >= 0,

i.e. the minimum-trace certificate proving the bound ``t0`` (in raw
objective units; the T-bound is ``t0 + objective_shift``).  Every matrix
entry is copied through verbatim as text, so no precision is lost.

The selector's yMat blocks (parse with sdpa_extract.parse_result) are the
canonical certificate; the final 1x1 block is the bound slack, and the
negated dual objective is the minimal trace.
"""

from __future__ import annotations

import argparse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("problem", help="input PROBLEM.dat-s")
    parser.add_argument("output", help="output SELECTOR.dat-s")
    parser.add_argument(
        "--bound",
        required=True,
        help=(
            "raw objective level t0 the certificate must reach "
            "(T-bound minus objective_shift, e.g. -0.001 - (-0.25) "
            "= 0.249 written as 2.49E-1)"
        ),
    )
    arguments = parser.parse_args()

    header: list[str] = []
    entries: list[str] = []
    with open(arguments.problem) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if len(header) < 4:
                header.append(line)
            else:
                entries.append(line)

    variable_count = int(header[0].split("=")[0].strip())
    block_count = int(header[1].split("=")[0].strip())
    sizes = [
        int(token)
        for token in header[2].split("=")[0].strip().strip("()").split(",")
    ]
    if len(sizes) != block_count:
        raise ValueError("Block structure does not match nBLOCK")
    objective = header[3].strip()
    if not (objective.startswith("{") and objective.endswith("}")):
        raise ValueError("Unexpected objective line")
    coefficients = objective[1:-1]

    slack_block = block_count + 1
    lines: list[str] = [
        f"{variable_count + 1} = mDIM",
        f"{block_count + 1} = nBLOCK",
        "("
        + ", ".join(str(size) for size in sizes)
        + ", 1) = bLOCKsTRUCT",
        "{" + coefficients + ", " + arguments.bound + "}",
    ]

    # New F0 = -Identity on the original blocks (dual objective -tr Y).
    for block_index, size in enumerate(sizes):
        for diagonal in range(size):
            lines.append(
                f"0 {block_index + 1} {diagonal + 1} {diagonal + 1} -1.0"
            )

    # Constraints 1..m: original F_i entries, untouched (verbatim text).
    # Constraint m+1: original F0 entries plus the -1 slack coupling.
    new_constraint = str(variable_count + 1)
    for line in entries:
        index, remainder = line.split(" ", 1)
        if index == "0":
            lines.append(f"{new_constraint} {remainder}")
        else:
            lines.append(line)
    lines.append(f"{new_constraint} {slack_block} 1 1 -1.0")

    with open(arguments.output, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print(
        f"wrote {arguments.output}: m={variable_count + 1}, "
        f"blocks={block_count + 1} (+1x1 slack), bound={arguments.bound}"
    )


if __name__ == "__main__":
    main()
