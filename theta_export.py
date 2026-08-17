#!/usr/bin/env python3
"""GMP-ready export: all-measures cone + gap scalar cuts + theta-atom cuts.

Produces an exactly rational SDPA sparse file for the epsilon = 0
feasibility problem over the proof-carrying (all-measures) degree-14
cone augmented with the theta-atom window cuts of
docs/THETA_ATOM_NOTE.md section 3:

    sum_{N' < |n|, within degree} q^{n^2} Q[Ghat_n m](y)  <=  T^f_q(N')

for every even-leaf-parity monomial multiplier m within the degree cap
(per-multiplier atoms tau_{f,q,m}) and every truncation N'.  On the
unprojected cone the lower (G) rows are implied by the full two-root
Gram blocks, so only the upper window cuts are new.

Pipeline: sos_search.py --export-sdpa writes the base problem plus its
exact .map.json (base point y0 and equality-kernel directions q_j with
normalizers N_j); the atom rows are then translated to the exported
z-coordinates exactly (Fractions end to end; no float round trip) and
appended as one diagonal block.  Every appended row is verified exactly
at the ONB measure before writing.

Usage:
  .venv/bin/python theta_export.py --out sdpa_runs/deg14_allm_gapcuts_theta.dat-s
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import theta_atoms
from sos_search import fraction_decimal, onb_label_value

BASE_TOGGLES = [
    "--degree", "14", "--no-pointwise-sos", "--harmonics",
    "--three-point-flags", "--four-point-flags", "--two-root-flags",
    "--rank-relations", "--gap-scalar-cuts",
]


def build_base(base_path: Path) -> dict:
    command = [
        sys.executable, "sos_search.py", "--export-sdpa", str(base_path)
    ] + BASE_TOGGLES
    completed = subprocess.run(
        command, capture_output=True, text=True, check=True
    )
    return json.loads(completed.stdout)


def atom_cut_rows(q: Fraction, degree: int = 14, localized: bool = False):
    """Exact label-space window cuts: (description, row, majorant).

    Constraint: sum_labels row[label] * y_label <= majorant.
    localized=True produces the h2-localized atom windows
    (h2 * Q[Ghat_n m]; h2 in [0, 1] so the same majorants apply).
    """
    cuts = []
    effective_degree = degree - (2 if localized else 0)
    for family in (2, 6):
        n_min, n_max = theta_atoms.family_window(family, effective_degree)
        leaves_by_multiplier: dict[tuple, dict[int, dict]] = {}
        for n in range(n_min, n_max + 1):
            for multiplier, leaf in theta_atoms.localized_leaves(
                family, n, effective_degree
            ):
                if sum(multiplier[:2]) % 2:
                    continue
                expansion = theta_atoms.expand_q_block(leaf)
                if not expansion:
                    continue
                if localized:
                    expansion = theta_atoms.h2_localize(expansion)
                leaves_by_multiplier.setdefault(multiplier, {})[n] = expansion
        for multiplier, by_order in sorted(leaves_by_multiplier.items()):
            orders = sorted(by_order)
            cap = max(max(orders), -min(orders))
            for low in range(0, cap):
                window = [n for n in orders if abs(n) > low]
                if not window:
                    continue
                row: dict[tuple, Fraction] = {}
                for n in window:
                    weight = q ** (n * n)
                    for label, value in by_order[n].items():
                        updated = row.get(label, Fraction(0)) + weight * value
                        if updated:
                            row[label] = updated
                        else:
                            row.pop(label, None)
                majorant = theta_atoms.tail_majorant(family, q, low)
                cuts.append(
                    (
                        {
                            "kind": (
                                "h2loc_window" if localized else "window"
                            ),
                            "family": family,
                            "q": str(q),
                            "multiplier": list(multiplier),
                            "window_low": low,
                            "window_orders": window,
                            "majorant": f"{majorant.numerator}/"
                                        f"{majorant.denominator}",
                        },
                        row,
                        majorant,
                    )
                )
    return cuts


def cap_cut_rows(degree: int = 14, localized: bool = False):
    """Pointwise sup-norm caps: Q[Ghat_n](y) <= (4/3 + c|n|)^2 (L7),
    and their h2-localized copies (h2 <= 1).  Valid at every measure;
    exact rational constants.  These are the per-n sharp instruments
    against sign-alternating escapes (docs/THETA_ATOM_NOTE.md sec. 8).
    """
    cuts = []
    effective_degree = degree - (2 if localized else 0)
    for family in (2, 6):
        n_min, n_max = theta_atoms.family_window(family, effective_degree)
        for n in range(n_min, n_max + 1):
            expansion = theta_atoms.expand_q_block(
                theta_atoms.hat_generator(family, n)
            )
            if localized:
                expansion = theta_atoms.h2_localize(expansion)
            bound = (
                theta_atoms.A_CONST
                + theta_atoms.SLOPE[family] * abs(n)
            ) ** 2
            cuts.append(
                (
                    {
                        "kind": "h2loc_cap" if localized else "cap",
                        "family": family,
                        "n": n,
                        "majorant": f"{bound.numerator}/"
                                    f"{bound.denominator}",
                    },
                    expansion,
                    bound,
                )
            )
    return cuts


def verify_at_onb(cuts) -> None:
    for description, row, majorant in cuts:
        value = sum(
            coefficient * onb_label_value(label)
            for label, coefficient in row.items()
        )
        if value > majorant:
            raise SystemExit(
                f"cut violated at the ONB measure (bug): {description} "
                f"value {value} > {majorant}"
            )
        if value < 0:
            raise SystemExit(
                f"negative window mass at the ONB (bug): {description} "
                f"value {value}"
            )


def parse_base(base_path: Path):
    lines = base_path.read_text().splitlines()
    header, entries = lines[:4], lines[4:]
    m_dim = int(header[0].split("=")[0].strip())
    n_block = int(header[1].split("=")[0].strip())
    block_struct = header[2].split("=")[0].strip()
    assert block_struct.startswith("(") and block_struct.rstrip().endswith(
        ")"
    ), block_struct
    sizes = [
        int(value) for value in block_struct.strip()[1:-1].split(",")
    ]
    objective_line = header[3]
    return m_dim, n_block, sizes, objective_line, entries


def load_map(map_path: Path):
    with open(map_path) as handle:
        data = json.load(handle)
    base_point = {
        label: Fraction(value)
        for label, value in data["base_point"].items()
    }
    directions = [
        (
            Fraction(entry["normalizer"]),
            {
                label: Fraction(value)
                for label, value in entry["coefficients"].items()
            },
        )
        for entry in data["directions"]
    ]
    return base_point, directions, data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--q", default="1/2")
    parser.add_argument("--digits", type=int, default=50)
    parser.add_argument(
        "--base",
        help="reuse an existing base export instead of rebuilding",
    )
    parser.add_argument(
        "--with-localized",
        action="store_true",
        help="also append the h2-localized atom windows",
    )
    parser.add_argument(
        "--with-caps",
        action="store_true",
        help="also append the per-n sup-norm caps (plain and h2-localized)",
    )
    args = parser.parse_args()

    out_path = Path(args.out)
    q = Fraction(args.q)

    if args.base:
        base_path = Path(args.base)
        base_info = {"export": str(base_path), "reused": True}
        print(json.dumps(base_info, indent=1))
    else:
        base_path = out_path.with_suffix(".base.dat-s")
        base_info = build_base(base_path)
        print(json.dumps(base_info, indent=1))
    map_path = base_path.with_name(base_path.name + ".map.json")

    print(f"computing exact theta-atom window cuts (q = {q}) ...")
    cuts = atom_cut_rows(q)
    if args.with_localized:
        cuts += atom_cut_rows(q, localized=True)
    if args.with_caps:
        cuts += cap_cut_rows()
        cuts += cap_cut_rows(localized=True)
    print(f"  {len(cuts)} cut rows")
    verify_at_onb(cuts)
    print("  all cuts verified exactly at the ONB measure")

    base_point, directions, map_data = load_map(map_path)
    label_universe = set(base_point)
    for _, coefficients in directions:
        label_universe.update(coefficients)
    missing = {
        str(label)
        for _, row, _ in cuts
        for label in row
        if str(label) not in label_universe
    }
    if missing:
        raise SystemExit(
            f"{len(missing)} cut labels missing from the base export "
            f"(examples: {sorted(missing)[:5]})"
        )

    def pair(row: dict, coefficients: dict) -> Fraction:
        return sum(
            (
                value * coefficients[str(label)]
                for label, value in row.items()
                if str(label) in coefficients
            ),
            Fraction(0),
        )

    m_dim, n_block, sizes, objective_line, entries = parse_base(base_path)
    if m_dim != len(directions):
        raise SystemExit(
            f"map/export mismatch: mDIM {m_dim} vs "
            f"{len(directions)} directions"
        )

    cut_block_index = n_block + 1
    cut_count = len(cuts)
    new_entries = list(entries)
    appended = 0
    # Constraint row i of the diagonal cut block encodes
    #   majorant - <row, y0> - sum_j z_j <row, q_j> / N_j >= 0,
    # i.e. F_0[i,i] = <row, y0> - majorant, F_j[i,i] = -<row, q_j>/N_j.
    for cut_index, (description, row, majorant) in enumerate(cuts):
        constant = pair(row, {str(k): v for k, v in base_point.items()})
        f0_value = constant - majorant
        if f0_value != 0:
            new_entries.append(
                f"0 {cut_block_index} {cut_index + 1} {cut_index + 1} "
                f"{fraction_decimal(f0_value, args.digits)}"
            )
            appended += 1
        for direction_index, (normalizer, coefficients) in enumerate(
            directions
        ):
            value = -pair(row, coefficients) / normalizer
            if value != 0:
                new_entries.append(
                    f"{direction_index + 1} {cut_block_index} "
                    f"{cut_index + 1} {cut_index + 1} "
                    f"{fraction_decimal(value, args.digits)}"
                )
                appended += 1

    header = [
        f"{m_dim} = mDIM",
        f"{n_block + 1} = nBLOCK",
        "("
        + ", ".join(str(size) for size in sizes + [-cut_count])
        + ") = bLOCKsTRUCT",
        objective_line,
    ]
    out_path.write_text("\n".join(header + new_entries) + "\n")

    sidecar = {
        "base_export": str(base_path),
        "base_map": str(map_path),
        "objective_shift": map_data["objective_shift"],
        "note": "bound = objValPrimal + objective_shift; cut block is "
                "the final diagonal block, one row per theta cut",
        "toggles": None if args.base else BASE_TOGGLES,
        "q": str(q),
        "cut_block_index": cut_block_index,
        "cuts": [description for description, _, _ in cuts],
    }
    sidecar_path = out_path.with_name(out_path.name + ".map.json")
    sidecar_path.write_text(json.dumps(sidecar, indent=1))
    print(
        f"wrote {out_path} ({appended} appended entries, "
        f"{cut_count} diagonal cut rows) and {sidecar_path}"
    )


if __name__ == "__main__":
    main()
