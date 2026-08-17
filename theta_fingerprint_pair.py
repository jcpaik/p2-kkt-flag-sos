#!/usr/bin/env python3
"""Pair the theta-atom window functionals against the weighted-escape D.

D (sdpa_runs/fingerprint_D_e3e4.json, docs/UNPROJECTED_ESCAPE_NOTE.md)
is the label expansion of the growth direction of the weighted
selector certificates between eps = 1e-3 and 1e-4.  This script pairs,
by plain dot product in label space as directed:

  plain atoms      l_n      = <expand Q[Ghat_n], D>
  localized atoms  l_n^loc  = <expand h2*Q[Ghat_n], D>

and the corresponding homogenized window functionals
-sum_{N'<|n|<=N} q^{n^2} l_n^(loc).  Positive window mass = the valid
upper cut pairs with the "cut sign" (kills the escape per the
recession-cone analysis of docs/THETA_ATOM_NOTE.md section 3).

Usage:
  .venv/bin/python theta_fingerprint_pair.py \
      --fingerprint sdpa_runs/fingerprint_D_e3e4.json
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction

import theta_atoms


def row_pairing(expansion, vector):
    hit_mass = 0.0
    miss_mass = 0.0
    pairing = 0.0
    for label, coefficient in expansion.items():
        key = str(label)
        if key in vector:
            pairing += float(coefficient) * vector[key]
            hit_mass += abs(float(coefficient))
        else:
            miss_mass += abs(float(coefficient))
    return pairing, hit_mass, miss_mass


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fingerprint", default="sdpa_runs/fingerprint_D_e3e4.json"
    )
    parser.add_argument("--key", default="D")
    parser.add_argument("--degree", type=int, default=14)
    parser.add_argument("--q-values", default="1/4,1/2,3/4,9/10")
    parser.add_argument("--out")
    args = parser.parse_args()

    with open(args.fingerprint) as handle:
        vector = json.load(handle)[args.key]
    q_values = [Fraction(text) for text in args.q_values.split(",")]

    constant_component = vector.get(str(("constant",)), 0.0)
    per_generator = []
    for localized in (False, True):
        # localized rows carry the extra p2 factor: degree budget - 2
        effective_degree = args.degree - (2 if localized else 0)
        for family in (2, 6):
            n_min, n_max = theta_atoms.family_window(
                family, effective_degree
            )
            for n in range(n_min, n_max + 1):
                expansion = theta_atoms.expand_q_block(
                    theta_atoms.hat_generator(family, n)
                )
                if localized:
                    expansion = theta_atoms.h2_localize(expansion)
                pairing, hit, miss = row_pairing(expansion, vector)
                cap = float(
                    (theta_atoms.A_CONST + theta_atoms.SLOPE[family]
                     * abs(n)) ** 2
                )
                per_generator.append(
                    {
                        "localized": localized,
                        "family": family,
                        "n": n,
                        "pairing": pairing,
                        "sup_cap": cap,
                        # valid cut cap - L_n(y) >= 0, paired against D
                        "cap_cut_pairing": cap * constant_component
                        - pairing,
                        "miss_mass_fraction": (
                            miss / (hit + miss) if hit + miss else 0.0
                        ),
                    }
                )

    print("per-generator pairings against D "
          "(l_n = <row, D>, plain dot product):")
    print(f"{'kind':>10} {'family':>6} {'n':>3} {'l_n(D)':>14} "
          f"{'cap-cut pair':>14} {'miss%':>6}")
    for row in per_generator:
        kind = "h2loc" if row["localized"] else "plain"
        print(
            f"{kind:>10} {row['family']:>6} {row['n']:>3} "
            f"{row['pairing']:>14.4f} {row['cap_cut_pairing']:>14.4f} "
            f"{100*row['miss_mass_fraction']:>5.1f}%"
        )
    killers = [
        row for row in per_generator if row["cap_cut_pairing"] < 0
    ]
    print(
        f"\nsup-cap cuts with negative (kill-sign) pairing: "
        f"{len(killers)} of {len(per_generator)}; strongest:"
    )
    for row in sorted(killers, key=lambda r: r["cap_cut_pairing"])[:6]:
        kind = "h2loc" if row["localized"] else "plain"
        print(
            f"  cap[{kind},{row['family']},n={row['n']}] "
            f"pairing {row['cap_cut_pairing']:.4g}"
        )

    windows = []
    for localized in (False, True):
        values = {
            (row["family"], row["n"]): row["pairing"]
            for row in per_generator
            if row["localized"] == localized
        }
        for family in (2, 6):
            orders = sorted(n for f, n in values if f == family)
            if not orders:
                continue
            cap = max(max(orders), -min(orders))
            for q in q_values:
                for low in range(0, cap):
                    window = [n for n in orders if abs(n) > low]
                    if not window:
                        continue
                    mass = sum(
                        float(q) ** (n * n) * values[(family, n)]
                        for n in window
                    )
                    windows.append(
                        {
                            "localized": localized,
                            "family": family,
                            "q": str(q),
                            "window_low": low,
                            "window_mass": mass,
                            "cut_pairing": -mass,
                        }
                    )

    print("\nwindow functionals (cut pairing = -window mass; "
          "negative pairing = cut sign = kills the escape):")
    print(f"{'kind':>10} {'family':>6} {'q':>5} {'low':>3} "
          f"{'window_mass':>14} {'cut_pairing':>14}")
    for row in windows:
        kind = "h2loc" if row["localized"] else "plain"
        print(
            f"{kind:>10} {row['family']:>6} {row['q']:>5} "
            f"{row['window_low']:>3} {row['window_mass']:>14.4f} "
            f"{row['cut_pairing']:>14.4f}"
        )

    out_path = args.out or (
        args.fingerprint.replace(".json", "") + "_theta_pairings.json"
    )
    with open(out_path, "w") as handle:
        json.dump(
            {
                "fingerprint": args.fingerprint,
                "key": args.key,
                "per_generator": per_generator,
                "windows": windows,
            },
            handle,
            indent=1,
        )
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
