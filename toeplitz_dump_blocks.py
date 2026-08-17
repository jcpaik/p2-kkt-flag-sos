#!/usr/bin/env python3
"""Write the block dump of the Jensen/Toeplitz-augmented problem.

Extends sdpa_runs/blocks_deg14_h2w_h2all.json (the --dump-blocks
payload of the base weighted problem) with the appended families from
the toeplitz_export.py capture, so fingerprint_blocks.py /
fingerprint_expand.py work unchanged on sel_toep_* results.

Usage:
  .venv/bin/python toeplitz_dump_blocks.py \
      [--out sdpa_runs/blocks_deg14_h2w_h2all_toep.json]
"""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

CAPTURE_PATH = Path(
    "/private/tmp/claude-501/-Users-jcpaik-Documents-research-"
    "p2-kkt-flag-sos/2d5da291-4f4d-44a8-bb7b-d3ca80702a32/scratchpad/"
    "toeplitz_capture_deg14_h2w.pkl"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base", default="sdpa_runs/blocks_deg14_h2w_h2all.json"
    )
    parser.add_argument(
        "--out", default="sdpa_runs/blocks_deg14_h2w_h2all_toep.json"
    )
    args = parser.parse_args()

    with open(args.base) as handle:
        payload = json.load(handle)
    with open(CAPTURE_PATH, "rb") as handle:
        capture = pickle.load(handle)

    labels = set(payload["labels"])
    for name, label_matrices in capture["extra_blocks"]:
        payload["blocks"][name] = {
            str(label): matrix.tolist()
            for label, matrix in label_matrices.items()
        }
        labels.update(
            str(label) for label in label_matrices
        )
    payload["labels"] = sorted(labels)
    Path(args.out).write_text(json.dumps(payload))
    print(
        f"wrote {args.out}: {len(payload['blocks'])} blocks, "
        f"{len(payload['labels'])} labels"
    )


if __name__ == "__main__":
    main()
