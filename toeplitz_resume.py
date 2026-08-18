#!/usr/bin/env python3
"""Resume a toeplitz_export.py run from its saved capture.

The capture (written after all families are built and exactly
verified, before the exact elimination) contains every argument of
sos_search.export_sdpa_problem plus the verified extra blocks; this
driver re-runs only the elimination + write phase, skipping the base
model rebuild and the family verification.

Usage:
  .venv/bin/python toeplitz_resume.py CAPTURE.pkl OUT.dat-s
"""

from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import sos_search


def main() -> None:
    capture_file, out_path = sys.argv[1], sys.argv[2]
    with open(capture_file, "rb") as handle:
        capture = pickle.load(handle)
    print(
        f"[resume] {len(capture['psd_blocks'])} base blocks + "
        f"{len(capture['extra_blocks'])} verified families, "
        f"{len(capture['ordered_labels'])} labels",
        file=sys.stderr,
    )
    result = sos_search.export_sdpa_problem(
        Path(out_path),
        50,
        capture["target"],
        capture["ordered_labels"],
        list(capture["psd_blocks"]) + list(capture["extra_blocks"]),
        capture["free_label_matrices"],
        capture["relations"],
    )
    print(json.dumps(result, indent=2, sort_keys=True, default=str))
    sidecar = {
        "note": (
            "resumed from capture (families pre-verified exactly); "
            "bound = objValPrimal + 2/3"
        ),
        "capture": capture_file,
        "families": [
            {
                "name": name,
                "size": next(iter(matrices.values())).shape[0],
                "labels": len(matrices),
            }
            for name, matrices in capture["extra_blocks"]
        ],
    }
    sidecar_path = Path(out_path + ".families.json")
    sidecar_path.write_text(json.dumps(sidecar, indent=1))
    print(f"wrote {sidecar_path}")


if __name__ == "__main__":
    main()
