"""Build one persistent native HAL segment index without retaining rows in RAM."""

from __future__ import annotations

import argparse
import json
import shutil

from hal_multisynteny.native_index import build_index


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hal", required=True)
    parser.add_argument("--child", required=True)
    parser.add_argument("--parent", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--executable", default=shutil.which("halAlignedExtract") or "halAlignedExtract")
    args = parser.parse_args()
    print(json.dumps(build_index(args.hal, args.child, args.parent, args.output, args.executable), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
