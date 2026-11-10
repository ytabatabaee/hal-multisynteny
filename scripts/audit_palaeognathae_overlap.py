"""Audit all available preserved Palaeognathae direct-edge PSL captures."""

import argparse

from hal_multisynteny.overlap_audit import audit_preserved_psl


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-dir", default="experiments/palaeognathae/audit")
    parser.add_argument("--output-prefix", default="experiments/psl-overlap/palaeognathae")
    args = parser.parse_args()
    result = audit_preserved_psl(args.audit_dir, args.output_prefix)
    print(
        f"Audited {result['input_block_count']} intervals across {len(result['edges'])} preserved edges"
    )


if __name__ == "__main__":
    main()
