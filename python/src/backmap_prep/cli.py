"""Command-line interface for backmap-prep."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .builder import build_system
from .schema import load_settings
from .table_converter import convert_tables
from .writers import write_lammps_data, write_lammps_input


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="backmap-prep",
        description="Generate LAMMPS input files for CG→AT backmapping",
    )
    parser.add_argument(
        "settings",
        type=Path,
        help="YAML settings file",
    )
    parser.add_argument(
        "--settings",
        dest="settings_flag",
        type=Path,
        default=None,
        help="Alternative: --settings settings.yaml (backward compatibility)",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default=None,
        help="Override output.prefix from YAML",
    )

    args = parser.parse_args(argv)
    settings_path = args.settings_flag or args.settings

    if not settings_path.exists():
        print(f"Error: settings file not found: {settings_path}", file=sys.stderr)
        return 1

    settings = load_settings(settings_path)

    prefix = args.output_prefix or settings.output.prefix
    out_dir = settings_path.parent

    system = build_system(settings, base_dir=out_dir)

    data_path = out_dir / f"{prefix}.data"
    write_lammps_data(system, data_path)
    print(f"Wrote {data_path}")

    input_path = out_dir / f"in.{prefix}"
    write_lammps_input(system, settings, input_path, data_filename=f"{prefix}.data")
    print(f"Wrote {input_path}")

    table_files = convert_tables(system, settings, out_dir)
    for tf in table_files:
        print(f"Wrote {tf}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
