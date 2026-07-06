"""Compare hybrid GROMACS topology sections (parity helper for rim135)."""

from __future__ import annotations

import sys

from backmap_prep.network.bakery import files_io


def _canon_angle_key(key: tuple[int, int, int]) -> tuple[int, int, int]:
    a, b, c = key
    return (a, b, c) if a <= c else (c, b, a)


def _section_keys_equal(section: str, left: dict, right: dict) -> bool:
    if section == "cross_angles":
        left_keys = {_canon_angle_key(k) for k in left}
        right_keys = {_canon_angle_key(k) for k in right}
        return left_keys == right_keys
    return set(left) == set(right)


def compare_topology_files(path_a: str, path_b: str) -> dict[str, bool]:
    top_a = files_io.GROMACSTopologyFile(path_a)
    top_a.read()
    top_b = files_io.GROMACSTopologyFile(path_b)
    top_b.read()

    sections = [
        "atoms",
        "dihedrals",
        "angles",
        "bonds",
        "cross_bonds",
        "cross_angles",
        "cross_dihedrals",
        "cross_pairs",
        "pairs",
    ]
    return {
        section: _section_keys_equal(section, getattr(top_a, section), getattr(top_b, section))
        for section in sections
    }


def main(argv: list[str] | None = None) -> int:
    args = argv if argv is not None else sys.argv[1:]
    if len(args) != 2:
        print("usage: compare-topology <ref.top> <generated.top>", file=sys.stderr)
        return 2
    results = compare_topology_files(args[0], args[1])
    for section, ok in results.items():
        print((section, ok))
    return 0 if all(results.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
