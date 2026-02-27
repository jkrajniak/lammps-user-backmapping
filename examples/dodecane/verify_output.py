# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Verify backmap-prep output for dodecane test case.

Checks:
  1. Data file structure and counts
  2. Input script keywords
  3. Bond type classification (static vs backmap/*)
  4. Atom type mapping (CG first, AT second)
  5. Table file conversion
"""

from __future__ import annotations

import sys
from pathlib import Path

PASS = 0
FAIL = 0


def check(condition: bool, msg: str) -> None:
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  PASS: {msg}")
    else:
        FAIL += 1
        print(f"  FAIL: {msg}")


def verify_data_file(path: Path) -> None:
    print(f"\n--- Verifying {path.name} ---")
    text = path.read_text()
    lines = text.splitlines()

    check("180 atoms" in text, "180 atoms (10 mol × 18 atoms/mol)")
    check("4 atom types" in text, "4 atom types (2 CG + 2 AT)")
    check("4 bond types" in text, "4 bond types")

    # Check masses section
    check("29.062000" in text, "CG type A mass = 29.062")
    check("28.054000" in text, "CG type B mass = 28.054")
    check("15.035000" in text, "AT type CH3 mass = 15.035")
    check("14.027000" in text, "AT type CH2 mass = 14.027")

    # Check box dimensions (5.75376 nm → 57.5376 Å)
    check("57.5376" in text, "Box dimension converted from nm to Angstrom")

    # Check Atoms section exists
    check("Atoms # full" in text, "Atoms section with 'full' style")
    check("Bonds" in text, "Bonds section exists")
    check("Angles" in text, "Angles section exists")

    # First atom should be CG (mol 1, type 1)
    atom_lines = []
    in_atoms = False
    for line in lines:
        if line.strip() == "Atoms # full":
            in_atoms = True
            continue
        if in_atoms and line.strip() == "":
            if atom_lines:
                break
            continue
        if in_atoms and line.strip():
            atom_lines.append(line.strip())

    if atom_lines:
        first = atom_lines[0].split()
        check(first[0] == "1", "First atom id = 1")
        check(first[1] == "1", "First atom mol_id = 1")
        check(first[2] in ("1", "2"), "First atom is CG type (1 or 2)")

    # Check that bond count is reasonable
    n_bonds = int(
        [ln for ln in lines if "bonds" in ln.lower() and "types" not in ln.lower()][
            0
        ].split()[0]
    )
    check(n_bonds > 0, f"Bond count = {n_bonds} (positive)")

    # Check angle count
    n_angles = int(
        [ln for ln in lines if "angles" in ln.lower() and "types" not in ln.lower()][
            0
        ].split()[0]
    )
    check(n_angles > 0, f"Angle count = {n_angles} (positive)")


def verify_input_script(path: Path) -> None:
    print(f"\n--- Verifying {path.name} ---")
    text = path.read_text()

    check("units real" in text, "LAMMPS real units")
    check("atom_style full" in text, "atom_style full")
    check("pair_style backmap" in text, "pair_style backmap present")
    check("bond_style hybrid" in text, "bond_style hybrid (static + backmap)")
    check("backmap/harmonic" in text, "backmap/harmonic sub-style")
    check("backmap/table" in text, "backmap/table sub-style")
    check("fix bm all backmap" in text, "fix backmap present")
    check("fix integrate at_atoms nve" in text, "NVE integration on AT atoms only")
    check("langevin" in text, "Langevin thermostat present")
    check("special_bonds" in text, "special_bonds for exclusions")
    check("group at_atoms type" in text, "AT atom group defined")
    check("group cg_atoms type" in text, "CG atom group defined")

    # Three-phase run
    check("Phase 1: CG equilibration" in text, "Phase 1 comment")
    check("Phase 2: Backmapping" in text, "Phase 2 comment")
    check("Phase 3: AT production" in text, "Phase 3 comment")
    check("fix_modify bm active no" in text, "Lambda frozen in phase 1")
    check("fix_modify bm active yes" in text, "Lambda ramp in phase 2")

    # Check pair_coeff routing
    check("atomistic" in text, "Atomistic pair coefficients present")
    check("pair_coeff" in text and "cg" in text, "CG pair coefficients present")
    check("none" in text, "'none' cross-type pair coefficients")

    # Check bond_coeff routing
    check(
        "bond_coeff" in text and "harmonic" in text, "Static harmonic bond coefficients"
    )
    check("bond_coeff" in text and "at" in text, "AT cross-bond coefficients")
    check("bond_coeff" in text and "cg" in text, "CG cross-bond coefficients")

    # Temperature
    check("298.0" in text, "Temperature = 298 K")

    # Alpha
    check("0.0001" in text, "Alpha = 0.0001")


def verify_table_file(path: Path) -> None:
    print(f"\n--- Verifying {path.name} ---")
    text = path.read_text()
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]

    check("ENTRY" in text, "Table keyword 'ENTRY' present")
    check(len(lines) > 10, f"Table has {len(lines)} data lines (>10)")

    # Check first data line has correct format: index r energy force
    data_start = False
    for line in text.splitlines():
        if line.strip().startswith("1 "):
            tokens = line.split()
            check(len(tokens) == 4, "Data line has 4 columns (idx, r, E, F)")
            r_val = float(tokens[1])
            check(r_val >= 0, f"Distance r={r_val:.4f} ≥ 0")
            data_start = True
            break
    check(data_start, "Found data entries in table")


def main() -> int:
    base = Path(__file__).parent

    data_path = base / "dodecane.data"
    input_path = base / "in.dodecane"
    table_path = base / "table_b1.table"

    if not data_path.exists():
        print(f"ERROR: {data_path} not found. Run backmap-prep first.")
        return 1
    if not input_path.exists():
        print(f"ERROR: {input_path} not found. Run backmap-prep first.")
        return 1

    verify_data_file(data_path)
    verify_input_script(input_path)

    if table_path.exists():
        verify_table_file(table_path)
    else:
        print(f"\nWARN: {table_path} not found, skipping table verification")

    print(f"\n=== Results: {PASS} passed, {FAIL} failed ===")
    return 1 if FAIL > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
