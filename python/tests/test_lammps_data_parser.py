"""Tests for backmap_prep.parsers.lammps_data_parser."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

import pytest

from backmap_prep.parsers.lammps_data_parser import (
    parse_at_fragment,
    parse_cg_system,
    parse_lammps_data,
)

SAMPLE_DATA = """\
LAMMPS data file — test

3 atoms
1 bonds
0 angles
0 dihedrals
0 impropers

2 atom types
1 bond types
0 angle types
0 dihedral types
0 improper types

0.0 50.000000 xlo xhi
0.0 50.000000 ylo yhi
0.0 50.000000 zlo zhi

Masses

1 29.062000 # A (CG)
2 14.027000 # CH2

Atoms # full

1 1 1 0.000000 10.500000 20.300000 30.100000 0 0 0
2 1 2 0.000000 11.000000 21.000000 31.000000 1 0 -1
3 1 2 0.000000 12.000000 22.000000 32.000000

Bonds

1 1 2 3
"""


class TestParseLammpsData:
    def test_parses_atoms(self, tmp_path: Path) -> None:
        p = tmp_path / "test.data"
        p.write_text(SAMPLE_DATA)
        result = parse_lammps_data(p)
        assert len(result.atoms) == 3
        assert result.atoms[0].x == 10.5
        assert result.atoms[0].mol_id == 1
        assert result.atoms[0].type_id == 1

    def test_parses_image_flags(self, tmp_path: Path) -> None:
        p = tmp_path / "test.data"
        p.write_text(SAMPLE_DATA)
        result = parse_lammps_data(p)
        assert result.atoms[1].ix == 1
        assert result.atoms[1].iy == 0
        assert result.atoms[1].iz == -1
        assert result.atoms[2].ix == 0  # default when not present

    def test_parses_bonds(self, tmp_path: Path) -> None:
        p = tmp_path / "test.data"
        p.write_text(SAMPLE_DATA)
        result = parse_lammps_data(p)
        assert len(result.bonds) == 1
        assert result.bonds[0].i == 2
        assert result.bonds[0].j == 3

    def test_parses_box(self, tmp_path: Path) -> None:
        p = tmp_path / "test.data"
        p.write_text(SAMPLE_DATA)
        result = parse_lammps_data(p)
        assert result.box[0] == 50.0
        assert result.box[1] == 50.0
        assert result.box[2] == 50.0

    def test_parses_type_counts(self, tmp_path: Path) -> None:
        p = tmp_path / "test.data"
        p.write_text(SAMPLE_DATA)
        result = parse_lammps_data(p)
        assert result.n_atom_types == 2
        assert result.n_bond_types == 1


CG_SYSTEM_DATA = """\
LAMMPS data file -- pure CG for equilibration

6 atoms
0 bonds
0 angles
0 dihedrals
0 impropers

2 atom types
0 bond types
0 angle types

0.0 57.537600 xlo xhi
0.0 57.537600 ylo yhi
0.0 57.537600 zlo zhi

Masses

1 29.062000
2 28.054000

Atoms # full

1 1 1 0.000000 1.0 2.0 3.0
2 1 2 0.500000 4.0 5.0 6.0
3 1 1 0.000000 7.0 8.0 9.0
4 2 1 0.000000 10.0 11.0 12.0
5 2 2 0.500000 13.0 14.0 15.0
6 2 1 0.000000 16.0 17.0 18.0
"""

CG_SYSTEM_DATA_NO_MASSES = CG_SYSTEM_DATA.replace("Masses\n\n1 29.062000\n2 28.054000\n\n", "")

CG_SYSTEM_DATA_NO_ATOMS = """\
LAMMPS data file -- no atoms section

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 29.062000
"""

CG_SYSTEM_DATA_UNEVEN = """\
LAMMPS data file -- 5 atoms, not a multiple of the 3-atom first block

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 29.062000

Atoms # full

1 1 1 0.0 1.0 1.0 1.0
2 1 1 0.0 2.0 2.0 2.0
3 1 1 0.0 3.0 3.0 3.0
4 2 1 0.0 4.0 4.0 4.0
5 2 1 0.0 5.0 5.0 5.0
"""


class TestParseCgSystem:
    def test_parses_box(self, tmp_path: Path) -> None:
        p = tmp_path / "cg_system.data"
        p.write_text(CG_SYSTEM_DATA)
        gro, _top = parse_cg_system(p)
        assert gro.box == (57.5376, 57.5376, 57.5376)

    def test_parses_full_atom_list_positionally(self, tmp_path: Path) -> None:
        p = tmp_path / "cg_system.data"
        p.write_text(CG_SYSTEM_DATA)
        gro, _top = parse_cg_system(p)
        assert len(gro.atoms) == 6
        assert (gro.atoms[0].x, gro.atoms[0].y, gro.atoms[0].z) == (1.0, 2.0, 3.0)
        assert (gro.atoms[5].x, gro.atoms[5].y, gro.atoms[5].z) == (16.0, 17.0, 18.0)

    def test_derives_molecule_template_and_count(self, tmp_path: Path) -> None:
        p = tmp_path / "cg_system.data"
        p.write_text(CG_SYSTEM_DATA)
        _gro, top = parse_cg_system(p)
        assert top.molecules == [("CG", 2)]
        mol = top.molecule_types["CG"]
        assert [a.type for a in mol.atoms] == ["1", "2", "1"]
        assert mol.atoms[1].charge == pytest.approx(0.5)

    def test_type_mass_from_masses_section(self, tmp_path: Path) -> None:
        p = tmp_path / "cg_system.data"
        p.write_text(CG_SYSTEM_DATA)
        _gro, top = parse_cg_system(p)
        mol = top.molecule_types["CG"]
        assert mol.atoms[0].mass == pytest.approx(29.062)
        assert mol.atoms[1].mass == pytest.approx(28.054)
        assert top.atom_types["1"].mass == pytest.approx(29.062)

    def test_missing_masses_section_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "cg_system.data"
        p.write_text(CG_SYSTEM_DATA_NO_MASSES)
        with pytest.raises(ValueError, match="Masses"):
            parse_cg_system(p)

    def test_missing_atoms_section_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "cg_system.data"
        p.write_text(CG_SYSTEM_DATA_NO_ATOMS)
        with pytest.raises(ValueError, match="Atoms"):
            parse_cg_system(p)

    def test_uneven_molecule_grouping_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "cg_system.data"
        p.write_text(CG_SYSTEM_DATA_UNEVEN)
        with pytest.raises(ValueError, match="not a multiple"):
            parse_cg_system(p)


AT_FRAGMENT_DATA = """\
LAMMPS data file -- AT fragment test (4-atom chain)

4 atoms
3 bonds
2 angles
1 dihedrals

2 atom types
1 bond types
1 angle types
1 dihedral types

0.0 50.0 xlo xhi
0.0 50.0 ylo yhi
0.0 50.0 zlo zhi

Masses

1 15.035
2 14.027

Atoms # full

1 1 1 0.0 0.0 0.0 0.0
2 1 2 0.0 1.5 0.0 0.0
3 1 2 0.0 3.0 0.0 0.0
4 1 1 0.0 4.5 0.0 0.0

Bonds

1 1 1 2
2 1 2 3
3 1 3 4

Angles

1 1 1 2 3
2 1 2 3 4

Dihedrals

1 1 1 2 3 4
"""

AT_FRAGMENT_SCRIPT = """\
units real
pair_style lj/cut 14.0
pair_coeff 1 1 0.207 3.748
pair_coeff 2 2 0.118 3.905
bond_style harmonic
bond_coeff 1 800.0 1.53
angle_style harmonic
angle_coeff 1 126.6 111.0
dihedral_style ryckaert
dihedral_coeff 1 9.27 12.15 -13.12 -3.06 26.24 -13.49
"""

AT_FRAGMENT_DATA_NO_MASSES = AT_FRAGMENT_DATA.replace("Masses\n\n1 15.035\n2 14.027\n\n", "")

AT_FRAGMENT_DATA_NO_ATOMS = """\
LAMMPS data file -- no atoms section

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 15.035
"""

AT_FRAGMENT_DATA_GAP = AT_FRAGMENT_DATA.replace(
    "4 1 1 0.0 4.5 0.0 0.0\n", "5 1 1 0.0 4.5 0.0 0.0\n"
)


class TestParseAtFragment:
    def test_parses_atom_count(self, tmp_path: Path) -> None:
        dp = tmp_path / "at.data"
        sp = tmp_path / "at.in"
        dp.write_text(AT_FRAGMENT_DATA)
        sp.write_text(AT_FRAGMENT_SCRIPT)
        gro, top = parse_at_fragment(dp, sp)
        assert len(gro.atoms) == 4
        mol = next(iter(top.molecule_types.values()))
        assert len(mol.atoms) == 4

    def test_atom_name_is_numeric_id(self, tmp_path: Path) -> None:
        dp = tmp_path / "at.data"
        sp = tmp_path / "at.in"
        dp.write_text(AT_FRAGMENT_DATA)
        sp.write_text(AT_FRAGMENT_SCRIPT)
        _gro, top = parse_at_fragment(dp, sp)
        mol = next(iter(top.molecule_types.values()))
        assert [a.name for a in mol.atoms] == ["1", "2", "3", "4"]

    def test_atom_types_namespaced_distinct_from_cg(self, tmp_path: Path) -> None:
        dp = tmp_path / "at.data"
        sp = tmp_path / "at.in"
        dp.write_text(AT_FRAGMENT_DATA)
        sp.write_text(AT_FRAGMENT_SCRIPT)
        _gro, top = parse_at_fragment(dp, sp)
        # Must not collide with parse_cg_system's bare "1"/"2" type-name convention.
        assert set(top.atom_types) == {"AT1", "AT2"}

    def test_bond_params_are_r0_k_order(self, tmp_path: Path) -> None:
        dp = tmp_path / "at.data"
        sp = tmp_path / "at.in"
        dp.write_text(AT_FRAGMENT_DATA)
        sp.write_text(AT_FRAGMENT_SCRIPT)
        _gro, top = parse_at_fragment(dp, sp)
        mol = next(iter(top.molecule_types.values()))
        bond = mol.bonds[0]
        assert bond.func == 1
        assert bond.params == pytest.approx([1.53, 800.0])

    def test_angle_params_are_theta0_k_order(self, tmp_path: Path) -> None:
        dp = tmp_path / "at.data"
        sp = tmp_path / "at.in"
        dp.write_text(AT_FRAGMENT_DATA)
        sp.write_text(AT_FRAGMENT_SCRIPT)
        _gro, top = parse_at_fragment(dp, sp)
        mol = next(iter(top.molecule_types.values()))
        angle = mol.angles[0]
        assert angle.func == 1
        assert angle.params == pytest.approx([111.0, 126.6])

    def test_dihedral_params_passed_through_unconverted(self, tmp_path: Path) -> None:
        dp = tmp_path / "at.data"
        sp = tmp_path / "at.in"
        dp.write_text(AT_FRAGMENT_DATA)
        sp.write_text(AT_FRAGMENT_SCRIPT)
        _gro, top = parse_at_fragment(dp, sp)
        mol = next(iter(top.molecule_types.values()))
        dihedral = mol.dihedrals[0]
        assert dihedral.func == 3
        assert dihedral.params == pytest.approx([9.27, 12.15, -13.12, -3.06, 26.24, -13.49])

    def test_pair_sigma_epsilon_from_diagonal_pair_coeff(self, tmp_path: Path) -> None:
        dp = tmp_path / "at.data"
        sp = tmp_path / "at.in"
        dp.write_text(AT_FRAGMENT_DATA)
        sp.write_text(AT_FRAGMENT_SCRIPT)
        _gro, top = parse_at_fragment(dp, sp)
        assert top.atom_types["AT1"].sigma == pytest.approx(3.748)
        assert top.atom_types["AT1"].epsilon == pytest.approx(0.207)

    def test_missing_masses_section_raises(self, tmp_path: Path) -> None:
        dp = tmp_path / "at.data"
        sp = tmp_path / "at.in"
        dp.write_text(AT_FRAGMENT_DATA_NO_MASSES)
        sp.write_text(AT_FRAGMENT_SCRIPT)
        with pytest.raises(ValueError, match="Masses"):
            parse_at_fragment(dp, sp)

    def test_missing_atoms_section_raises(self, tmp_path: Path) -> None:
        dp = tmp_path / "at.data"
        sp = tmp_path / "at.in"
        dp.write_text(AT_FRAGMENT_DATA_NO_ATOMS)
        sp.write_text(AT_FRAGMENT_SCRIPT)
        with pytest.raises(ValueError, match="Atoms"):
            parse_at_fragment(dp, sp)

    def test_non_contiguous_atom_ids_raise(self, tmp_path: Path) -> None:
        dp = tmp_path / "at.data"
        sp = tmp_path / "at.in"
        dp.write_text(AT_FRAGMENT_DATA_GAP)
        sp.write_text(AT_FRAGMENT_SCRIPT)
        with pytest.raises(ValueError, match="contiguous"):
            parse_at_fragment(dp, sp)

    def test_missing_bond_coeff_raises(self, tmp_path: Path) -> None:
        dp = tmp_path / "at.data"
        sp = tmp_path / "at.in"
        dp.write_text(AT_FRAGMENT_DATA)
        sp.write_text("units real\n")
        with pytest.raises(ValueError, match="bond type 1"):
            parse_at_fragment(dp, sp)

    def test_missing_pair_coeff_raises(self, tmp_path: Path) -> None:
        dp = tmp_path / "at.data"
        sp = tmp_path / "at.in"
        dp.write_text(AT_FRAGMENT_DATA)
        sp.write_text(
            "units real\n"
            "bond_style harmonic\nbond_coeff 1 800.0 1.53\n"
            "angle_style harmonic\nangle_coeff 1 126.6 111.0\n"
            "dihedral_style ryckaert\ndihedral_coeff 1 9.27 12.15 -13.12 -3.06 26.24 -13.49\n"
        )
        with pytest.raises(ValueError, match="pair type"):
            parse_at_fragment(dp, sp)
