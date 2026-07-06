"""Tests for Phase 3 network hybrid preparation."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from backmap_prep.network.api import build_hybrid_gromacs, build_network_lammps
from backmap_prep.network.compare_topology import compare_topology_files
from backmap_prep.network.pbc import max_bond_length, max_euclidean_bond_length
from backmap_prep.schema import load_settings, resolve_data_dir, resolve_tables_dir
from backmap_prep.table_converter import convert_tables
from backmap_prep.writers import write_cross_pairs_file, write_lammps_data, write_lammps_input

RIM135_ENV = "BACKMAP_RIM135_DIR"
_DEFAULT_RIM135 = Path(__file__).resolve().parents[4] / "bakery" / "tests" / "rim135"


def rim135_dir() -> Path | None:
    """Return rim135 test data directory if available."""
    override = os.environ.get(RIM135_ENV)
    if override:
        path = Path(override)
        return path if path.is_dir() else None
    if _DEFAULT_RIM135.is_dir() and (_DEFAULT_RIM135 / "settings.xml").is_file():
        return _DEFAULT_RIM135
    return None


@pytest.mark.integration
def test_rim135_build_hybrid_smoke() -> None:
    """Run network engine on rim135; check atom count and output files."""
    data_dir = rim135_dir()
    if data_dir is None:
        pytest.skip(f"rim135 fixtures not found (set {RIM135_ENV} or clone bakery beside hub)")

    for name in ("hyb_conf.gro", "hyb_topol.top", "missing_definitions.txt"):
        path = data_dir / name
        if path.exists():
            path.unlink()

    result = build_hybrid_gromacs(
        data_dir / "settings.xml",
        base_dir=data_dir,
        chain_rng_seed=42,
    )

    assert result.n_atoms == 13653
    assert result.coordinates_path.is_file()
    assert result.topology_path.is_file()
    assert result.coordinates_path.stat().st_size > 0
    assert result.topology_path.stat().st_size > 0


@pytest.mark.integration
def test_rim135_build_hybrid_v2_parity() -> None:
    """Native v2 YAML produces the same hybrid outputs as bakery XML."""
    data_dir = rim135_dir()
    if data_dir is None:
        pytest.skip(f"rim135 fixtures not found (set {RIM135_ENV} or clone bakery beside hub)")

    ref_top = data_dir / "ref_hyb_topol.top"
    ref_gro = data_dir / "ref_hyb_conf.gro"
    if not ref_top.is_file() or not ref_gro.is_file():
        pytest.skip("rim135 reference outputs missing")

    v2_yaml = Path(__file__).resolve().parents[2] / "examples" / "epoxy" / "settings.v2.yaml"
    if not v2_yaml.is_file():
        pytest.skip("settings.v2.yaml example missing")

    from backmap_prep.schema import load_settings, resolve_data_dir

    settings = load_settings(v2_yaml)
    work_dir = resolve_data_dir(v2_yaml, settings)

    for name in ("hyb_conf.gro", "hyb_topol.top", "missing_definitions.txt"):
        path = work_dir / name
        if path.exists():
            path.unlink()

    result = build_hybrid_gromacs(
        settings,
        base_dir=work_dir,
        chain_rng_seed=42,
    )

    assert result.n_atoms == 13653
    assert result.coordinates_path.read_bytes() == ref_gro.read_bytes()

    parity = compare_topology_files(str(ref_top), str(result.topology_path))
    assert all(parity.values()), parity


@pytest.mark.integration
def test_rim135_topology_parity() -> None:
    """Generated hybrid topology matches bakery ref (section key sets)."""
    data_dir = rim135_dir()
    if data_dir is None:
        pytest.skip(f"rim135 fixtures not found (set {RIM135_ENV} or clone bakery beside hub)")

    ref_top = data_dir / "ref_hyb_topol.top"
    ref_gro = data_dir / "ref_hyb_conf.gro"
    if not ref_top.is_file() or not ref_gro.is_file():
        pytest.skip("rim135 reference outputs missing")

    for name in ("hyb_conf.gro", "hyb_topol.top", "missing_definitions.txt"):
        path = data_dir / name
        if path.exists():
            path.unlink()

    result = build_hybrid_gromacs(
        data_dir / "settings.xml",
        base_dir=data_dir,
        chain_rng_seed=42,
    )

    assert result.coordinates_path.read_bytes() == ref_gro.read_bytes()

    parity = compare_topology_files(str(ref_top), str(result.topology_path))
    assert all(parity.values()), parity


@pytest.mark.integration
def test_rim135_build_v2_lammps_smoke() -> None:
    """Settings v2 network path builds a LAMMPS data/input pair for rim135."""
    v2_yaml = Path(__file__).resolve().parents[2] / "examples" / "epoxy" / "settings.v2.yaml"
    if not v2_yaml.is_file():
        pytest.skip("settings.v2.yaml example missing")

    settings = load_settings(v2_yaml)
    work_dir = resolve_data_dir(v2_yaml, settings)
    if not work_dir.is_dir():
        pytest.skip("rim135 fixtures not available for settings.v2.yaml")

    prefix = settings.output.prefix
    data_path = work_dir / f"{prefix}.data"
    input_path = work_dir / f"in.{prefix}"
    for output_path in (data_path, input_path):
        if output_path.exists():
            output_path.unlink()

    result = build_network_lammps(settings, v2_yaml)
    write_lammps_data(result.system, data_path)
    if result.system.has_cross_pairs:
        pairs_path = work_dir / result.system.cross_pairs_file
        write_cross_pairs_file(result.system, pairs_path)
    write_lammps_input(result.system, settings, input_path, data_filename=data_path.name)

    table_search = [work_dir]
    tables_dir = resolve_tables_dir(v2_yaml, settings)
    if tables_dir is not None:
        table_search.append(tables_dir)
    convert_tables(result.system, settings, work_dir, extra_dirs=table_search[1:])

    assert data_path.is_file()
    assert input_path.is_file()
    assert data_path.stat().st_size > 0
    assert input_path.stat().st_size > 0

    atom_line = next(
        (line for line in data_path.read_text().splitlines() if line.strip().endswith(" atoms")),
        "",
    )
    atom_count = int(atom_line.split()[0]) if atom_line else 0
    assert 13_000 <= atom_count <= 14_000

    max_bond = max_bond_length(result.system.atoms, result.system.bonds, result.system.box)
    assert max_bond < 20.0, f"longest min-image bond {max_bond:.1f} Å after PBC prep"

    input_text = input_path.read_text()
    assert "angle_style hybrid backmap/harmonic backmap/table linear 1000" in input_text
    assert "angle_coeff" in input_text and "backmap/table cg table_a1.table" in input_text
    assert "angle_coeff" in input_text and "backmap/table cg table_a2.table" in input_text
    assert "cg 0.000000 0.0000" not in input_text
    assert (work_dir / "table_a1.table").is_file()
    assert (work_dir / "table_a2.table").is_file()
    table_angle_types = {
        line.split()[1]
        for line in input_text.splitlines()
        if line.startswith("angle_coeff") and "backmap/table cg" in line
    }
    assert len(table_angle_types) == 2

    dihedral_line = next(
        (
            line
            for line in data_path.read_text().splitlines()
            if line.strip().endswith(" dihedrals")
        ),
        "",
    )
    dihedral_count = int(dihedral_line.split()[0]) if dihedral_line else 0
    assert 33_000 <= dihedral_count <= 34_000, dihedral_line
    assert "dihedral_style hybrid" in input_text
    assert "backmap/ryckaert" in input_text
    assert "dihedral_coeff" in input_text

    comm_line = next(
        (line for line in input_text.splitlines() if line.startswith("comm_modify cutoff")),
        "",
    )
    assert comm_line
    comm_cutoff = float(comm_line.split()[-1])
    max_euclidean = max_euclidean_bond_length(result.system.atoms, result.system.bonds)
    assert comm_cutoff >= max(max_bond, max_euclidean), comm_line
    assert "reset_atoms image all" in input_text

    assert result.system.has_cross_pairs
    pairs_path = work_dir / result.system.cross_pairs_file
    assert pairs_path.is_file()
    pair_lines = pairs_path.read_text().splitlines()
    pair_count = int(pair_lines[0])
    assert 3_500 <= pair_count <= 3_700, pair_count
    assert len(pair_lines) == pair_count + 1
    assert "fix pairs all backmap/pairs at file pairs.dat" in input_text

    # Tier B bakery protocol (PR4): velocity init, cap_force, gamma=15, 1 fs dt, 10k ramp
    assert "velocity all create" in input_text
    assert "fix cap all backmap/capforce" in input_text
    assert "langevin 298.0" in input_text
    langevin_line = next(
        (line for line in input_text.splitlines() if "fix thermo at_atoms langevin" in line),
        "",
    )
    assert langevin_line
    damp = float(langevin_line.split()[-2])
    assert 66.0 <= damp <= 67.0, langevin_line
    assert "timestep 1.00" in input_text
    assert "run 10000" in input_text
    assert "fix_modify bm active yes" in input_text
    assert "fix_modify bm active no" not in input_text


@pytest.mark.integration
def test_rim135_build_cg_only_smoke() -> None:
    """CG-only network build: 700 beads, plain table FF, ~15 Å comm."""
    v2_yaml = Path(__file__).resolve().parents[2] / "examples" / "epoxy" / "settings.v2.yaml"
    if not v2_yaml.is_file():
        pytest.skip("settings.v2.yaml example missing")

    settings = load_settings(v2_yaml)
    work_dir = resolve_data_dir(v2_yaml, settings)
    if not work_dir.is_dir():
        pytest.skip("rim135 fixtures not available for settings.v2.yaml")

    from backmap_prep.network.lammps_builder import build_system_from_cg

    tables_dir = resolve_tables_dir(v2_yaml, settings)
    table_search = [tables_dir] if tables_dir is not None else []
    system = build_system_from_cg(settings, v2_yaml, table_search_dirs=table_search)

    assert len(system.atoms) == 700
    assert all(atom.is_cg for atom in system.atoms)
    assert system.bond_types
    assert all(bt.style == "table" for bt in system.bond_types if bt.table_file)
    assert not system.has_cross_pairs

    max_bond = max_bond_length(system.atoms, system.bonds, system.box)
    assert max_bond < 20.0, f"longest bonded distance {max_bond:.1f} Å after PBC prep"


@pytest.mark.integration
def test_rim135_finalize_cg_from_equil_smoke() -> None:
    """Finalize equilibrated CG frame: image flags + unwrapped GRO."""
    v2_yaml = Path(__file__).resolve().parents[2] / "examples" / "epoxy" / "settings.v2.yaml"
    if not v2_yaml.is_file():
        pytest.skip("settings.v2.yaml example missing")

    settings = load_settings(v2_yaml)
    work_dir = resolve_data_dir(v2_yaml, settings)
    equil_path = work_dir / "rim135_cg_equil.data"
    if not equil_path.is_file():
        pytest.skip("rim135_cg_equil.data missing (run CG equil first)")

    from backmap_prep.network.finalize_cg import finalize_cg_from_equil, write_cg_gro

    tables_dir = resolve_tables_dir(v2_yaml, settings)
    table_search = [tables_dir] if tables_dir is not None else []
    result = finalize_cg_from_equil(settings, v2_yaml, equil_path, table_search_dirs=table_search)

    assert result.max_bond_ang < 20.0
    assert any(atom.ix or atom.iy or atom.iz for atom in result.system.atoms)

    final_path = work_dir / "test_rim135_cg_final.data"
    gro_path = work_dir / "test_rim135_cg_equil.gro"
    try:
        from backmap_prep.writers import write_lammps_data

        write_lammps_data(result.system, final_path)
        atom_line = next(
            line for line in final_path.read_text().splitlines() if line.strip().endswith(" atoms")
        )
        assert int(atom_line.split()[0]) == 700
        atoms_section = final_path.read_text().split("Atoms # full", 1)[1].split("Bonds", 1)[0]
        assert any(len(line.split()) >= 10 for line in atoms_section.splitlines() if line.strip())

        template_gro = work_dir / settings.cg_system.coordinates
        write_cg_gro(result.system, template_gro, gro_path, unwrapped=True)
        assert gro_path.is_file()
        assert int(gro_path.read_text().splitlines()[1].strip()) == 700
    finally:
        if final_path.exists():
            final_path.unlink()
        if gro_path.exists():
            gro_path.unlink()


def test_rim135_rebuild_from_finalized_cg() -> None:
    """Rebuild hybrid from equilibrated CG: sane min-image PBC, network comm cutoff."""
    v2_yaml = Path(__file__).resolve().parents[2] / "examples" / "epoxy" / "settings.v2.yaml"
    if not v2_yaml.is_file():
        pytest.skip("settings.v2.yaml example missing")

    settings = load_settings(v2_yaml)
    work_dir = resolve_data_dir(v2_yaml, settings)
    cg_final = work_dir / "rim135_cg_final.data"
    if not cg_final.is_file():
        pytest.skip("rim135_cg_final.data missing (run finalize-cg first)")

    from backmap_prep.network.rebuild import rebuild_network_lammps
    from backmap_prep.writers import _compute_params

    result = rebuild_network_lammps(settings, v2_yaml, cg_final)
    system = result.system

    assert len(system.atoms) == 13653
    assert system.write_image_flags
    max_bond = max_bond_length(system.atoms, system.bonds, system.box)
    assert max_bond < 20.0, f"longest min-image bond {max_bond:.1f} Å after rebuild"
    assert sum(1 for atom in system.atoms if atom.ix or atom.iy or atom.iz) > 0

    max_euclidean = max_euclidean_bond_length(system.atoms, system.bonds)
    assert max_euclidean > max_bond, "network crosslinks should span the primary cell"

    params = _compute_params(system, settings)
    assert params["comm_cutoff_ang"] >= max_euclidean + 0.9, (
        f"comm cutoff {params['comm_cutoff_ang']:.1f} Å must cover box-spanning bonds "
        f"(max euclidean bond {max_euclidean:.1f} Å)"
    )
