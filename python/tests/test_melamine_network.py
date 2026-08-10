"""Tier A parity: crosslinked MF network hybrid topology generation.

Migrated (Task 2b) from bakery `settings.xml` passthrough mode to native
Settings v2 YAML (`settings.yaml`): the bead mapping, active sites, and
remove rules are now expressed directly as `molecules[].beads[].mapping`
(a faithful translation of bakery's vendored `settings.xml`, which is kept
around only as a passive, unread provenance artifact -- see git history at
commit 88f0f61 for the pristine reference), and the crosslink-site bond/
angle/dihedral parameters are supplied directly as real numeric values in
`cross_interactions` instead of bare atom-name patterns that relied on a
downstream GROMACS grompp step (which the LAMMPS-conversion path never
actually ran, silently producing zero-force-constant crosslink terms).
All builds below go through `load_settings()` + the native
`build_hybrid_gromacs(settings, ...)` path (`has_native_network_config`
true), matching how `uv run backmap-prep build-hybrid settings.yaml` and
the other v2 examples (PET, epoxy) already work.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

import pytest

from backmap_prep.network.api import build_hybrid_gromacs, build_network_lammps
from backmap_prep.schema import load_settings
from backmap_prep.table_converter import convert_tables
from backmap_prep.units import distance, gromacs_rb_to_lammps, spring_angle, spring_bond
from backmap_prep.writers import write_cross_pairs_file, write_lammps_data, write_lammps_input

MF_NETWORK_DIR = Path(__file__).resolve().parents[2] / "examples" / "melamine_network" / "large"
SETTINGS_YAML = MF_NETWORK_DIR / "settings.yaml"


def _build():
    """Build the hybrid GROMACS system from the native v2 settings.yaml."""
    settings = load_settings(SETTINGS_YAML)
    return build_hybrid_gromacs(settings, base_dir=MF_NETWORK_DIR, chain_rng_seed=42)


def test_mf_network_assets_present() -> None:
    """Vendored bakery assets exist before attempting a build.

    settings.xml is kept as a passive, unread provenance artifact -- the
    pipeline itself is now driven entirely by settings.yaml.
    """
    for name in (
        "settings.yaml",
        "settings.xml",
        "at_topol.top",
        "single_mf.gro",
        "cg_conf.gro",
        "cg_topol.top",
        "hyb_topol.top",
    ):
        assert (MF_NETWORK_DIR / name).is_file(), f"missing vendored asset: {name}"


def test_mf_network_build_hybrid_smoke() -> None:
    """Hybrid build succeeds and produces the crosslinked atom count."""
    result = _build()

    assert result.coordinates_path.is_file()
    assert result.topology_path.is_file()
    assert result.coordinates_path.stat().st_size > 0
    assert result.topology_path.stat().st_size > 0
    # 500 molecules x 27 AT atoms (unreacted baseline) = 13500 atoms.
    # settings.yaml's per-active-site `remove` lists delete 1 or 2
    # leaving-group atoms (the reacted arm's OH hydrogen, and -- for the
    # arm whose carbon is the active site -- its O too) at each of the 675
    # crosslink sites; the net effect nets out to exactly 525 atoms removed
    # system-wide, confirmed from this run's real output (13500 - 525 =
    # 12975). Identical to Task 3's bakery_xml-passthrough result -- the
    # underlying chemistry is unchanged, only how it's expressed changed.
    assert result.n_atoms == 12975


def test_mf_network_charge_conserved() -> None:
    """Total system charge is conserved within 1e-6 e."""
    result = _build()

    total_charge = 0.0
    in_atoms = False
    for line in result.topology_path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("[ atoms ]"):
            in_atoms = True
            continue
        if stripped.startswith("["):
            in_atoms = False
            continue
        if in_atoms and stripped and not stripped.startswith(";"):
            parts = stripped.split()
            if len(parts) >= 7:
                total_charge += float(parts[6])
    assert abs(total_charge) < 1e-6, f"charge not conserved: {total_charge}"


def test_mf_network_crosslink_bond_count() -> None:
    """Generated hybrid topology has exactly 675 inter-molecular AT bonds --
    one real covalent bond per CG-level crosslink (spec:
    example-melamine-network, 'Crosslink count matches bakery's reference').

    Detection method changed from Task 3's bakery_xml-passthrough test: that
    version searched for a literal "; AT cross bonds" comment tag, which was
    how bakery's XML declared the crosslink bond pattern (no inline numeric
    params -- delegated to a downstream grompp resolution step that this
    project's LAMMPS-conversion path never actually ran, i.e. the exact bug
    this migration fixes). The native v2 `cross_interactions.bonds` entry
    below carries real params (`1 0.143 267776.0`, sourced from
    at_topol.top's own intramolecular C-O bond) directly, so the generated
    `[ cross_bonds ]` lines have no comment tag at all -- they are now
    identified by their literal, untagged params string instead.
    """
    result = _build()

    text = result.topology_path.read_text()

    in_bonds = False
    n_bonds = 0
    in_cross_bonds = False
    n_cross_at_bonds = 0
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("[ bonds ]"):
            in_bonds = True
            continue
        if stripped.startswith("[ cross_bonds ]"):
            in_cross_bonds = True
            continue
        if stripped.startswith("["):
            in_bonds = False
            in_cross_bonds = False
            continue
        if in_bonds and stripped and not stripped.startswith(";"):
            n_bonds += 1
        if in_cross_bonds and stripped.endswith("1 0.143 267776.0"):
            n_cross_at_bonds += 1
    print(f"total AT bonds in [ bonds ] (intramolecular only): {n_bonds}")
    print(
        f"real AT crosslink bonds in [ cross_bonds ] (params '1 0.143 267776.0'): "
        f"{n_cross_at_bonds}"
    )

    # Byte-for-byte confirmed (during Task 2b diagnosis) that this build's
    # [ bonds ] section is identical to Task 3's bakery_xml-passthrough
    # result -- crosslinks were never going to show up there; they live
    # exclusively in [ cross_bonds ], matching the CG-level crosslink count
    # independently verified in Task 1 (`cg_topol.top`: chem=675).
    assert n_cross_at_bonds == 675


def test_mf_network_charges_match_bakery_charge_map() -> None:
    """Every per-atom charge value appearing in the crosslinked network's
    generated topology also appears among the values listed in the vendored
    settings.xml's own <charge_map> elements -- i.e. arms (reacted or not)
    are built from bakery's real, CG-mapping-consistent charge scheme, not
    at_topol.top's plain AA charges and not some separately perturbed set
    (spec: example-melamine-network, 'Unreacted arms match bakery's own
    reference charge scheme').

    Re-verified empirically for the native v2 pipeline during Task 2b: this
    is NOT automatic just because the pipeline changed source format --
    charge_map is only applied by structures.py (BackmapperSettings2) when
    the XML's <beads> element actually carries a <charge_map> child, and
    that child is only emitted (via v2_loader.py's `_emit_bead_extras`) if
    the corresponding `mapping` entry in settings.yaml explicitly sets
    `charge_map`. This migration's settings.yaml deliberately does carry
    `charge_map` on every degree=2/degree=3 mapping entry (copied verbatim
    from bakery's own settings.xml values), so bakery's real charge scheme
    is preserved. Confirmed by inspecting the freshly-built topology's
    `[ atoms ]` section directly: charges there are e.g. -0.407377 /
    0.682493 / -0.299731 / 0.024669 (bakery's charge_map for the A1 arm),
    NOT at_topol.top's native AA charges for the same atoms (e.g. N11's
    native charge is -0.901187, not -0.407377). Had settings.yaml omitted
    charge_map from the mapping entries, this test's premise would be false
    and the pipeline would use at_topol.top's native AA charges instead --
    this is a property of the authored YAML, not an automatic pipeline
    behavior, which is why it needed re-verification rather than being
    assumed unchanged.
    """
    xml_path = MF_NETWORK_DIR / "settings.xml"

    # Parse the ground truth directly from settings.xml (still a valid
    # reference for the expected *values* even though the pipeline no
    # longer reads this file -- settings.yaml's charge_map entries were
    # copied verbatim from it).
    charge_map_values: set[str] = set()
    root = ET.parse(xml_path).getroot()
    charge_map_elements = root.findall(".//charge_map")
    assert charge_map_elements, "settings.xml has no <charge_map> elements to check against"
    for charge_map in charge_map_elements:
        for token in (charge_map.text or "").split():
            charge_map_values.add(f"{float(token):.6f}")

    result = _build()
    generated_charges: set[str] = set()
    in_atoms = False
    for line in result.topology_path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("[ atoms ]"):
            in_atoms = True
            continue
        if stripped.startswith("["):
            in_atoms = False
            continue
        if in_atoms and stripped and not stripped.startswith(";"):
            parts = stripped.split()
            if len(parts) >= 7:
                generated_charges.add(f"{float(parts[6]):.6f}")

    unexpected = generated_charges - charge_map_values
    assert not unexpected, (
        f"generated topology has {len(unexpected)} charge values not present "
        f"in settings.xml's <charge_map> elements -- arms may have been built "
        f"from an unexpected charge source: {sorted(unexpected)[:10]}"
    )


def test_mf_network_no_missing_definitions() -> None:
    """Crosslink-site angle/dihedral terms are fully parameterized.

    Migrated (Task 2b) from Task 4's two-part fix (bare MF:-atom-name
    patterns registered in bakery's settings.xml + generic OPLS-AA
    angletype/dihedraltype entries in ffbonded.itp for GROMACS/grompp to
    resolve at preprocessing time) to native v2 `cross_interactions`
    entries that carry real numeric params directly -- no pattern-to-value
    resolution step is needed at all for this path, so
    missing_definitions.txt is empty for a structurally different reason
    than before (nothing was ever a "pattern" requiring resolution; the 9
    bond / 36 angle / 63 dihedral MF:-atom-name signatures all resolve
    immediately from settings.yaml's own literal params).
    """
    result = _build()
    if result.missing_definitions_path and result.missing_definitions_path.is_file():
        content = result.missing_definitions_path.read_text().strip()
        assert content == "", f"unresolved force-field gaps:\n{content}"


def test_mf_network_crosslink_angles_dihedrals_nonzero() -> None:
    """Crosslink-site angle/dihedral terms carry real, nonzero force
    constants in the generated topology -- the specific defect (Task 5)
    that motivated this migration away from bakery_xml passthrough: that
    mode combined with the LAMMPS-conversion step produced crosslink
    bonds/angles/dihedrals with zero force constants, so crosslinks had no
    physical effect in the simulation at all.
    """
    result = _build()
    text = result.topology_path.read_text()

    sections: dict[str, list[str]] = {"cross_angles": [], "cross_dihedrals": []}
    current: str | None = None
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("[ cross_angles ]"):
            current = "cross_angles"
            continue
        if stripped.startswith("[ cross_dihedrals ]"):
            current = "cross_dihedrals"
            continue
        if stripped.startswith("["):
            current = None
            continue
        if current and stripped:
            sections[current].append(stripped)

    # The 3 new crosslink angle param sets sourced in Task 2b/4.
    expected_angle_params = {
        "1 109.500 502.080",  # CT-OH-CT
        "1 109.500 292.880",  # HC-CT-OH
        "1 109.500 418.400",  # NT-CT-OH
    }
    found_angle_params = {
        " ".join(line.split()[3:]) for line in sections["cross_angles"] if len(line.split()) >= 6
    }
    missing_angles = expected_angle_params - found_angle_params
    assert not missing_angles, f"expected crosslink angle params not found: {missing_angles}"

    # The 4 distinct new crosslink dihedral param sets sourced in Task 2b/4
    # (CT-OH-CT-HC and CT-OH-CT-NT share the same numeric RB coefficients).
    expected_dihedral_params = {
        "3 1.58992 4.76976 0.00000 -6.35968 0.00000 0.00000",
        "3 1.78866 3.49154 0.53555 -5.81576 0.00000 0.00000",
        "3 -1.26775 3.02085 1.74473 -3.49782 0.00000 0.00000",
    }
    found_dihedral_params = {
        " ".join(line.split()[4:11])
        for line in sections["cross_dihedrals"]
        if len(line.split()) >= 11
    }
    missing_dihedrals = expected_dihedral_params - found_dihedral_params
    assert not missing_dihedrals, (
        f"expected crosslink dihedral params not found: {missing_dihedrals}"
    )


def _section_lines(topology_text: str, section: str) -> list[str]:
    """Non-empty lines of a `[ section ]` block in a GROMACS topology."""
    lines: list[str] = []
    on = False
    for line in topology_text.splitlines():
        stripped = line.strip()
        if stripped.startswith(f"[ {section} ]"):
            on = True
            continue
        if stripped.startswith("["):
            on = False
            continue
        if on and stripped:
            lines.append(stripped)
    return lines


def _typed_bonds(data_text: str) -> dict[frozenset[int], int]:
    """`Bonds` section of a LAMMPS `.data` file: {frozenset(atom pair): type}."""
    mapping: dict[frozenset[int], int] = {}
    on = False
    for line in data_text.splitlines():
        stripped = line.strip()
        if stripped == "Bonds":
            on = True
            continue
        if not on:
            continue
        if stripped == "":
            if mapping:
                break
            continue
        parts = stripped.split()
        if len(parts) != 4:
            break
        _, btype, a1, a2 = (int(x) for x in parts)
        mapping[frozenset((a1, a2))] = btype
    return mapping


def _typed_angles(data_text: str) -> dict[tuple[int, int, int], int]:
    """`Angles` section: {(min(i,k), j, max(i,k)): type} (i-j-k == k-j-i)."""
    mapping: dict[tuple[int, int, int], int] = {}
    on = False
    for line in data_text.splitlines():
        stripped = line.strip()
        if stripped == "Angles":
            on = True
            continue
        if not on:
            continue
        if stripped == "":
            if mapping:
                break
            continue
        parts = stripped.split()
        if len(parts) != 5:
            break
        _, atype, a1, a2, a3 = (int(x) for x in parts)
        mapping[(min(a1, a3), a2, max(a1, a3))] = atype
    return mapping


def _typed_dihedrals(data_text: str) -> dict[tuple[int, int, int, int], int]:
    """`Dihedrals` section: {min(fwd, reversed(fwd)): type} (i-j-k-l == l-k-j-i)."""
    mapping: dict[tuple[int, int, int, int], int] = {}
    on = False
    for line in data_text.splitlines():
        stripped = line.strip()
        if stripped == "Dihedrals":
            on = True
            continue
        if not on:
            continue
        if stripped == "":
            if mapping:
                break
            continue
        parts = stripped.split()
        if len(parts) != 6:
            break
        _, dtype, a1, a2, a3, a4 = (int(x) for x in parts)
        fwd = (a1, a2, a3, a4)
        mapping[min(fwd, tuple(reversed(fwd)))] = dtype
    return mapping


def test_mf_network_lammps_crosslink_types_nonzero(tmp_path: Path) -> None:
    """LAMMPS-level regression test (Task 5b): the crosslink bond/angle/
    dihedral LAMMPS types carry real, nonzero coefficients, converted from
    Task 2b's `cross_interactions` GROMACS params via the existing,
    already-tested `units.py` conversion (`spring_bond`/`spring_angle`/
    `gromacs_rb_to_lammps`) -- not re-derived here.

    This is the LAMMPS-layer counterpart to
    `test_mf_network_crosslink_angles_dihedrals_nonzero` /
    `test_mf_network_crosslink_bond_count` above, which only checked the
    GROMACS-level hybrid topology. Task 5's original attempt found that a
    correct-looking GROMACS topology could still be converted into a LAMMPS
    `.data`/`in.*` pair with `K=0` crosslink bond/angle/dihedral types (the
    bug this whole migration fixes), so this test regenerates the LAMMPS
    output end to end (`build_network_lammps` + the `write_lammps_*`
    writers, matching the CLI `build` command) and inspects the actual
    `bond_coeff`/`angle_coeff`/`dihedral_coeff` lines.

    LAMMPS type IDs are identified by exact atom-ID pair/triple/quadruple
    cross-reference between `hyb_topol.top`'s tagged `cross_bonds`/
    `cross_angles`/`cross_dihedrals` sections and the generated `.data`
    file's `Bonds`/`Angles`/`Dihedrals` sections -- NOT by eyeballing
    coefficient values. Task 5's original attempt initially misidentified
    an unrelated, same-K-looking intramolecular bond type this way before
    checking atom pairs/counts directly (see task-5-report.md); the same
    risk applies here since the crosslink ether C-O bond now shares its
    LAMMPS bond type with the pre-existing intramolecular ether C-O bond
    (same chemistry, same params -- see settings.yaml's cross_interactions
    comment), so the type is expected to carry more instances (1500) than
    just the 675 crosslink ones, and that is not itself a bug.
    """
    settings = load_settings(SETTINGS_YAML)
    result = build_network_lammps(settings, SETTINGS_YAML)

    data_path = tmp_path / f"{settings.output.prefix}.data"
    input_path = tmp_path / f"in.{settings.output.prefix}"
    write_lammps_data(result.system, data_path)
    if result.system.has_cross_pairs:
        write_cross_pairs_file(result.system, tmp_path / result.system.cross_pairs_file)
    write_lammps_input(result.system, settings, input_path, data_filename=data_path.name)
    convert_tables(result.system, settings, tmp_path, extra_dirs=[MF_NETWORK_DIR])

    topology_text = result.topology_path.read_text()
    data_text = data_path.read_text()
    input_text = input_path.read_text()

    # --- Step 1: crosslink atom pairs/triples/quadruples tagged by Task 2b's
    # exact cross_interactions param strings, read straight from the fresh
    # hyb_topol.top (GROMACS level -- already independently verified by
    # Task 2b and its reviewer; re-extracted here only as the join key). ---
    bond_param = "1 0.143 267776.0"
    angle_params = {
        "CT-OH-CT": "1 109.500 502.080",
        "HC-CT-OH": "1 109.500 292.880",
        "NT-CT-OH": "1 109.500 418.400",
    }
    dihedral_params = {
        "CT-OH-CT-HC_or_NT": "3 1.58992 4.76976 0.00000 -6.35968 0.00000 0.00000",
        "OH-CT-NT-CA": "3 1.78866 3.49154 0.53555 -5.81576 0.00000 0.00000",
        "OH-CT-NT-H": "3 -1.26775 3.02085 1.74473 -3.49782 0.00000 0.00000",
    }

    cross_bond_pairs: list[frozenset[int]] = []
    for line in _section_lines(topology_text, "cross_bonds"):
        if line.endswith(bond_param):
            i, j = (int(x) for x in line.split()[:2])
            cross_bond_pairs.append(frozenset((i, j)))
    assert len(cross_bond_pairs) == 675, len(cross_bond_pairs)
    assert len(set(cross_bond_pairs)) == 675, "duplicate crosslink bond atom pairs"

    cross_angle_triples: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for line in _section_lines(topology_text, "cross_angles"):
        for label, params in angle_params.items():
            if line.endswith(params):
                i, j, k = (int(x) for x in line.split()[:3])
                cross_angle_triples[label].append((i, j, k))
    n_cross_angles = sum(len(v) for v in cross_angle_triples.values())
    assert n_cross_angles == 2700, n_cross_angles

    cross_dihedral_quads: dict[str, list[tuple[int, int, int, int]]] = defaultdict(list)
    for line in _section_lines(topology_text, "cross_dihedrals"):
        core = line.split(";")[0].strip()
        for label, params in dihedral_params.items():
            if core.endswith(params):
                i, j, k, m = (int(x) for x in core.split()[:4])
                cross_dihedral_quads[label].append((i, j, k, m))
    n_cross_dihedrals = sum(len(v) for v in cross_dihedral_quads.values())
    assert n_cross_dihedrals == 5400, n_cross_dihedrals

    # --- Step 2: cross-reference those atom-ID pairs/triples/quadruples
    # against the generated LAMMPS .data file's Bonds/Angles/Dihedrals
    # sections by exact match + instance count, to find the LAMMPS type
    # ID(s) that actually carry the crosslink interactions. ---
    bond_type_by_pair = _typed_bonds(data_text)
    matched_bond_types = Counter(
        bond_type_by_pair[p] for p in cross_bond_pairs if p in bond_type_by_pair
    )
    assert sum(matched_bond_types.values()) == 675, (
        f"not all 675 crosslink bond atom pairs were found in the .data Bonds "
        f"section: {matched_bond_types}"
    )
    assert len(matched_bond_types) == 1, (
        f"crosslink bond atom pairs map to more than one LAMMPS bond type "
        f"(expected exactly one): {matched_bond_types}"
    )
    crosslink_bond_type = next(iter(matched_bond_types))

    angle_type_by_triple = _typed_angles(data_text)
    crosslink_angle_types: dict[str, int] = {}
    for label, triples in cross_angle_triples.items():
        keys = [(min(i, k), j, max(i, k)) for (i, j, k) in triples]
        matched = Counter(angle_type_by_triple[key] for key in keys if key in angle_type_by_triple)
        assert sum(matched.values()) == len(triples), (
            f"{label}: not all {len(triples)} crosslink angle triples were "
            f"found in the .data Angles section: {matched}"
        )
        assert len(matched) == 1, (
            f"{label}: crosslink angle triples map to more than one LAMMPS angle type: {matched}"
        )
        crosslink_angle_types[label] = next(iter(matched))

    dihedral_type_by_quad = _typed_dihedrals(data_text)
    crosslink_dihedral_types: dict[str, int] = {}
    for label, quads in cross_dihedral_quads.items():
        keys = [min(q, tuple(reversed(q))) for q in quads]
        matched = Counter(
            dihedral_type_by_quad[key] for key in keys if key in dihedral_type_by_quad
        )
        assert sum(matched.values()) == len(quads), (
            f"{label}: not all {len(quads)} crosslink dihedral quadruples "
            f"were found in the .data Dihedrals section: {matched}"
        )
        assert len(matched) == 1, (
            f"{label}: crosslink dihedral quadruples map to more than one "
            f"LAMMPS dihedral type: {matched}"
        )
        crosslink_dihedral_types[label] = next(iter(matched))

    # --- Step 3: the identified types' bond_coeff/angle_coeff/dihedral_coeff
    # lines in in.melamine_network are real and nonzero, and match Task 2b's
    # declared GROMACS params converted to LAMMPS real units via units.py
    # (spring_bond halves+converts kJ/(mol*nm^2) -> kcal/(mol*Angstrom^2);
    # spring_angle halves+converts kJ/(mol*rad^2) -> kcal/(mol*rad^2);
    # gromacs_rb_to_lammps alternates sign and converts kJ/mol -> kcal/mol
    # per RB coefficient). ---
    bond_coeff_re = re.compile(
        r"^bond_coeff\s+(\d+)\s+backmap/harmonic\s+at\s+([-\d.]+)\s+([-\d.]+)\s*$",
        re.MULTILINE,
    )
    bond_coeffs = {
        int(m.group(1)): (float(m.group(2)), float(m.group(3)))
        for m in bond_coeff_re.finditer(input_text)
    }
    got_bond_k, got_bond_r0 = bond_coeffs[crosslink_bond_type]
    assert got_bond_k != 0.0, (
        f"crosslink bond type {crosslink_bond_type} has K=0 -- the Task 5 defect"
    )
    assert got_bond_k == pytest.approx(spring_bond(267776.0), rel=1e-5), (
        crosslink_bond_type,
        got_bond_k,
    )
    assert got_bond_r0 == pytest.approx(distance(0.143), rel=1e-5), (
        crosslink_bond_type,
        got_bond_r0,
    )

    angle_coeff_re = re.compile(
        r"^angle_coeff\s+(\d+)\s+at\s+([-\d.]+)\s+([-\d.]+)\s*$", re.MULTILINE
    )
    angle_coeffs = {
        int(m.group(1)): (float(m.group(2)), float(m.group(3)))
        for m in angle_coeff_re.finditer(input_text)
    }
    expected_angle_k = {
        "CT-OH-CT": spring_angle(502.080),
        "HC-CT-OH": spring_angle(292.880),
        "NT-CT-OH": spring_angle(418.400),
    }
    for label, type_id in crosslink_angle_types.items():
        got_k, got_theta0 = angle_coeffs[type_id]
        assert got_k != 0.0, f"{label} (angle type {type_id}) has K=0"
        assert got_k == pytest.approx(expected_angle_k[label], rel=1e-5), (
            label,
            type_id,
            got_k,
        )
        assert got_theta0 == pytest.approx(109.5, abs=1e-3), (label, type_id, got_theta0)

    dihedral_coeff_re = re.compile(
        r"^dihedral_coeff\s+(\d+)\s+backmap/ryckaert\s+at\s+"
        r"([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$",
        re.MULTILINE,
    )
    dihedral_coeffs = {
        int(m.group(1)): [float(m.group(i)) for i in range(2, 8)]
        for m in dihedral_coeff_re.finditer(input_text)
    }
    expected_dihedral_c = {
        "CT-OH-CT-HC_or_NT": gromacs_rb_to_lammps(
            [1.58992, 4.76976, 0.00000, -6.35968, 0.00000, 0.00000]
        ),
        "OH-CT-NT-CA": gromacs_rb_to_lammps(
            [1.78866, 3.49154, 0.53555, -5.81576, 0.00000, 0.00000]
        ),
        "OH-CT-NT-H": gromacs_rb_to_lammps(
            [-1.26775, 3.02085, 1.74473, -3.49782, 0.00000, 0.00000]
        ),
    }
    for label, type_id in crosslink_dihedral_types.items():
        got = dihedral_coeffs[type_id]
        assert any(abs(c) > 1e-6 for c in got), (
            f"{label} (dihedral type {type_id}) is all-zero: {got}"
        )
        expected = expected_dihedral_c[label]
        for got_c, exp_c in zip(got, expected, strict=False):
            assert got_c == pytest.approx(exp_c, abs=1e-3), (label, type_id, got, expected)
