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

import xml.etree.ElementTree as ET
from pathlib import Path

from backmap_prep.network.api import build_hybrid_gromacs
from backmap_prep.schema import load_settings

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
