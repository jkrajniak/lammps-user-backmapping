"""Convert Settings format v2 YAML into bakery XML for the network engine."""

from __future__ import annotations

import xml.etree.ElementTree as etree

from backmap_prep.schema import (
    BeadDef,
    BeadMappingEntry,
    CrossAngle,
    CrossBond,
    DegreeSourceFile,
    MoleculeDef,
    Settings,
)


def _sub(parent: etree.Element, tag: str, text: str | None = None) -> etree.Element:
    elem = etree.SubElement(parent, tag)
    if text is not None:
        elem.text = text
    return elem


def _flatten_atoms(atoms: list[str]) -> list[str]:
    """Normalize YAML atom lists (handles inline comma groups)."""
    flat: list[str] = []
    for item in atoms:
        for part in item.replace(",", " ").split():
            if part:
                flat.append(part)
    return flat


def _prefix_atoms(atoms: list[str], ident: str) -> list[str]:
    return [f"1:{ident}:{atom}" for atom in _flatten_atoms(atoms)]


def _resolve_mapping(
    entries: list[BeadMappingEntry],
    ident: str,
) -> list[tuple[BeadMappingEntry, list[str]]]:
    """Resolve delta mapping entries to prefixed atom lists."""
    by_degree: dict[int, list[str]] = {}
    resolved: list[tuple[BeadMappingEntry, list[str]]] = []
    for entry in sorted(entries, key=lambda item: item.degree):
        if entry.atoms is not None:
            atoms = _flatten_atoms(list(entry.atoms))
        elif entry.base is not None:
            if entry.base not in by_degree:
                raise ValueError(
                    f"mapping degree {entry.degree} references missing base degree {entry.base}"
                )
            atoms = list(by_degree[entry.base])
            for atom in entry.remove:
                if atom in atoms:
                    atoms.remove(atom)
            atoms.extend(entry.add)
        else:
            raise ValueError(f"mapping degree {entry.degree} must define atoms or base")
        by_degree[entry.degree] = atoms
        resolved.append((entry, _prefix_atoms(atoms, ident)))
    return resolved


def _active_site_attr(entry: BeadMappingEntry, ident: str) -> str | None:
    if not entry.active_sites:
        return None
    parts = [f"{ident}:{site.atom}:{site.max_degree}" for site in entry.active_sites]
    return " ".join(parts)


def _format_beads_text(atoms: list[str]) -> str:
    """Match bakery XML wrapping (~10 atoms per line)."""
    lines: list[str] = []
    for idx in range(0, len(atoms), 10):
        lines.append(" ".join(atoms[idx : idx + 10]))
    return "\n                ".join(lines)


def _append_source_files(
    parent: etree.Element,
    tag: str,
    source: str | list[DegreeSourceFile],
) -> None:
    container = _sub(parent, tag)
    if isinstance(source, str):
        container.text = source
        return
    for entry in source:
        attrs: dict[str, str] = {"molecule_degree": str(entry.molecule_degree)}
        if entry.when:
            attrs["when"] = entry.when
        file_elem = etree.SubElement(container, "file", attrib=attrs)
        file_elem.text = entry.file


def _emit_cg_bead(parent: etree.Element, bead: BeadDef, ident: str) -> None:
    cg_bead = _sub(parent, "cg_bead")
    _sub(cg_bead, "name", bead.name)
    _sub(cg_bead, "type", bead.type)

    if bead.mapping:
        carried_active: str | None = None
        for entry, atoms in _resolve_mapping(bead.mapping, ident):
            attrs: dict[str, str] = {"degree": str(entry.degree)}
            if entry.molecule_degree:
                attrs["molecule_degree"] = entry.molecule_degree
            active_site = _active_site_attr(entry, ident) or carried_active
            if entry.active_sites:
                carried_active = active_site
            if active_site:
                attrs["active_site"] = active_site
            beads_elem = etree.SubElement(cg_bead, "beads", attrib=attrs)
            beads_elem.text = f"\n                {_format_beads_text(atoms)}\n            "
        return

    if bead.atoms:
        beads_elem = etree.SubElement(cg_bead, "beads", attrib={"degree": "*"})
        beads_elem.text = " ".join(_prefix_atoms(bead.atoms, ident))

    if bead.atoms_by_degree:
        for abd in bead.atoms_by_degree:
            attrs: dict[str, str] = {"degree": str(abd.degree)}
            if abd.molecule_degree:
                attrs["molecule_degree"] = abd.molecule_degree
            if abd.active_site:
                attrs["active_site"] = abd.active_site
            beads_elem = etree.SubElement(cg_bead, "beads", attrib=attrs)
            beads_elem.text = " ".join(abd.atoms)


def _emit_cg_molecule(parent: etree.Element, mol: MoleculeDef) -> None:
    ident = mol.ident or mol.name
    cg_mol = _sub(parent, "cg_molecule")
    _sub(cg_mol, "name", mol.name)
    _sub(cg_mol, "ident", ident)
    _append_source_files(cg_mol, "source_coordinate", mol.source.coordinates)
    _append_source_files(cg_mol, "source_topology", mol.source.topology)
    for bead in mol.beads:
        _emit_cg_bead(cg_mol, bead, ident)


def _emit_cross_bonds(parent: etree.Element, bonds: list[CrossBond]) -> None:
    for bond in bonds:
        bonds_elem = etree.SubElement(parent, "bonds", attrib={"params": bond.params or ";"})
        lines = [f"{pair[0]} {pair[1]}" for pair in bond.pairs]
        bonds_elem.text = "\n            " + "\n            ".join(lines) + "\n        "


def _emit_cross_angles(parent: etree.Element, angles: list[CrossAngle]) -> None:
    for angle in angles:
        angles_elem = etree.SubElement(parent, "angles", attrib={"params": angle.params})
        lines = [f"{triple[0]} {triple[1]} {triple[2]}" for triple in angle.triples]
        angles_elem.text = "\n" + "\n".join(lines) + "\n        "


def settings_to_xml_root(settings: Settings) -> etree.Element:
    """Build a bakery-compatible <settings> element tree from v2 Settings."""
    if settings.cg_system is None or settings.hybrid is None:
        raise ValueError("network v2 settings require cg_system and hybrid sections")

    root = etree.Element("settings")

    for mol in settings.molecules:
        _emit_cg_molecule(root, mol)

    cg_conf = _sub(root, "cg_configuration")
    _sub(cg_conf, "format", "GROMACS")
    _sub(cg_conf, "file", settings.cg_system.coordinates)
    _sub(cg_conf, "topology", settings.cg_system.topology)
    if settings.cg_system.predefined_active_sites:
        _sub(cg_conf, "predefined_active_sites", settings.cg_system.predefined_active_sites)

    hyb_conf = _sub(root, "hybrid_configuration")
    _sub(hyb_conf, "file", settings.hybrid.coordinates)
    _sub(hyb_conf, "format", "GROMACS")

    hyb_top = _sub(root, "hybrid_topology")
    _sub(hyb_top, "file", settings.hybrid.topology)
    _emit_cross_bonds(hyb_top, settings.cross_interactions.bonds)
    _emit_cross_angles(hyb_top, settings.cross_interactions.angles)

    mol_type = _sub(hyb_top, "molecule_type")
    _sub(mol_type, "name", settings.hybrid.molecule_type_name)
    _sub(mol_type, "exclusion", str(settings.hybrid.exclusion))

    if settings.hybrid.topology_includes:
        include_text = " ".join(
            f'"{inc}"' if " " in inc else inc for inc in settings.hybrid.topology_includes
        )
        _sub(hyb_top, "include", include_text)

    _sub(hyb_top, "system", settings.hybrid.system_name)
    return root


def has_native_network_config(settings: Settings) -> bool:
    """True when settings carry a full v2 network definition (no bakery_xml)."""
    return (
        settings.prep.engine == "network"
        and settings.prep.bakery_xml is None
        and bool(settings.molecules)
        and settings.cg_system is not None
        and settings.hybrid is not None
    )
