"""Command-line interface for backmap-prep."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable

from .builder import CGPositionOverride, System, build_system
from .network.api import build_hybrid_gromacs, build_network_lammps
from .network.compare_topology import main as compare_topology_main
from .network.finalize_cg import finalize_cg_from_equil, write_cg_gro
from .network.lammps_builder import build_system_from_cg
from .network.rebuild import rebuild_network_lammps
from .parsers import parse_lammps_data
from .schema import Settings, load_settings, resolve_data_dir, resolve_tables_dir
from .table_converter import convert_tables
from .writers import write_cross_pairs_file, write_lammps_data, write_lammps_input


def _cmd_build(args: argparse.Namespace) -> int:
    """Default build: generate hybrid data + input from GROMACS sources."""
    settings = load_settings(args.settings)
    prefix = args.output_prefix or settings.output.prefix
    out_dir = args.settings.parent
    if settings.prep.engine == "network":
        out_dir = resolve_data_dir(args.settings, settings)
        result = build_network_lammps(settings, args.settings)
        system = result.system
    else:
        system = build_system(settings, base_dir=out_dir, settings_path=args.settings)

    data_path = out_dir / f"{prefix}.data"
    write_lammps_data(system, data_path)
    print(f"Wrote {data_path}")

    if system.has_cross_pairs:
        pairs_path = out_dir / system.cross_pairs_file
        write_cross_pairs_file(system, pairs_path)
        print(f"Wrote {pairs_path}")

    input_path = out_dir / f"in.{prefix}"
    write_lammps_input(system, settings, input_path, data_filename=f"{prefix}.data")
    print(f"Wrote {input_path}")

    table_search: list[Path] = [out_dir]
    tables_dir = resolve_tables_dir(args.settings, settings)
    if tables_dir is not None:
        table_search.append(tables_dir)
    table_files = convert_tables(system, settings, out_dir, extra_dirs=table_search[1:])
    for tf in table_files:
        print(f"Wrote {tf}")

    return 0


def _cmd_rebuild(args: argparse.Namespace) -> int:
    """Rebuild hybrid data from equilibrated CG coordinates (.data or unwrapped .gro)."""
    settings = load_settings(args.settings)
    prefix = args.output_prefix or settings.output.prefix
    out_dir = args.settings.parent
    table_search: list[Path] = [out_dir]

    if settings.prep.engine == "network":
        out_dir = resolve_data_dir(args.settings, settings)
        table_search = [out_dir]
        tables_dir = resolve_tables_dir(args.settings, settings)
        if tables_dir is not None:
            table_search.append(tables_dir)
        result = rebuild_network_lammps(settings, args.settings, args.cg_frame.resolve())
        system = result.system
    else:
        cg_frame = parse_lammps_data(args.cg_frame)
        positions = [(a.x, a.y, a.z) for a in sorted(cg_frame.atoms, key=lambda a: a.atom_id)]
        box = (float(cg_frame.box[0]), float(cg_frame.box[1]), float(cg_frame.box[2]))
        override = CGPositionOverride(positions=positions, box=box)
        system = build_system(
            settings,
            base_dir=out_dir,
            cg_override=override,
            settings_path=args.settings,
        )

    data_path = out_dir / f"{prefix}.data"
    write_lammps_data(system, data_path)
    print(f"Wrote {data_path}")

    if system.has_cross_pairs:
        pairs_path = out_dir / system.cross_pairs_file
        write_cross_pairs_file(system, pairs_path)
        print(f"Wrote {pairs_path}")

    input_path = out_dir / f"in.{prefix}"
    write_lammps_input(system, settings, input_path, data_filename=f"{prefix}.data")
    print(f"Wrote {input_path}")

    table_files = convert_tables(system, settings, out_dir, extra_dirs=table_search[1:])
    for tf in table_files:
        print(f"Wrote {tf}")

    return 0


def _cmd_cg_only(args: argparse.Namespace) -> int:
    """Extract CG-only system for pre-equilibration."""
    settings = load_settings(args.settings)
    prefix = args.output_prefix or settings.output.prefix
    out_dir = args.settings.parent
    table_search: list[Path] = []

    if settings.prep.engine == "network":
        if settings.cg_system is None:
            print("network cg-only requires cg_system in settings", file=sys.stderr)
            return 1
        out_dir = resolve_data_dir(args.settings, settings)
        table_search = [out_dir]
        tables_dir = resolve_tables_dir(args.settings, settings)
        if tables_dir is not None:
            table_search.append(tables_dir)
        cg_system = build_system_from_cg(
            settings,
            args.settings,
            table_search_dirs=table_search[1:],
        )
    else:
        system = build_system(settings, base_dir=out_dir, settings_path=args.settings)

        cg_atoms = [a for a in system.atoms if a.is_cg]
        cg_atom_ids = {a.atom_id for a in cg_atoms}
        cg_bonds = [b for b in system.bonds if b.i in cg_atom_ids and b.j in cg_atom_ids]
        cg_angles = [
            a
            for a in system.angles
            if a.i in cg_atom_ids and a.j in cg_atom_ids and a.k in cg_atom_ids
        ]

        cg_system = System(
            atoms=cg_atoms,
            bonds=cg_bonds,
            angles=cg_angles,
            atom_types=[at for at in system.atom_types if at.is_cg],
            bond_types=[
                bt
                for bt in system.bond_types
                if bt.keyword == "cg" or bt.style in {"harmonic", "table"}
            ],
            angle_types=[
                at
                for at in system.angle_types
                if at.keyword == "cg" or at.style in {"harmonic", "table"}
            ],
            pair_types=[pt for pt in system.pair_types if pt.kind == "cg"],
            box=system.box,
            table_files=system.table_files,
            angle_table_files=system.angle_table_files,
            pair_table_files=system.pair_table_files,
        )

    data_path = out_dir / f"{prefix}_cg.data"
    write_lammps_data(cg_system, data_path)
    print(f"Wrote {data_path}")

    equil_path = out_dir / f"in.{prefix}_cg_equil"
    _write_cg_equilibration(cg_system, settings, equil_path, f"{prefix}_cg.data")
    print(f"Wrote {equil_path}")

    table_files = convert_tables(
        cg_system,
        settings,
        out_dir,
        extra_dirs=table_search[1:] if table_search else None,
    )
    for tf in table_files:
        print(f"Wrote {tf}")

    return 0


def _write_cg_equilibration(
    system: System,
    settings: Settings,
    path: Path,
    data_filename: str,
) -> None:
    """Write a LAMMPS input script for CG-only equilibration."""
    from . import units

    sim = settings.simulation
    temp = sim.temperature
    tdamp_fs = units.time(sim.thermostat_tdamp) if sim.thermostat_tdamp > 0 else 100.0
    lj_cut_ang = units.distance(sim.lj_cutoff)
    cg_cut_ang = units.distance(sim.cg_cutoff)
    comm_cutoff_ang = max(lj_cut_ang, cg_cut_ang) + 1.0

    omp_threads = sim.omp_threads

    with open(path, "w") as f:
        f.write("# CG equilibration — generated by backmap-prep cg-only\n\n")
        if omp_threads > 0:
            f.write(f"package omp {omp_threads}\n")
            f.write("suffix omp\n\n")
        f.write("units real\n")
        f.write("atom_style full\n")
        f.write("boundary p p p\n\n")
        f.write(f"read_data {data_filename}\n\n")
        f.write(f"comm_modify cutoff {comm_cutoff_ang:.2f}\n\n")

        f.write("pair_style table linear 1000\n")
        for pt in system.pair_types:
            if pt.table_file:
                f.write(f"pair_coeff {pt.itype} {pt.jtype} {pt.table_file} {pt.table_keyword}\n")
            else:
                f.write(f"pair_coeff {pt.itype} {pt.jtype} 0.0 0.0\n")
        f.write("\n")

        has_bond_table = any(bt.table_file for bt in system.bond_types)
        has_bond_harmonic = any(bt.style == "harmonic" and bt.params for bt in system.bond_types)
        if has_bond_table and has_bond_harmonic:
            f.write("bond_style hybrid harmonic table linear 1000\n")
        elif has_bond_table:
            f.write("bond_style table linear 1000\n")
        elif has_bond_harmonic:
            f.write("bond_style harmonic\n")

        for bt in system.bond_types:
            if bt.table_file:
                if has_bond_table and has_bond_harmonic:
                    f.write(f"bond_coeff {bt.type_id} table {bt.table_file} {bt.table_keyword}\n")
                else:
                    f.write(f"bond_coeff {bt.type_id} {bt.table_file} {bt.table_keyword}\n")
            elif bt.style == "harmonic" and bt.params:
                if has_bond_table and has_bond_harmonic:
                    f.write(
                        f"bond_coeff {bt.type_id} harmonic {bt.params[0]:.6f} {bt.params[1]:.6f}\n"
                    )
                else:
                    f.write(f"bond_coeff {bt.type_id} {bt.params[0]:.6f} {bt.params[1]:.6f}\n")
        f.write("\n")

        has_angle_table = any(at.table_file for at in system.angle_types)
        has_angle_harmonic = any(at.style == "harmonic" and at.params for at in system.angle_types)
        if has_angle_table or has_angle_harmonic:
            if has_angle_table and has_angle_harmonic:
                f.write("angle_style hybrid harmonic table linear 1000\n")
            elif has_angle_table:
                f.write("angle_style table linear 1000\n")
            else:
                f.write("angle_style harmonic\n")

            for at in system.angle_types:
                if at.table_file:
                    if has_angle_table and has_angle_harmonic:
                        f.write(
                            f"angle_coeff {at.type_id} table {at.table_file} {at.table_keyword}\n"
                        )
                    else:
                        f.write(f"angle_coeff {at.type_id} {at.table_file} {at.table_keyword}\n")
                elif at.style == "harmonic" and at.params:
                    if has_angle_table and has_angle_harmonic:
                        f.write(
                            f"angle_coeff {at.type_id} harmonic "
                            f"{at.params[0]:.6f} {at.params[1]:.4f}\n"
                        )
                    else:
                        f.write(f"angle_coeff {at.type_id} {at.params[0]:.6f} {at.params[1]:.4f}\n")
            f.write("\n")

        f.write("special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0\n")
        f.write("neigh_modify delay 0 every 1 check yes\n")
        f.write("reset_atoms image all\n\n")

        if sim.ensemble == "npt":
            p_atm = units.pressure(sim.pressure)
            pdamp = units.time(sim.barostat_pdamp)
            f.write(
                f"fix equil all npt temp {temp:.1f} {temp:.1f} {tdamp_fs:.1f} "
                f"iso {p_atm:.1f} {p_atm:.1f} {pdamp:.1f}\n\n"
            )
        else:
            f.write(f"fix equil all nvt temp {temp:.1f} {temp:.1f} {tdamp_fs:.1f}\n\n")

        ts_fs = units.time(sim.timestep)
        f.write(f"velocity all create {temp:.1f} {48279} mom yes rot yes dist gaussian\n\n")
        f.write(f"timestep {ts_fs:.2f}\n")
        f.write(f"thermo {sim.energy_interval}\n")
        f.write("thermo_style custom step temp pe ke etotal press\n\n")

        equil_steps = sim.equilibration_steps if sim.equilibration_steps > 0 else 10000
        f.write(f"run {equil_steps}\n\n")
        f.write("reset_atoms image all\n")
        f.write(f"write_data {data_filename.replace('.data', '_equil.data')}\n")


def _cmd_finalize_cg(args: argparse.Namespace) -> int:
    """Fix PBC/image flags on an equilibrated CG LAMMPS frame."""
    settings = load_settings(args.settings)
    prefix = args.output_prefix or settings.output.prefix
    out_dir = resolve_data_dir(args.settings, settings)

    table_search: list[Path] = [out_dir]
    tables_dir = resolve_tables_dir(args.settings, settings)
    if tables_dir is not None:
        table_search.append(tables_dir)

    result = finalize_cg_from_equil(
        settings,
        args.settings,
        args.equil.resolve(),
        table_search_dirs=table_search[1:],
    )

    data_path = out_dir / f"{prefix}_cg_final.data"
    write_lammps_data(result.system, data_path)
    print(f"Wrote {data_path}")
    print(
        f"PBC ok: max bonded distance {result.max_bond_ang:.2f} Å, "
        f"{result.nonzero_image_flags} atoms with non-zero image flags"
    )

    if not args.no_gro:
        if settings.cg_system is None:
            print("Warning: cg_system missing; skipped GRO export", file=sys.stderr)
        else:
            template_gro = (out_dir / settings.cg_system.coordinates).resolve()
            gro_name = (
                args.gro_name or f"{settings.cg_system.coordinates.replace('.gro', '_equil.gro')}"
            )
            gro_path = out_dir / gro_name
            write_cg_gro(result.system, template_gro, gro_path, unwrapped=not args.folded_gro)
            label = "folded" if args.folded_gro else "unwrapped"
            print(f"Wrote {gro_path} ({label} coordinates, nm)")

    return 0


def _cmd_migrate_settings(args: argparse.Namespace) -> int:
    """Convert legacy bakery XML settings to v2 YAML (stub)."""
    print(
        f"migrate-settings is not implemented yet ({args.input_xml} → {args.output_yaml})",
        file=sys.stderr,
    )
    return 2


_KNOWN_COMMANDS = {
    "build",
    "rebuild",
    "cg-only",
    "finalize-cg",
    "build-hybrid",
    "compare-topology",
    "migrate-settings",
}


def _resolve_bakery_xml(settings: Settings, settings_path: Path) -> Path:
    """Resolve bakery XML path from settings or fall back to a sibling .xml file."""
    if settings.prep.bakery_xml:
        xml_path = (settings_path.parent / settings.prep.bakery_xml).resolve()
        if not xml_path.is_file():
            raise FileNotFoundError(f"bakery XML not found: {xml_path}")
        return xml_path
    sibling = settings_path.with_suffix(".xml")
    if sibling.is_file():
        return sibling.resolve()
    raise ValueError(
        f"network prep requires prep.bakery_xml in settings.yaml or a sibling file {sibling.name}"
    )


def _cmd_build_hybrid(args: argparse.Namespace) -> int:
    """Build GROMACS hybrid coordinate + topology (Phase 3 network prep)."""
    settings_path: Path = args.settings

    if settings_path.suffix.lower() in {".yaml", ".yml"}:
        settings = load_settings(settings_path)
        allow_no_bonds = settings.prep.allow_no_bonds
        chain_rng_seed = settings.prep.chain_rng_seed
        work_dir = resolve_data_dir(settings_path, settings)
        if settings.prep.bakery_xml:
            xml_path = _resolve_bakery_xml(settings, settings_path)
            result = build_hybrid_gromacs(
                xml_path,
                base_dir=xml_path.parent,
                allow_no_bonds=allow_no_bonds,
                chain_rng_seed=chain_rng_seed,
            )
        else:
            result = build_hybrid_gromacs(
                settings,
                base_dir=work_dir,
                allow_no_bonds=allow_no_bonds,
                chain_rng_seed=chain_rng_seed,
            )
    elif settings_path.suffix.lower() == ".xml":
        xml_path = settings_path.resolve()
        allow_no_bonds = args.allow_no_bonds
        chain_rng_seed = args.chain_rng_seed
        result = build_hybrid_gromacs(
            xml_path,
            base_dir=xml_path.parent,
            allow_no_bonds=allow_no_bonds,
            chain_rng_seed=chain_rng_seed,
        )
    else:
        print(
            f"Error: unsupported settings format: {settings_path} (use .yaml or bakery .xml)",
            file=sys.stderr,
        )
        return 1

    print(f"Wrote {result.coordinates_path}")
    print(f"Wrote {result.topology_path}")
    print(f"Hybrid atoms: {result.n_atoms}")
    if result.missing_definitions_path is not None:
        print(f"Wrote {result.missing_definitions_path} (review missing cross terms)")
    return 0


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("settings", type=Path, help="YAML settings file")
    parser.add_argument("--output-prefix", type=str, default=None, help="Override output.prefix")


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    # Backward compatibility: if first arg is not a subcommand, assume "build"
    if argv and argv[0] not in _KNOWN_COMMANDS and not argv[0].startswith("-"):
        argv = ["build", *argv]

    parser = argparse.ArgumentParser(
        prog="backmap-prep",
        description="Generate LAMMPS input files for CG→AT backmapping",
    )
    subparsers = parser.add_subparsers(dest="command")

    build_parser = subparsers.add_parser(
        "build", help="Generate hybrid CG+AT data and input files (default)"
    )
    _add_common_args(build_parser)

    rebuild_parser = subparsers.add_parser(
        "rebuild",
        help="Rebuild hybrid data from an equilibrated CG LAMMPS data file",
    )
    _add_common_args(rebuild_parser)
    rebuild_parser.add_argument(
        "--cg-frame",
        type=Path,
        required=True,
        help="Equilibrated CG frame: rim135_cg_final.data or unwrapped cg_cl_conf_equil.gro",
    )

    cg_parser = subparsers.add_parser(
        "cg-only",
        help="Extract CG-only system for pre-equilibration",
    )
    _add_common_args(cg_parser)

    finalize_parser = subparsers.add_parser(
        "finalize-cg",
        help="Fix image flags on equilibrated CG LAMMPS data and export GRO",
    )
    _add_common_args(finalize_parser)
    finalize_parser.add_argument(
        "--equil",
        type=Path,
        required=True,
        help="Equilibrated CG LAMMPS data file (e.g. rim135_cg_equil.data)",
    )
    finalize_parser.add_argument(
        "--gro-name",
        type=str,
        default=None,
        help="Output GRO filename (default: cg_system coordinates stem + _equil.gro)",
    )
    finalize_parser.add_argument(
        "--no-gro",
        action="store_true",
        help="Skip unwrapped GRO export",
    )
    finalize_parser.add_argument(
        "--folded-gro",
        action="store_true",
        help="Write primary-cell (folded) coords to GRO instead of unwrapped",
    )

    hybrid_parser = subparsers.add_parser(
        "build-hybrid",
        help="Build GROMACS hybrid GRO/TOP for reactive network systems (Phase 3)",
    )
    hybrid_parser.add_argument(
        "settings",
        type=Path,
        help="YAML settings (with prep.bakery_xml) or bakery settings.xml",
    )
    hybrid_parser.add_argument(
        "--allow-no-bonds",
        action="store_true",
        help="Allow skipping AT bonds (XML-only mode)",
    )
    hybrid_parser.add_argument(
        "--chain-rng-seed",
        type=int,
        default=None,
        help="RNG seed for multi-chain AT source selection",
    )

    compare_parser = subparsers.add_parser(
        "compare-topology",
        help="Compare hybrid topology sections (rim135 parity helper)",
    )
    compare_parser.add_argument("reference", type=Path, help="Reference .top file")
    compare_parser.add_argument("candidate", type=Path, help="Generated .top file")

    migrate_parser = subparsers.add_parser(
        "migrate-settings",
        help="Convert legacy bakery settings.xml to v2 YAML (not yet implemented)",
    )
    migrate_parser.add_argument(
        "input_xml",
        type=Path,
        help="Source bakery settings.xml",
    )
    migrate_parser.add_argument(
        "-o",
        "--output",
        dest="output_yaml",
        type=Path,
        required=True,
        help="Destination settings YAML path",
    )

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 1

    if hasattr(args, "settings") and not args.settings.exists():
        print(f"Error: settings file not found: {args.settings}", file=sys.stderr)
        return 1

    commands: dict[str, Callable[[argparse.Namespace], int]] = {
        "build": _cmd_build,
        "rebuild": _cmd_rebuild,
        "cg-only": _cmd_cg_only,
        "finalize-cg": _cmd_finalize_cg,
        "build-hybrid": _cmd_build_hybrid,
        "compare-topology": lambda a: compare_topology_main([str(a.reference), str(a.candidate)]),
        "migrate-settings": _cmd_migrate_settings,
    }
    return commands[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
