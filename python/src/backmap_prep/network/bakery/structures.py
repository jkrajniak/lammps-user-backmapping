"""
Copyright (C) 2015-2019 Jakub Krajniak <jkrajniak@gmail.com>

This file is part of bakery.

This file is distributed under free software licence:
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import collections
import copy
import logging
import os
import pickle as cPickle
import random
import sys
import warnings
import xml.etree.ElementTree as etree
from functools import cmp_to_key

import networkx
import numpy as np

from . import files_io, tools

logger = logging.getLogger(__name__)

__doc__ = "Data structures."

BeadID = collections.namedtuple("BeadID", ["name", "degree"])

CGMolecule = collections.namedtuple(
    "CGMolecule", ["name", "ident", "fragment_name", "source_coordinate", "source_topology"]
)

BondedParams = collections.namedtuple("BondedParams", ["params", "typeid"])

logger = logging.getLogger(__name__)


def _chain_has_atoms(chain_atoms: dict, bead_list: list[str]) -> bool:
    """Return True when every atom in bead_list exists in chain_atoms."""
    for bead in bead_list:
        atom_name = bead.split(":")[2]
        if atom_name not in chain_atoms:
            return False
    return True


def _pick_chain_atoms(chains: dict, bead_list: list[str]) -> dict | None:
    """Pick a chain conformation that contains all atoms referenced by bead_list."""
    for chain_idx in sorted(chains.keys()):
        chain_atoms = chains[chain_idx]
        if _chain_has_atoms(chain_atoms, bead_list):
            return chain_atoms
    return None


def _fragment_atoms_available(
    coordinate: files_io.CoordinateFile,
    fragment_name: str | None,
    bead_list: list[str],
) -> bool:
    """Return True when some chain in the coordinate file contains all bead atoms."""
    if not coordinate.chains:
        return False
    if fragment_name:
        chain_name = fragment_name
    elif len(coordinate.chains) == 1:
        chain_name = next(iter(coordinate.chains))
    else:
        return True
    return _pick_chain_atoms(coordinate.chains.get(chain_name, {}), bead_list) is not None


class CGFragment:
    """Complex struct, with the construct that does a bit of processing."""

    def __init__(
        self,
        cg_molecule,
        fragment_list,
        active_sites=None,
        charge_map=None,
        as_remove=None,
        equilibrate_charges=False,
        type_map=None,
    ):
        self.topology = cg_molecule.source_topology
        self.coordinate = cg_molecule.source_coordinate
        self.fragment_list = fragment_list

        self.cg_molecule = cg_molecule

        self.charge_map = charge_map
        self.type_map = type_map
        self.active_sites_remove_map = as_remove

        # Select the set of atomistic coordinates.
        fragment_from_coordinate = False
        if cg_molecule.fragment_name:
            fragment_name = cg_molecule.fragment_name
        else:
            if len(self.coordinate.chains.values()) > 1:
                raise RuntimeError(
                    f"Please specify fragment_name for {cg_molecule.name} as it contains more than one chain type"
                )
            fragment_name = self.coordinate.chains.keys().__iter__().__next__()
            fragment_from_coordinate = True
        chains = self.coordinate.chains[fragment_name]
        # Select random chain from the ensemble of chains (bakery behaviour).
        chain_keys = list(chains.keys())
        atoms = None
        for chain_idx in random.sample(chain_keys, len(chain_keys)):
            candidate = chains[chain_idx]
            if _chain_has_atoms(candidate, self.fragment_list):
                atoms = candidate
                break
        if atoms is None:
            atoms = chains[random.sample(chain_keys, 1)[0]]
        self.atomparams = {}
        if fragment_name not in self.topology.chain_atom_names:
            logger.debug(80 * "!")
            logger.debug("\nError while preparing cg fragments")
            logger.debug(f"Fragment topology: {self.topology.file_name}")
            logger.debug(f"Fragment coordinates: {self.coordinate.file_name}\n")
            logger.debug(
                f"Fragment {fragment_name} not found in the list of residues {self.topology.chain_atom_names}"
            )
            if fragment_from_coordinate:
                logger.debug(
                    f"\nPlease check your coordinate file {self.coordinate.file_name}. Very likely the name of molecules ({fragment_name}) does not cooresponds"
                    " to the names defined in the settings file and the topology file"
                    f" ({self.topology.chains.keys().__iter__().__next__()})."
                )
            logger.debug(80 * "!")
            raise RuntimeError("Wrong fragment")
        for k, v in self.topology.chain_atom_names[fragment_name].items():
            self.atomparams[k] = v[0]

        self.active_sites = {}

        self.com = np.zeros(3)
        total_mass = 0.0
        tmp_atom_list = []
        at_charges = []
        num_at_change = 0
        for bid, bead in enumerate(self.fragment_list):
            atom_name = bead.split(":")[2]
            atom = atoms[atom_name]._replace()
            atom_mass = self.atomparams[atom_name].mass
            if charge_map:
                if charge_map[bid] == "*":
                    at_charges.append(self.atomparams[atom_name].charge)
                    num_at_change += 1
                else:
                    at_charges.append(float(charge_map[bid]))
            else:
                at_charges.append(self.atomparams[atom_name].charge)
                num_at_change += 1
            self.com += atom_mass * atom.position
            total_mass += atom_mass
            tmp_atom_list.append(atom)

        total_charge = sum(at_charges)
        if equilibrate_charges and total_charge != 0.0:
            num_atoms = len(self.fragment_list)
            dcharge = total_charge / num_at_change
            if charge_map:
                for bid, c in enumerate(charge_map):
                    if c == "*":
                        self.charge_map[bid] = at_charges[bid] - dcharge
                    else:
                        self.charge_map[bid] = float(self.charge_map[bid])
            else:
                self.charge_map = []
                for c in at_charges:
                    self.charge_map.append(c - dcharge)
            logger.info(f"Equilibrate charge, total: {total_charge} -> {sum(self.charge_map)}")

        self.com /= total_mass
        self.cg_mass = total_mass
        # Move the atoms in fragments to the origin, always by substracting the com of the fragment.
        self.atom_in_fragments = [
            files_io.Atom(x[0], x[1], x[2], x[3], x[4] - self.com) for x in tmp_atom_list
        ]
        if active_sites:
            for at_as in active_sites:
                t = at_as.split(":")  # Format: atom name -> maximum degree
                self.active_sites[t[1]] = int(t[2])


def _atom_names_from_param_key(key: tuple) -> tuple:
    """Strip chain prefixes from bakery topology param keys (e.g. EPO:C1 -> C1)."""
    names = []
    for part in key:
        if isinstance(part, str) and ":" in part:
            names.append(part.split(":", 1)[1])
        else:
            names.append(part)
    return tuple(names)


def _chain_names_from_param_key(key: tuple) -> tuple:
    chains = []
    for part in key:
        if isinstance(part, str) and ":" in part:
            chains.append(part.split(":", 1)[0])
        else:
            chains.append(part)
    return tuple(chains)


def _get_params(input_dict, key_list, raise_exception=True):
    param = None
    for key in key_list:
        param = input_dict.get(tuple(key))
        if param:
            return param
    # Bakery falls back to atom-name tuples when chain:name keys miss, but only
    # when exactly one template shares the same atom-name and chain multiset.
    if len(key_list) >= 2 and len(key_list[0]) == len(key_list[1]):
        chain_key = tuple(key_list[0])
        name_key = tuple(key_list[1])
        if all(isinstance(x, str) for x in name_key):
            triplet_chains = _chain_names_from_param_key(chain_key)
            matches = []
            for template_key, template_param in input_dict.items():
                if len(template_key) != len(name_key):
                    continue
                tmpl_names = _atom_names_from_param_key(template_key)
                if tmpl_names != name_key and tmpl_names != tuple(reversed(name_key)):
                    continue
                if sorted(_chain_names_from_param_key(template_key)) != sorted(triplet_chains):
                    continue
                matches.append(template_param)
            if len(matches) == 1:
                return matches[0]
    if not param:
        if raise_exception:
            raise RuntimeError(f"Missing definition for {key_list}")
        logger.error(f"Missing definition for {key_list}")
        return False
    return param


class BackmapperSettings2:
    def __init__(self, input_xml, allow_no_bonds=False, generate_only_graph=False):
        self.res2atom = collections.defaultdict(list)
        self.cg_active_sites = collections.defaultdict(list)
        self.atom_id2fragment = {}
        self.atom_ids = []  # List of ids of atomistic particles
        self.cg_old_new_id = {}
        self.cross_prefix = "cross_"
        self.cg_new_id_old = {}
        self.atom2cg = {}
        self.fragments = {}  # key-> molecule name, val-> dict
        self.mol_atomid_map = collections.defaultdict(
            dict
        )  # Map between old id and new id, divided into molecules.
        self.mol_atomname_map = collections.defaultdict(
            dict
        )  # Map between old id and new id, divided into molecules.
        self.cg2atom = collections.defaultdict(list)

        self.allow_no_bonds = allow_no_bonds
        self.generate_only_graph = generate_only_graph

        self.charge_transfer = {}  # Map with charge transfer.

        self.global_graph = networkx.Graph()

        # Parse XML file or accept an in-memory Element root (v2 loader).
        if isinstance(input_xml, etree.Element):
            self.root = input_xml
        else:
            tree = etree.parse(input_xml)
            self.root = tree.getroot()
        self._parse()

        # AT cross terms
        self.at_cross_bonds = {}
        self.at_cross_angles = {}
        self.at_cross_dihedrals = {}

    def _parse(self):
        cg_configuration = self.root.find("cg_configuration")
        self.cg_topology = files_io.GROMACSTopologyFile(
            cg_configuration.find("topology").text.strip()
        )
        self.cg_topology.read()

        # Fix the res_id;
        last_chain_idx = 0
        ch_idx = 0
        for at_id in sorted(self.cg_topology.atoms):
            if self.cg_topology.atoms[at_id].chain_idx != last_chain_idx:
                last_chain_idx = self.cg_topology.atoms[at_id].chain_idx
                ch_idx += 1
            self.cg_topology.atoms[at_id].chain_idx = ch_idx
        # Replicate if the topology is multimolecule.
        self.cg_topology.replicate()

        self.cg_graph = self.cg_topology.get_graph()
        coordinate_file = None
        if cg_configuration.find("file") is not None:
            coordinate_file = cg_configuration.find("file").text.strip()
            warnings.warn(
                "cg_configuration.file is deprected, use cg_configuration.coordinate",
                DeprecationWarning,
            )
        elif cg_configuration.find("coordinate") is not None:
            coordinate_file = cg_configuration.find("coordinate").text.strip()
        else:
            raise RuntimeError("Missing tag cg_configuration.coordinate")
        self.cg_coordinate = files_io.read_coordinates(coordinate_file)
        logger.debug(
            f"Read coordinate file {coordinate_file} (num: {len(self.cg_coordinate.atoms)})"
        )

        if len(self.cg_coordinate.atoms) != len(self.cg_topology.atoms):
            raise RuntimeError(
                f"Number of atoms defined in the coordinate file {coordinate_file} (num: {len(self.cg_coordinate.atoms)}) is not the same as in topology (num: {len(self.cg_topology.atoms)})"
            )

        # Predefined active sites map bid1 bid2 label_active_site1 label_active_site2
        self.predefined_active_sites = {}
        if cg_configuration.find("predefined_active_sites") is not None:
            predef_fname = cg_configuration.find("predefined_active_sites").text.strip()
            if os.path.exists(predef_fname):
                for l in open(predef_fname):
                    t = l.split()
                    self.predefined_active_sites[tuple(map(int, t[:2]))] = t[2:]
            else:
                logger.warn(f"File {predef_fname} with predefined active sites not found")
        logger.info(f"Predefined active sites: {len(self.predefined_active_sites)}")

        hyb_configuration = self.root.find("hybrid_configuration")
        self.hybrid_configuration = {
            "file": files_io.GROFile(hyb_configuration.find("file").text.strip())
        }

        hyb_topology = self.root.find("hybrid_topology")
        self.hyb_topology = files_io.GROMACSTopologyFile(hyb_topology.find("file").text.strip())
        self.hyb_topology.init()
        self.hyb_topology.new_data.update(
            {f"{self.cross_prefix}{k}": {} for k in self.hyb_topology.new_data}
        )
        self.hyb_topology.header_section.append(
            "; GROMACS like hybrid topology file.\n; ===================== \n\n"
        )
        if hyb_topology.find("include") is not None:
            for include_entry in hyb_topology.find("include").text.split():
                ie = include_entry.strip()
                if ie.startswith(";"):
                    self.hyb_topology.header_section.append(
                        ";#include {}\n".format(ie.replace(";", ""))
                    )
                else:
                    self.hyb_topology.header_section.append(f"#include {ie}\n")
            self.hyb_topology.header_section.append("\n\n")

        self.hyb_topology.atomtypes = {}
        # Change atomtypes to 'V' for CG.
        for k, d in self.cg_topology.atomtypes.items():
            d["type"] = "V"
            self.hyb_topology.atomtypes[k] = d

        # Separate exclusion for AT and CG; if exclusion present then
        # excl_AT=excl_AT = exclusion
        #
        excl_at = excl_cg = excl = None
        if hyb_topology.find("molecule_type").find("exclusion_at") is not None:
            excl_at = hyb_topology.find("molecule_type").find("exclusion_at").text.strip()
        if hyb_topology.find("molecule_type").find("exclusion_cg") is not None:
            excl_cg = hyb_topology.find("molecule_type").find("exclusion_cg").text.strip()
        if hyb_topology.find("molecule_type").find("exclusion") is not None:
            excl = hyb_topology.find("molecule_type").find("exclusion").text.strip()
        else:
            if excl_at:
                excl = excl_at
            elif excl_cg:
                excl = excl_cg
            else:
                raise RuntimeError("exclusion not specified")
        if excl and not excl_at:
            excl_at = excl
        if excl and not excl_cg:
            excl_cg = excl

        self.hyb_topology.moleculetype = {
            "name": hyb_topology.find("molecule_type").find("name").text.strip(),
            "nrexcl": excl,
            "excl_at": excl_at,
            "excl_cg": excl_cg,
        }
        self.hyb_topology.molecules = {"name": hyb_topology.find("system").text.strip(), "mol": 1}
        self.hyb_topology.bondtypes = self.cg_topology.bondtypes
        self.hyb_topology.angletypes = self.cg_topology.angletypes
        self.hyb_topology.dihedraltypes = self.cg_topology.dihedraltypes

        # Allowed hybrid bond patterns.
        self.restricted_cross_bond_patterns = None
        if hyb_topology.find("restricted_cross_bond_patterns") is not None:
            self.restricted_cross_bond_patterns = set()
            for x in hyb_topology.find("restricted_cross_bond_patterns").text.strip().split():
                self.restricted_cross_bond_patterns.add(tuple(x.split("-")))
                self.restricted_cross_bond_patterns.add(tuple(reversed(x.split("-"))))
            logger.info(
                f"Found restrictions on the cross bonds: {self.restricted_cross_bond_patterns}"
            )

        # Reads parameters of bonded terms.
        self.bond_params = collections.defaultdict(dict)
        self.angle_params = collections.defaultdict(dict)
        self.dihedral_params = collections.defaultdict(dict)
        bonded_terms = hyb_topology.findall("bonds")
        for bond_term in bonded_terms:
            params = bond_term.attrib["params"].split()
            typeid = bond_term.attrib.get("typeid")
            bead_list = bond_term.text.strip().split()
            if len(bead_list) % 2 != 0:
                raise RuntimeError(f"Wrong pairs in <bonds> section, found {bead_list}")
            pair_list = zip(bead_list[::2], bead_list[1::2], strict=False)
            for p1, p2 in pair_list:
                if (p1, p2) in self.bond_params:
                    raise RuntimeError(f"Parameters for pair: {(p1, p2)} alread defined {params}")
                self.bond_params[(p1, p2)] = BondedParams(params, typeid)
                self.bond_params[(p2, p1)] = BondedParams(params, typeid)
        angle_terms = hyb_topology.findall("angles")
        for angle_term in angle_terms:
            params = angle_term.attrib["params"].split()
            typeid = angle_term.attrib.get("typeid")
            bead_list = angle_term.text.strip().split()
            if len(bead_list) % 3 != 0:
                raise RuntimeError(f"Wrong triplets in <angles> section, found {bead_list}")
            for idx in range(0, len(bead_list), 3):
                p1, p2, p3 = bead_list[idx], bead_list[idx + 1], bead_list[idx + 2]
                if (p1, p2, p3) in self.angle_params:
                    raise RuntimeError(
                        f"Parameters for triplet: {(p1, p2, p3)} already defined {params}"
                    )
                self.angle_params[(p1, p2, p3)] = BondedParams(params, typeid)
                self.angle_params[(p3, p2, p1)] = BondedParams(params, typeid)
                logger.info(f"Read angle params for {p1}-{p2}-{p3}: {params}")
        dihedral_terms = hyb_topology.findall("dihedrals")
        for dihedral_term in dihedral_terms:
            params = dihedral_term.attrib["params"].split()
            typeid = dihedral_term.attrib.get("typeid")
            bead_list = dihedral_term.text.strip().split()
            if len(bead_list) % 4 != 0:
                raise RuntimeError(f"Wrong quadruplets in <dihedrals> section, found {bead_list}")
            for idx in range(0, len(bead_list), 4):
                p1, p2, p3, p4 = (
                    bead_list[idx],
                    bead_list[idx + 1],
                    bead_list[idx + 2],
                    bead_list[idx + 3],
                )
                if (p1, p2, p3, p4) in self.dihedral_params:
                    raise RuntimeError(
                        f"Parameters for triplet: {(p1, p2, p3, p4)} already defined {params}"
                    )
                self.dihedral_params[(p1, p2, p3, p4)] = BondedParams(params, typeid)
                self.dihedral_params[(p4, p3, p2, p1)] = BondedParams(params, typeid)
                logger.info(f"Read dihedral params for {p1}-{p2}-{p3}-{p4}: {params}")

        self._parse_cg_molecules()

    def _parse_cg_molecules(self):
        """Parse the cg_molecule section and prepares data structures."""
        for r in self.root.findall("cg_molecule"):
            equilibrate_charges = r.attrib.get("equilibrate_charges", "0") == "1"
            cg_mol_name = r.find("name").text
            cg_mol_ident = r.find("ident").text
            fragment_name = r.find("fragment_name")
            if fragment_name:
                fragment_name = fragment_name.text
            # For every residue degree the molecule can have different topology/coordinate file
            cg_mol_source_topologies = {}
            for x in r.find("source_topology").findall("file"):
                mol_degree = x.attrib.get("molecule_degree", "*")
                mol_when = x.attrib.get("when", "*")
                file_name = x.text
                cg_mol_source_topologies[(mol_degree, mol_when)] = file_name

            if len(cg_mol_source_topologies) != len(r.find("source_topology").findall("file")):
                raise RuntimeError(
                    "The degree atritbute in source_topology.file tag has to be unique"
                )

            if not cg_mol_source_topologies:
                cg_mol_source_topologies[("*", "*")] = r.find("source_topology").text

            cg_mol_source_coordinate = {}
            for x in r.find("source_coordinate").findall("file"):
                mol_degree = x.attrib.get("molecule_degree", "*")
                mol_when = x.attrib.get("when", "*")
                file_name = x.text
                cg_mol_source_coordinate[(mol_degree, mol_when)] = file_name

            if len(cg_mol_source_coordinate) != len(r.find("source_coordinate").findall("file")):
                raise RuntimeError(
                    "The degree atritbute in source_coordinate.file tag has to be unique"
                )
            if not cg_mol_source_coordinate:
                cg_mol_source_coordinate[("*", "*")] = r.find("source_coordinate").text

            if not all(
                map(
                    lambda x: x[0] == x[1],
                    zip(cg_mol_source_coordinate, cg_mol_source_topologies, strict=False),
                )
            ):
                raise RuntimeError(
                    "The set of degree for source_coordinate and source_topology has to be the same"
                )

            related_bead_names = [x[1] for x in cg_mol_source_coordinate]

            for mol_deg_bead_name in cg_mol_source_coordinate:
                mol_deg, related_bead_name = mol_deg_bead_name
                cg_molecule = CGMolecule(
                    cg_mol_name,
                    cg_mol_ident,
                    fragment_name,
                    files_io.read_coordinates(cg_mol_source_coordinate[mol_deg_bead_name]),
                    files_io.GROMACSTopologyFile(cg_mol_source_topologies[mol_deg_bead_name]),
                )
                cg_molecule.source_topology.read()
                cg_molecule.source_coordinate.read()

                if (mol_deg, cg_molecule.name, related_bead_name) in self.fragments:
                    raise RuntimeError(
                        f"Fragment ({mol_deg}, {cg_molecule.name}, {related_bead_name}) already defined"
                    )

                self.fragments[(mol_deg, cg_molecule.name, related_bead_name)] = {
                    "cg_molecule": cg_molecule
                }
                # Generate the beads.
                for cg_bead in r.findall("cg_bead"):
                    name = cg_bead.find("name").text  # bead name
                    self.fragments[(mol_deg, cg_molecule.name, related_bead_name)][name] = {}
                    for beads in cg_bead.findall("beads"):
                        degree = beads.attrib.get("degree", "*")
                        mol_degrees = beads.attrib.get("molecule_degree")
                        # if mol_degrees and mol_deg not in mol_degrees.split(','):
                        #    continue
                        active_sites = beads.attrib.get("active_site")
                        as_remove = {}
                        if active_sites:
                            active_sites = list(map(str.strip, active_sites.split()))
                            tmp_as = [":".join(x.split(":")[:2]) for x in active_sites]
                            # Remove with active sites. So whenever the active site will be used then this list
                            # of atoms will be removed.
                            removes = beads.findall("remove")
                            for asr in removes:
                                as_name = asr.attrib["active_site"]
                                if as_name not in tmp_as:
                                    raise RuntimeError(
                                        f"Defined <remove> entry for active_site {as_name} but it's not on the list"
                                    )
                                beads_to_remove = asr.text
                                if not beads_to_remove:
                                    raise RuntimeError(
                                        f"No beads to remove in active_site {as_name}"
                                    )
                                as_remove[as_name] = list(map(str.strip, asr.text.split()))
                        bead_list = beads.text.split()

                        cm = beads.find("charge_map")
                        charge_map = None
                        if cm is not None:
                            charge_map = cm.text.split()
                            if len(charge_map) != len(bead_list):
                                raise RuntimeError(
                                    f"Number of entries in charge_map {len(charge_map)} does not match number of beads {len(bead_list)}"
                                )
                        tm = beads.find("type_map")
                        type_map = None
                        if tm is not None:
                            type_map = tm.text.split()
                            if len(type_map) != len(bead_list):
                                raise RuntimeError(
                                    f"Number of entries in type_map {len(type_map)} does not match number of beads {len(bead_list)}"
                                )

                        if not _fragment_atoms_available(
                            cg_molecule.source_coordinate, fragment_name, bead_list
                        ):
                            continue

                        cg_fragment = CGFragment(
                            cg_molecule,
                            bead_list,
                            active_sites,
                            charge_map,
                            as_remove,
                            equilibrate_charges,
                            type_map,
                        )
                        try:
                            self.fragments[(mol_deg, cg_molecule.name, related_bead_name)][name][
                                degree
                            ] = cg_fragment
                        except KeyError:
                            continue

            logger.debug("Found following fragment definitions")
            for k, v in self.fragments.items():
                logger.debug(f"==== {k} ====")
                for kk, vv in v.items():
                    if kk == "cg_molecule":
                        continue
                    logger.debug(f"++ {kk} degs={vv.keys()}")

            # Read charge transfer.
            # <charge_transfer on="IPD:N1:2" from="IPD:H8" to="EPO:C23#H25,EPO:C41#H66,HDD:C21#H43,HDD:C32#H44" />
            # <charge_transfer on="when" from="atom label" to="set of atom labels" />
            charge_transfers = r.find("charge_transfers")
            if charge_transfers is not None:
                for ct in charge_transfers.findall("charge_transfer"):
                    t = ct.attrib["when"].split(":")
                    on_key = f"{t[0]}:{t[1]}"
                    on_deg = int(t[2])
                    if on_key not in self.charge_transfer:
                        self.charge_transfer[on_key] = {}
                    if on_deg in self.charge_transfer[on_key]:
                        raise RuntimeError(f"Mistake in charge transfer, key {t} already defined")
                    transfer_from = ct.attrib["from_atom"]
                    transfer_to_map = {}
                    self.charge_transfer[on_key][on_deg] = {
                        "from_atom": transfer_from.split(":"),
                        "when_to_atom": transfer_to_map,
                    }
                    transfer_to = ct.attrib["when_to_atom"].split(",")
                    for tt in transfer_to:
                        tt_to_on, tt_to = tt.split("#")
                        transfer_to_map[tt_to_on] = tt_to

    def prepare_hybrid(self):
        """Creates hybrid files."""
        outfile = self.hybrid_configuration["file"]
        new_at_id = 1

        cg_ids = sorted(self.cg_graph.nodes())
        cg_atomtypes = []

        # Residue graph for getting the residue degree.
        residue_graph = networkx.MultiGraph()

        for cg_id in self.cg_graph.nodes():
            cg_bead = self.cg_graph.nodes[cg_id]
            res_id = cg_bead["res_id"]
            if res_id not in residue_graph:
                residue_graph.add_node(
                    cg_bead["res_id"], chain_name=cg_bead["chain_name"], cg_nodes=[]
                )
            residue_graph.nodes[res_id]["cg_nodes"].append(cg_id)

        for cg_bonds in self.cg_graph.edges():
            cg_nodes = list(map(self.cg_graph.nodes.get, cg_bonds))
            if cg_nodes[0]["res_id"] != cg_nodes[1]["res_id"]:  # Ignore self-loops
                residue_graph.add_edge(cg_nodes[0]["res_id"], cg_nodes[1]["res_id"])

        for res_id, deg in residue_graph.degree():
            residue_graph.nodes[res_id]["degree"] = str(deg)
            residue_graph.nodes[res_id]["fragment_key"] = None

        for res_id in sorted(residue_graph.nodes()):
            cg_nodes = residue_graph.nodes[res_id]["cg_nodes"]
            residue_degree = residue_graph.nodes[res_id]["degree"]
            residue_name = residue_graph.nodes[res_id]["chain_name"]
            possible_fragments = [
                x
                for x in self.fragments.keys()
                if x[1] == residue_name and (x[0] == "*" or x[0] == residue_degree)
            ]

            selected_fragment = None
            wrong_cg_nodes = []
            for possible_fragment in possible_fragments:
                f = self.fragments[possible_fragment]
                selected_fragment = possible_fragment
                # Iterates over all cg nodes for given residue.
                for cg_id in cg_nodes:
                    cg_node = self.cg_graph.nodes[cg_id]
                    if (
                        str(cg_node["degree"]) not in f[cg_node["name"]]
                        and "*" not in f[cg_node["name"]]
                    ):
                        selected_fragment = None
                        wrong_cg_nodes.append((cg_node, f[cg_node["name"]]))
                if selected_fragment is not None:
                    break
            if selected_fragment is None:
                logger.error("Last processed molecule:")
                logger.error(f"Residue id: {res_id}")
                logger.error(f"Bead ids: {cg_nodes}")
                logger.error(f"Residue degree: {residue_degree}")
                logger.error(f"Residue name: {residue_name}")

                logger.error(
                    "\nIt is very likely that your .xml file does not"
                    f" contain definition of a fragment for residue {residue_name}"
                )

                logger.error("====== List of problematic CG nodes ======")
                for wcn in wrong_cg_nodes:
                    logger.error(wcn)

                logger.error(40 * "=")

                raise RuntimeError(
                    f"Problem with the option file, could not find correct fragment for molecule {res_id}"
                )

            for cg_id in cg_nodes:
                cg_atom = self.cg_coordinate.atoms[cg_id]
                cg_bead = self.cg_graph.nodes[cg_id]
                cg_bead_id = new_at_id

                fragments = self.fragments[selected_fragment][cg_bead["name"]]

                cg_fragment = fragments.get(str(cg_bead["degree"]), fragments.get("*"))
                if not cg_fragment:
                    raise RuntimeError("Problem with getting atomistic fragments")
                outfile.atoms[cg_bead_id] = files_io.Atom(
                    atom_id=cg_bead_id,
                    name=cg_bead["name"],
                    chain_name=cg_fragment.cg_molecule.ident,
                    # Change chain name to one from <ident> tag in the cg_molecule
                    chain_idx=res_id,
                    position=cg_atom.position,
                )
                self.cg_new_id_old[cg_bead_id] = cg_id
                self.cg_old_new_id[cg_id] = cg_bead_id

                self.atom_id2fragment[new_at_id] = cg_fragment

                self.global_graph.add_node(
                    cg_bead_id, atom_id=new_at_id, bead_type="CG", cg_bead_id=cg_bead_id, **cg_bead
                )

                # Change the mass of CG bead
                chain_name = cg_fragment.cg_molecule.ident
                topol_atom = copy.copy(self.cg_topology.atoms[cg_id])
                topol_atom.mass = cg_fragment.cg_mass
                topol_atom.atom_id = cg_bead_id
                topol_atom.chain_idx = res_id
                topol_atom.chain_name = (
                    cg_fragment.cg_molecule.ident
                )  # Change chain name to one from <ident> tag in the cg_molecule
                topol_atom.cgnr = cg_bead_id
                self.hyb_topology.atoms[new_at_id] = topol_atom
                if cg_fragment.cg_molecule.ident not in self.hyb_topology.chains:
                    self.hyb_topology.chains[cg_fragment.cg_molecule.ident] = {}
                if res_id not in self.hyb_topology.chains[cg_fragment.cg_molecule.ident]:
                    self.hyb_topology.chains[cg_fragment.cg_molecule.ident][res_id] = {}
                if (
                    topol_atom.name
                    in self.hyb_topology.chains[cg_fragment.cg_molecule.ident][res_id]
                ):
                    raise RuntimeError(
                        f"{topol_atom.name} already defined, please make sure that atom names are unique"
                    )
                cg_atomtypes.append(topol_atom.atom_type)
                self.hyb_topology.chains[cg_fragment.cg_molecule.ident][res_id][topol_atom.name] = (
                    topol_atom
                )

                # Set the atomistic coordinates for this fragment, put it in the
                # output coordinate file.
                cg_com = cg_atom.position
                new_at_id += 1
                for idx, at in enumerate(cg_fragment.atom_in_fragments):
                    new_at_atom = files_io.Atom(
                        new_at_id,
                        at.name,
                        cg_fragment.cg_molecule.ident,  # Change chain name to one from <ident> tag in the cg_molecule
                        res_id,
                        at.position + cg_com,
                    )
                    outfile.atoms[new_at_id] = new_at_atom
                    self.global_graph.add_node(
                        new_at_id,
                        atom_id=new_at_id,
                        name=at.name,
                        bead_type="AT",
                        res_id=res_id,
                        cg_bead_id=cg_bead_id,
                        orig_cg_bead_id=cg_id,
                        position=at.position + cg_com,
                        chain_name=chain_name,
                    )
                    if "atom_ids" not in self.global_graph.nodes[cg_bead_id]:
                        self.global_graph.nodes[cg_bead_id]["atom_ids"] = []
                    self.global_graph.nodes[cg_bead_id]["atom_ids"].append(new_at_id)

                    # Set topology atom
                    topol_atom = copy.copy(
                        cg_fragment.topology.chain_atom_names[at.chain_name][at.name][0]
                    )
                    if res_id not in self.mol_atomid_map[chain_name]:
                        self.mol_atomid_map[chain_name][res_id] = {}
                        self.mol_atomname_map[chain_name][res_id] = {}
                    self.mol_atomid_map[chain_name][res_id][topol_atom.atom_id] = new_at_id
                    self.mol_atomname_map[chain_name][res_id][topol_atom.name] = new_at_id

                    if cg_fragment.charge_map:
                        topol_atom.charge = cg_fragment.charge_map[idx]
                    if cg_fragment.type_map and cg_fragment.type_map[idx] != "*":
                        topol_atom.atom_type = cg_fragment.type_map[idx]
                    topol_atom.atom_id = new_at_id
                    topol_atom.chain_idx = res_id
                    topol_atom.chain_name = cg_fragment.cg_molecule.ident
                    topol_atom.cgnr = new_at_id
                    self.hyb_topology.atoms[new_at_id] = topol_atom
                    if cg_fragment.cg_molecule.ident not in self.hyb_topology.chains:
                        self.hyb_topology.chains[cg_fragment.cg_molecule.ident] = {}
                    if res_id not in self.hyb_topology.chains[cg_fragment.cg_molecule.ident]:
                        self.hyb_topology.chains[cg_fragment.cg_molecule.ident][res_id] = {}
                    if (
                        topol_atom.name
                        in self.hyb_topology.chains[cg_fragment.cg_molecule.ident][res_id]
                    ):
                        raise RuntimeError(
                            f"{topol_atom.name} already defined, please make sure that atom names are unique"
                        )
                    self.hyb_topology.chains[cg_fragment.cg_molecule.ident][res_id][
                        topol_atom.name
                    ] = topol_atom

                    # Set active sites.
                    if at.name in cg_fragment.active_sites:
                        self.cg_active_sites[cg_bead_id].append(
                            (new_at_atom, cg_fragment.active_sites[at.name])
                        )

                    self.atom_id2fragment[new_at_id] = cg_fragment
                    self.atom_ids.append(new_at_id)
                    self.atom2cg[new_at_id] = cg_bead_id
                    self.cg2atom[cg_bead_id].append(new_at_id)
                    self.res2atom[res_id].append(new_at_id)

                    new_at_id += 1

        outfile.box = self.cg_coordinate.box

        # Rebuild hybrid topology.
        self.rebuild_hybrid_topology()
        self.hyb_topology.write()
        # Write the hybrid coordinate file.
        outfile.write(force=True)
        # Write only atomistic topology
        hyb_top = files_io.GROMACSTopologyFile(self.hyb_topology.file_name)
        hyb_top.read()
        at_topol = tools.get_atomistic_topology(hyb_top, virtual_atomtypes=cg_atomtypes)
        at_topol.write(f"at_{self.hyb_topology.file_name}")

        # Write the list of bonds, angles and dihedrals to separate files.
        out_cross_bonds = "cross_bonds_{}".format(
            self.hyb_topology.file_name.replace(".top", ".dat")
        )
        out_cross_angles = "cross_angles_{}".format(
            self.hyb_topology.file_name.replace(".top", ".dat")
        )
        out_cross_dihedrals = "cross_dihedrals_{}".format(
            self.hyb_topology.file_name.replace(".top", ".dat")
        )
        with open(out_cross_bonds, "w") as outbond:
            outl = []
            for k, p in self.at_cross_bonds.items():
                new_k = list(map(at_topol.old2new_ids.get, k))
                if p.typeid:
                    outl.append([int(p.typeid)] + new_k)
            outl.sort(key=lambda x: x[0])
            outbond.write("\n".join([" ".join(map(str, p)) for p in outl]))
        logger.info(f"Saved {out_cross_bonds}")

        with open(out_cross_angles, "w") as outbond:
            outl = []
            for k, p in self.at_cross_angles.items():
                new_k = list(map(at_topol.old2new_ids.get, k))
                if p.typeid:
                    outl.append([int(p.typeid)] + new_k)
            outl.sort(key=lambda x: x[0])
            outbond.write("\n".join([" ".join(map(str, p)) for p in outl]))
        logger.info(f"Saved {out_cross_angles}")

        with open(out_cross_dihedrals, "w") as outbond:
            outl = []
            for k, p in self.at_cross_dihedrals.items():
                new_k = list(map(at_topol.old2new_ids.get, k))
                if p.typeid:
                    outl.append([int(p.typeid)] + new_k)
            outl.sort(key=lambda x: x[0])
            outbond.write("\n".join([" ".join(map(str, p)) for p in outl]))
        logger.info(f"Saved {out_cross_dihedrals}")
        # Generate exclusion list.
        self._generate_exclusion_lists()

    def rebuild_hybrid_topology(self):
        """Regenerate the hybrid topology based on the new particle ids."""

        # First build coarse-grained topology.
        def generate_cg_b_list(old_list):
            if not old_list:
                return {}
            return {
                tuple(map(self.cg_old_new_id.get, k)): v + [" ; cg_bonded"]
                for k, v in old_list.items()
            }

        def generate_at_b_list(at):
            """Generates atomistic list of bonds"""
            topology = self.atom_id2fragment[at.atom_id].topology
            bonded_lists = [
                (2, topology.bonds, "bonds"),
                (3, topology.angles, "angles"),
                (4, topology.dihedrals, "dihedrals"),
                (4, topology.improper_dihedrals, "dihedrals"),
                (2, topology.pairs, "pairs"),
            ]
            for ternary, top_list, output_name in bonded_lists:
                for p, params in top_list.items():
                    old2new_id = self.mol_atomid_map[at.chain_name][at.chain_idx]
                    new_tuple = tuple(map(old2new_id.get, p))

                    # Skip the terms that involves the missing atoms. The atoms can be missing because of
                    # the degree dependent atomistic fragments.
                    if None not in new_tuple and at.atom_id in new_tuple:  # correct pair
                        cg_beads = list(map(self.atom2cg.get, new_tuple))
                        if cg_beads.count(cg_beads[0]) == len(
                            cg_beads
                        ):  # Atoms in the same CG bead
                            out_name = output_name
                        else:
                            out_name = f"{self.cross_prefix}{output_name}"
                        self.hyb_topology.new_data[out_name][new_tuple] = params

        logger.info("Renumering atomistic and coarse-grained bonds")
        self.hyb_topology.new_data[f"{self.cross_prefix}bonds"] = generate_cg_b_list(
            self.cg_topology.bonds
        )
        self.hyb_topology.new_data[f"{self.cross_prefix}angles"] = generate_cg_b_list(
            self.cg_topology.angles
        )
        self.hyb_topology.new_data[f"{self.cross_prefix}dihedrals"] = generate_cg_b_list(
            self.cg_topology.dihedrals
        )
        self.hyb_topology.new_data[f"{self.cross_prefix}dihedrals"].update(
            generate_cg_b_list(self.cg_topology.improper_dihedrals)
        )

        # Create the atomistic topology.
        for atid in self.atom_ids:
            generate_at_b_list(self.hyb_topology.atoms[atid])

        for (b1, b2), params in self.hyb_topology.new_data[f"{self.cross_prefix}bonds"].items():
            self.global_graph.add_edge(b1, b2)
        for (b1, b2), params in self.hyb_topology.new_data["bonds"].items():
            self.global_graph.add_edge(b1, b2)

        # Generate list of CG bonds that are not defined at AT level
        cg_cross_bonds = set()
        for b1, b2 in self.global_graph.edges():
            if b1 in self.cg_new_id_old and b2 in self.cg_new_id_old:  # Look only on CG bonds
                ats_b1 = self.cg2atom[b1]
                ats_b2 = set(self.cg2atom[b2])
                is_connected = False
                for at1 in ats_b1:
                    edges = set(self.global_graph[at1])
                    intersect = ats_b2.intersection(edges)
                    if intersect:
                        is_connected = True
                        break
                if not is_connected:
                    cg_cross_bonds.add(tuple(sorted([b1, b2])))

        logger.info(f"Found {len(cg_cross_bonds)} bonds between coarse-grained beads")
        logger.info(
            "Generating atomistic cross-bonds between coarse-grained beads; It will take a while..."
        )

        with open("graph_before_cross_bonds.pck", "wb+") as outf:
            cPickle.dump(self.global_graph, outf, protocol=2)
        logger.info("Saved graph before cross-bonds in graph_before_cross_bonds.pck")

        if self.generate_only_graph:
            logger.info("That's all for now - graph_before_cross_bonds.pck was generated.")
            sys.exit(0)

        # Create the atomistic bonds across the coarse-grained beads.
        # We iterate over bonds in cg_graph and then generate the bonds.
        # At the level of topology, CG bonds are already defined.
        at_cross_bonds = []
        atoms_to_remove = []
        charge_to_transfer = []
        progress_indc = 0.0
        progress_indc_total = len(cg_cross_bonds)
        global_degree = dict(self.global_graph.degree())
        cg_cross_bonds = sorted(cg_cross_bonds)
        for b1, b2 in cg_cross_bonds:
            n1 = self.global_graph.nodes[b1]
            n2 = self.global_graph.nodes[b2]
            if n1["res_id"] != n2["res_id"]:  # Cross bond between beads in different chains.
                # Look for active sites on both CG molecules.
                if self.predefined_active_sites:
                    ats1, ats2, global_degree = self._get_predefined_active_sites(
                        b1, b2, global_degree, atoms_to_remove
                    )
                else:
                    ats1, ats2, global_degree = self._search_active_sites(
                        b1, b2, global_degree, atoms_to_remove
                    )

                # End for at1, for at2 loops
                if ats1 is not None and ats2 is not None:
                    # Selected proper active sites.
                    at_cross_bonds.append((ats1.atom_id, ats2.atom_id))
                    self.global_graph.add_edge(ats1.atom_id, ats2.atom_id)

                    self.global_graph.nodes[ats1.atom_id]["cg_bead_id"] = b1
                    self.global_graph.nodes[ats1.atom_id]["cg_res_id"] = n1["res_id"]
                    self.global_graph.nodes[ats2.atom_id]["cg_bead_id"] = b2
                    self.global_graph.nodes[ats2.atom_id]["cg_res_id"] = n2["res_id"]

                    global_degree[ats1.atom_id] += 1
                    global_degree[ats2.atom_id] += 1
                    # Look for charge to transfer.
                    deg1 = global_degree[ats1.atom_id]
                    deg2 = global_degree[ats2.atom_id]
                    b_key1 = f"{ats1.chain_name}:{ats1.name}"
                    b_key2 = f"{ats2.chain_name}:{ats2.name}"

                    charge_transfer_fragment = None
                    charge_transfer_cfg = None
                    target_atom = None
                    target_atom_key = None
                    if b_key1 in self.charge_transfer:
                        charge_transfer_cfg = self.charge_transfer[b_key1][deg1]
                        target_atom_key = b_key2
                        target_atom = ats2
                        charge_transfer_fragment = self.atom_id2fragment[ats1.atom_id]
                    elif b_key2 in self.charge_transfer:
                        charge_transfer_cfg = self.charge_transfer[b_key2][deg2]
                        target_atom_key = b_key1
                        target_atom = ats1
                        charge_transfer_fragment = self.atom_id2fragment[ats2.atom_id]

                    if charge_transfer_cfg is not None:
                        at_to_transfer = charge_transfer_cfg["when_to_atom"][target_atom_key]
                        # Gets the charge from the source topology.
                        from_chain_name, from_atom_name = charge_transfer_cfg["from_atom"]
                        at_from_transfer = (
                            charge_transfer_fragment.cg_molecule.source_topology.chain_atom_names[
                                from_chain_name
                            ][from_atom_name][0]
                        )
                        # Now update directly hybrid topology with new charge.
                        self.hyb_topology.chains[target_atom.chain_name][target_atom.chain_idx][
                            at_to_transfer
                        ].charge = at_from_transfer.charge
                else:
                    logger.error(
                        f"We could not find any active sites for the CG bonds between beads {b1}-{b2}"
                    )
                    if self.allow_no_bonds:
                        logger.info(f"We skip creation of AT bond between CG beads {b1}-{b2}")
                        with open("missing_cross_bonds.txt", "a+") as missing_cross_f:
                            missing_cross_f.write(f"{b1} - {b2}\n")
                    else:
                        logger.error(f"node 1: {n1} ats 1 {ats1}")
                        logger.error(f"node 2: {n2} ats 2 {ats2}")
                        with open("graph_error.pck", "wb+") as outf:
                            cPickle.dump(self.global_graph, outf, protocol=2)
                        raise RuntimeError("Something is really wrong!")
            sys.stdout.write(f"{100.0 * (progress_indc / progress_indc_total)} %\r")
            progress_indc += 1.0

        # Generate entries for AT cross bonds.
        logger.info(f"Found {len(at_cross_bonds)} atomistic cross bonds")
        self._generate_atomistic_bonds(at_cross_bonds)
        self._remove_atomistic_particles(set(atoms_to_remove))

    def _get_predefined_active_sites(self, b1, b2, global_degree, atoms_to_remove):
        """Look for correct active sites for CG bonds b1-b2."""
        ats1, ats2 = None, None
        old_cg_id1 = self.cg_new_id_old[b1]
        old_cg_id2 = self.cg_new_id_old[b2]
        predefined_as1, predefined_as2 = self.predefined_active_sites.get(
            (old_cg_id1, old_cg_id2), self.predefined_active_sites.get((old_cg_id2, old_cg_id1))
        )
        if not predefined_as1 or not predefined_as2:
            raise RuntimeError(
                f"Cannot read predefined active sites for CG bond {old_cg_id1}-{old_cg_id2}"
            )
        as1 = [s for s in self.cg_active_sites[b1] if s[0].name == predefined_as1]
        if len(as1) > 1:
            raise RuntimeError(f"Found multiple candidates for active site of CG1 {b1}")
        if not as1:
            raise RuntimeError(f"Cannot find active site for CG bond {old_cg_id1}-{old_cg_id2}")
        as2 = [s for s in self.cg_active_sites[b2] if s[0].name == predefined_as2]
        if len(as2) > 1:
            raise RuntimeError(f"Found multiple candidates for active site of CG2 {b2}")
        if not as2:
            raise RuntimeError(f"Cannot find active site for CG bond {old_cg_id1}-{old_cg_id2}")
        at1, max_d1 = as1[0]
        at2, max_d2 = as2[0]

        # for at1, max_d1 in self.cg_active_sites[b1]:
        #    for at6, max_d2 in self.cg_active_sites[b2]:
        b1_key = f"{at1.chain_name}:{at1.name}"
        b2_key = f"{at2.chain_name}:{at2.name}"
        test_bond = self.bond_params.get((b1_key, b2_key))
        if not test_bond:
            test_bond = self.bond_params.get((at1.name, at2.name))
        if (
            test_bond
            and self.global_graph.has_node(at1.atom_id)
            and self.global_graph.has_node(at2.atom_id)
        ):
            at1_deg = global_degree[at1.atom_id]
            at2_deg = global_degree[at2.atom_id]
            at_remove1 = self.atom_id2fragment[b1].active_sites_remove_map.get(b1_key)
            at_remove2 = self.atom_id2fragment[b2].active_sites_remove_map.get(b2_key)
            tmp_atoms_to_remove = []
            valid = True
            # Before we check the degree, let's first try to remove atoms (if they are defined to
            # remove and after that compare the degree.
            if at_remove1:
                for atr1 in at_remove1:
                    atr_chain_name, atr_name = atr1.split(":")[1], atr1.split(":")[2]
                    atr1_id = self.mol_atomname_map[atr_chain_name][at1.chain_idx][atr_name]
                    if not self.global_graph.has_node(atr1_id):
                        valid = False
                    tmp_atoms_to_remove.append(atr1_id)
                    if self.global_graph.has_edge(atr1_id, at1.atom_id):
                        at1_deg -= 1
            if at_remove2 and valid:
                for atr2 in at_remove2:
                    atr_chain_name, atr_name = atr2.split(":")[1], atr2.split(":")[2]
                    atr2_id = self.mol_atomname_map[atr_chain_name][at2.chain_idx][atr_name]
                    if not self.global_graph.has_node(atr2_id):
                        valid = False
                    tmp_atoms_to_remove.append(atr2_id)
                    if self.global_graph.has_edge(atr2_id, at2.atom_id):
                        at2_deg -= 1

            # Check the degree after update with virtual removing of atoms.
            if at1_deg < max_d1 and at2_deg < max_d2 and valid:
                # Found correct pair of active sites. Remove the atoms that were
                # connected to this atom if in settings the set of atoms to remove were defined.
                ats1, ats2 = at1, at2  # Pair of selected atoms.
                atoms_to_remove.extend(tmp_atoms_to_remove)
                for ai in tmp_atoms_to_remove:
                    self.global_graph.remove_node(ai)
                global_degree = dict(self.global_graph.degree())
            else:
                logger.debug(
                    f"{b1_key}({b1})-{b2_key}({b2}) deg1:{at1_deg} < {max_d1} deg2:{at2_deg} < {max_d2} valid: {valid}"
                )
        else:
            logger.debug(f"Params for {b1_key}-{b2_key} not found")
        if ats1 is None or ats2 is None:
            # Found active sites, break the for at1, for at2 loops
            raise RuntimeError("Active sites not found")

        return ats1, ats2, global_degree

    def _search_active_sites(self, b1, b2, global_degree, atoms_to_remove):
        """Look for correct active sites for CG bonds b1-b2."""
        ats1, ats2 = None, None

        def _sort_active_sites(a, b):
            if a[0].atom_id not in global_degree or b[0].atom_id not in global_degree:
                return 0

            at1_deg = global_degree[a[0].atom_id]
            at2_deg = global_degree[b[0].atom_id]
            if at1_deg > at2_deg:
                return 1
            if at1_deg < at2_deg:
                return -1
            return 0

        for at1, max_d1 in sorted(self.cg_active_sites[b1], key=cmp_to_key(_sort_active_sites)):
            for at2, max_d2 in sorted(self.cg_active_sites[b2], key=cmp_to_key(_sort_active_sites)):
                b1_key = f"{at1.chain_name}:{at1.name}"
                b2_key = f"{at2.chain_name}:{at2.name}"
                if self.restricted_cross_bond_patterns and (
                    (b1_key, b2_key) in self.restricted_cross_bond_patterns
                    or (at1.name, at2.name) in self.restricted_cross_bond_patterns
                ):
                    logger.info(f"Bond {b1_key}-{b2_key} not used due to the restriction")
                    continue

                test_bond = self.bond_params.get((b1_key, b2_key))
                if not test_bond:
                    test_bond = self.bond_params.get((at1.name, at2.name))
                if not test_bond:
                    logger.info(f"params for {b1_key} {b2_key} not found - skip")
                    continue

                if (
                    test_bond
                    and self.global_graph.has_node(at1.atom_id)
                    and self.global_graph.has_node(at2.atom_id)
                ):
                    at1_deg = global_degree[at1.atom_id]
                    at2_deg = global_degree[at2.atom_id]
                    at_remove1 = self.atom_id2fragment[b1].active_sites_remove_map.get(b1_key)
                    at_remove2 = self.atom_id2fragment[b2].active_sites_remove_map.get(b2_key)
                    tmp_atoms_to_remove = []
                    valid = True
                    # Before we check the degree, let's first try to remove atoms (if they are defined to
                    # remove and after that compare the degree.
                    if at_remove1:
                        for atr1 in at_remove1:
                            atr_chain_name, atr_name = atr1.split(":")[1], atr1.split(":")[2]
                            atr1_id = self.mol_atomname_map[atr_chain_name][at1.chain_idx][atr_name]
                            if not self.global_graph.has_node(atr1_id):
                                valid = False
                                break
                            tmp_atoms_to_remove.append(atr1_id)
                            if self.global_graph.has_edge(atr1_id, at1.atom_id):
                                at1_deg -= 1
                    if at_remove2 and valid:
                        for atr2 in at_remove2:
                            atr_chain_name, atr_name = atr2.split(":")[1], atr2.split(":")[2]
                            atr2_id = self.mol_atomname_map[atr_chain_name][at2.chain_idx][atr_name]
                            if not self.global_graph.has_node(atr2_id):
                                valid = False
                                break
                            tmp_atoms_to_remove.append(atr2_id)
                            if self.global_graph.has_edge(atr2_id, at2.atom_id):
                                at2_deg -= 1

                    # Check the degree after update with virtual removing of atoms.
                    if at1_deg < max_d1 and at2_deg < max_d2 and valid:
                        # Found correct pair of active sites. Remove the atoms that were
                        # connected to this atom if in settings the set of atoms to remove were defined.
                        ats1, ats2 = at1, at2  # Pair of selected atoms.
                        atoms_to_remove.extend(tmp_atoms_to_remove)
                        for ai in tmp_atoms_to_remove:
                            self.global_graph.remove_node(ai)
                        global_degree = dict(self.global_graph.degree())
                        break
                    valid = False
                    logger.debug(
                        f"{b1_key}({self.cg_new_id_old[b1]})-{b2_key}({self.cg_new_id_old[b2]}) deg1:{at1_deg} < {max_d1} deg2:{at2_deg} < {max_d2} valid: {valid}"
                    )
                    logger.debug(f"active site b1: {self.cg_active_sites[b1]}")
                    logger.debug(f"active site b2: {self.cg_active_sites[b2]}")
                    if at1_deg >= max_d1:
                        logger.debug(f"AT1 has degree {at1_deg} >= {max_d1} ({at1})")
                        nbs_1 = self.global_graph[at1.atom_id]
                        logger.debug(f"Found neighbours of {at1}: ")
                        for nb in nbs_1:
                            logger.debug(self.global_graph.nodes[nb])
                    if at2_deg >= max_d2:
                        logger.debug(f"AT2 has degree {at2_deg} >= {max_d2} ({at2})")
                        nbs_2 = self.global_graph[at2.atom_id]
                        logger.debug(f"Found neighbours of {at2}: ")
                        for nb in nbs_2:
                            logger.debug(self.global_graph.nodes[nb])
                else:
                    logger.info(f"params for {b1_key}-{b2_key} not found")
            if ats1 is not None and ats2 is not None:
                # Found active sites, break the for at1, for at2 loops
                break

        return ats1, ats2, global_degree

    def _generate_atomistic_bonds(self, at_cross_bonds):
        """Generates parameters for atomistic bonds."""
        fout_filename = "missing_definitions.txt"
        fout = open(fout_filename, "w")
        missing_definitions = set()

        bond_key = f"{self.cross_prefix}bonds"
        angle_key = f"{self.cross_prefix}angles"
        dihedral_key = f"{self.cross_prefix}dihedrals"

        for b1, b2 in at_cross_bonds:
            n1 = self.global_graph.nodes[b1]
            n2 = self.global_graph.nodes[b2]
            b1_key = n1_key = "{}:{}".format(n1["chain_name"], n1["name"])
            b2_key = n2_key = "{}:{}".format(n2["chain_name"], n2["name"])
            s1_key, s2_key = n1["name"], n2["name"]
            if (b1, b2) in self.hyb_topology.new_data[bond_key] or (
                b2,
                b1,
            ) in self.hyb_topology.new_data[bond_key]:
                raise RuntimeError(f"Bond {b1}-{b2} already defined in the topology, ??")

            param = _get_params(
                self.bond_params, [(n1_key, n2_key), (s1_key, s2_key)], raise_exception=True
            )
            self.hyb_topology.new_data[bond_key][(b1, b2)] = param.params
            self.at_cross_bonds[(b1, b2)] = param
            # Generate angles.
            triplets = tools.gen_bonded_tuples(self.global_graph, 3, (b1, b2))
            quadruplets = tools.gen_bonded_tuples(self.global_graph, 4, (b1, b2))
            for triplet in triplets:
                n1_key, n2_key, n3_key = [
                    "{}:{}".format(x["chain_name"], x["name"])
                    for x in map(self.global_graph.nodes.get, triplet)
                ]
                s1_key, s2_key, s3_key = [
                    x["name"] for x in map(self.global_graph.nodes.get, triplet)
                ]

                triplet_defined = triplet in self.hyb_topology.new_data[angle_key]
                reverse_triplet_defined = (
                    tuple(reversed(triplet)) in self.hyb_topology.new_data[angle_key]
                )
                # Here the exception can be raised with missing definition
                if not (triplet_defined or reverse_triplet_defined):
                    param = _get_params(
                        self.angle_params,
                        [(n1_key, n2_key, n3_key), (s1_key, s2_key, s3_key)],
                        raise_exception=False,
                    )
                    if param:
                        self.hyb_topology.new_data[angle_key][triplet] = param.params
                        self.at_cross_angles[triplet] = param
                    else:
                        logger.info(
                            f"Missing definition of angle: {n1_key}-{n2_key}-{n3_key} (bond {b1_key}-{b2_key})"
                        )
                        missing_definitions.add(tuple([n1_key, n2_key, n3_key] + list(triplet)))
            # Generate dihedrals.
            for quadruplet in quadruplets:
                n1_key, n2_key, n3_key, n4_key = [
                    "{}:{}".format(x["chain_name"], x["name"])
                    for x in map(self.global_graph.nodes.get, quadruplet)
                ]
                s_keys = [x["name"] for x in map(self.global_graph.nodes.get, quadruplet)]
                if not (
                    quadruplet in self.hyb_topology.new_data[dihedral_key]
                    or tuple(reversed(quadruplet)) in self.hyb_topology.new_data[dihedral_key]
                ):
                    param = _get_params(
                        self.dihedral_params,
                        [(n1_key, n2_key, n3_key, n4_key), s_keys],
                        raise_exception=False,
                    )
                    if param:
                        self.hyb_topology.new_data[dihedral_key][quadruplet] = param.params
                        self.at_cross_dihedrals[quadruplet] = param
                    else:
                        logger.info(
                            f"Missing definition of dihedral: {n1_key}-{n2_key}-{n3_key}-{n4_key} (bond {b1_key}-{b2_key})"
                        )
                        missing_definitions.add(
                            tuple([n1_key, n2_key, n3_key, n4_key] + list(quadruplet))
                        )

        fout.writelines(
            "\n".join(
                [" ".join(map(str, x)) for x in sorted(missing_definitions, key=lambda l: len(l))]
            )
        )
        logger.info(f"Wrote missing definitions in {fout_filename}")
        fout.close()

    def _remove_atomistic_particles(self, atoms_to_remove):
        """Update coordinate and topology file by removing atoms and renumbering"""
        if atoms_to_remove:
            logger.info(
                f"Clean up atomistic particles after creating bonds, atoms to remove: {len(atoms_to_remove)}"
            )
            self.hyb_topology.remove_atoms(atoms_to_remove, renumber=False)
            self.hybrid_configuration["file"].remove_atoms(atoms_to_remove, renumber=False)
            # Renumber data files
            self.hyb_topology.renumber()
            self.hybrid_configuration["file"].renumber()

    def _generate_exclusion_lists(self):
        # Collect all bonds and greate global_graph again...
        at_g = networkx.Graph()
        cg_g = networkx.Graph()

        bonds = [
            b
            for b in self.hyb_topology.bonds.keys()
            if b[0] in self.atom_ids and b[1] in self.atom_ids
        ]
        for k in self.hyb_topology.new_data:
            if "bonds" in k:
                bonds.extend(
                    [
                        b
                        for b in self.hyb_topology.new_data[k]
                        if b[0] in self.atom_ids and b[1] in self.atom_ids
                    ]
                )
        at_g.add_edges_from(bonds)

        bonds = [
            b
            for b in self.hyb_topology.bonds.keys()
            if b[0] in self.cg2atom and b[1] in self.cg2atom
        ]
        for k in self.hyb_topology.new_data:
            if "bonds" in k:
                bonds.extend(
                    [
                        b
                        for b in self.hyb_topology.new_data[k]
                        if b[0] in self.cg2atom and b[1] in self.cg2atom
                    ]
                )
        cg_g.add_edges_from(bonds)

        logger.info(
            "Generating exclusion lists AT, nrexcl={}".format(
                self.hyb_topology.moleculetype["excl_at"]
            )
        )
        paths = dict(
            networkx.all_pairs_shortest_path(at_g, int(self.hyb_topology.moleculetype["excl_at"]))
        )
        exclusions = set()
        for l in paths.values():
            for p in l.values():
                if len(p) > 1:
                    exclusions.add(tuple(sorted([p[0], p[-1]])))

        logger.info(
            "Generating exclusion lists CG, nrexcl={}".format(
                self.hyb_topology.moleculetype["excl_cg"]
            )
        )
        paths = dict(
            networkx.all_pairs_shortest_path(cg_g, int(self.hyb_topology.moleculetype["excl_cg"]))
        )
        for l in paths.values():
            for p in l.values():
                if len(p) > 1:
                    exclusions.add(tuple(sorted([p[0], p[-1]])))

        output_filename = "exclusion_{}.list".format(self.hyb_topology.file_name.split(".")[0])
        out_file = open(output_filename, "w")
        out_file.writelines("\n".join(["{} {}".format(*d) for d in sorted(exclusions)]))
        out_file.close()
        logger.info(f"Generated {len(exclusions)} exclusions, writen to {output_filename}")
