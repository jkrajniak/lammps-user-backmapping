"""Unit tests for Settings v2 → bakery XML conversion."""

from __future__ import annotations

from pathlib import Path

import yaml

from backmap_prep.network.v2_loader import settings_to_xml_root
from backmap_prep.schema import load_settings


def test_load_settings_v2_unified_source(tmp_path: Path) -> None:
    raw = {
        "version": 2,
        "prep": {"engine": "network"},
        "molecules": [
            {
                "name": "MOL",
                "source": [
                    {
                        "molecule_degree": 0,
                        "coordinates": "a.gro",
                        "topology": "a.itp",
                    }
                ],
                "beads": [
                    {
                        "name": "B1",
                        "type": "B",
                        "mapping": [{"degree": 1, "atoms": ["C1", "C2"]}],
                    }
                ],
            }
        ],
        "cg_system": {"coordinates": "cg.gro", "topology": "cg.top"},
        "hybrid": {"coordinates": "hyb.gro", "topology": "hyb.top"},
    }
    path = tmp_path / "settings.yaml"
    path.write_text(yaml.dump(raw))
    settings = load_settings(path)
    root = settings_to_xml_root(settings)
    bead = root.find(".//beads[@degree='1']")
    assert bead is not None
    assert "1:MOL:C1" in (bead.text or "")


def test_cross_bond_params_default_from_comment() -> None:
    from backmap_prep.schema import Settings

    settings = Settings(
        prep={"engine": "network"},
        molecules=[
            {
                "name": "EPO",
                "source": {"coordinates": "a.gro", "topology": "a.itp"},
                "beads": [
                    {"name": "A1", "type": "A", "mapping": [{"degree": 1, "atoms": ["C1"]}]},
                ],
            }
        ],
        cg_system={"coordinates": "cg.gro", "topology": "cg.top"},
        hybrid={"coordinates": "hyb.gro", "topology": "hyb.top"},
        cross_interactions={
            "bonds": [{"comment": "EPO-IPD", "pairs": [["EPO:C1", "IPD:N1"]]}],
        },
    )
    root = settings_to_xml_root(settings)
    bonds = root.find(".//bonds")
    assert bonds is not None
    assert bonds.attrib["params"] == "; EPO-IPD"


def test_charge_map_and_equilibrate_emitted() -> None:
    from backmap_prep.schema import Settings

    settings = Settings(
        prep={"engine": "network"},
        molecules=[
            {
                "name": "DIO",
                "ident": "DIO",
                "source": {"coordinates": "a.gro", "topology": "a.itp"},
                "charge_management": {"equilibrate": True},
                "beads": [
                    {
                        "name": "E1",
                        "type": "E",
                        "mapping": [
                            {
                                "degree": 1,
                                "atoms": ["C1", "C2", "O2", "H2", "H3", "H4", "H5", "H6"],
                            },
                            {
                                "degree": 2,
                                "atoms": ["C1", "C2", "H2", "H3", "H4", "H5"],
                                "charge_map": [0.0, 0.0, "*", "*", "*", "*"],
                            },
                        ],
                    }
                ],
            }
        ],
        cg_system={"coordinates": "cg.gro", "topology": "cg.top"},
        hybrid={"coordinates": "hyb.gro", "topology": "hyb.top"},
    )
    root = settings_to_xml_root(settings)
    mol = root.find(".//cg_molecule")
    assert mol is not None
    assert mol.attrib.get("equilibrate_charges") == "1"
    beads_deg2 = root.find(".//beads[@degree='2']")
    assert beads_deg2 is not None
    cm = beads_deg2.find("charge_map")
    assert cm is not None
    assert cm.text == "0.0 0.0 * * * *"
    # degree-1 beads must NOT carry a charge_map
    beads_deg1 = root.find(".//beads[@degree='1']")
    assert beads_deg1 is not None
    assert beads_deg1.find("charge_map") is None


def test_cross_dihedrals_emitted() -> None:
    from backmap_prep.schema import Settings

    settings = Settings(
        prep={"engine": "network"},
        molecules=[
            {
                "name": "TER",
                "source": {"coordinates": "a.gro", "topology": "a.itp"},
                "beads": [
                    {"name": "A1", "type": "A", "mapping": [{"degree": 1, "atoms": ["O1"]}]},
                    {"name": "B1", "type": "B", "atoms": ["C2"]},
                    {"name": "A2", "type": "A", "mapping": [{"degree": 1, "atoms": ["O3"]}]},
                    {
                        "name": "Q1",
                        "type": "C",
                        "mapping": [{"degree": 2, "atoms": ["O1"]}, {"degree": 2, "atoms": ["O1"]}],
                    },
                    {
                        "name": "Q2",
                        "type": "C",
                        "mapping": [{"degree": 2, "atoms": ["O4"]}, {"degree": 2, "atoms": ["O4"]}],
                    },
                ],
            }
        ],
        cg_system={"coordinates": "cg.gro", "topology": "cg.top"},
        hybrid={"coordinates": "hyb.gro", "topology": "hyb.top"},
        cross_interactions={
            "dihedrals": [
                {
                    "params": "8 0 1.0",
                    "quadruples": [
                        ["TER:A1", "TER:B1", "TER:Q1", "DIO:E1"],
                        ["DIO:E1", "TER:Q1", "TER:B1", "TER:A1"],
                    ],
                }
            ]
        },
    )
    root = settings_to_xml_root(settings)
    dih = root.find(".//dihedrals")
    assert dih is not None
    assert dih.attrib["params"] == "8 0 1.0"
    tokens = (dih.text or "").split()
    assert "TER:A1" in tokens
    assert len(tokens) % 4 == 0
