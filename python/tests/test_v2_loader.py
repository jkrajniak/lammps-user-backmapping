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
