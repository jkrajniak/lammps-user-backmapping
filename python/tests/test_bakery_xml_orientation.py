"""`build` passes fragment_orientation on the bakery_xml route, as `build-hybrid` does."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import pytest

from backmap_prep.network import api
from backmap_prep.schema import load_settings

if TYPE_CHECKING:
    from pathlib import Path

SETTINGS = """\
prep:
  engine: network
  bakery_xml: settings.xml
  fragment_orientation: align
hybrid:
  coordinates: hyb_conf.gro
  topology: hyb_topol.top
"""


class _StopError(Exception):
    pass


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> list[dict[str, Any]]:
    """Record the arguments of build_hybrid_gromacs and stop there."""
    (tmp_path / "settings.xml").write_text("<settings/>\n")
    (tmp_path / "settings.yaml").write_text(SETTINGS)
    seen: list[dict[str, Any]] = []

    def fake(*args: Any, **kwargs: Any) -> None:
        seen.append(kwargs)
        raise _StopError

    monkeypatch.setattr(api, "build_hybrid_gromacs", fake)
    return seen


def test_build_passes_fragment_orientation(calls: list[dict[str, Any]], tmp_path: Path) -> None:
    path = tmp_path / "settings.yaml"
    with pytest.raises(_StopError):
        api.build_network_lammps(load_settings(path), path, output_dir=tmp_path)
    assert calls[-1]["fragment_orientation"] == "align"
