# Crosslinked Melamine (MF) Network Example Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `examples/melamine_network/` — the same 500-molecule melamine-formaldehyde
(MF) system as `examples/melamine/`, but crosslinked to match bakery's actual reference
network (675 real inter-molecular bonds, confirmed by direct inspection of
`cg_topol.top`), built through the already-proven `network-backmap-prep` engine, run
through LAMMPS on the VM, and RDF-compared against the same reference `.xvg` files —
for the first time a fair, like-for-like comparison.

**Architecture:** Reuse the existing `prep.bakery_xml` passthrough (network engine),
which already vendors bakery's own Python topology-building code
(`backmap_prep/network/bakery/structures.py`) and already computes each CG bead's
degree (crosslinked vs. not) directly from the actual bond graph of whatever CG
topology file is supplied. Vendor bakery's MF asset files (settings.xml, at_topol.top,
cg_topol.top with its real 675 crosslink bonds, cg_conf.gro, hyb_topol.top, CG tables)
into `examples/melamine_network/large/` so the example is self-contained (not
dependent on the `bakery_full_unused` VM directory it was sourced from — that
directory's name signals it isn't a stable long-term dependency). No new Python
generator code is required: `_resolve_mapping()` in `network/v2_loader.py` already
supports the `base`/`remove`/`active_sites` delta pattern bakery's own
`degree="2"`/`degree="3"` beads use, but since we're using the `bakery_xml`
passthrough (not v2-native YAML), that loader code path isn't even exercised — bakery's
`settings.xml` is consumed directly by `structures.py`.

**Tech Stack:** Python 3.10+ (backmap-prep, pydantic, networkx), GROMACS-format
topology/coordinate files, LAMMPS (C++ package already built on the VM), pytest.

## Global Constraints

- All LAMMPS compilation/execution happens on the VM (`sc@<vm_host>`), in `tmux` —
  never locally. (Established project practice, followed throughout this plan.)
- `examples/melamine/` (the uncrosslinked example) and
  `openspec/specs/example-melamine/spec.md` are never modified by this work.
- No live/reactive bond formation — the 675-bond network is imported as a fixed,
  static topology exactly as bakery's `cg_topol.top` already defines it.
- No ESPResSo++ build — confirmed not built/importable on the VM; only the
  already-vendored pure-Python `backmap_prep/network/bakery/structures.py` is used.
- Success criterion is a stable, structurally faithful run — RDF results are computed
  and reported honestly, with no committed target pass-count to hit before calling a
  task "done."
- Commit after every task (see each task's final step). Do not batch commits across
  tasks.

---

### Task 1: Vendor bakery's MF reference assets into the new example directory

**Files:**
- Create: `examples/melamine_network/large/` (new directory)
- Create: `examples/melamine_network/large/settings.xml` (copied from bakery, byte-identical)
- Create: `examples/melamine_network/large/at_topol.top` (copied from bakery)
- Create: `examples/melamine_network/large/single_mf.gro` (copied from bakery)
- Create: `examples/melamine_network/large/cg_conf.gro` (copied from bakery — has the real crosslinked network's coordinates)
- Create: `examples/melamine_network/large/cg_topol.top` (copied from bakery — has the real 675 crosslink bonds)
- Create: `examples/melamine_network/large/hyb_topol.top` (copied from bakery — needed by `structures.py` as an output-naming/include template)
- Create: `examples/melamine_network/large/table_A_A.xvg`, `table_b1.xvg` (copied from bakery's CG tabulated potentials)
- Create: `examples/melamine_network/large/README.md`

**Interfaces:**
- Produces: the on-disk asset set that Task 2's `settings.yaml` references via
  relative paths (`prep.bakery_xml: settings.xml`, `prep.data_dir: .`).

- [ ] **Step 1: Locate and pin the bakery MF reference data path on the VM**

Run on the VM (`sc@<vm_host>`):

```bash
ssh sc@<vm_host> "ls /home/sc/sc/e++/bakery_full_unused/examples/network_backmapping/mf/backmapping/"
```

Expected: the file list confirmed during design (`settings.xml`, `at_topol.top`,
`single_mf.gro`, `cg_conf.gro`, `cg_topol.top`, `hyb_topol.top`, `table_A_A.xvg`,
`table_b1.xvg`, `table_A_A.pot`, `table_b1.pot`, plus GROMACS `.mdp`/`.ndx` files we
don't need). If any of the 8 files listed under "Files" above is missing, stop and
report — do not substitute a different source.

- [ ] **Step 2: Copy the needed files from VM to local**

```bash
mkdir -p examples/melamine_network/large
scp sc@<vm_host>:/home/sc/sc/e++/bakery_full_unused/examples/network_backmapping/mf/backmapping/settings.xml \
    sc@<vm_host>:/home/sc/sc/e++/bakery_full_unused/examples/network_backmapping/mf/backmapping/at_topol.top \
    sc@<vm_host>:/home/sc/sc/e++/bakery_full_unused/examples/network_backmapping/mf/backmapping/single_mf.gro \
    sc@<vm_host>:/home/sc/sc/e++/bakery_full_unused/examples/network_backmapping/mf/backmapping/cg_conf.gro \
    sc@<vm_host>:/home/sc/sc/e++/bakery_full_unused/examples/network_backmapping/mf/backmapping/cg_topol.top \
    sc@<vm_host>:/home/sc/sc/e++/bakery_full_unused/examples/network_backmapping/mf/backmapping/hyb_topol.top \
    sc@<vm_host>:/home/sc/sc/e++/bakery_full_unused/examples/network_backmapping/mf/backmapping/table_A_A.xvg \
    sc@<vm_host>:/home/sc/sc/e++/bakery_full_unused/examples/network_backmapping/mf/backmapping/table_b1.xvg \
    examples/melamine_network/large/
```

- [ ] **Step 3: Verify the crosslink count survived the copy**

```bash
python3 -c "
lines = open('examples/melamine_network/large/cg_topol.top').read()
static = lines.count('; static')
chem = lines.count('; chem')
print(f'static={static} chem={chem}')
assert static == 1500, f'expected 1500 static bonds, got {static}'
assert chem == 675, f'expected 675 chem (crosslink) bonds, got {chem}'
print('OK: crosslink network intact')
"
```

Expected: `static=1500 chem=675`, `OK: crosslink network intact`.

- [ ] **Step 4: Write `examples/melamine_network/large/README.md`**

```markdown
# Melamine-formaldehyde (MF) — crosslinked network example

Same 500-molecule melamine-formaldehyde system as `examples/melamine/`, but
crosslinked to match bakery's own reference network exactly: 675 real
inter-molecular covalent bonds across 1500 CG beads (~45% of each molecule's
3 pendant hydroxymethyl arms condensed into an ether bridge with another
molecule), vendored from
`bakery_full_unused/examples/network_backmapping/mf/backmapping/` (see
`research/experiments/20260802_melamine-bakery-rerun.md`, 2026-08-06 update,
for how this crosslinking mismatch was found).

Unlike `examples/melamine/`, this network is imported as a fixed, static
topology (bakery's `cg_topol.top`) — not generated by live reaction dynamics
in this pipeline.

## Files

- `settings.xml`, `at_topol.top`, `single_mf.gro`, `cg_conf.gro`,
  `cg_topol.top`, `hyb_topol.top`, `table_A_A.xvg`, `table_b1.xvg` — vendored
  unmodified from bakery's reference (see provenance above).
- `settings.yaml` — backmap-prep entry point (`prep.bakery_xml` passthrough).
- Generated: `hyb_conf.gro`, `hyb_topol.top` (Tier A), `melamine_network.data`,
  `in.melamine_network` (LAMMPS).

## Running

```
uv run backmap-prep build-hybrid examples/melamine_network/large/settings.yaml
```

then (VM only): `lmp -in in.melamine_network -var eq_steps ... -var ramp_steps ... -var prod_steps ...`
```

- [ ] **Step 5: Commit**

```bash
git add examples/melamine_network/large/
git commit -m "feat(examples): vendor bakery's crosslinked MF reference assets"
```

---

### Task 2: `settings.yaml` using the `bakery_xml` passthrough

**Files:**
- Create: `examples/melamine_network/large/settings.yaml`

**Interfaces:**
- Consumes: `PrepConfig` fields `engine`, `bakery_xml`, `data_dir`, `tables_dir`,
  `forcefield_dir` (`python/src/backmap_prep/schema.py:279-304`); CLI resolution in
  `_resolve_bakery_xml` (`python/src/backmap_prep/cli.py:360`).
- Produces: a `Settings` object with `prep.bakery_xml` set, consumed by Task 3's
  `build_hybrid_gromacs()` call.

- [ ] **Step 1: Write the settings file**

```yaml
# examples/melamine_network/large/settings.yaml
# Crosslinked melamine (MF) network -- bakery_xml passthrough.
# Assets vendored from bakery's own network_backmapping/mf/backmapping/
# reference (see README.md for provenance). No new molecule/bead YAML is
# authored here -- bakery's settings.xml already defines the degree=2
# (unreacted) / degree=3 (reacted) bead templates, and structures.py computes
# each CG bead's actual degree from cg_topol.top's real bond graph.

version: 2

prep:
  engine: network
  bakery_xml: settings.xml
  data_dir: .
  tables_dir: .
  forcefield_dir: ../../epoxy/forcefield
  chain_rng_seed: 42

output:
  prefix: melamine_network
  format: lammps
  units: real
```

- [ ] **Step 2: Verify it loads without error**

```bash
uv run python3 -c "
from pathlib import Path
from backmap_prep.schema import load_settings
s = load_settings(Path('examples/melamine_network/large/settings.yaml'))
print('engine:', s.prep.engine)
print('bakery_xml:', s.prep.bakery_xml)
assert s.prep.bakery_xml == 'settings.xml'
print('OK')
"
```

Expected: `engine: network`, `bakery_xml: settings.xml`, `OK`.

- [ ] **Step 3: Commit**

```bash
git add examples/melamine_network/large/settings.yaml
git commit -m "feat(examples): add settings.yaml for crosslinked MF network"
```

---

### Task 3: Tier A parity check (hybrid topology generation)

**Files:**
- Create: `python/tests/test_melamine_network.py`

**Interfaces:**
- Consumes: `build_hybrid_gromacs(xml_path, base_dir, chain_rng_seed) -> HybridBuildResult`
  (`python/src/backmap_prep/network/api.py`), same signature already exercised by
  `python/tests/test_network.py::test_rim135_build_hybrid_smoke`.
- `HybridBuildResult` fields used: `.n_atoms`, `.coordinates_path`, `.topology_path`,
  `.missing_definitions_path`.

- [ ] **Step 1: Write the failing test**

```python
# python/tests/test_melamine_network.py
"""Tier A parity: crosslinked MF network hybrid topology generation."""

from __future__ import annotations

from pathlib import Path

import pytest

from backmap_prep.network.api import build_hybrid_gromacs

MF_NETWORK_DIR = Path(__file__).resolve().parents[2] / "examples" / "melamine_network" / "large"


def test_mf_network_assets_present() -> None:
    """Vendored bakery assets exist before attempting a build."""
    for name in (
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
    xml_path = MF_NETWORK_DIR / "settings.xml"
    result = build_hybrid_gromacs(xml_path, base_dir=MF_NETWORK_DIR, chain_rng_seed=42)

    assert result.coordinates_path.is_file()
    assert result.topology_path.is_file()
    assert result.coordinates_path.stat().st_size > 0
    assert result.topology_path.stat().st_size > 0
    # 500 molecules x 27 AT atoms, minus 2 atoms removed per crosslink (H2O
    # leaving group is 1 O + varies -- exact count confirmed by this test,
    # not assumed; see Step 2 for how to read the actual value first.
    assert result.n_atoms > 0


def test_mf_network_charge_conserved() -> None:
    """Total system charge is conserved within 1e-6 e."""
    xml_path = MF_NETWORK_DIR / "settings.xml"
    result = build_hybrid_gromacs(xml_path, base_dir=MF_NETWORK_DIR, chain_rng_seed=42)

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
    example-melamine-network, 'Crosslink count matches bakery's reference')."""
    xml_path = MF_NETWORK_DIR / "settings.xml"
    result = build_hybrid_gromacs(xml_path, base_dir=MF_NETWORK_DIR, chain_rng_seed=42)

    in_bonds = False
    n_bonds = 0
    for line in result.topology_path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("[ bonds ]"):
            in_bonds = True
            continue
        if stripped.startswith("["):
            in_bonds = False
            continue
        if in_bonds and stripped and not stripped.startswith(";"):
            n_bonds += 1
    # Total AT bonds = (per-molecule intramolecular bonds x 500) + 675 crosslink
    # bonds. Per-molecule bond count comes from examples/melamine/large/topol_aa.top
    # (already known: 26 bonds/molecule for the unreacted case -- but crosslinked
    # molecules lose bonds to their removed leaving-group atoms too). Read the
    # actual total here first (informational), then hand-verify the crosslink
    # component specifically:
    print(f"total AT bonds in generated topology: {n_bonds}")
    # The 675 crosslink bonds are additional to whatever the per-molecule count
    # is; assert the total is large enough to plausibly include them as a first
    # check, then tighten in Step 3 below once the exact per-molecule bond count
    # for reacted vs. unreacted beads is confirmed from this run's own output.
    assert n_bonds > 500 * 26, "bond count too low to plausibly include the 675 crosslinks"


def test_mf_network_charges_subset_of_uncrosslinked_example() -> None:
    """Every per-atom charge value appearing in the crosslinked network's
    generated topology also appears in topol_aa.top's per-atom charge list
    (the same 27-atom template examples/melamine/ uses) -- i.e. unreacted
    arms are built from unmodified template charges, not a separately
    perturbed set (spec: example-melamine-network, 'Unreacted arms match the
    uncrosslinked example')."""
    template_charges: set[str] = set()
    template_top = MF_NETWORK_DIR / "at_topol.top"
    in_atoms = False
    for line in template_top.read_text().splitlines():
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
                template_charges.add(f"{float(parts[6]):.6f}")

    xml_path = MF_NETWORK_DIR / "settings.xml"
    result = build_hybrid_gromacs(xml_path, base_dir=MF_NETWORK_DIR, chain_rng_seed=42)
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

    unexpected = generated_charges - template_charges
    assert not unexpected, (
        f"generated topology has {len(unexpected)} charge values not present "
        f"in the unreacted template -- unreacted arms may have been perturbed: "
        f"{sorted(unexpected)[:10]}"
    )
```

- [ ] **Step 2: Run all 5 tests, record the real atom count and bond count (some
      will fail informatively — that's expected on this first run, not an error)**

```bash
uv run pytest python/tests/test_melamine_network.py -v -s 2>&1 | tail -60
```

Expected: `test_mf_network_assets_present` PASSES (if Task 1 completed correctly).
`test_mf_network_build_hybrid_smoke` PASSES (its only assertion is `n_atoms > 0`)
— read the actual `result.n_atoms` value from output (add a temporary
`print(result.n_atoms)` if not visible) and record it.
`test_mf_network_charge_conserved` PASSES or FAILS informatively (a real charge
bug, if so — investigate before proceeding). `test_mf_network_crosslink_bond_count`
PASSES (coarse check) and prints the real total AT bond count — record it.
`test_mf_network_charges_subset_of_uncrosslinked_example` PASSES (if it fails,
that's a real finding — unreacted-arm charges were perturbed somewhere in the
generation path, investigate before proceeding, don't loosen the assertion).

- [ ] **Step 3: Tighten the atom-count assertion to the real value**

Replace `assert result.n_atoms > 0` in `test_mf_network_build_hybrid_smoke` with
the exact value recorded in Step 2:

```python
    assert result.n_atoms == <exact_value_from_step_2>
```

- [ ] **Step 4: Tighten `test_mf_network_crosslink_bond_count` to the exact value**

Using the real total AT bond count recorded in Step 2, and a one-off `grep -c`
on `at_topol.top`'s own `[ bonds ]` section for the per-molecule template bond
count, confirm by hand that
`generated_total - 500 * <template bond count> == 675` (accounting for any bonds
lost at reacted sites where a leaving-group atom's own bonds are removed too,
e.g. the removed O-H bond disappears along with the removed H — so the
per-crosslink delta may not be exactly +1 bond; compute and record the actual
delta before asserting on it). Replace the coarse
`assert n_bonds > 500 * 26` with the exact equality assertion once confirmed.

- [ ] **Step 5: Check `missing_definitions_path` — this is Task 4's input, not a failure here**

```bash
uv run python3 -c "
from pathlib import Path
from backmap_prep.network.api import build_hybrid_gromacs
r = build_hybrid_gromacs(
    Path('examples/melamine_network/large/settings.xml'),
    base_dir=Path('examples/melamine_network/large'),
    chain_rng_seed=42,
)
if r.missing_definitions_path and r.missing_definitions_path.is_file():
    print(r.missing_definitions_path.read_text())
else:
    print('no missing definitions')
"
```

Expected: either "no missing definitions" (meaning `forcefield_dir`'s OPLS-AA table
already resolved everything — go straight to Task 5 after this), or a list of
missing angle/dihedral types at crosslink sites (expected, matching bakery's own
`missing_definitions.txt` gap) — hand this output to Task 4.

- [ ] **Step 6: Run full test file, confirm all 5 tests pass**

```bash
uv run pytest python/tests/test_melamine_network.py -v
```

Expected: all 5 tests PASS (charge conservation, tightened atom count, tightened
crosslink bond count, and unreacted-charge-subset check all included).

- [ ] **Step 7: Commit**

```bash
git add python/tests/test_melamine_network.py
git commit -m "test(network): add Tier A parity tests for crosslinked MF example"
```

---

### Task 4: Resolve crosslink-site force-field parameters

**Files:**
- Modify: `examples/epoxy/forcefield/oplsaa.ff/ffbonded.itp` (only if the check
  below finds a genuine gap — this file is shared with epoxy and melamine, so any
  edit must be additive, never removing/changing existing entries)
- Modify: `python/tests/test_melamine_network.py` (extend, don't replace)

**Interfaces:**
- Consumes: Task 3's `missing_definitions_path` output.

- [ ] **Step 1: Check whether the shared OPLS-AA table already covers the gap**

Using the specific missing angle/dihedral atom-type triples/quadruples from Task 3
Step 4's output (e.g. if it lists `opls_157 opls_154 opls_157` for a C-O-C angle),
grep the shared force field directly:

```bash
grep -A2 "\[ angletypes \]" examples/epoxy/forcefield/oplsaa.ff/ffbonded.itp | head -5
grep "opls_157.*opls_154\|opls_154.*opls_157" examples/epoxy/forcefield/oplsaa.ff/ffbonded.itp
grep "opls_157.*opls_154.*opls_157\|opls_157.*opls_157" examples/epoxy/forcefield/oplsaa.ff/ffbonded.itp
```

(Substitute the actual missing atom types from Task 3 Step 4's output if they differ
from `opls_157`/`opls_154` — those are the bridge-C/hydroxyl-O types identified
earlier in `topol_aa.top`, but confirm against the real
`missing_definitions.txt`-equivalent output, don't assume.)

- [ ] **Step 2a (if Step 1 finds existing entries):** No force-field edit needed —
      re-run `build_hybrid_gromacs` and confirm `missing_definitions_path` is now
      empty or absent:

```bash
uv run python3 -c "
from pathlib import Path
from backmap_prep.network.api import build_hybrid_gromacs
r = build_hybrid_gromacs(
    Path('examples/melamine_network/large/settings.xml'),
    base_dir=Path('examples/melamine_network/large'),
    chain_rng_seed=42,
)
print('missing:', r.missing_definitions_path)
"
```

Expected: `missing: None` or a file that, when read, is empty.

- [ ] **Step 2b (if Step 1 finds a genuine gap):** Add the missing generic
      angletype/dihedraltype entries to `ffbonded.itp`, sourced from the standard
      published OPLS-AA parameter set for the matching chemical environment
      (aliphatic ether C-O-C angle / C-C-O-C dihedral — the same functional group
      as any dialkyl ether, not melamine-specific). Do not invent numeric values —
      if the standard OPLS-AA table entry for this exact type combination cannot be
      confirmed, use the closest documented generic aliphatic-ether entry, comment
      the source inline, and flag it explicitly in the commit message and in
      `research/decisions/` (see Step 4) as an approximation pending review — this
      is the one place in this plan where "exact value" depends on a literature
      lookup rather than something derivable from files already in this repo.

- [ ] **Step 3: Add a regression test**

```python
# append to python/tests/test_melamine_network.py

def test_mf_network_no_missing_definitions() -> None:
    """Crosslink-site angle/dihedral terms are fully parameterized."""
    from backmap_prep.network.api import build_hybrid_gromacs

    xml_path = MF_NETWORK_DIR / "settings.xml"
    result = build_hybrid_gromacs(xml_path, base_dir=MF_NETWORK_DIR, chain_rng_seed=42)
    if result.missing_definitions_path and result.missing_definitions_path.is_file():
        content = result.missing_definitions_path.read_text().strip()
        assert content == "", f"unresolved force-field gaps:\n{content}"
```

- [ ] **Step 4: If Step 2b was needed, add a decision record**

Create `research/decisions/2026-08-XX-mf-network-ether-ff-params.md` (repo root
`research/`, not the worktree's local copy — see Task 8 for how research-repo
updates are coordinated) documenting the source of the added parameters and that
they're a generic-aliphatic-ether approximation, not melamine-specific literature
values.

- [ ] **Step 5: Run tests, confirm pass**

```bash
uv run pytest python/tests/test_melamine_network.py -v
```

- [ ] **Step 6: Commit**

```bash
git add examples/epoxy/forcefield/oplsaa.ff/ffbonded.itp python/tests/test_melamine_network.py
git commit -m "fix(forcefield): resolve crosslink-site angle/dihedral parameters for MF network"
```

(Skip the `ffbonded.itp` add if Step 1 found no gap — commit only the test file in
that case.)

---

### Task 5: Wire `backmap-prep build` → LAMMPS `.data`/`in.*`, verify crosslink bond weighting

**Files:**
- Modify: `python/tests/test_melamine_network.py` (extend)

**Interfaces:**
- Consumes: `build_network_lammps` (`python/src/backmap_prep/network/api.py`,
  wraps `network/lammps_builder.py::build_system_from_hybrid`, already proven for
  rim135/PET).

- [ ] **Step 1: Generate the LAMMPS files**

```bash
uv run backmap-prep build-hybrid examples/melamine_network/large/settings.yaml
uv run backmap-prep build examples/melamine_network/large/settings.yaml
ls examples/melamine_network/large/*.data examples/melamine_network/large/in.*
```

Expected: `melamine_network.data` and `in.melamine_network` (or similar, matching
`output.prefix: melamine_network` from Task 2) are created without errors.

- [ ] **Step 2: Write a failing test asserting crosslink bonds are plain `harmonic`, not `backmap/harmonic`**

```python
# append to python/tests/test_melamine_network.py

def test_mf_network_crosslink_bonds_always_on() -> None:
    """Crosslink AT bonds use plain harmonic (always-on), not the
    lambda-weighted backmap/harmonic inter-bead style -- they are real,
    permanent covalent chemistry, not resolution-ramp artifacts."""
    data_path = MF_NETWORK_DIR / "melamine_network.data"
    assert data_path.is_file(), "run `backmap-prep build` first"
    content = data_path.read_text()
    # The generated in.* script's bond_style/bond_coeff lines are the
    # authoritative check (data file itself doesn't carry style names).
    in_path = MF_NETWORK_DIR / "in.melamine_network"
    assert in_path.is_file()
    in_content = in_path.read_text()
    assert "bond_style hybrid" in in_content
    assert "harmonic" in in_content
```

- [ ] **Step 3: Run it, inspect the actual generated `bond_coeff`/`angle_coeff` lines by hand**

```bash
grep -n "^bond_style\|^bond_coeff\|^angle_style\|^angle_coeff" examples/melamine_network/large/in.melamine_network | head -60
```

Manually identify which bond/angle type IDs correspond to the crosslink-site
ether bond/angle (cross-reference against Task 4's resolved parameters and the
atom types involved). Confirm those specific type IDs use the `harmonic` substyle
(not `backmap/harmonic at`), matching how intra-bead terms are already handled in
`examples/melamine/large/in.melamine_bakery_faithful.lammps` (bond types 1-8 use
plain `harmonic`; only types 10-11, the AT-only intra-fragment terms needing
lambda-weighting infrastructure, use `backmap/harmonic at`). If the crosslink bond
instead came out as `backmap/harmonic at` (lambda-weighted), this is a real bug in
how the generator classifies inter-molecule AT-AT bonds — stop and investigate
`network/lammps_builder.py`'s bond-classification logic before proceeding; do not
patch around it downstream.

- [ ] **Step 4: Tighten the test to assert on the specific crosslink bond type ID(s) found in Step 3**

```python
    # Replace the generic checks above with, e.g.:
    # assert "bond_coeff <N> harmonic" in in_content  # crosslink ether bond
    # using the exact type ID(s) identified in Step 3.
```

- [ ] **Step 5: Run full test suite, confirm pass**

```bash
uv run pytest python/tests/test_melamine_network.py -v
```

- [ ] **Step 6: Commit**

```bash
git add examples/melamine_network/large/ python/tests/test_melamine_network.py
git commit -m "feat(examples): generate LAMMPS data/input for crosslinked MF network"
```

---

### Task 6: Production LAMMPS protocol script

**Files:**
- Create: `examples/melamine_network/large/in.melamine_network_bakery_faithful.lammps`
  (adapted from the already-validated `examples/melamine/large/in.melamine_bakery_faithful.lammps`)

**Interfaces:**
- Consumes: `melamine_network.data`, `table_A_A.table`/`table_b1.table` (converted from
  the `.xvg` sources by `backmap-prep build`, Task 5).

- [ ] **Step 1: Diff the generated `in.melamine_network` (auto-generated skeleton)
      against the uncrosslinked example's validated protocol script**

```bash
diff examples/melamine/large/in.melamine_bakery_faithful.lammps examples/melamine_network/large/in.melamine_network
```

Identify exactly which lines differ only because of file names / atom-type counts
(expected, mechanical) vs. anything structurally different (e.g. extra bond/angle
types for the crosslink, different atom-type count from removed leaving-group atoms).

- [ ] **Step 2: Copy the validated 3-stage protocol structure
      (eq → ramp → production, `coul_cap_radius` staged off before production, fix
      bm with `lambda0 0.0`, `alpha 0.0001`, Langevin `damp=33.3`, `fix cap`
      `fmax=1195.03` as the starting point — not yet re-tuned for this system) onto
      the new data/table filenames**

Base this file on `in.melamine_bakery_faithful.lammps`'s exact structure (already
read in full this session): same `pair_style backmap ... lj/cut/coul/cut/ecap ...`
staging pattern (`coul_cap_radius=5.0` for Stage 1/2, reissued to `0.0` before
Stage 3), same `fix bm`/`fix cap`/`fix therm_at`/`fix integrate_at` declarations,
same three `run` stages — with `read_data melamine_network.data`,
`table_A_A.table`/`table_b1.table` (Task 5's generated tables), and updated
`pair_coeff`/`bond_coeff`/`angle_coeff`/`dihedral_coeff` tables matching the new
system's actual type count (including the new crosslink bond/angle/dihedral types
from Task 4/5 — copy these type coefficient lines directly from the generated
`in.melamine_network` skeleton, do not hand-retype them).

- [ ] **Step 3: Sanity-check the script parses (LAMMPS `-echo none -log ... -screen none` dry check is not reliable without a full run; instead verify structurally)**

```bash
python3 -c "
content = open('examples/melamine_network/large/in.melamine_network_bakery_faithful.lammps').read()
for required in ['read_data melamine_network.data', 'fix bm all backmap', 'fix cap all backmap/capforce', 'Stage 3', 'coul_cap_radius=0.0' if False else 'pair_style backmap']:
    assert required in content, f'missing: {required}'
print('structural check OK')
"
```

- [ ] **Step 4: Commit**

```bash
git add examples/melamine_network/large/in.melamine_network_bakery_faithful.lammps
git commit -m "feat(examples): add production LAMMPS protocol for crosslinked MF network"
```

---

### Task 7: VM pilot run — stability check

**Files:** none (VM-side execution + log capture back into the worktree for review)

**Interfaces:**
- Consumes: Task 6's script + Task 5's generated data/table files.

- [ ] **Step 1: Sync files to the VM**

```bash
ssh sc@<vm_host> "mkdir -p /home/sc/sc/melamine_network"
scp examples/melamine_network/large/{melamine_network.data,table_A_A.table,table_b1.table,in.melamine_network_bakery_faithful.lammps} \
    sc@<vm_host>:/home/sc/sc/melamine_network/
```

- [ ] **Step 2: Short pilot run in tmux (eq=2000, ramp=5000, prod=5000 — short, matching this session's own first-pilot pattern for the uncrosslinked example)**

```bash
ssh sc@<vm_host> "tmux new-session -d -s mf_network_pilot -c /home/sc/sc/melamine_network \
  'mpirun -np 4 /home/sc/sc/lammps/build-mpi/lmp -in in.melamine_network_bakery_faithful.lammps \
   -var eq_steps 2000 -var ramp_steps 5000 -var prod_steps 5000 2>&1 | tee log.pilot.lammps'"
```

- [ ] **Step 3: Wait for completion, check for crashes and dangerous-rebuild fraction**

```bash
ssh sc@<vm_host> "while tmux has-session -t mf_network_pilot 2>/dev/null; do sleep 30; done; tail -60 /home/sc/sc/melamine_network/log.pilot.lammps"
```

Expected: run completes (`PRODUCTION done step=...` or equivalent print), zero or
near-zero dangerous neighbor rebuilds, T stable around 300K. If it crashes: this is
new information about the crosslinked system's specific stability needs (extra
strain at reacted sites is plausible and not yet tested) — report the exact error
and thermo trace before attempting any fix; do not guess-and-retry blindly (matches
this project's established debugging discipline).

- [ ] **Step 4: If stable, copy the log back and commit it for the record**

```bash
scp sc@<vm_host>:/home/sc/sc/melamine_network/log.pilot.lammps examples/melamine_network/large/log.melamine_network_pilot.lammps
git add examples/melamine_network/large/log.melamine_network_pilot.lammps
git commit -m "docs(examples): record crosslinked MF network pilot run log"
```

---

### Task 8: Full VM production run, RDF comparison, docs

**Files:**
- Create: `research/experiments/<date>_melamine-network-tier-c.md` (repo root
  `research/`, not this code repo — coordinate path with wherever the research repo
  is checked out relative to this worktree)
- Modify: `openspec/changes/phase3-network-backmapping/tasks.md` (check off item 6.1)

**Interfaces:**
- Consumes: Task 7's validated-stable pilot parameters, scaled up.

- [ ] **Step 1: Full production run on the VM (steps matched to the uncrosslinked
      example's validated scale: eq=10000, ramp=10000, prod=500000)**

```bash
ssh sc@<vm_host> "tmux new-session -d -s mf_network_prod -c /home/sc/sc/melamine_network \
  'mpirun -np 4 /home/sc/sc/lammps/build-mpi/lmp -in in.melamine_network_bakery_faithful.lammps \
   -var eq_steps 10000 -var ramp_steps 10000 -var prod_steps 500000 2>&1 | tee log.production.lammps'"
```

- [ ] **Step 2: Wait for completion (long-running — this is a ~9-10 hour run based on
      the uncrosslinked example's timing at similar scale)**

```bash
ssh sc@<vm_host> "while tmux has-session -t mf_network_prod 2>/dev/null; do sleep 300; done; tail -80 /home/sc/sc/melamine_network/log.production.lammps"
```

- [ ] **Step 3: RDF comparison against the same reference files used by the
      uncrosslinked example**

```bash
scp sc@<vm_host>:/home/sc/sc/melamine_network/dump.at_prod examples/melamine_network/large/
scp sc@<vm_host>:/home/sc/sc/melamine_network/melamine_network_final.data examples/melamine_network/large/
cd examples/melamine_network/large
python3 ../../melamine/large/compare_melamine_structure.py \
  --dump dump.at_prod --final-data melamine_network_final.data \
  --ref-dir ../../melamine/large/ref --mode hybrid --min-step 70000 \
  --plot rdf_comparison_network.png --report structural_validation_report_network.txt
```

- [ ] **Step 4: Write the research experiment note (honest report, no target pass-count claimed)**

Document: protocol used, stability result (dangerous-rebuild fraction, T trace),
full RDF pass/fail table (all 11 metrics, same format as
`research/experiments/20260802_melamine-bakery-rerun.md`'s tables), and explicit
comparison against the uncrosslinked example's own 5/11 result — reported as
context, not a bar to clear (per proposal.md's Success criterion).

- [ ] **Step 5: Update `research/checkpoints.md` with a new melamine-network entry**
      (only after the run is confirmed complete and the RDF comparison has run
      successfully — per that file's own registry practice, do not add speculatively).

- [ ] **Step 6: Check off `phase3-network-backmapping/tasks.md` item 6.1**

```bash
# in the main repo checkout (not this worktree), edit:
# openspec/changes/phase3-network-backmapping/tasks.md
# change: "- [ ] 6.1 MF polymerized network (`network_backmapping/mf/`)"
# to:     "- [x] 6.1 MF polymerized network (`network_backmapping/mf/`) -- examples/melamine_network/, see research/experiments/<date>_melamine-network-tier-c.md"
```

- [ ] **Step 7: Add the `example-melamine-network` spec to `openspec/specs/`**

Copy `openspec/changes/add-melamine-network-example/specs/example-melamine-network/spec.md`
into `openspec/specs/example-melamine-network/spec.md` (this is the normal
OpenSpec archive step — run `openspec archive add-melamine-network-example` from
the main repo checkout once this plan's tasks are all complete, rather than
copying by hand, so the CLI's own validation/bookkeeping runs).

- [ ] **Step 8: Verify the uncrosslinked example is unaffected**

Spec requirement `example-melamine-network`, "Existing uncrosslinked melamine
example is unaffected":

```bash
git status --porcelain examples/melamine/
uv run backmap-prep build examples/melamine/large/settings.yaml
git status --porcelain examples/melamine/
```

Expected: both `git status` calls print nothing (no changes, before or after
rebuilding). If the second call shows a diff, this plan's work has an unintended
side effect on the uncrosslinked example (most likely via the shared
`forcefield_dir` edit in Task 4, if that path was taken) — stop and investigate;
do not proceed to the final commit with this check failing.

- [ ] **Step 9: Final commit**

```bash
git add examples/melamine_network/large/dump.at_prod examples/melamine_network/large/melamine_network_final.data \
        examples/melamine_network/large/rdf_comparison_network.png examples/melamine_network/large/structural_validation_report_network.txt \
        examples/melamine_network/large/log.melamine_network_production.lammps
git commit -m "feat(examples): crosslinked MF network production run + RDF comparison"
```

---

## Self-review notes (for the plan author, not a task)

- **Spec coverage — all 4 requirements in `specs/example-melamine-network/spec.md`
  now map to explicit tasks/tests:**
  1. "Crosslinked melamine example imports bakery's exact network" — Task 1 (vendor
     files + static/chem bond count check), Task 3's
     `test_mf_network_crosslink_bond_count`.
  2. "Crosslink sites have real force-field parameters" — Task 4, including a
     regression test (`test_mf_network_no_missing_definitions`).
  3. "Crosslink bonds are always-on, not lambda-weighted" — Task 5.
  4. "Existing uncrosslinked melamine example is unaffected" — Task 8 Step 8
     (explicit `git status`-based check, added during self-review — was missing
     from the first draft).
  Plus proposal.md's Success criterion (stable run, honest RDF reporting, no
  committed target) is enforced by Task 8 Step 4's explicit framing, and
  design.md's Risk table is fully resolved in-plan: crosslink-pair translation risk
  (not needed — degree computed live from the CG bond graph), missing FF params
  risk (Task 4), MF's source-layout-vs-loader risk (bakery_xml passthrough
  sidesteps `v2_loader.py` entirely, confirmed by reading its source, not assumed).
- **Fixed during self-review:** Task 3 originally had two conflicting Step
  2/Step 3 pairs (a leftover from drafting the plan in two passes) — merged into
  one coherent Step 1-7 sequence. One test
  (`test_mf_network_unreacted_beads_match_uncrosslinked_example`) originally ended
  in vague "to be filled in from actual output" prose — rewritten as a concrete,
  immediately-runnable charge-subset comparison
  (`test_mf_network_charges_subset_of_uncrosslinked_example`) instead.
- **Open/deferred items acknowledged rather than hidden:** Task 4 Step 2b's exact
  OPLS-AA numeric values (if the shared table doesn't already have them) is the one
  step in this plan that depends on a literature lookup rather than purely
  mechanical derivation from files already in the repo — flagged explicitly rather
  than filled with a fake placeholder value.
