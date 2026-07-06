## ADDED Requirements

### Requirement: Degree-dependent bead mapping

The `backmap-prep` generator SHALL support bead definitions where the atom list
depends on the CG bead's topological degree in the cured network, matching bakery
`atoms_by_degree` semantics.

#### Scenario: EPO degree 0 bead

- **WHEN** an EPO CG bead has degree 0 in the cured topology
- **THEN** the generator SHALL map it to the degree-0 AT fragment from settings

#### Scenario: EPO degree 2 bead

- **WHEN** an EPO CG bead has degree 2
- **THEN** the generator SHALL map it to the degree-2 AT fragment (cross-linked)

---

### Requirement: Multi-molecule hybrid systems

The generator SHALL produce a single hybrid GRO/TOP/data file containing multiple
molecule types (EPO, IPD, HDD) with correct molecule IDs and cross-type exclusions.

#### Scenario: RIM135 three-component system

- **WHEN** epoxy settings define EPO, IPD, and HDD molecules
- **THEN** generated LAMMPS data SHALL have correct atom counts per species

---

### Requirement: Charge conservation

Charge transfer rules SHALL conserve total system charge within 1e-6 e.

#### Scenario: Epoxy charge sum

- **WHEN** charge transfers are applied for cured epoxy
- **THEN** sum of atomic charges before and after SHALL differ by less than 1e-6 e

---

### Requirement: Epoxy Tier A parity

Generated hybrid files SHALL match bakery `tests/rim135/ref_hyb_topol.top` and
`ref_hyb_conf.gro` within migration-parity tolerances.

#### Scenario: Rim135 generator test

- **WHEN** `backmap-prep examples/epoxy/settings.yaml` is run
- **THEN** topology section counts SHALL match bakery `compare_top.py` expectations
