# LAMMPS Backmapping Package

Adaptive-resolution backmapping (reverse mapping) from coarse-grained to
atomistic resolution, implemented as a LAMMPS user package.

The method ramps a per-atom resolution parameter **&lambda;** from 0 (pure CG)
to 1 (pure AT) uniformly across the simulation box, gradually restoring
atomistic detail while maintaining thermodynamic consistency.

!!! cite "References"
    Krajniak et al., "Generic Adaptive Resolution Method for Reverse Mapping of
    Polymers from Coarse-Grained to Atomistic Descriptions",
    *J. Chem. Theory Comput.* 2016, 12, 5549--5562.
    [DOI: 10.1021/acs.jctc.6b00595](https://doi.org/10.1021/acs.jctc.6b00595)

    Krajniak, Zhang, Pandiyan, Nies, Samaey, "Reverse Mapping Method for Complex
    Polymer Systems", *J. Comput. Chem.* **2018**, *39*, 648--664 (online
    **Dec. 2017**) — networks (epoxy, melamine), hyperbranched polymer, PET.
    [DOI: 10.1002/jcc.25129](https://doi.org/10.1002/jcc.25129)

??? cite "Applications (optional)"
    Zhang, Wheatle, Krajniak, Keith, Ganesan, "Ion Mobilities, Transference
    Numbers, and Inverse Haven Ratios of Polymeric Ionic Liquids",
    *ACS Macro Lett.* **2020**, *9*, 84--89 — charged polymers (polymeric
    ionic liquids of different chain lengths).
    [DOI: 10.1021/acsmacrolett.9b00908](https://doi.org/10.1021/acsmacrolett.9b00908)

    Zhang, Krajniak, Keith, Ganesan, "Mechanisms of Ion Transport in Block
    Copolymeric Polymerized Ionic Liquids", *ACS Macro Lett.* **2019**, *8*,
    1096--1101 — microphase-separated block copolymers.
    [DOI: 10.1021/acsmacrolett.9b00478](https://doi.org/10.1021/acsmacrolett.9b00478)

    Zhang, Krajniak, Ganesan, "A Multiscale Simulation Study of Influence of
    Morphology on Ion Transport in Block Copolymeric Ionic Liquids",
    *Macromolecules* **2021**, *54*, 4997--5010 — microphase-separated block
    copolymers.
    [DOI: 10.1021/acs.macromol.1c00025](https://doi.org/10.1021/acs.macromol.1c00025)

    Zhang, Sass, Krajniak, Ganesan, "Ion Correlations and Partial Ionicities
    in the Lamellar Phases of Block Copolymeric Ionic Liquids",
    *ACS Macro Lett.* **2022**, *11*, 1265--1271 — microphase-separated block
    copolymers (lamellar phase).
    [DOI: 10.1021/acsmacrolett.2c00401](https://doi.org/10.1021/acsmacrolett.2c00401)

    Zhang, Krajniak, Sass, Sachar, Marioni, Duncan, Ganesan, "Conductivity and
    Transference Numbers in Lithium Salt-Doped Block Copolymeric Ionic Liquid
    Electrolytes", *Macromolecules* **2023**, *56*, 9750--9765 —
    microphase-separated block copolymers with salt.
    [DOI: 10.1021/acs.macromol.3c01791](https://doi.org/10.1021/acs.macromol.3c01791)

    Zhang, Krajniak, Samaey, Nies, "A Parallel Multiscale Simulation Framework
    for Complex Polymerization: AB2-Type Monomer Hyperbranched Polymerization
    as an Example", *Adv. Theory Simul.* **2019**, *2*, 1800102 — topologies
    formed by polymerization at the CG level (hyperbranched polymers).
    [DOI: 10.1002/adts.201800102](https://doi.org/10.1002/adts.201800102)

## Features

- **Smooth resolution transition** -- lambda ramp with configurable rate drives
  CG &rarr; AT conversion over the course of the simulation
- **Lambda-weighted interactions** -- pair, bond, and angle styles that
  automatically weight forces by the current resolution
- **CG force distribution** -- CG-level forces are redistributed to AT atoms
  proportional to mass
- **Automated input generation** -- the `backmap-prep` CLI reads a YAML
  settings file and produces LAMMPS data files, input scripts, and interaction
  tables from GROMACS topologies (the CG system and/or each molecule's AT
  fragment can alternatively be supplied as native LAMMPS input -- see
  [Settings Reference](settings-reference.md#cg_system))

## Components

| Component | Description |
|-----------|-------------|
| [`fix backmap`](components/fix-backmap.md) | Lambda ramp, CG-AT mapping, COM tracking, CG force distribution |
| [`pair_style backmap`](components/pair-backmap.md) | Lambda-weighted non-bonded pair forces |
| [`bond_style backmap/harmonic`](components/bond-styles.md) | Lambda-weighted harmonic bond forces |
| [`bond_style backmap/table`](components/bond-styles.md#backmap-table) | Lambda-weighted tabulated bond forces |
| [`angle_style backmap/harmonic`](components/angle-styles.md) | Lambda-weighted harmonic angle forces |
| [`backmap-prep`](cli/backmap-prep.md) | Python CLI for generating LAMMPS input files |

## Quick Links

- [Getting Started](getting-started.md) -- install and run your first backmapping
  simulation
- [Large-scale examples](large-scale-examples.md) -- production-scale variants (75 chains, 500 molecules) and how to run them
- [Settings Reference](settings-reference.md) -- complete YAML settings documentation
- [Tutorial: Setting Up a New System](tutorial-new-system.md) -- step-by-step guide for
  your own molecule
- [Theory](theory.md) -- how the backmapping method works

## License

GPL-3.0-or-later
