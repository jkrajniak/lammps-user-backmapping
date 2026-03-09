# LAMMPS Backmapping Package

Time-dependent backmapping (reverse mapping) from coarse-grained to atomistic
resolution, implemented as a LAMMPS user package.

The method ramps a per-atom resolution parameter **&lambda;** from 0 (pure CG)
to 1 (pure AT) uniformly across the simulation box, gradually restoring
atomistic detail while maintaining thermodynamic consistency.

!!! cite "References"
    Krajniak et al., "Generic Adaptive Resolution Method for Reverse Mapping of
    Polymers from Coarse-Grained to Atomistic Descriptions",
    *J. Chem. Theory Comput.* 2016, 12, 5549--5562.
    [DOI: 10.1021/acs.jctc.6b00595](https://doi.org/10.1021/acs.jctc.6b00595)

    Krajniak, Zhang et al., "Reverse Mapping Method for Complex Polymer Systems",
    *J. Comput. Chem.* 2018.
    [DOI: 10.1002/jcc.25129](https://doi.org/10.1002/jcc.25129)

## Features

- **Smooth resolution transition** -- lambda ramp with configurable rate drives
  CG &rarr; AT conversion over the course of the simulation
- **Lambda-weighted interactions** -- pair, bond, and angle styles that
  automatically weight forces by the current resolution
- **CG force distribution** -- CG-level forces are redistributed to AT atoms
  proportional to mass
- **Automated input generation** -- the `backmap-prep` CLI reads a YAML
  settings file and produces LAMMPS data files, input scripts, and interaction
  tables from GROMACS topologies
- **Restart support** -- per-atom lambda values are saved and restored across
  LAMMPS restarts

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
