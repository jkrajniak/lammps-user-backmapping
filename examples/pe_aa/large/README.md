# Large-scale PE-AA example

Large-scale inputs are sourced from the [bakery](https://github.com/bakery-cg2at/bakery) project.

## How to obtain large-scale inputs

1. Clone or have the bakery repository available.
2. Copy the required files from `bakery/examples/pe_aa/` (e.g. from `prepare/` or `backmapping/` subdirs) into this directory: CG/AT topologies, coordinates, `pe_single.gro`, and table files.
3. Run backmap-prep from repository root: `uv run backmap-prep examples/pe_aa/large/settings.yaml`.
4. Run LAMMPS with the generated input script when ready.
