# Large-scale melamine example (~500 molecules)

Large-scale inputs are sourced from the [bakery](https://github.com/bakery-cg2at/bakery) project.

## How to obtain large-scale inputs

1. Clone or have the bakery repository available.
2. Copy the required files from `bakery/examples/melamine/` (or bakery tests) into this directory: `cg_conf.gro` (500 molecules), `topol_cg.top`, `topol_aa.top`, `single_mf.gro`, and `table_*.xvg` / `table_*.table` as needed.
3. Run backmap-prep from repository root: `uv run backmap-prep examples/melamine/large/settings.yaml`.
4. Run LAMMPS with the generated input script when ready.
