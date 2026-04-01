# Large-scale PE-10 example (75 chains)

Large-scale inputs are sourced from the [bakery](https://github.com/bakery-cg2at/bakery) project.

## How to obtain large-scale inputs

1. Clone or have the bakery repository available.
2. Copy the required files from `bakery/examples/pe_10/` (or `bakery/tests/pe_10/`) into this directory: `cg_conf.gro`, `topol_cg.top`, `topol_aa.top`, `pe_single.gro`, and all `table_*.xvg` files.
3. Run backmap-prep from repository root: `uv run backmap-prep examples/pe_10/large/settings.yaml`.
4. Run LAMMPS with the generated input script when ready.
