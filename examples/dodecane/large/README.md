# Large-scale dodecane example

Large-scale inputs for this example are sourced from the [bakery](https://github.com/bakery-cg2at/bakery) project.

## How to obtain large-scale inputs

1. Clone or have the bakery repository available.
2. Copy the required files from `bakery/examples/dodecane/` into this directory:
   - `cg_conf.gro` — CG coordinates (many chains)
   - `topol_cg.top`, `topol_aa.top`
   - `dodecane_single.gro`
   - `table_*.xvg` (table_b1.xvg, table_a1.xvg, table_a2.xvg, table_A_A.xvg, table_A_B.xvg, table_B_B.xvg).
   Angle tables `table_a1.xvg` and `table_a2.xvg` are large; if missing here, copy them from the parent directory `examples/dodecane/`.

   Example (from this repo root):
   ```bash
   cp /path/to/bakery/examples/dodecane/{cg_conf.gro,topol_cg.top,topol_aa.top,dodecane_single.gro,table_*.xvg} examples/dodecane/large/
   ```
   Or copy angle tables from the small-scale example: `cp ../table_a1.xvg ../table_a2.xvg .`

3. Run backmap-prep to generate LAMMPS input and data files:
   ```bash
   cd python && uv run backmap-prep ../examples/dodecane/large/settings.yaml
   ```
   Output will be written to `examples/dodecane/large/` (or the current directory; check `output.prefix` in settings).

4. Run LAMMPS (optional): `lmp -in in.dodecane` (from the directory containing the generated files).
