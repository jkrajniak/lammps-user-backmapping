# Large-scale examples

The repository includes **small-scale** variants of each example (reduced system sizes) for quick testing. For production-scale validation and regression testing, **large-scale** variants are available. Their inputs are sourced from the [bakery](https://github.com/bakery-cg2at/bakery) project.

## Layout

Each example directory has a `large/` subdirectory:

| Example   | Small-scale (root)     | Large-scale (`large/`)        |
|----------|------------------------|--------------------------------|
| Dodecane | `examples/dodecane/`   | `examples/dodecane/large/`     |
| PE       | `examples/pe/`         | `examples/pe/large/`           |
| PE4      | `examples/pe4/`        | `examples/pe4/large/`          |
| PE-10    | `examples/pe_10/`      | `examples/pe_10/large/`        |
| PE-AA    | `examples/pe_aa/`      | `examples/pe_aa/large/`        |
| Melamine | `examples/melamine/`   | `examples/melamine/large/`     |

The `large/` directory holds (or documents how to obtain) bakery-scale inputs and generated LAMMPS files.

## Source of large-scale inputs

Large-scale coordinates, topologies, and tabulated potentials come from **bakery**:

- **Bakery repository**: [bakery-cg2at/bakery](https://github.com/bakery-cg2at/bakery)
- **Paths**: `bakery/examples/<name>/` (and `bakery/tests/<name>/` where needed)
- **Typical sizes**: e.g. 75 chains for PE systems, ~500 molecules for melamine

Each `examples/<name>/large/README.md` describes which files to copy and how to refresh them from bakery.

## How to run large-scale examples

1. **Obtain inputs**
   Copy the required files from bakery into `examples/<name>/large/` as described in that directory’s README.

2. **Generate LAMMPS input and data**
   From the repository root:
   ```bash
   cd python && uv run backmap-prep ../examples/<name>/large/settings.yaml
   ```
   Output (`.data` and `in.<name>`) is written according to the paths in `settings.yaml` (often into the same `large/` directory).

3. **Run LAMMPS** (optional)
   From the directory that contains the generated files:
   ```bash
   lmp -in in.<name>
   ```

Large-scale runs can be slow; use short test runs or reduced output frequency if you only need to verify the pipeline.

## Optional validation script

A script is provided to verify that `backmap-prep` runs successfully on at least one large-scale example (file generation only; it does not run LAMMPS). From the repository root:

```bash
./scripts/validate-large-scale-prep.sh dodecane
```

You can pass any example name (`dodecane`, `pe`, `pe4`, `pe_10`, `pe_aa`, `melamine`) to validate that example's large-scale settings. This is suitable for CI or manual checks.
