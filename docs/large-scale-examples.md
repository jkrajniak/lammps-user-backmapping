# Large-scale examples

The repository includes **small-scale** variants of each example (reduced system sizes) for quick testing. For production-scale validation and regression testing, **large-scale** variants are available. Their inputs are sourced from the [bakery](https://github.com/bakery-cg2at/bakery) project.

## Layout

Each example directory has a `large/` subdirectory:

| Example   | Small-scale (root)     | Large-scale (`large/`)        |
|----------|------------------------|--------------------------------|
| Dodecane | `examples/dodecane/`   | `examples/dodecane/large/`     |
| Dodecane (250 mol) | — | `examples/dodecane/n250/` (medium; CG from bakery `large/cg_conf.gro`) |
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
   uv run backmap-prep examples/<name>/large/settings.yaml
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

## Paper-grade RDF validation (dodecane large)

`examples/dodecane/large/` ships a publication-quality validation workflow that
compares the backmapped melt against an independent all-atom reference
simulation at the same density and temperature.

**Inputs (already in the directory):**

- `in.dodecane_at_long` — reads the extracted post-backmap AT system
  (`dodecane_at.data`), runs a multi-phase equilibration (minimise → NPT 100→298 K
  → NPT 298 K → NVT 298 K) followed by 1 ns NVT production. RDF averaged in
  five 200 ps blocks with `fix ave/time 100 2000 200000`.
- `in.dodecane_at_ref_long` — independent all-atom reference for 500 dodecane
  molecules (NVT high-T → `fix deform` compression → NPT anneal → NVT 298 K),
  followed by 1 ns NVT production with the same RDF averaging settings.
- `compare_rdf_blocks.py` — parses the multi-block `fix ave/time` output of
  both runs, computes per-pair mean and standard error across blocks, plots
  mean ± SEM bands, and reports PASS/FAIL on first-peak position, height, and
  L2(g_bm − g_ref) tolerances.

**Run:**

```bash
lmp -in in.dodecane_at_ref_long -log log.at_ref_long.lammps   # ~3.5 h serial
lmp -in in.dodecane_at_long    -log log.at_long.lammps        # ~4.5 h serial

uv run compare_rdf_blocks.py \
    --backmap rdf_backmap_long.dat \
    --reference rdf_reference_long.dat \
    --plot rdf_comparison_long.png \
    --report rdf_comparison_long.txt \
    --skip-first 1
```

The committed `rdf_comparison_long.{png,txt}` are the reference outputs and
should be regenerated whenever the C++ styles or the dodecane example change.
