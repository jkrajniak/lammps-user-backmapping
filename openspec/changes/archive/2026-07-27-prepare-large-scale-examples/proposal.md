## Why

The current examples (dodecane, pe, pe4, pe_10, pe_aa, melamine) use reduced system sizes (e.g. 10 chains, 50 molecules) for quick testing. Bakery provides the same systems at production scale (e.g. 75 chains, 500 molecules). To validate the LAMMPS backmapping package under realistic load, support regression testing, and give users reference setups for larger runs, we need larger-scale example variants that mirror bakery’s sizes and are fully migrated to the `backmap-prep` + LAMMPS workflow.

## What Changes

- **Add large-scale example variants** for dodecane, pe, pe4, pe_10, pe_aa, and melamine, using system sizes aligned with bakery (e.g. 75 chains for PE systems, 500 molecules for melamine).
- **Source inputs** (GRO, TOP, XVG tables) from bakery where applicable; document any scaling or regeneration steps.
- **Per-example layout**: each example keeps its current small-scale version; large-scale data lives in a dedicated subdir (e.g. `examples/pe/large/` or `examples/pe_75chains/`) or is documented so users can generate it from bakery inputs.
- **Documentation**: README updates and/or a “Large-scale examples” doc describing how to obtain/run the large variants and how they map to bakery.
- **CI/validation**: Optional integration test or script that verifies `backmap-prep` and LAMMPS run successfully on at least one large-scale example (e.g. dodecane or pe).

## Capabilities

### New Capabilities

- `large-scale-examples`: Specification for larger-scale variants of existing examples (dodecane, pe, pe4, pe_10, pe_aa, melamine), including directory layout, source from bakery, settings.yaml and LAMMPS inputs, and documentation/run instructions.

### Modified Capabilities

- None. Existing example specs (example-pe-systems, example-melamine) describe the small-scale examples; the new spec defines the large-scale extension only.

## Impact

- **`examples/`**: New or extended directories for large-scale variants; possible duplication of topology/table files or symlinks/scripts to pull from bakery.
- **`docs/`**: New or updated page for “Large-scale examples” and updates to getting-started if needed.
- **`README.md` / `CHANGELOG.md`**: Updated to mention large-scale examples and how to run them.
- **CI**: Optional job or script to run one large-scale example; runtime and resource limits must be considered.
