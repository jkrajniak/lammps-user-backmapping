## 1. Layout and documentation skeleton

- [x] 1.1 Add `large/` subdirectory under each example (dodecane, pe, pe4, pe_10, pe_aa, melamine) with a README or placeholder stating that large-scale inputs come from bakery and how to obtain them
- [x] 1.2 Create docs page or section "Large-scale examples" describing layout, bakery as source, and how to run backmap-prep and LAMMPS for large variants

## 2. Dodecane large-scale (pilot)

- [x] 2.1 Copy or document bakery dodecane large-scale inputs (e.g. from `bakery/examples/dodecane/`) into `examples/dodecane/large/`
- [x] 2.2 Add `examples/dodecane/large/settings.yaml` (or adapt paths from small-scale) so backmap-prep uses large-scale CG/AT files
- [x] 2.3 Run backmap-prep for dodecane large-scale and verify generation of `.data` and `in.dodecane`
- [x] 2.4 Document in `examples/dodecane/README.md` (or `large/README.md`) how to obtain and run the large-scale dodecane example

## 3. PE systems large-scale (pe, pe4, pe_10, pe_aa)

- [x] 3.1 Add large-scale inputs and settings for `examples/pe/large/` from bakery (75 chains); add or adapt `settings.yaml` and document in README
- [x] 3.2 Add large-scale inputs and settings for `examples/pe4/large/` from bakery (75 chains); add or adapt `settings.yaml` and document in README
- [x] 3.3 Add large-scale inputs and settings for `examples/pe_10/large/` from bakery (75 chains); add or adapt `settings.yaml` and document in README
- [x] 3.4 Add large-scale inputs and settings for `examples/pe_aa/large/` from bakery; add or adapt `settings.yaml` and document in README

## 4. Melamine large-scale

- [x] 4.1 Add large-scale inputs and settings for `examples/melamine/large/` from bakery (500 molecules); add or adapt `settings.yaml` and document in README

## 5. Main README and docs

- [x] 5.1 Update top-level README to mention large-scale examples and link to the Large-scale examples doc and per-example instructions
- [x] 5.2 Update CHANGELOG.md under [Unreleased] with Added entry for large-scale example variants

## 6. Optional CI or validation script

- [x] 6.1 Add script or CI job that runs backmap-prep and LAMMPS (short run) on at least one large-scale example (e.g. dodecane or pe), marked optional or on-demand so default CI stays fast
