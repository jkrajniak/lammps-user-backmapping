## Context

The LAMMPS backmapping package ships with small-scale examples (dodecane, pe, pe4, pe_10, pe_aa, melamine) migrated from bakery, using reduced sizes (e.g. 10 chains, 50 molecules) for fast testing. Bakery’s reference systems use production-scale sizes (e.g. 75 PE chains, 500 melamine molecules). Users and CI need larger-scale variants to validate performance, regression, and real-world workflows without changing the existing small examples.

## Goals / Non-Goals

**Goals:**

- Provide larger-scale example variants for all six systems (dodecane, pe, pe4, pe_10, pe_aa, melamine) with sizes aligned to bakery where applicable.
- Use bakery as the source for GRO/TOP/XVG and single-molecule references; document how to obtain or regenerate large-scale data.
- Keep the same `backmap-prep` + LAMMPS workflow; large-scale examples are additional data/config, not new tooling.
- Document layout (e.g. `large/` subdirs or separate dirs) and how to run large-scale examples; update docs and README.

**Non-Goals:**

- Changing `backmap-prep` schema or CLI for large-scale only; any schema changes are out of scope unless required by existing examples.
- Automatically generating large systems from small ones (no “scale-up” generator); we use or document bakery-derived inputs.
- Guaranteeing CI runs of all large examples on every commit (optional; at least one large example run is desirable).

## Decisions

1. **Layout for large-scale data**
   Use a `large/` subdirectory under each existing example (e.g. `examples/pe/large/`) to hold large-scale coordinates, topology, and generated LAMMPS files. Small-scale remains at `examples/<name>/` root.
   *Rationale*: Keeps one example per system; avoids duplicating settings/docs at the top level; clear separation of “quick” vs “large” runs.
   *Alternative*: Separate top-level dirs (e.g. `examples/pe_75chains/`) — rejected to avoid proliferation of example dirs.

2. **Source of large-scale inputs**
   Prefer copying or symlinking from bakery’s `examples/` (and `tests/` where needed) into `examples/<name>/large/`. Document in each example’s README that large inputs come from bakery and how to refresh them.
   *Rationale*: Single source of truth in bakery; no custom generators.
   *Alternative*: Scripts to download/generate — deferred; manual copy or documented path is sufficient for this change.

3. **Settings and LAMMPS inputs**
   Reuse the same `settings.yaml` schema as the small-scale example where possible. Large-scale may use a copy (e.g. `large/settings.yaml`) only if file paths or parameters differ (e.g. different `cg_conf.gro`). Generated `in.<name>` and `<name>.data` live under `large/` when paths differ.
   *Rationale*: Minimal duplication; same tool, same schema.

4. **CI and validation**
   Add an optional integration step (script or CI job) that runs `backmap-prep` and LAMMPS on at least one large-scale example (e.g. dodecane or pe). Mark as optional or time-limited so CI stays fast.
   *Rationale*: Proves the pipeline works at scale without blocking every PR on long runs.

## Risks / Trade-offs

- **Large binaries in repo**
  [Risk] Large GRO/DATA files may bloat the repo.
  [Mitigation] Prefer documenting “copy from bakery” or using Git LFS if we commit large files; start with one or two large examples committed and document the rest.

- **Divergence from bakery**
  [Risk] Bakery updates may make our copied inputs stale.
  [Mitigation] README and docs state that large-scale inputs are from bakery and give the path/version; users can refresh manually.

- **Runtime of large examples**
  [Risk] Full 75-chain or 500-molecule runs may be too slow for casual testing.
  [Mitigation] Document expected runtimes; CI runs a single short test (e.g. few steps) or skips long runs by default.

## Migration Plan

1. Add `large/` layout and docs (README, large-scale doc) without adding all large files yet.
2. Migrate one system (e.g. dodecane or pe) end-to-end: bakery inputs → `large/` → `backmap-prep` → LAMMPS run; document steps.
3. Repeat for remaining systems (pe4, pe_10, pe_aa, melamine); add optional CI script/job.
4. Update main README and docs to point to “Large-scale examples” and per-example READMEs.

No rollback beyond reverting commits; no database or API changes.

## Open Questions

- Whether to commit any large GRO/DATA files or rely entirely on “copy from bakery” instructions.
- Exact CI job name and trigger (e.g. `large-scale-integration` on schedule or manual only).
