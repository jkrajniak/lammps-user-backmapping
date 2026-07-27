# PET / Dacron network example (3 species: TER + DIO + H2O)

Migrates the polyethylene terephthalate (Dacron) polyester network from the
legacy bakery `pete` example to `backmap-prep` + LAMMPS `fix backmap`.

Reference: Krajniak, Zhang et al., *J. Comput. Chem.* 2018, 39, 648-664
(10.1002/jcc.25129). Source inputs: JCC 2018 published reproducibility
archive, `paper-reverse-mapping-polymer-networks/preparation/dacron/backmapping/`.

## System

- 3 CG molecule species: TER (terephthalate), DIO (ethylene glycol), H2O
- 7 CG bead types: A, B, C, D, E, W (Z transitional, unused post-reaction)
- 5757 CG atoms / 33757 AT (hybrid) atoms, 7.13 nm cubic box
- Degree-dependent bead maps, `charge_map` + `equilibrate_charges` on DIO/TER

## Files in this directory

| File | Committed? | Notes |
|------|-----------|-------|
| `settings.v2.yaml` | yes | backmap-prep v2 settings; `data_dir`/`tables_dir` point to the JCC 2018 paper-data archive |
| `cross_angles_pete.yaml`, `cross_dihedrals_pete.yaml` | yes | placeholders (`[]`) for the 54 ester-spanning cross-bonded terms bakery could not auto-derive (see `missing_definitions.txt` in the paper-data archive) |
| `report_dacron_ref_peaks.py` | yes | Tier C reference peak extractor (published RDFs) |
| `in.pet` | yes | LAMMPS input, robust 3-phase protocol (minimize -> nve/limit -> lambda ramp with `capforce` -> NVT), **hand-patched**: `comm_modify cutoff 75`, `reset_atoms image all` (see Known Issues; the previously-applied 2602-bond exclusion-list workaround was removed 2026-07-13 — see notebook `2026-07-13_pet-exclusion-list-verification.md`) |
| `table_A_*.table` ... `table_W_W.table` (21 pair tables) | yes | small (<40 KB each), converted IBI pair potentials |
| `pairs.dat` | yes (483 KB) | 1-4 pair list for `fix backmap/pairs` |
| `pet.data` | yes (4.8 MB) | LAMMPS data file (system topology + coordinates); regenerate via `backmap-prep build` if it ever needs refreshing (see below) |
| `table_a0.table`..`table_a3.table`, `table_d0.table`, `table_d1.table` | yes (1.8-3.7 MB each) | large angle/dihedral tables — regenerate the same way |
| `log.pet_*.lammps`, `pet_tierb.out` | no | VM run logs; regenerate on demand, don't keep stale ones around (gitignored) |

## How to regenerate the large files

1. Have the JCC 2018 paper-data archive available at
   `../../../../paper-reverse-mapping-polymer-networks/preparation/dacron/backmapping/`
   relative to this directory (already the workspace layout).
2. From the repo root:
   ```bash
   uv run backmap-prep build examples/pet/large/settings.v2.yaml
   ```
   This writes `pet.data`, `in.pet`, `pairs.dat`, and all `table_*.table` files
   into the **paper-data archive directory** (backmap-prep's `data_dir`), not
   here. Copy them into this directory, then apply the hand-patches below.
3. **Hand-patches required** (backmap-prep does not yet do these automatically
   for PET — see Known Issues):
   - `comm_modify cutoff 75.00` (network bonds span the box; auto-detected
     cutoff of 15 Å is wrong for this system)
   - `reset_atoms image all` (missing after `read_data`)
   - Do **not** add synthetic zero-strength exclusion bonds. A prior
     debugging session (2026-07-11/12) added 2602 fake bond-type-20
     "exclusions" to paper over close AT-AT contacts; re-verifying on
     2026-07-13 showed the real exclusion generator (bakery's
     `_generate_exclusion_lists`) already matches LAMMPS's bond-graph-derived
     special-neighbor list **exactly** (0 mismatches) on a clean build — the
     synthetic bonds were based on a stale build and actually introduced
     18,392 incorrect over-exclusions. See
     `research/notebook/2026-07-13_pet-exclusion-list-verification.md`.

## Tier status (2026-07-11/12)

| Tier | Status |
|------|--------|
| A (generator parity) | **PASS** — `hyb_conf.gro` byte-identical to bakery, `[ bonds ]` line-identical |
| B (LAMMPS dynamics) | **blocked** — see Known Issues |
| C (RDF vs published) | not started (blocked on B) |

## Known issues (2026-07-11/12 debugging)

1. **AT fragment overlap across CG-bead boundaries.** Shift-and-lift places
   each bead's AT fragment independently, centered on its own CG bead
   position, with correct *internal* geometry but *uncoordinated* orientation
   relative to neighboring beads. For beads only ~1.5-2.6 Å apart (bonded CG
   neighbors), this can place AT atoms from *different* beads within
   0.1-0.2 Å of each other, even when those atoms are far apart in the
   covalent bond graph (e.g. found a real "1-6" pair, `TER` ring atom 16840 to
   `TER` ester atom 16844 within the *same* TER molecule, 0.116 Å apart, 5
   bonds distant — see `research/notebook/2026-07-11_pet-overlap-diagnosis.md`).
2. **`special_bonds lj 0 0 0` (nrexcl=3) does not exclude these pairs — and
   that is correct, not a bug.** Re-verified 2026-07-13: bakery's
   `exclusion_hyb_topol.list` and LAMMPS's bond-graph-derived special-neighbor
   list match **exactly** (125,683 pairs, 0 differing members) on a clean
   build. The AT-AT overlaps described in point 1 are genuinely *not*
   bonded-graph-adjacent (they are spatially close only because of the
   placement bug) — the correct topological exclusion list rightly does not
   exclude them, so they are computed at full non-bonded strength, which is
   the actual source of instability once AT-AT intra-bead terms are always-on
   (see `docs/theory.md`).
3. **`fix backmap/capforce`** (added to the repo 2026-07-12, was previously
   VM-only) correctly bounds the *force* magnitude (confirmed via source
   review: `f[i] *= fmax/|f|` in `POST_FORCE`), evidenced by stable
   Temp/KinEng even while `PotEng` is enormous (~6e8) — but does not resolve
   the underlying overlap fast enough. The system reaches a persistent
   high-PotEng plateau, and `fix backmap`'s COM-tracking eventually drags a
   CG bead into a position that violates a **different, unrelated CG-CG
   table's inner cutoff** (`ERROR: Pair distance < table inner cutoff`,
   `src/pair_table.cpp:991`), currently around lambda ~0.15.
4. **Best next step (2026-07-13):** the exclusion-list question is closed
   (point 2, above) — it is not the blocker. The real blocker is the rigid,
   translate-only AT-fragment placement (point 1): TER's ring-split geometry
   spans CG-bead boundaries, and shift-and-lift has no mechanism to
   coordinate each bead's fragment orientation with its neighbors', so
   inter-bead angles/contacts can land arbitrarily close. Fixing this
   requires an orientation-aware placement step in `backmap-prep`'s network
   builder (`network/bakery/structures.py::prepare_hybrid`), not further
   LAMMPS-input-script tuning (ramp schedule, thermostat, or anneal cycles
   were all tried in the 2026-07-13 session without success — see
   `research/notebook/2026-07-13_at-intrabead-always-on-fix.md`).
