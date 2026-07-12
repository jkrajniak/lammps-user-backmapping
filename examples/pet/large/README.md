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
| `in.pet` | yes | LAMMPS input, robust 3-phase protocol (minimize -> nve/limit -> lambda ramp with `capforce` -> NVT), **hand-patched**: `comm_modify cutoff 75`, `reset_atoms image all`, 3048 zero-strength exclusion bonds (type 20) appended to `pet.data` (see Known Issues) |
| `table_A_*.table` ... `table_W_W.table` (21 pair tables) | yes | small (<40 KB each), converted IBI pair potentials |
| `pairs.dat` | yes (483 KB, under 500 KB pre-commit limit) | 1-4 pair list for `fix backmap/pairs` |
| `pet.data` | **no** (4.8 MB) | generated LAMMPS data file — regenerate via `backmap-prep build` (see below) |
| `table_a0.table`..`table_a3.table`, `table_d0.table`, `table_d1.table` | **no** (1.8-3.7 MB each) | large angle/dihedral tables — regenerate the same way |
| `log.pet_tierb.lammps`, `pet_tierb.out` | no | most recent VM run log, kept locally for debugging reference |

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
   - Append the 3048 exclusion bonds (bond type 20, `at 0.0 1.0`) — see
     `add_exclusion_bonds.py` (TODO: not yet committed; currently done by an
     ad-hoc script during debugging, needs to be turned into a proper CLI
     step or fixed upstream in `backmap-prep`)

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
2. **`special_bonds lj 0 0 0` (nrexcl=3) does not exclude these pairs** — only
   1-2/1-3/1-4. The bakery reference `exclusion_hyb_topol.list` has 125683
   pairs total (same count as LAMMPS's computed special-neighbor list) but
   **differs by 3048 members** — some close pairs bakery excludes are not
   caught by LAMMPS's bond-graph walk, and vice versa. Root cause of the
   count-matches-but-members-differ mismatch not yet identified.
3. **`fix backmap/capforce`** (added to the repo 2026-07-12, was previously
   VM-only) correctly bounds the *force* magnitude (confirmed via source
   review: `f[i] *= fmax/|f|` in `POST_FORCE`), evidenced by stable
   Temp/KinEng even while `PotEng` is enormous (~6e8) — but does not resolve
   the underlying overlap fast enough. The system reaches a persistent
   high-PotEng plateau, and `fix backmap`'s COM-tracking eventually drags a
   CG bead into a position that violates a **different, unrelated CG-CG
   table's inner cutoff** (`ERROR: Pair distance < table inner cutoff`,
   `src/pair_table.cpp:991`), currently around lambda ~0.15.
4. **Best next step (undecided as of 2026-07-12):** extend the CG-frozen AT
   relaxation phase (Phase 0b) much longer before starting the lambda ramp,
   so overlaps resolve before COM-tracking begins moving CG beads; or
   investigate why the exclusion-list members differ from bakery's.
