# PE Examples: Findings and Required Fixes

**Date:** 2025-05-15
**Branch:** `test/pe-examples`
**Status:** Investigation complete, two issues identified

---

## Summary

All four PE examples (`pe`, `pe4`, `pe_10`, `pe_aa`) were tested end-to-end
with `backmap-prep` generation and LAMMPS simulation. Two distinct issues were
found:

1. **Non-uniform atoms-per-bead** — `pe_10` and `pe_aa` fail at startup because
   all-atom PE has different atom counts per bead type (CH3 vs CH2). **Fixed** by
   adding the `apb` keyword to `fix backmap`.

2. **Production phase blowup** — All four examples crash during the backmapping
   ramp (Phase 2) with "Bond atoms missing" errors. This is a simulation
   stability issue, not a code bug.

---

## Issue 1: Non-Uniform Atoms-Per-Bead (FIXED)

### Problem

`fix backmap` required `n_at % n_cg == 0` — all CG beads must map to the same
number of AT atoms. This is impossible for all-atom PE with terminated chains:

| Example | CG beads | AT atoms/mol | Bead A (CH3-end) | Bead B (CH2-interior) |
|---------|----------|-------------|------------------|-----------------------|
| pe      | 50       | 100 (UA)    | 2                | 2                     |
| pe4     | 25       | 100 (UA)    | 4                | 4                     |
| pe_10   | 10       | 302 (AA)    | 31               | 30                    |
| pe_aa   | 50       | 302 (AA)    | 7                | 6                     |

United-atom examples (pe, pe4) have uniform mapping because CH3 and CH2 are
each one atom. All-atom examples have one extra hydrogen per terminal CH3 group
(3H vs 2H), making end beads (type A) contain one more atom than interior beads
(type B).

### Root Cause

`fix_backmap.cpp:build_bead_map()` divided AT atoms uniformly:
`bead i → AT atoms [i*apb, (i+1)*apb)` where `apb = n_at / n_cg`.

### Fix

Added `apb` keyword to `fix backmap`:

```
fix bm all backmap cg_type 1 2 alpha 0.0001 lambda0 0.0 apb 1:7 2:6
```

Syntax: `apb T1:N1 T2:N2 ...` — maps CG type `T` to `N` AT atoms per bead.
When absent, legacy uniform mode is used (backward compatible).

**Changes:**

| File | Change |
|------|--------|
| `src/fix_backmap.h` | Added `std::map<int,int> apb_map_` member |
| `src/fix_backmap.cpp` | Parse `apb` keyword; per-type slicing in `build_bead_map()` |
| `python/src/backmap_prep/builder.py` | Compute `apb_by_cg_type` from bead definitions |
| `python/src/backmap_prep/writers.py` | Emit `apb` keyword when CG types differ |
| `python/tests/test_builder.py` | Tests for uniform and non-uniform apb computation |
| `python/tests/test_writers.py` | Tests for apb keyword emission |

**Commits:** `ccdd1ac`, `b713a49`, `a4c2d03`

**Verification:** After the fix, all four PE examples pass `backmap-prep`
generation and the backmapping phase (Phase 1/2) completes with stable energies.
168 Python unit tests pass.

---

## Issue 2: Backmapping Ramp Blowup (OPEN)

### Problem

All four PE examples crash approximately 1000 steps into the backmapping ramp
(Phase 2) with "Bond atoms missing" errors. The energy increases by two orders
of magnitude within ~1000 steps, then atoms escape the simulation box.

### Observed Behaviour

| Example | Crash step | Energy before ramp | Energy at crash |
|---------|-----------|-------------------|-----------------|
| pe      | 11021     | ~1800             | ~111000         |
| pe4     | 11002     | ~1800             | ~21000          |
| pe_10   | 11002     | ~7400             | ~13800          |
| pe_aa   | 11002     | ~7100             | ~85000          |

The script structure is:

```
# Phase 1: λ frozen, CG equilibration (steps 0–10000)
fix_modify bm active no
run 10000

# Phase 2: Backmapping ramp (steps 10000–20000)
fix_modify bm active yes
run 10000

# Phase 3: AT production (steps 20000–30000)
run 10000
```

Phase 1 completes normally with stable energies. The crash occurs in Phase 2
after `fix_modify bm active yes` enables the lambda ramp.

### Analysis

When the ramp activates:

1. `end_of_step()` increments `lambda[i] += alpha` (alpha=0.0001), so after
   1000 steps lambda ≈ 0.1.
2. `initial_integrate()` moves CG beads toward AT centre-of-mass:
   `x_cg += lambda * (COM_AT - x_cg)`.
3. `post_force()` distributes CG forces to AT atoms proportional to mass ratio.

The blowup at lambda ≈ 0.1 suggests the CG→AT force distribution creates
instability. Possible causes:

- **CG forces too large for AT timestep.** The CG pair/bond potentials may
  produce forces appropriate for CG dynamics but too strong when distributed to
  AT atoms at 1 fs timestep.
- **No energy minimization before ramp.** Phase 1 runs NVE+Langevin but AT
  atoms may still have high-energy contacts after random placement within CG
  beads. A `minimize` step before Phase 2 could relieve these.
- **Missing CG force scaling by lambda.** The `post_force` distributes the full
  CG force `f_cg` to AT atoms regardless of lambda. The reference method
  (Krajniak et al., JCTC 2016) scales CG forces by `(1 - lambda)` as the
  resolution transitions. Currently `post_force` applies full CG forces at all
  times.
- **Thermostat coupling too weak.** Langevin damping of 30 (gamma_inv = 33.3 fs)
  may not dissipate energy fast enough during the ramp.

### Recommended Investigation

1. **Check force scaling.** In the reference paper (Eq. 5–7), CG potentials are
   scaled by `(1 - λ)` and AT potentials by `λ`. The current implementation
   does NOT scale CG forces — `post_force` applies full CG force. This is
   likely the primary cause.

2. **Add minimize before ramp.** Insert `minimize 1.0e-4 1.0e-6 1000 10000`
   between Phase 1 and Phase 2 to relieve high-energy AT contacts.

3. **Reduce backmapping timestep.** Try 0.1 fs during Phase 2.

4. **Increase thermostat coupling.** Try gamma=100 or higher during the ramp.

### Priority

**High.** Without this fix, no PE example produces a complete backmapped
configuration. The `(1 - λ)` force scaling is the most likely root cause and
should be implemented first.

---

## Test Matrix

| Example | backmap-prep | Phase 1 (equil) | Phase 2 (ramp) | Phase 3 (prod) |
|---------|-------------|-----------------|----------------|----------------|
| pe      | PASS        | PASS            | FAIL (blowup)  | —              |
| pe4     | PASS        | PASS            | FAIL (blowup)  | —              |
| pe_10   | PASS (apb)  | PASS            | FAIL (blowup)  | —              |
| pe_aa   | PASS (apb)  | PASS            | FAIL (blowup)  | —              |
| Unit tests (168) | PASS | — | — | — |
