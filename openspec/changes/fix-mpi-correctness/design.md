# Design — MPI correctness for fix backmap

## Current (MPI-incorrect) data flow

`build_bead_map` builds `bead_map`: one entry per *local* CG bead, holding the
local+ghost indices of its mapped AT atoms. `initial_integrate` and
`post_force` loop over `bead_map` and operate on the CG-owner rank only.

- COM: sum over `bm.at_local` (local+ghost AT). Missing AT atoms on non-ghost
  ranks.
- Force: `f[at] += ratio * f_cg` for `at` in `bm.at_local` incl. ghosts. Ghost
  `f` writes are lost; AT atoms with a remote non-ghosted CG get nothing.

## New (MPI-correct) data flow

### Per-atom arrays (grown with `nmax`)

- `int *atom2cg` (size `maxatom`): for every AT atom, the local-or-ghost index
  of its mapped CG bead, `-1` for CG atoms and unmapped atoms. Built in
  `build_bead_map` over `nlocal + nghost`.
- `double *com_buf` (size `maxatom * 4`): scratch accumulator per CG atom
  `(mass, mdx, mdy, mdz)`. Indexed by CG atom local-or-ghost index.
- `double *cg_fwd` (size `maxatom * 4`): forward-comm scratch for ghost CG
  atoms `(fx, fy, fz, mass)`. Local CG atoms read `f[cg]` and `mass[type[cg]]`
  directly; ghost CG atoms read `cg_fwd`.

### Communication sizes

- `comm_forward = 5`: `lambda` + `(cg_fx, cg_fy, cg_fz, cg_mass)` (zeros for
  non-CG atoms; only read for CG atoms).
- `comm_reverse = 4`: COM `(mass, mdx, mdy, mdz)`.

### `initial_integrate` (COM update)

1. Zero `com_buf` for all CG atoms (local+ghost).
2. For each **local** AT atom `i`: `cg = atom2cg[i]`; if `cg < 0` skip. Compute
   minimum-image `dx = x[i] - x[cg]`. `com_buf[cg] += (mass[type[i]],
   mass[type[i]]*dx)`.
3. `comm->reverse_comm(this, 4)` — sums ghost-CG `com_buf` into local-CG owners.
4. For each **local** CG atom `cg`: read `com_buf[cg] = (m_tot, mdx, mdy, mdz)`;
   if `m_tot > 0`: `x[cg] += lambda[cg] * (mdx, mdy, mdz) / m_tot`.
5. Zero CG velocities/forces on local CG atoms (unchanged).

Each AT atom is local on exactly one rank, so every AT atom contributes once;
`reverse_comm` routes contributions to the CG owner. Correctness no longer
depends on the CG bead being the "owner" rank of the loop.

### `post_force` (force distribution)

1. `comm->forward_comm(this)` — pack `lambda` + CG `(fx, fy, fz, mass)`; ghost
   CG atoms now carry their owner's full CG force.
2. For each **local** AT atom `i`: `cg = atom2cg[i]`; if `cg < 0` skip. Read
   `(fx_cg, fy_cg, fz_cg, m_cg)` from `f[cg]`+`mass` if `cg` local, else from
   `cg_fwd[cg]`. `f[i] += (mass[type[i]] / m_cg) * (fx_cg, fy_cg, fz_cg)`.
3. Zero `f` on **local** CG atoms (after distribution).

All force additions target local AT atoms → no reverse-comm needed. The CG
force is complete on its owner rank (LAMMPS `newton on` gives full forces on
local atoms) and is forwarded to ghost images where remote AT atoms read it.

### `build_bead_map`

- Iterate `nlocal + nghost`, group by molecule, deduplicate by tag (prefer
  local), as today.
- For each molecule with a CG bead: assign each AT atom (by tag-sorted index)
  to its CG bead (by tag-sorted index, `apb` AT per CG).
- Fill `atom2cg[at_local_idx] = cg_local_idx` for every AT atom (local or
  ghost) whose CG bead is present (local or ghost) on this rank.
- Keep `bead_map` for the setup-time count/log and `validate_masses`.

### Ghost-cutoff requirement

For `atom2cg[i]` to be non-negative, each AT atom's CG bead must be local or
ghost on the AT atom's rank. This requires `comm_modify cutoff >=` maximum
CG-AT distance within a bead. For backmapped fragments this is a few Å; the
existing examples already set `comm_modify cutoff 15`. The fix documents this
requirement and warns if `atom2cg` is `-1` for any local AT atom (indicates
insufficient cutoff).

## Serial behaviour

On 1 rank, `nghost` atoms still exist (within `comm_modify cutoff`).
`reverse_comm` and `forward_comm` are no-ops on 1 rank. The atom-centric flow
produces identical results to the current CG-owner-centric flow because every
AT atom is local and every CG bead is local. No numerical change to serial
validation.

## Test plan

1. **Serial regression**: rebuild, run `examples/dodecane/large/in.dodecane`
   for 1000 steps; compare final `pe`/`etotal`/positions to pre-fix baseline
   within `1e-10` (bit-identical expected).
2. **MPI correctness** (VM, 4 ranks): small dodecane, `mpirun -np 4 lmp -in
   in.dodecane_mpi`; compare final data file to serial run within tolerance
   (positions `1e-4` Å, energy `1e-6`). Script
   `examples/dodecane/large/test_mpi_serial_vs_4rank.sh`.
3. Unit: `python/tests` unchanged; add a C++ unit if a `fix_backmap` harness
   exists (none currently — deferred).
