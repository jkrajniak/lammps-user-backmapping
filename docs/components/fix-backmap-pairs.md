# fix backmap/pairs

## Syntax

```
fix ID group backmap/pairs at|cg file pairs.dat [cut Rc]
```

- **at|cg** -- lambda weighting of the listed pairs (`at`: atomistic 1-4
  pairs; `cg`: CG pairs)
- **file** -- pairs file (format below)
- **cut** -- optional cutoff; pairs farther apart are skipped (default: no
  cutoff, as for any bonded term)

## Description

Applies the 1-4 interactions of the atomistic force field. `special_bonds`
excludes every 1-4 pair from the pair style, so each 1-4 pair of the source
topology (within a bead and across beads) is listed here with its own
parameters:

\[
E_{ij} = w \left[ 4\epsilon_{ij}\left( \left(\frac{\sigma_{ij}}{r}\right)^{12}
 - \left(\frac{\sigma_{ij}}{r}\right)^{6} \right)
 + s_{ij}\, C \frac{q_i q_j}{r} \right]
\]

with \( C \) = `force->qqrd2e`, \( s_{ij} \) the 1-4 Coulomb scale (GROMACS
`fudgeQQ`) and \( w \) the lambda weight of `fix backmap` (1 for a pair
within one bead; \( \lambda \) for an `at` pair across beads).

`backmap-prep` writes `pairs.dat` from the hybrid topology's `[ pairs ]` and
`[ cross_pairs ]` sections.

## Pairs file

```
N
id1 id2 sigma epsilon [qq_scale]
...
```

`sigma` in distance units, `epsilon` in energy units. The optional fifth
column is the 1-4 Coulomb scale; a four-column line is LJ only. Charges come
from the atoms, so a nonzero scale needs an atom style with charges.

## Output

- Global scalar: total 1-4 energy (LJ + Coulomb).
- Global vector (2): 1-4 LJ energy, 1-4 Coulomb energy.

The energy is part of the potential energy (`pe`) and the virial is part of
the pressure by default; `fix_modify ID energy no` / `virial no` remove them.
Forces and energy are included from the first force evaluation of a run
(`run 0`) and in energy minimization.

## Parallel runs

The partner of every pair must be present as an owned atom or a ghost on each
rank that owns the other atom; otherwise the run stops with an error. Set
`comm_modify cutoff` above the longest 1-4 distance (`backmap-prep` does this).

Each MPI rank evaluates every pair with at least one of its own atoms and
applies the force to its own atoms only, so the result does not depend on
the domain decomposition or on `newton_pair`.

## Related

- [fix backmap](fix-backmap.md)
- [pair_style backmap](pair-backmap.md)
