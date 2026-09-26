# Sources of the third-party files

| Files | Origin | Version | Licence |
|-------|--------|---------|---------|
| `martini_v3.0.0.itp`, `martini_v3.0.0_phospholipids_v1.itp`, `martini_v3.0.0_solvents_v1.itp` (downloaded by `fetch_sources.sh`) | [marrink-lab/martini-forcefields](https://github.com/marrink-lab/martini-forcefields) `martini_forcefields/regular/v3.0.0/gmx_files` | commit `784591ebdc91d762ed4df986c4650546c938f776` | Apache-2.0 |
| `amber99sb-ildn_slipids.ff/{forcefield,ffnonbonded,ffbonded,gbsa}.itp` (downloaded by `fetch_sources.sh`) | [owenvickery/cg2at](https://github.com/owenvickery/cg2at) `database/forcefields/amber99sb-ildn_slipids.ff` | commit `5fd704204067782e2ed01bec9c36491f088c12c4` | GPL-3.0 |
| `POPC.itp`, `POPC.pdb` (committed) | cg2at `database/fragments/martini_3-0_slipids/non_protein/POPC` | same commit | GPL-3.0 |
| `w4_single.gro`, `topol_w4.top` (committed, derived) | cg2at `database/fragments/martini_3-0_slipids/solvent/W/{TIP3P.pdb,tip3p.itp}` | same commit | GPL-3.0 |

`fetch_sources.sh` checks the SHA-256 of every download.

Derived files and what changed (all committed text files have trailing whitespace
stripped by the repository's pre-commit hooks; no other change to `POPC.itp` / `POPC.pdb`):

- `popc_single.gro`: `POPC.pdb` coordinates (Å -> nm) in the atom order of `POPC.itp`.
- `topol_popc.top`: `POPC.itp` inlined (bakery does not follow `#include`), position-restraint include dropped.
- `w4_single.gro`, `topol_w4.top`: the four TIP3P waters of one W bead as one
  residue `W` (bakery fragments are single residues) with unique atom names
  (`OW1 HW11 HW12` ... `OW4 HW41 HW42`); flexible TIP3P bonds and angle
  (the `FLEXIBLE` branch of `tip3p.itp`), nrexcl 2.
- `topol_cg.top`: written by `make_cg_topology.py` from the MARTINI files.
- `cg_conf.gro`: last frame of a 1 µs MARTINI 3 NPT run (310 K) of an
  `insane` bilayer (128 POPC, 1752 W); run record in the research notebook
  (`experiments/20260925_martini-popc-feasibility.md`).
