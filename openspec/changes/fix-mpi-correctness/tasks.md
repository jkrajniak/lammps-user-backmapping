# Tasks — MPI correctness for fix backmap

- [x] 1.1 Add per-atom `atom2cg`, `com_buf`, `cg_fwd` arrays to `fix_backmap.h/cpp`
  with `grow_arrays` / destructor management; set `comm_forward = 5`,
  `comm_reverse = 4`.
- [x] 1.2 Rework `build_bead_map` to populate `atom2cg` (local+ghost) while
  keeping `bead_map` for logging/`validate_masses`.
- [x] 1.3 Rewrite `initial_integrate` COM update to the atom-centric +
  `reverse_comm(4)` flow.
- [x] 1.4 Rewrite `post_force` force distribution to the `forward_comm` +
  local-AT-only flow.
- [x] 1.5 Update `pack_forward_comm` / `unpack_forward_comm` to 5 doubles; add
  `pack_reverse_comm` / `unpack_reverse_comm`.
- [x] 1.6 Warn if any local AT atom has `atom2cg == -1` (insufficient
  `comm_modify cutoff`).
- [ ] 2.1 Build serial locally; run dodecane 1000-step regression; confirm
  bit-identical to pre-fix. *(blocked: VM down; local build is MPI STUBS)*
- [x] 2.2 Create `examples/dodecane/large/in.dodecane_mpi` and
  `test_mpi_serial_vs_4rank.sh`.
- [ ] 2.3 On VM: `mpirun -np 4` vs serial; compare positions/energy; record
  evidence in `research/`. *(blocked: VM down)*
- [x] 3.1 Update `CHANGELOG.md`, `docs/components/fix-backmap.md` (MPI section).
- [ ] 3.2 Update `STATUS.md`, `paper/main.tex` limitations after 4-rank PASS.
