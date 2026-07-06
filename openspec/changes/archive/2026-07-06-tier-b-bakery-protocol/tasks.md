# Tasks: tier-b-bakery-protocol (PR4a)

- [x] OpenSpec proposal, design, delta specs
- [x] Update `examples/epoxy/settings.v2.yaml` (alpha, dt, gamma, cap_force, T)
- [x] Add `_write_initial_velocities` and `_write_cap_force` in `writers.py`
- [x] Unit tests for velocity, cap_force, langevin damp
- [x] Extend `test_rim135_build_v2_lammps_smoke` for PR4 protocol gates
- [x] PR4b: LAMMPS C++ `fix backmap/capforce` implementation
- [x] Local LAMMPS smoke100 (PR4c step 1)
- [x] Local 8k ramp (PR4c step 2 — full 8k+8k PASS, log.rim135.pr4, Tier B gate)
- [x] Local full 10k ramp (PR4c step 3 — optional; ramp-only done in log.pr4.full)
- [x] VM smoke when local passes (`log.rim135.ref.full.out`, Jul 2026)
- [x] Archive change when Tier B gate passes (local + VM)
