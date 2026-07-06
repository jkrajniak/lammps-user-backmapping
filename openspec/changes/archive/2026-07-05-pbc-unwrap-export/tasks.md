# Tasks — molecule-aware PBC unwrap

## 1. OpenSpec

- [x] 1.1 Create change artifacts and delta specs

## 2. Python

- [x] 2.1 Per-`mol_id` intra-fragment unwrap in `pbc.py`
- [x] 2.2 Spanning-tree + iterative inter-mol molecule translation
- [x] 2.3 Bond-tree image flags; min-image gate in `validate_bond_geometry`
- [x] 2.4 Network comm cutoff + `reset_atoms image all` in writers

## 3. Validation

- [x] 3.1 Unit tests (cycle, cross-mol bond)
- [x] 3.2 Rim135 integration: min-image max bond < 20 Å
- [x] 3.3 Rebuild rim135 data; comm cutoff ≥ 115 Å in input
