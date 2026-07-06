## ADDED Requirements

### Requirement: Hybrid TOP cross-pair export

The network builder SHALL parse `[ cross_pairs ]`, resolve func-1 LJ parameters, and
emit `pairs.dat` plus `fix backmap/pairs`.

#### Scenario: Rim135 cross-pair count

- **WHEN** `backmap-prep build examples/epoxy/settings.v2.yaml` runs
- **THEN** `pairs.dat` SHALL contain approximately 3,600 pairs
