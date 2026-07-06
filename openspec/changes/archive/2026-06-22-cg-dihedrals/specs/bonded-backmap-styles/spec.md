## ADDED Requirements

### Requirement: Dihedral style `ryckaert` (static)

The package SHALL provide `dihedral_style ryckaert` for intra-bead RB dihedrals without
lambda weighting:

```
dihedral_style ryckaert
dihedral_coeff N C0 C1 C2 C3 C4 C5
```

#### Scenario: Intra-bead RB dihedral

- **WHEN** all four atoms of a dihedral belong to the same CG bead
- **THEN** the generator SHALL assign a static `ryckaert` type with converted C0..C5
