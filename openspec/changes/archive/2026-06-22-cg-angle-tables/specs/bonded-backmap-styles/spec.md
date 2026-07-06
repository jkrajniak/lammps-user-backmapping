## ADDED Requirements

### Requirement: Angle style `backmap/table`

The `angle_style backmap/table` SHALL compute tabulated angle forces with lambda
weighting:

```
F = w × F_table(θ)
E = w × E_table(θ)
```

The weight `w` is computed from the lambda values of the first and last atoms of
the angle (atoms i and k in the i-j-k triple), using the same AT/CG rules as
`backmap/harmonic`.

Command syntax:

```
angle_style backmap/table linear N
angle_coeff M at/cg filename keyword
```

Used within `angle_style hybrid`:

```
angle_style hybrid backmap/harmonic backmap/table linear 1000
angle_coeff 2 backmap/table cg table_a1.table ENTRY
```

#### Scenario: CG tabulated angle at lambda=0

- **WHEN** `angle_coeff 2 backmap/table cg table_a1.table ENTRY` and λ_i = λ_k = 0
- **THEN** the angle force SHALL be at full tabulated strength (weight = 1.0)

#### Scenario: CG tabulated angle at lambda=0.5

- **WHEN** `angle_coeff 2 backmap/table cg table_a1.table ENTRY` and λ_i = λ_k = 0.5
- **THEN** the angle energy SHALL be scaled by 0.75 (weight = 1 − 0.25)

#### Scenario: Missing fix backmap

- **WHEN** a `backmap/table` angle is defined but no `fix backmap` exists
- **THEN** the style SHALL abort with: "angle_style backmap/table requires fix backmap"
