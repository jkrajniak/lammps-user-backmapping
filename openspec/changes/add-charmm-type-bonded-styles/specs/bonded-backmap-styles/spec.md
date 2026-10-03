## ADDED Requirements

### Requirement: Urey-Bradley angle backmap style

The package SHALL provide `angle_style backmap/charmm` with coefficients
`at|cg K theta0 K_ub r_ub` and energy
`E = w [K (theta - theta0)^2 + K_ub (r13 - r_ub)^2]`, where w is the
`compute_weight3` weight of the angle.

#### Scenario: Equals the stock style at full weight
- **WHEN** an `at` angle spans two CG beads and lambda_global = 1
- **THEN** its energy and forces SHALL equal those of LAMMPS `angle_style charmm`
  with the same K, theta0, K_ub, r_ub

#### Scenario: Weighted across beads
- **WHEN** lambda_global = 0.5 and the angle is `at` across beads
- **THEN** energy, forces and virial SHALL be half the stock values

### Requirement: Multi-term periodic dihedral backmap style

The package SHALL provide `dihedral_style backmap/fourier` with coefficients
`at|cg m K1 n1 d1 ... Km nm dm` and energy
`E = w sum_i K_i [1 + cos(n_i phi - d_i)]` with arbitrary phase d_i (degrees).

#### Scenario: Equals the stock style at full weight
- **WHEN** an intra-bead dihedral is evaluated at any lambda
- **THEN** its energy and forces SHALL equal those of LAMMPS
  `dihedral_style fourier` with the same terms

### Requirement: Harmonic improper backmap style

The package SHALL provide `improper_style backmap/harmonic` with coefficients
`at|cg K chi` and energy `E = w K (chi - chi0)^2`.

#### Scenario: Weighted CG improper
- **WHEN** a `cg` improper is evaluated at lambda_global = 0.25
- **THEN** its energy SHALL be 0.75 times the LAMMPS `improper_style harmonic` energy
