# Tasks

## 1. C++
- [ ] 1.1 `dihedral_style backmap/opls` (subclass of backmap/ryckaert; coeff
      conversion; write_data and restart keep K1..K4)
- [ ] 1.2 Regression test vs stock `dihedral_style opls` (intra-bead, across
      beads at lambda 1 and 0.5, cg at lambda 0.25), energy and forces
- [ ] 1.3 Finite-difference force entry

## 2. backmap-prep
- [ ] 2.1 Native fragment parser accepts `dihedral_style opls`, converts to RB
- [ ] 2.2 Unit tests: conversion equals the OPLS energy at several angles;
      ryckaert path unchanged

## 3. Docs
- [ ] 3.1 `docs/components/dihedral-styles.md`, README, CHANGELOG
