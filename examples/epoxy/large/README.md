# rim135 Tier B/C production script

`in.rim135` is the LAMMPS input for the full 3-phase backmapping protocol
(minimize + AT relax with CG frozen -> lambda ramp with thermostat -> NVT
production with dt ramp-up). See the header comment in the file for phase
details and the rationale for the thermostat during the ramp.

## Regenerating inputs

The `.data` / `.table` / `pairs.dat` files this script reads are not
committed (large, regenerable from bakery reference data):

```bash
export BACKMAP_RIM135=/path/to/bakery/tests/rim135
cd "$BACKMAP_RIM135"
uv run --directory ../../../../../lammps-user-backmapping backmap-prep build-hybrid \
  ../../../../../lammps-user-backmapping/examples/epoxy/settings.v2.yaml
```

Then copy `in.rim135` next to the generated `rim135.data` / `table_*.table` /
`pairs.dat` and run:

```bash
mpirun -np 4 lmp -in in.rim135
```

## Status

VM regression (2026-07-12, np=4): **PASS**. T stable at ~296-303 K throughout
Phase 1 (lambda ramp) and Phase 2 (production up to dt=1.0 fs). See
`research/notebook/2026-07-12_lambda-weighting-vm-regression.md`.
