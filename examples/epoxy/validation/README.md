# rim135 structural validation artefacts

Validation of the VM backmap run against GROMACS reference RDFs from
[paper-reverse-mapping-polymer-networks](https://github.com/jkrajniak/paper-reverse-mapping-polymer-networks)
(`paper/rim135_small/`).

## Paper-facing files (git-tracked)

| File | Role |
|------|------|
| `structural_validation_report.txt` | Full PASS/FAIL summary from validator |
| `rdf_comparison.pdf` | C–O and C–N overlay vs paper AA reference — **vector, promoted to the manuscript** |
| `rdf_comparison.png` | Same plot, raster; kept for quick viewing |
| `../compare_rim135_structure.py` | Reproducible validator |

## Large inputs (local / VM only)

| File | Source |
|------|--------|
| `dump.backmap` | VM `/home/<vm_user>/sc/rim135/dump.backmap` |
| `rim135_final.data` | VM `/home/<vm_user>/sc/rim135/rim135_final.data` |

Re-fetch from VM:

```bash
scp -i ~/.ssh/<vm_key> \
  <vm_user>@<vm_host>:/home/<vm_user>/sc/rim135/{dump.backmap,rim135_final.data} \
  .
```

## Re-run

```bash
cd lammps-user-backmapping
REF=../paper-reverse-mapping-polymer-networks
uv run examples/epoxy/compare_rim135_structure.py \
  --dump examples/epoxy/validation/dump.backmap \
  --gro "$REF/preparation/rim135_small/backmapping/aa/hyb_conf.gro" \
  --final-data examples/epoxy/validation/rim135_final.data \
  --paper-dir "$REF/paper/rim135_small" \
  --plot examples/epoxy/validation/rdf_comparison.pdf \
  | tee examples/epoxy/validation/structural_validation_report.txt
```

The validator exits non-zero whenever any check fails, which is the expected
state here (the L2 and total-charge checks fail by design — see
*Interpretation* below). A non-zero exit is not a regression on its own;
compare the report against the committed one.

**Use the `aa/` `hyb_conf.gro`, not `ua/` or `pure/`.** `--gro` supplies the
atom-id → element map that selects which atoms enter each RDF, so the wrong
file silently produces plausible-looking but wrong curves rather than an
error. Only `aa/` has 13,653 atoms, matching `dump.backmap` and
`rim135_final.data`; `ua/` has 6,720 and `pure/` has 13,680. Checking that
atom count is the quickest way to confirm the right input. Passing `ua/`
shifts the computed C–O peak from 0.241 nm to 0.217 nm and turns three
passing metrics into failures, with no other symptom.

## Interpretation for manuscript

**Claimed:** first-shell C–O / C–N peak positions and heights vs paper AA reference
(4/4 PASS, Jul 2026).

**Not claimed:** full-curve L2, free volume, or long-trajectory statistical match —
see `research/notebook/2026-07-06_rim135-structural-validation-vs-jcc2017.md`.
