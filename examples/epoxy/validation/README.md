# rim135 structural validation artefacts

Validation of the VM backmap run against GROMACS reference RDFs from
[paper-reverse-mapping-polymer-networks](https://github.com/jkrajniak/paper-reverse-mapping-polymer-networks)
(`paper/rim135_small/`).

## Paper-facing files (git-tracked)

| File | Role |
|------|------|
| `structural_validation_report.txt` | Full PASS/FAIL summary from validator |
| `rdf_comparison.png` | C–O and C–N overlay vs paper AA reference |
| `../compare_rim135_structure.py` | Reproducible validator |

## Large inputs (local / VM only)

| File | Source |
|------|--------|
| `dump.backmap` | VM `/home/azureuser/sc/rim135/dump.backmap` |
| `rim135_final.data` | VM `/home/azureuser/sc/rim135/rim135_final.data` |

Re-fetch from VM:

```bash
scp -i ~/.ssh/ai_vm_key \
  azureuser@10.110.0.4:/home/azureuser/sc/rim135/{dump.backmap,rim135_final.data} \
  .
```

## Re-run

```bash
cd lammps-user-backmapping
uv run examples/epoxy/compare_rim135_structure.py \
  --dump examples/epoxy/validation/dump.backmap \
  --gro "$BACKMAP_RIM135/hyb_conf.gro" \
  --final-data examples/epoxy/validation/rim135_final.data \
  --paper-dir ../paper-reverse-mapping-polymer-networks/paper/rim135_small \
  --plot examples/epoxy/validation/rdf_comparison.png \
  | tee examples/epoxy/validation/structural_validation_report.txt
```

## Interpretation for manuscript

**Claimed:** first-shell C–O / C–N peak positions and heights vs paper AA reference
(4/4 PASS, Jul 2026).

**Not claimed:** full-curve L2, free volume, or long-trajectory statistical match —
see `research/notebook/2026-07-06_rim135-structural-validation-vs-jcc2017.md`.
