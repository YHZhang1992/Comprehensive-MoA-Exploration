# Comprehensive mechanism-of-action exploration

Runnable reference workflow for prioritizing signals across model characterization, treatment reversal, network importance, and mouse-to-human translation.

## Quick start

```bash
python workflow.py --input examples/input.csv --config config/workflow.json --output output
python -m unittest discover -s tests -v
```

Outputs include gene-level rankings, pathway summaries, and a checksum-bearing run manifest. The example is synthetic and dependency-free so repository health can be verified on any Python 3.10+ installation.

## Intended analysis stages

1. Characterize disease-model signals at gene, pathway, and network levels.
2. Measure whether ASO/siRNA treatment reverses model-associated effects.
3. Evaluate concordance with human disease evidence using explicit orthology and annotation releases.
4. Prioritize results for experimental review; do not interpret the composite score as causal proof.

See `docs/RUNBOOK.md` for required provenance and review steps.
