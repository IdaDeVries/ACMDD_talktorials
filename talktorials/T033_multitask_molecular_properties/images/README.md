# T033 · Molecular properties: multitask modelling

## Images

This folder stores images used in the Jupyter notebook. All of them are generated
by `../scripts/make_figures.py` and can be regenerated with

```bash
python scripts/make_figures.py
```

- `multitask_concept.png`: single-task models vs. one shared multitask model.
- `multitask_strategies.png`: the three ways of obtaining several outputs
  (wrapper, native multi-output, boosted multi-output trees).
- `sparse_label_matrix.png`: the sparse label matrix and the two ways of dealing
  with missing measurements (complete cases vs. per-task masking).
- `splitting_methods.png`: random, scaffold, cluster and stratified random splits.
