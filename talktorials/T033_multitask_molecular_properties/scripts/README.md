# T033 · Molecular properties: multitask modelling

## Scripts

This folder stores scripts used to prepare the material for the Jupyter notebook.
Students do **not** need to run them: their outputs are committed to the repository.

- `prepare_ochem_data.py`: turns the raw OCHEM property exports in `../data/raw`
  into the single tidy file the notebook reads,
  `../data/ochem_properties.csv.gz`. See `../data/README.md` for which property
  sets were exported and what the cleaning steps do.
- `make_figures.py`: generates the schematic figures in `../images`, so that they
  can be regenerated and tweaked together with the notebook.
