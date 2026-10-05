# T033 · Molecular properties: multitask modelling

## Data

This folder stores the input data for the Jupyter notebook.

- `ochem_properties.csv.gz`: the tidy, long-format property table the notebook reads.
- `raw/`: the raw OCHEM exports the table was built from. Not committed — see below.

### Where the data come from

The measurements come from the [OCHEM database](https://ochem.eu) (the *Online
CHEmical Modeling environment*), a curated collection of experimental
physicochemical, ADME and toxicological data compiled from the literature
([*J. Comput. Aided Mol. Des.*, **25**, 533–554 (2011)](https://doi.org/10.1007/s10822-011-9440-2)).
Data contributed to OCHEM are released under a Creative Commons CC-BY 4.0 licence.

**OCHEM requires a (free) user account and offers no anonymous programmatic
interface**, so — unlike the ChEMBL, PubChem and Papyrus talktorials — this
notebook cannot download its data while it runs. The exports were therefore
downloaded once through the OCHEM web interface, preprocessed with
`../scripts/prepare_ochem_data.py`, and the small result committed here. This
mirrors how Talktorial T032 ships a pre-computed ClustalO alignment rather than
calling the web service on every run.

### The exported property sets

Each set was exported from the OCHEM *properties browser*
(`https://ochem.eu/properties/show.do`) by selecting the property and exporting
all of its records. Melting point exceeds the size of a single download and was
exported in two parts, which the preprocessing script concatenates and
de-duplicates.

| Property in OCHEM | Raw file | Records | Unit of the converted column | Task in the notebook | Final unit |
|---|---|---|---|---|---|
| Melting Point | `meltingpoint.csv`, `meltingpoint_2.csv` | 403,333 | Celsius | `melting_point` | °C |
| logPow | `logpow.csv` | 155,831 | Log unit | `logp` | log unit |
| Water solubility | `watersolubility.csv` | 53,487 | **mol/L** | `logs` | **log(mol/L)** |
| Boiling Point | `boilingpoint.csv` | 29,056 | Celsius | `boiling_point` | °C |
| Vapor Pressure | `vaporpressure.csv` | 13,977 | **log(Pa)** | `vapour_pressure` | **log(Pa)** |

These endpoints were chosen because they are physically coupled — the General
Solubility Equation ties melting point, logP and aqueous solubility together, and
boiling point and vapour pressure are two views of the same intermolecular
forces — so a multitask model has something real to exploit. On this dataset the
GSE reaches R² = 0.65 (Pearson r = 0.84) on the 5,525 compounds that carry all
three, with no fitting at all, which is both the justification for the task set
and the baseline the solubility model has to beat.

They also span a wide range of difficulty, from logP (easy, nearly additive) to
melting point (hard, governed by crystal packing that molecular descriptors
cannot see).

### Reading an OCHEM export: which value column to use

Every export carries **two** value columns per property:

```
"<Property>"                        + "UNIT {<Property>}"      # as the source reported it
"<Property> {measured, converted}"  + "UNIT {<Property>}.1"    # after OCHEM's conversion
```

The first column's unit is a *mixture* — melting point arrives as `°C`, `K` and
`°F`; vapour pressure as `log(mmHg)`, `log(Pa)`, `mm Hg`, `miliPa`, `Pa` and
`ln(Pa)`; water solubility as `log(mol/L)`, `-log(M)`, `log(mmol/L)`, `mg/L` and
more. Using it without converting would silently mix scales.

The preprocessing script therefore uses **only the `{measured, converted}` column
together with its own `UNIT {...}.1` column**, and asserts that the unit is the
expected one before using the numbers, so that a different export setting fails
loudly instead of changing the results quietly.

OCHEM's conversions were verified against the original values and are correct:
re-deriving the converted number from the source value and unit agrees to
~3 × 10⁻⁵ log units for water solubility (including the `mg/L` records, which
need the molecular weight) and to ~2 × 10⁻⁴ for vapour pressure. The one
exception is Kelvin, which OCHEM converts with 273 rather than 273.15, giving a
constant **0.15 °C** bias on the 4,971 melting- and boiling-point records that
were reported in K — far below the experimental noise on those endpoints, and
left uncorrected.

### Two derived units

Two endpoints do not arrive in the unit the notebook models them in:

* **Water solubility** is converted by OCHEM to *linear* mol/L. The modelling
  literature, and the General Solubility Equation, work in log units, so the
  script takes the base-10 logarithm. Round-tripping the 38,231 records whose
  source unit was already `log(mol/L)` agrees to ~3 × 10⁻⁵ log units, so this
  costs no precision. Two records with a converted value of exactly `0` are
  dropped, since their logarithm is undefined.
* **Vapour pressure** is converted to `log(Pa)` rather than the `log(mmHg)` of
  the OCHEM catalogue listing. The SI unit is kept as-is, so no conversion is
  applied; the notebook's volatility threshold is expressed in `log(Pa)`
  accordingly.

### How the raw exports were processed

`../scripts/prepare_ochem_data.py` performs, per property:

1. Reads every raw export belonging to the property, concatenates them, and drops
   records duplicated where split downloads overlapped (4,201 for melting point).
2. Checks the converted unit, then applies the derived transform above.
3. Standardises every structure with RDKit: largest organic fragment, neutralised,
   isotopes and stereochemistry removed; the InChIKey of the result is the
   compound identity used to join properties together.
4. Discards records that cannot be used: unparsable SMILES, inorganics, molecules
   outside 3–70 heavy atoms, and values outside a physically sensible range for
   the property.
5. Aggregates the surviving literature records per compound into a **median**,
   keeping the number of records aggregated and their standard deviation.
6. Optionally caps how many compounds carrying only a *single* measurement are
   kept, which is what keeps the committed file small. Compounds measured for two
   or more properties are never dropped, so the multitask structure of the data is
   preserved exactly.

Re-running it requires the raw exports in `raw/`:

```bash
python scripts/prepare_ochem_data.py --max-single-task 25000   # what was committed
python scripts/prepare_ochem_data.py --inspect                 # report files, columns and units
```

The `raw/` directory is not committed: the exports total ~110 MB, and
redistributing them wholesale is neither necessary nor courteous to the source.
Anyone with an OCHEM account can reproduce them from the table above.

### What ended up in the committed file

45,276 compounds and 79,654 compound–property pairs:

| Task | Compounds | Unit |
|---|---|---|
| `melting_point` | 37,586 | °C |
| `logp` | 20,685 | log unit |
| `logs` | 9,349 | log(mol/L) |
| `boiling_point` | 8,110 | °C |
| `vapour_pressure` | 3,924 | log(Pa) |

The label matrix is sparse, which is the point: 25,000 compounds carry one
measurement, 11,327 carry two, 5,286 three, 2,173 four and 1,490 all five. The
notebook models the melting point / logP / logS trio, for which 5,525 compounds
are complete.

### Schema of `ochem_properties.csv.gz`

One row per (compound, property) pair — the long format, which represents the
sparsity of the data honestly. The notebook pivots it into the wide label matrix
that a multitask model needs.

| Column | Description |
|---|---|
| `inchikey` | Standard InChIKey of the standardised structure; the compound identity |
| `smiles` | Canonical SMILES of the standardised structure |
| `casrn` | CAS registry number, where the export carried one |
| `property` | Task name (see the tables above) |
| `value` | Aggregated (median) measurement, in the final unit |
| `unit` | Unit of `value` |
| `n_records` | Number of independent literature records aggregated into `value` |
| `value_std` | Standard deviation across those records; empty when `n_records` is 1 |

A note on `value_std`: its median is close to zero for most endpoints, because
literature compilations copy from one another and "several independent records"
is often one measurement propagated several times. The 90th percentile is the
more honest estimate of the experimental noise floor.

### Citation

If you use these data, please cite OCHEM and the original literature sources
recorded with each measurement:

> Sushko, I. *et al.* Online chemical modeling environment (OCHEM): web platform
> for data storage, model development and publishing of chemical information.
> *J. Comput. Aided Mol. Des.* **25**, 533–554 (2011).
