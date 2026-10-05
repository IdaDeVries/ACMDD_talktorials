# T033 · Molecular properties: multitask modelling

**Note:** This talktorial is a part of TeachOpenCADD, a platform that aims to teach domain-specific skills and to provide pipeline templates as starting points for research projects.

Authors:

- Adapted for ACMDD (Advanced Computational Methods in Drug Discovery), MSc course, Leiden University, 2026.

This talktorial is the direct continuation of **Talktorial T032** by Marina Gorostiola González, Olivier J. M. Béquignon and Willem Jespers (Computational Drug Discovery, Drug Discovery & Safety Leiden University, The Netherlands). It reuses T032's modelling scaffolding and replaces its endpoint.


## Aim of this talktorial

Talktorial T032 covered several protein targets with one model by *describing* the target: proteochemometrics puts the task on the **input** side, as a block of protein descriptors. This talktorial covers several prediction problems by putting them on the **output** side instead — **multitask modelling** — and swaps the endpoint from binding affinity to the physicochemical properties of the molecules themselves, taken from the [OCHEM](https://ochem.eu) database.

Along the way it works through the three ways of obtaining several outputs from a scikit-learn-style model — estimators with native multi-output support, the `MultiOutputRegressor` / `MultiOutputClassifier` wrappers, and XGBoost's `multi_strategy='multi_output_tree'` — and it takes data splitting seriously, comparing every model under four different splits with repeated seeds.


### Contents in *Theory*

* From proteochemometrics to multitask modelling
* Data preparation
    * The OCHEM database
    * The properties we model, and why they belong together
    * Molecule encoding: molecular descriptors
* Multitask machine learning
    * Why sharing a model between tasks can help — and when it hurts
    * Three implementations of multi-output learning
    * Missing labels: the sparse label matrix
* Data splitting methods
    * Random, scaffold, cluster and stratified splits
    * Model selection vs. model assessment
    * Repeating a split: how much of a difference is a real difference?
* Evaluation metrics for multi-output models


### Contents in *Practical*

* Load the OCHEM property data
* Data preparation
    * Exploring the endpoints
    * From long to wide: the label matrix
    * Task overlap and task correlation
    * Calculate compound descriptors
    * Building the feature matrix
    * How many tasks can we afford?
* Data splitting
    * Grouping compounds: scaffolds and similarity clusters
    * Looking at what each split actually did
    * Quantifying the leakage: nearest-neighbour similarity between test and training set
* Multitask regression
    * Single-task baselines, native multi-output models, wrapped models, XGBoost
    * The full comparison: models × splits × repeats
    * Cross-validation that matches the split
    * A learning-curve experiment: does multitask learning actually help?
    * Sharing structure, or keeping the data?
* Multitask classification
    * `MultiOutputClassifier`, native multi-label forests, and XGBoost
    * Metrics for multi-label problems
* Discussion
* Quiz


### References

* OCHEM database: [*J. Comput. Aided Mol. Des.*, **25**, 533–554 (2011)](https://doi.org/10.1007/s10822-011-9440-2), [ochem.eu](https://ochem.eu)
* Multitask learning: Caruana, [*Machine Learning*, **28**, 41–75 (1997)](https://doi.org/10.1023/A:1007379606734)
* Multitask networks for QSAR: [*J. Chem. Inf. Model.*, **55**, 263–274 (2015)](https://doi.org/10.1021/ci500747n)
* Multi-output decision trees: [*Machine Learning*, **73**, 185–214 (2008)](https://doi.org/10.1007/s10994-008-5077-3)
* scikit-learn, multiclass and multioutput algorithms [(Docs)](https://scikit-learn.org/stable/modules/multiclass.html)
* XGBoost, multiple outputs [(Docs)](https://xgboost.readthedocs.io/en/stable/tutorials/multioutput.html)
* Bemis–Murcko scaffolds: [*J. Med. Chem.*, **39**, 2887–2893 (1996)](https://doi.org/10.1021/jm9602928)
* On splitting and over-optimism in QSAR: [*J. Chem. Inf. Model.*, **58**, 916–932 (2018)](https://doi.org/10.1021/acs.jcim.7b00403)
* The General Solubility Equation: [*J. Pharm. Sci.*, **69**, 912–922 (1980)](https://doi.org/10.1002/jps.2600690814)
* Evaluation metrics [(scikit-learn Docs)](https://scikit-learn.org/stable/modules/model_evaluation.html)
