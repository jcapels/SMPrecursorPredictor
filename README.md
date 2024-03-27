# SMPrecursorPredictor
A ML pipeline for the prediction of specialised metabolites starting substances.

### Table of contents:

- [Installation](#installation)
    - [Manually](#Manually)
    - [Pypi](#pypi)
- [Making predictions](#making-predictions)
- [Methods](#methods)
    - [Data](#data)
    - [AutoML](#automl)
    - [Analysis of the results](#analysis-of-the-results)
- [License](#licensing)

## Installation

### Manually

1. Clone the repository and move into the directory:

```bash
git clone
cd SMPrecursorPredictor
```

2. Create a conda environment and activate it:

```bash
conda create -n sm_precursor_predictor python=3.10
conda activate sm_precursor_predictor
```

3. Install the dependencies:

```bash
pip install -r requirements.txt
```

4. Install the package:

```bash
pip install .
```

### Pypi

1. Create a conda environment and activate it:

```bash
conda create -n sm_precursor_predictor python=3.10
conda activate sm_precursor_predictor
pip install SMPrecursorPrediction
```

## Making predictions

```python
from sm_precursor_predictor import predict_precursors
precursors = predict_precursors(
            ["[H][C@]89CN(CCc1c([nH]c2ccccc12)[C@@](C(=O)OC)(c3cc4c(cc3OC)N(C)[C@@]5([H])[C@@]"
             "(O)(C(=O)OC)[C@H](OC(C)=O)[C@]7(CC)C=CCN6CC[C@]45[C@@]67[H])C8)C[C@](O)(CC)C9",
             "COC1=C(C=CC(=C1)C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O[C@H]4[C@@H]([C@H]([C@H]([C@H](O4)CO)O)O)O)O"])
print(precursors)
```

or

read a csv file with a column of SMILES and a column of IDs and save the predictions in a csv file:

```python
from sm_precursor_predictor import predict_from_csv
predictions = predict_from_csv("path_to_csv", 
                               smiles_field="SMILES", 
                               ids_field="ID")
predictions.to_csv("path_to_save_predictions.csv")
```

## Making and explaining predictions

This is only possible with one model: *Morgan FP + Ridge Classifier*.

Example with linalool:

```python
from sm_precursor_predictor import get_prediction_and_explanation

prediction, images, plots = get_prediction_and_explanation(smiles="CC(=CCCC(C)(C=C)O)C", threshold=0.20)
```
![feature_importance](feature_importance.png)

```
prediction
```


```
['Geranyl diphosphate']
```

```
images[0]
```
![Linalool](molecule_Geranyl_diphosphate.png)



## Methods

### Data

The final dataset can be found in **[final_dataset.csv]**(models_and_datasets/final_dataset/final_dataset.csv). The LotusDB compounds predictions can be found at **[predictions]**(models_and_datasets/predictions_lotusdb/LOTUSDB_predictions_plants.csv).

The exploration of the dataset can be found at **[dataset_analysis.ipynb]**(models_and_datasets/analysis/dataset_analysis.ipynb).

### AutoML 

The AutoML was run using docker. To run using docker you must consider the following files:

- **[Dockerfile](models_and_datasets/analysis/results_for_new_dataset/Dockerfile)**
- **[run.sh](models_and_datasets/analysis/results_for_new_dataset/run.sh)**

Alternatively, if you're rather interested in running the AutoML with a python script, consider the following:

- **[train_models.py](models_and_datasets/analysis/results_for_new_dataset/train_models.py)**

### Analysis of the results

For the analysis of the results refer to the following files:

- **[Main results](models_and_datasets/analysis/results_for_new_dataset/results_analysis.ipynb)**
- **[Alkaloids dataset - Eguchi et al. 2019](models_and_datasets/analysis/results_for_alkaloids_dataset/pipeline_for_alkaloids_data.ipynb)**
- **[Challenging datasets](models_and_datasets/analysis/results_for_challenging_datasets/assess_for_challenging_datasets.ipynb)**
- For checking the model interpretability - **[Monoterpenoid indole alkaloids](models_and_datasets/analysis/results_for_challenging_compounds/alkaloids_derived_from_tryptophan.ipynb)** and **[others](models_and_datasets/analysis/case_studies/predict_for_case_studies.ipynb)**.










