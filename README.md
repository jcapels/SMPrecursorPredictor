# SMPrecursorPredictor
A ML pipeline for the prediction of specialised metabolites precursors.

### Table of contents:

- [Installation](#installation)
    - [Manually](#manually)
- [Methods](#methods)
    - [Problem setup](#problem-setup)
- [License](#licensing)

## Installation

### Manually


## Methods

### Problem setup

Data integration:
- Alkaloids: from Eguchi et al 2019;
- Terpenoids, phenols and gluconates: curated data from KEGG

Multi-label classification problem: 

![molecular_starters.png](imgs/molecular_starters.png)

### Data integration results

![data_integration.png](imgs/data_integration.png)

### Splitted dataset

![split_results_2.png](imgs%2Fsplit_results_2.png)

![split_results.png](imgs%2Fsplit_results.png)

### Model results

With a few lines of code we tested:

- 5 molecular fingerprints alone and combinations;
- 3 different standardization methods
- 7 models from sklearn for multilabel classification.
- 6 different optimization methods: an evolutionary algorithm, MOTPE, random search, TPE, CMAES and quasi monte carlo;
- 
In total, we tested 3000 combinations, 500 for each method of optimization;

![models_results.png](imgs%2Fmodels_results.png)

### Best trials stats - best fingerprints

![best_fingerprints.png](imgs%2Fbest_fingerprints.png)

### Best trials stats - best models

![best_models.png](imgs%2Fbest_models.png)

### Best pipelines

- Standardizer: ChEMBLStandardizer
- Fingerprints: 
- - Layered fingerprints (size: 2048, minimum path: 3, maximum path: 8)
- - AtomPair fingerprints (size: 2048, minimum length: 2, maximum length: 40, does not include information on chirality)
- Model: 
- - RidgeClassifier(alpha=7.338782054460601, fit_intercept=False,solver='sparse_cg')






