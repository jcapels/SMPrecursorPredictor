import os

import numpy as np
import pandas as pd
from deepmol.datasets import SmilesDataset
from deepmol.loaders import CSVLoader
from deepmol.pipeline import Pipeline


def convert_predictions_into_names_model_1(predictions):
    import numpy as np
    labels = np.array(['C00073', 'C00078', 'C00079', 'C00082', 'C00235', 'C00341',
                       'C00353', 'C00448', 'C01789', 'C03506',
                       'C00047', 'C00108', 'C00187', 'C00148', 'C00041',
                       'C00129', 'C00062', 'C01852', 'C00049', 'C00135'])

    labels_ = {
        'C00341': 'Geranyl diphosphate',
        'C01789': 'Campesterol',
        'C00078': 'Tryptophan',
        'C00049': 'L-Aspartate',
        'C00183': 'L-Valine',
        'C03506': 'Indoleglycerol phosphate',
        'C00187': 'Cholesterol',
        'C00079': 'L-Phenylalanine',
        'C00047': 'L-Lysine',
        'C01852': 'Secologanin',
        'C00407': 'L-Isoleucine',
        'C00129': 'Isopentenyl diphosphate',
        'C00235': 'Dimethylallyl diphosphate',
        'C00062': 'L-Arginine',
        'C00353': 'Geranylgeranyl diphosphate',
        'C00148': 'L-Proline',
        'C00073': 'L-Methionine',
        'C00108': 'Anthranilate',
        'C00123': 'L-Leucine',
        'C00135': 'L-Histidine',
        'C00448': 'Farnesyl diphosphate',
        'C00082': 'L-Tyrosine',
        'C00041': 'L-Alanine'
    }

    labels_names = np.array([labels_[label] for label in labels])
    ones = predictions == 1
    labels_all = []
    for i, prediction in enumerate(ones):
        labels_all.append(";".join(labels_names[prediction]))
    return labels_all


def convert_predictions_into_names_model_2(predictions):
    labels_ = {
        'C00341': 'Geranyl diphosphate',
        'C01789': 'Campesterol',
        'C00078': 'Tryptophan',
        'C00049': 'L-Aspartate',
        'C00183': 'L-Valine',
        'C03506': 'Indoleglycerol phosphate',
        'C00187': 'Cholesterol',
        'C00079': 'L-Phenylalanine',
        'C00047': 'L-Lysine',
        'C01852': 'Secologanin',
        'C00407': 'L-Isoleucine',
        'C00129': 'Isopentenyl diphosphate',
        'C00235': 'Dimethylallyl diphosphate',
        'C00062': 'L-Arginine',
        'C00353': 'Geranylgeranyl diphosphate',
        'C00148': 'L-Proline',
        'C00073': 'L-Methionine',
        'C00108': 'Anthranilate',
        'C00123': 'L-Leucine',
        'C00135': 'L-Histidine',
        'C00448': 'Farnesyl diphosphate',
        'C00082': 'L-Tyrosine',
        'C00041': 'L-Alanine',
        'C00540': 'Cinnamoyl-CoA',
        'C01477': 'Apigenin',
        'C05903': 'Kaempferol',
        'C05904': 'Pelargonin',
        'C05905': 'Cyanidin',
        'C05908': 'Delphinidin',
        'C00389': 'Quercetin',
        'C01514': 'Luteolin',
        'C09762': "Liquiritigenin",
        'C00509': 'Naringenin',
        'C00223': 'p-Coumaroyl-CoA'
    }

    labels = ['C00073', 'C00078', 'C00079', 'C00082', 'C00235', 'C00341', 'C00353',
              'C00448', 'C01789', 'C03506', 'C00047', 'C00108', 'C00187', 'C00148',
              'C00041', 'C00129', 'C00062', 'C01852', 'C00049', 'C00135', 'C00223',
              'C00509', 'C00540', 'C01477', 'C05903', 'C05904', 'C05905', 'C05908',
              'C09762']
    labels_names = np.array([labels_[label] for label in labels])
    ones = predictions == 1
    labels_all = []
    for i, prediction in enumerate(ones):
        labels_all.append(";".join(labels_names[prediction]))
    return labels_all


def predict_from_dataset(dataset):
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    pretrain_model_path = os.path.join(BASE_DIR,
                                       "prediction",
                                       "sm_predictor_pipeline_tpe",
                                       "trial_289")
    best_pipeline = Pipeline.load(pretrain_model_path)

    additional_model = os.path.join(BASE_DIR,
                                    "prediction",
                                    "sm_predictor_pipeline_motpe",
                                    "trial_311")
    additional_model = Pipeline.load(additional_model)

    predictions = best_pipeline.predict(dataset)
    predictions = convert_predictions_into_names_model_1(predictions)

    additional_predictions = additional_model.predict(dataset)
    additional_predictions = convert_predictions_into_names_model_2(additional_predictions)

    for i, prediction in enumerate(predictions):
        if "l-phenylalanine" in prediction.lower():
            predictions[i] = additional_predictions[i]
    return predictions


def predict_precursors(smiles: list):
    """Predicts the precursor of a given SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        list: List of SMILES strings of the predicted precursors.
    """

    dataset = SmilesDataset(smiles=smiles)
    return predict_from_dataset(dataset)


def predict_from_csv(csv_path, smiles_field, ids_field=None, **kwargs):
    """Predicts the precursor of a given SMILES string.

    Args:
        csv_path (str): Path to the csv file.
        smiles_field (str): Name of the column containing the SMILES strings.
        ids_field (str): Name of the column containing the IDs.

    Returns:
        list: List of SMILES strings of the predicted precursors.
    """
    dataset = CSVLoader(csv_path, smiles_field, ids_field).create_dataset(**kwargs)
    predictions = predict_from_dataset(dataset)
    dataset = pd.read_csv(csv_path, **kwargs)
    dataset["predicted_precursors"] = predictions
    return dataset
