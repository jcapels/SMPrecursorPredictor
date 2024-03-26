import os

import numpy as np
import pandas as pd
from deepmol.datasets import SmilesDataset
from deepmol.loaders import CSVLoader
from deepmol.pipeline import Pipeline


models = {

    "Layered FP + Low Variance FS + Ridge Classifier": "ridge_classifier_layered_fingerprints_variance_selector",
    "Morgan FP + Ridge Classifier": "ridge_classifier_morgan_fp"
}

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

kegg_labels = ['C00073', 'C00078', 'C00079', 'C00082', 'C00235', 'C00341', 'C00353',
              'C00448', 'C01789', 'C03506', 'C00047', 'C00108', 'C00187', 'C00148',
              'C00041', 'C00129', 'C00062', 'C01852', 'C00049', 'C00135', 'C00223',
              'C00509', 'C00540', 'C01477', 'C05903', 'C05904', 'C05905', 'C05908',
              'C09762']


def convert_predictions_into_names_model(predictions):

    labels_names = np.array([labels_[label] for label in kegg_labels])
    ones = predictions == 1
    labels_all = []
    for i, prediction in enumerate(ones):
        labels_all.append(";".join(labels_names[prediction]))
    return labels_all

import os
from rdkit.Chem import rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Draw

def _draw_bits(on_bits, mol, molecule):
    
    info = {}
    #best_pipeline.steps[0][1].draw_bit(grape_vine_dataset.mols[molecule], on_bits[2], file_path=f"bit_{bit}.png")
    rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2,
                                                   1024, bitInfo=info)
    
    atomsToUse = []
    atoms_ids = []
    highlightAtomColors = {}
    for bit in on_bits:
        aid, rad = info[bit][0]
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, 2, aid)
        atoms_ids.append(aid)
        highlightAtomColors[aid] = (0.3, 0.3, 1)
        for b in env:
            atomsToUse.append(mol.GetBondWithIdx(b).GetBeginAtomIdx())
            atomsToUse.append(mol.GetBondWithIdx(b).GetEndAtomIdx())
    atomsToUse = list(set(atomsToUse))
    drawer = Draw.rdMolDraw2D.MolDraw2DCairo(1000, 700)
    drawer.drawOptions().addAtomIndices = True
    drawer.drawOptions().annotationFontScale = 1
    drawer.drawOptions().setHighlightColour(((0.0, 0.0, 1.0, 0.1)))
    drawer.DrawMolecule(mol, highlightAtoms=atomsToUse, highlightAtomColors=highlightAtomColors)
    drawer.FinishDrawing()
    drawer.WriteDrawingText(f"molecule_{molecule}.png")   
    return drawer.GetDrawingText()

import re
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def draw_feature_importance_plot(coefs, on_bits, title="", threshold=0.15):
    coefs = pd.Series(coefs)
    bits = []
    coefs_to_plot = []
    for i, coef in enumerate(coefs):
        if coef >= threshold or coef <= -threshold or i in on_bits:
            bits.append("bit_{}".format(i))
            coefs_to_plot.append(coef)
    figure = plt.figure(figsize=(5, 5))
    bar_plot = sns.barplot(x=coefs_to_plot, y=bits, orient='h')
    plt.xticks(rotation=90)
    # paint bars above 0 with green and below 0 with red
    for i, coef in enumerate(coefs_to_plot):
        if int(re.search("\d+", bits[i]).group()) in on_bits:
            plt.barh(i, coef, color='blue')
        elif coef > 0:
            plt.barh(i, coef, color='green')
        else:
            plt.barh(i, coef, color='red')
    
    plt.title(f"Precursor: {title}", pad =20, fontsize=15)
    plt.xlabel('Ridge classifier coefficient', fontsize=15)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    return figure

from PIL import Image
import io

def draw_important_bits(smiles, predicted_labels, threshold=0.15, molecule_name="molecule"):

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    pretrain_model_path = os.path.join(BASE_DIR,
                                       "prediction",
                                       models["Morgan FP + Ridge Classifier"])
    pipeline_morgan_fp = Pipeline.load(pretrain_model_path)
    coefs = pipeline_morgan_fp.steps[-1][1].model.coef_

    predicted_labels = predicted_labels.split(";")
    for i, label in enumerate(predicted_labels):
        for _label in labels_:
            if labels_[_label] == label:
                predicted_labels[i] = _label

    numbers_of_predicted_label = [kegg_labels.index(label) for label in predicted_labels]

    images = []
    plots = []
    for j, number_predicted_label in enumerate(numbers_of_predicted_label):
        coefs_ = coefs[number_predicted_label,:]
        # get on bits from the morgan fingerprints
        dataset = SmilesDataset(smiles=[smiles])
        dataset = pipeline_morgan_fp.transform(dataset)
        bits = dataset.X[0, :]
        on_bits = []
        coefs_list = []
        for i, bit in enumerate(bits):
            if bit == 1 and coefs_[i] > threshold:
                on_bits.append(i)
                coefs_list.append(coefs_[i])
        
        image = _draw_bits(on_bits, dataset.mols[0], molecule_name + "_" + str(labels_[predicted_labels[j]]))
        image = Image.open(io.BytesIO(image))
        images.append(image)
        plots.append(draw_feature_importance_plot(coefs_, on_bits, title=str(labels_[predicted_labels[j]]), threshold=threshold))
    
    return images, plots

def predict_from_dataset(dataset, model):
    if model not in models:
        raise ValueError(f"Model {model} not found. Available models are: {models.keys()}")
    else:
        model = models[model]

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    pretrain_model_path = os.path.join(BASE_DIR,
                                       "prediction",
                                       model)
    best_pipeline = Pipeline.load(pretrain_model_path)

    predictions = best_pipeline.predict(dataset)
    predictions = convert_predictions_into_names_model(predictions)

    return predictions

def get_prediction_and_explanation(smiles):
    dataset = SmilesDataset(smiles=[smiles])
    model="Morgan FP + Ridge Classifier"
    predictions = predict_from_dataset(dataset, model)

    images, plots = draw_important_bits(smiles, predictions[0])

    return predictions, images, plots


def predict_precursors(smiles: list, model: str="Layered FP + Low Variance FS + Ridge Classifier"):
    """Predicts the precursor of a given SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.
        model (str): Model to use for prediction.

    Returns:
        list: List of SMILES strings of the predicted precursors.
    """

    dataset = SmilesDataset(smiles=smiles)
    return predict_from_dataset(dataset, model)


def predict_from_csv(csv_path, smiles_field, ids_field=None, model = "Layered FP + Low Variance FS + Ridge Classifier", 
                     **kwargs):
    """Predicts the precursor of a given SMILES string.

    Args:
        csv_path (str): Path to the csv file.
        smiles_field (str): Name of the column containing the SMILES strings.
        ids_field (str): Name of the column containing the IDs.

    Returns:
        list: List of SMILES strings of the predicted precursors.
    """
    dataset = CSVLoader(csv_path, smiles_field, ids_field).create_dataset(**kwargs)
    predictions = predict_from_dataset(dataset, model)
    dataset = pd.read_csv(csv_path, **kwargs)
    dataset["predicted_precursors"] = predictions
    return dataset
