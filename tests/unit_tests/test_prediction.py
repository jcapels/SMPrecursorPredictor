import os
from unittest import TestCase

from sm_precursor_predictor import predict_precursors, predict_from_csv, get_prediction_and_explanation
from tests import TEST_DIR


class TestPrediction(TestCase):

    def test_prediction(self):
        precursors = predict_precursors(
            ["[H][C@]89CN(CCc1c([nH]c2ccccc12)[C@@](C(=O)OC)(c3cc4c(cc3OC)N(C)[C@@]5([H])[C@@]"
             "(O)(C(=O)OC)[C@H](OC(C)=O)[C@]7(CC)C=CCN6CC[C@]45[C@@]67[H])C8)C[C@](O)(CC)C9",
             "COC1=C(C=CC(=C1)C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O[C@H]4[C@@H]([C@H]([C@H]([C@H](O4)CO)O)O)O)O"])

        self.assertTrue("Tryptophan" in precursors[0])
        self.assertTrue("Secologanin" in precursors[0])

        self.assertTrue("L-Phenylalanine;Kaempferol" in precursors[1])

    def test_from_dataset(self):
        predictions = predict_from_csv(os.path.join(TEST_DIR, "data", "data.csv"), "SMILES", "ID")

        self.assertTrue("Cholesterol" in predictions.loc[0, "predicted_precursors"])

    def test_from_dataset_morgan_fp(self):

        predictions = predict_from_csv(os.path.join(TEST_DIR, "data", "data.csv"), "SMILES", "ID", 
                                       model = "Morgan FP + Ridge Classifier")
        self.assertTrue("Cholesterol" in predictions.loc[0, "predicted_precursors"])

    def test_explainability(self):

        get_prediction_and_explanation(smiles="COC1=C(C=CC(=C1)C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O[C@H]4[C@@H]([C@H]([C@H]([C@H](O4)CO)O)O)O)O")

