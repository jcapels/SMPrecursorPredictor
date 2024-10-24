from unittest import TestCase, skip

import pandas as pd

from sm_precursor_predictor.data_integration.kegg_api import KeggApi

@skip
class TestKeggApi(TestCase):

    def test_get_list(self):
        list_of_pathways = KeggApi.get_list("pathway")
        self.assertEqual(list_of_pathways.shape, (571, 2))

    def test_get_info(self):
        list_of_pathways = KeggApi.get_list("brite")
        self.assertEqual(list_of_pathways.shape, (146, 2))
        self.assertIn("br08011", list_of_pathways[0].values)

    def test_get_metabolites_in_pathway(self):
        metabolites_in_pathway = KeggApi.get_links("compound", "path:map00902")
        self.assertEqual(metabolites_in_pathway.shape, (61, 2))

    def test_get_metabolite(self):
        metabolite = KeggApi.get("C00200")
        self.assertIn("C00200", metabolite)

    @skip
    def test_get_secondary_metabolites(self):
        secondary_metabolites = pd.read_html("https://www.genome.jp/kegg/compound/br08011.html")
        self.assertEqual(secondary_metabolites[0].shape, (154, 6))
