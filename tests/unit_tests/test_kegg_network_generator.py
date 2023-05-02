from unittest import TestCase

from sm_precursor_predictor.data_integration.generate_kegg_networks import KeggNetworkGenerator


class TestKeggNetworkGenerator(TestCase):

    def test_kegg_network_generator(self):
        KeggNetworkGenerator.get_kegg_network("map00902")