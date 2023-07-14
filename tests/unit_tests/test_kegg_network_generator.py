from unittest import TestCase, skip

from sm_precursor_predictor.data_integration.generate_kegg_networks import KeggNetworkGenerator


class TestKeggNetworkGenerator(TestCase):

    @skip("Not implemented")
    def test_kegg_network_generator(self):
        KeggNetworkGenerator.get_kegg_network("map00960")