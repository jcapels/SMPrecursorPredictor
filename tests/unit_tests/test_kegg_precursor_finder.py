from unittest import TestCase
import pandas as pd
import networkx as nx
from sm_precursor_predictor.data_integration.kegg_precursor_finder import KEGGPrecursorFinder


class TestKeggPrecursorFinder(TestCase):

    def test_get_precursors_in_pathway(self):
        path = "/home/joao/Downloads/precursors_map_curated_3.csv"
        precursors_data = pd.read_csv(path)

        finder = KEGGPrecursorFinder(map_id="map00940", precursors_data=precursors_data, create_graph=False)

        precursors = finder.get_precursors_in_pathway()

        expected_precursors = set(["C00079", "C00082", "C00315"])
        self.assertEqual(precursors, expected_precursors)

    def test_find_path_from_source_to_target(self):
        graph = nx.DiGraph()
        graph.add_edges_from([
            ("C00341", "C12345"),
            ("C12345", "C54321"),
            ("C54321", "C98765")
        ])

        finder = KEGGPrecursorFinder(map_id="map00902", precursors_data=None)
        finder.graph = graph

        source_compound = "C00341"
        target_compound = "C98765"

        paths = finder.find_path_from_source_to_target(source_compound, target_compound)

        expected_paths = ["C00341", "C12345", "C54321", "C98765"]
        self.assertEqual(paths, expected_paths)

    def test_check_path_from_compound_to_precursor2(self):
        graph = nx.DiGraph()
        graph.add_edges_from([
            ("C00341", "C12345"),
            ("C12345", "C54321"),
            ("C54321", "C98765")
        ])

        precursors_data = pd.DataFrame({
            "pathway": ["map00902", "map00902"],
            "precursors": ["C98765;C54321", "C54321"]
        })

        finder = KEGGPrecursorFinder(map_id="map00902", precursors_data=precursors_data, create_graph=False)
        finder.graph = graph

        result = finder.check_path_from_compound_to_precursor()

        self.assertTrue(result, "Path from compound to precursor check failed")

    def test_get_compounds_for_precursors(self):
        precursors_data = pd.DataFrame({
            "pathway": ["map00902", "map00909"],
            "precursors": ["C00341", "C00448"]
        })

        expected_mapping = {
            'C00341': ['C03190', 'C03161', 'C06074', 'C11673', 'C11636', 'C00808', 'C01767', 'C00843', 'C20944',
                       'C06099', 'C04433', 'C18025', 'C21203', 'C22225', 'C22014', 'C01123', 'C02869', 'C09804',
                       'C00553', 'C09893', 'C20789', 'C06308', 'C00521', 'C11409', 'C11951', 'C01500', 'C06070',
                       'C20943', 'C02485', 'C03092', 'C11388', 'C00341', 'C18027', 'C17621', 'C00080', 'C00964',
                       'C02452', 'C01433', 'C02462', 'C02745', 'C06305', 'C17622', 'C01957', 'C06066', 'C11383',
                       'C03024', 'C00400', 'C11382', 'C01512', 'C02344', 'C09782', 'C06071', 'C09844', 'C11389',
                       'C20790', 'C04718', 'C01852', 'C06307', 'C11952', 'C11672', 'C06304', 'C20221', 'C11393',
                       'C01765', 'C09769', 'C00848'],
            'C00448': ['C09665', 'C20200', 'C20193', 'C01841', 'C19829', 'C03161', 'C22150', 'C00751', 'C09704',
                       'C16775', 'C20181', 'C19832', 'C19749', 'C19819', 'C03428', 'C19708', 'C09672', 'C19752',
                       'C07667', 'C09382', 'C16141', 'C17277', 'C19746', 'C20195', 'C09723', 'C06394', 'C19834',
                       'C19833', 'C20479', 'C20478', 'C01860', 'C19732', 'C19828', 'C19748', 'C19973', 'C19742',
                       'C08637', 'C00207', 'C06080', 'C16814', 'C17955', 'C19753', 'C19676', 'C18245', 'C00139',
                       'C01126', 'C20187', 'C16829', 'C19678', 'C09737', 'C19918', 'C19711', 'C19908', 'C03461',
                       'C00080', 'C19747', 'C03220', 'C12142', 'C06309', 'C16269', 'C19820', 'C00138', 'C22151',
                       'C16286', 'C06079', 'C09666', 'C17954', 'C03024', 'C17953', 'C16142', 'C19939', 'C20724',
                       'C16028', 'C20191', 'C08616', 'C06083', 'C19835', 'C20192', 'C09699', 'C09629', 'C19801',
                       'C02004', 'C01054', 'C19755', 'C19974', 'C20159', 'C16143', 'C09627', 'C20163', 'C00448',
                       'C08615', 'C08628', 'C09684', 'C19863', 'C06310', 'C19740']
        }

        finder = KEGGPrecursorFinder(map_id="map00942", precursors_data=precursors_data)
        precursor_compound_dict = finder.get_compounds_for_precursors()

        self.assertEqual(precursor_compound_dict, expected_mapping)
