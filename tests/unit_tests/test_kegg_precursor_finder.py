from unittest import TestCase
from sm_precursor_predictor.data_integration.kegg_precursor_finder import KEGGPrecursorFinder
import pandas as pd
import networkx as nx 

class TestKeggPrecursorFinder(TestCase):

    def test_get_precursors_in_pathway(self):
        path = "/home/joao/Downloads/precursors_map_curated_3.csv"
        precursors_data = pd.read_csv(path)

        finder = KEGGPrecursorFinder(map_id="map00902",precursors_data=precursors_data)

        precursors = finder.get_precursors_in_pathway()

        expected_precursors = set(["C00341", "C00342", "C00343"])
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


    def test_check_path_from_compound_to_precursor(self):
        
        graph = nx.DiGraph()
        graph.add_edges_from([
            ("C00341", "C12345"),
            ("C12345", "C54321"),
            ("C54321", "C98765")
        ])

        precursors_data = pd.DataFrame({
            "pathway": ["map00902", "map00902", "map00902"],
            "precursors": ["C00341;C00342", "C12345", "C54321"]
        })
        
        finder = KEGGPrecursorFinder(map_id="map00902", precursors_data=precursors_data, 
                                     create_graph=False)
        finder.graph = graph
        
        finder.check_path_from_compound_to_precursor()      