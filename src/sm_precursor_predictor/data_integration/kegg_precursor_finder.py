import networkx as nx
import re
from sm_precursor_predictor.data_integration.generate_kegg_networks import KeggNetworkGenerator
from sm_precursor_predictor.data_integration.kegg_api import KeggApi

class KEGGPrecursorFinder:

    def __init__(self, map_id, precursors_data, create_graph=True) -> None:
        self.map_id = map_id
        self.precursors_data = precursors_data
        self.graph=None
        if create_graph:
            self.graph = KeggNetworkGenerator.get_kegg_network(self.map_id)
        
    
    def get_precursors_in_pathway(self):
        """
        Get all precursors in a pathway based on the map ID.
        """
        precur = set()
        pathway_data = self.precursors_data[self.precursors_data["pathway"] == self.map_id]
        for _, row in pathway_data.iterrows():
            precursors = row["precursors"].split(";")
            precur.update(precursors)
        return precur
    
    def find_path_from_source_to_target(self, source_compound, target_compound):
        """
        Find a path from the source compound to the target compound in the graph.
        """
        paths = nx.shortest_path(self.graph, source_compound, target_compound)
        
        return paths
    
    
    def get_allcompounds_of_pathway(self):
    
        pattern = r'^C\d{5}$'  
        allcomp = [node for node in self.graph if re.search(pattern, node)]
        return allcomp
    
    def check_path_from_compound_to_precursor(self):
        """
        Check if there is a path from a compound to a precursor in the graph.
        """
        compounds = self.get_allcompounds_of_pathway()
        precursors = self.get_precursors_in_pathway()
        
        path_found = False

        for compound in compounds:
            for precursor in precursors:
                try:
                    path = KEGGPrecursorFinder.find_path_from_source_to_target(self.graph, compound, precursor)
                    if path:
                        path_str = " -> ".join(path)
                        print(f"Path found: {path_str}")
                        path_found = True
                except nx.NodeNotFound:
                    print(f"No path found from compound {compound} to precursor {precursor}.")

        if not path_found:
            print("No path found between any compound and precursor.")
