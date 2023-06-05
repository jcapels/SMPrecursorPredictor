import networkx as nx
import re
from sm_precursor_predictor.data_integration.generate_kegg_networks import KeggNetworkGenerator


class KEGGPrecursorFinder:

    def __init__(self, map_id, precursors_data, create_graph=True) -> None:
        self.map_id = map_id
        self.precursors_data = precursors_data
        self.graph = None
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

    def get_all_compounds_of_pathway(self):

        pattern = r'^C\d{5}$'
        allcomp = [node for node in self.graph if re.search(pattern, node)]
        return allcomp

    def check_path_from_compound_to_precursor(self):
        """
        Check if there is a path from a compound to a precursor in the graph.
        """
        compounds = self.get_all_compounds_of_pathway()
        precursors = self.get_precursors_in_pathway()

        for compound in compounds:
            for precursor in precursors:
                try:
                    path = self.find_path_from_source_to_target(compound, precursor)
                    if path:
                        return True
                except nx.NodeNotFound:
                    pass

        return False

    def get_compounds_for_precursors(self):
        """
        Get a dictionary mapping precursors to their associated compounds for all map IDs.
        """
        precursor_compound_dict = {}

        map_ids = self.iloc[:, 0].unique().tolist()

        for map_id in map_ids:

            graph = KeggNetworkGenerator.get_kegg_network(map_id)
            precursors = self.get_precursors_in_pathway()

            for precursor in precursors:
                precursor_compound_dict.setdefault(precursor, [])

            compounds = self.get_all_compounds_of_pathway()

            for precursor in precursors:
                for compound in compounds:
                    try:
                        if compound in graph and precursor in graph:
                            path = self.find_path_from_source_to_target(compound, precursor)
                            if path:
                                if precursor not in precursor_compound_dict:
                                    precursor_compound_dict[precursor] = []
                                precursor_compound_dict[precursor] = list(
                                    set(precursor_compound_dict[precursor]) | {compound})
                    except (nx.NodeNotFound, nx.NetworkXNoPath):
                        pass

        return precursor_compound_dict
