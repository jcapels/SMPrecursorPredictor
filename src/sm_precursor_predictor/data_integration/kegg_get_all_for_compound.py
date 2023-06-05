from sm_precursor_predictor.data_integration.generate_kegg_networks import KeggNetworkGenerator
from sm_precursor_predictor.data_integration.kegg_precursor_finder import KEGGPrecursorFinder
import networkx as nx

class KEGG_all_compound_info:

    
    @staticmethod
    def get_precursors_for_compound(graph, map_ids, data):

        """
        Get a dictionary mapping compounds to their associated precursor for all map IDs.
        """
        compound_prec_dict = {}
        map_ids = data.iloc[:, 0].unique().tolist()

        for map_id in map_ids:
            graph = KeggNetworkGenerator.get_kegg_network(map_id)
            compounds = KEGGPrecursorFinder.get_allcompounds_of_pathway(graph)

            precursors = KEGGPrecursorFinder.get_precursors_in_pathway(map_id, data)

            for compound in compounds:
                if compound not in precursors:
                    compound_prec_dict.setdefault(compound, [])

                    for precursor in precursors:
                        try:
                            if precursor in graph and compound in graph:
                                path = KEGGPrecursorFinder.find_path_from_source_to_target(graph, precursor, compound)
                                if path:
                                    compound_prec_dict[compound] = list(set(compound_prec_dict[compound]) | set([precursor]))
                        except (nx.NodeNotFound, nx.NetworkXNoPath):
                            pass

        return compound_prec_dict

    @staticmethod
    def get_all_maps_for_compound(graph, map_ids, data):

        """
        Get a dictionary mapping compounds to their associated maps for all map IDs.
        """
        compound_map_dict = {}
        map_ids = data.iloc[:, 0].unique().tolist()

        for map_id in map_ids:
            graph = KeggNetworkGenerator.get_kegg_network(map_id)
            compounds = KEGGPrecursorFinder.get_allcompounds_of_pathway(graph)

            precursors = KEGGPrecursorFinder.get_precursors_in_pathway(map_id, data)

            for compound in compounds:
                if compound not in precursors:
                    if compound not in compound_map_dict:
                        compound_map_dict[compound] = []
                    if map_id not in compound_map_dict[compound]:
                        compound_map_dict[compound].append(map_id)

        return compound_map_dict