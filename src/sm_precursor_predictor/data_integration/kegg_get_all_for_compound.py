from sm_precursor_predictor.data_integration.kegg_precursor_finder import KEGGPrecursorFinder
import networkx as nx
from rdkit import Chem


class KEGGAllCompoundInfoExtractor:

    @staticmethod
    def get_compounds_for_precursors(graphs, map_ids, data):
        """
        Get a dictionary mapping precursors to their associated compounds for all map IDs.
        """
        precursor_compound_dict = {}

        for map_id in map_ids:
            kegg_precursors_finder = KEGGPrecursorFinder(map_id, data, create_graph=False)
            precursors = kegg_precursors_finder.get_precursors_in_pathway()

            if len(precursors) == 1:
                graph = graphs.get(map_id)
                kegg_precursors_finder.graph = graph
                compounds = kegg_precursors_finder.get_all_compounds_of_pathway()
                precursor_compound_dict[list(precursors)[0]] = compounds
            else:
                for precursor in precursors:
                    precursor_compound_dict.setdefault(precursor, [])

                graph = graphs.get(map_id)
                kegg_precursors_finder.graph = graph

                if graph is not None:
                    compounds = kegg_precursors_finder.get_all_compounds_of_pathway()

                    for precursor in precursors:
                        for compound in compounds:
                            try:
                                if compound in graph.nodes and precursor in graph.nodes:
                                    path = kegg_precursors_finder.find_path_from_source_to_target(compound, precursor)
                                    if path:
                                        if precursor not in precursor_compound_dict:
                                            precursor_compound_dict[precursor] = []
                                        precursor_compound_dict[precursor] = list(
                                            set(precursor_compound_dict[precursor]) | {compound})
                            except (nx.NodeNotFound, nx.NetworkXNoPath):
                                pass

        return precursor_compound_dict

    @staticmethod
    def get_precursors_for_compound(graphs, map_ids, data):
        """
        Get a dictionary mapping compounds to their associated precursor for all map IDs.
        """
        compound_prec_dict = {}

        for map_id in map_ids:
            kegg_precursors_finder = KEGGPrecursorFinder(map_id, data, create_graph=False)
            kegg_precursors_finder.graph = graphs.get(map_id)
            compounds = kegg_precursors_finder.get_all_compounds_of_pathway()

            precursors = kegg_precursors_finder.get_precursors_in_pathway()

            for compound in compounds:
                if compound not in precursors:
                    compound_prec_dict.setdefault(compound, [])

                    for precursor in precursors:
                        try:
                            if precursor in graphs[map_id] and compound in graphs[map_id]:
                                path = KEGGPrecursorFinder.find_path_from_source_to_target(graphs[map_id], precursor,
                                                                                           compound)
                                if path:
                                    compound_prec_dict[compound] = list(
                                        set(compound_prec_dict[compound]) | {precursor})
                        except (nx.NodeNotFound, nx.NetworkXNoPath):
                            pass

        return compound_prec_dict

    @staticmethod
    def get_all_maps_for_compound(graphs, map_ids, data):
        """
        Get a dictionary mapping compounds to their associated maps for all map IDs.
        """
        compound_map_dict = {}

        for map_id in map_ids:
            kegg_precursors_finder = KEGGPrecursorFinder(map_id, data, create_graph=False)

            kegg_precursors_finder.graph = graphs.get(map_id)

            compounds = kegg_precursors_finder.get_all_compounds_of_pathway()

            precursors = kegg_precursors_finder.get_precursors_in_pathway()

            for compound in compounds:
                if compound not in precursors:
                    if compound not in compound_map_dict:
                        compound_map_dict[compound] = []
                    if map_id not in compound_map_dict[compound]:
                        compound_map_dict[compound].append(map_id)

        return compound_map_dict

    @staticmethod
    def generate_sdf(compound_precursor_dict, compound_map_dict, graphs):

        """
        generates an SDF by mapping compounds to their associated precursors and map IDs.
        """

        writer = Chem.SDWriter('output.sdf')

        for compound, precursors in compound_precursor_dict.items():
            compound_structure = None
            compound_map_ids = compound_map_dict.get(compound, [])

            for map_id, graph in graphs.items():
                if compound in graph.nodes:
                    compound_structure = graph.nodes[compound]["mol"]
                    break

            if compound_structure:

                mol = Chem.MolFromMolBlock(compound_structure)
                mol.SetProp("Precursors", ";".join(precursors))
                mol.SetProp("Map_IDs", ";".join(compound_map_ids))

                writer.write(mol)
            else:
                print(f"Failed to retrieve structure for compound {compound}")

        writer.close()
