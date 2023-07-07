from sm_precursor_predictor.data_integration.generate_kegg_networks import KeggNetworkGenerator
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
                precursor_compound_dict[list(precursors)[0]] = compounds[1:]

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

            if len(precursors) == 1:
                graph = graphs.get(map_id)
                kegg_precursors_finder.graph = graph

                for compound in compounds:
                    compound_prec_dict[compound] = [list(precursors)[0]]
                    compound_prec_dict = {compound: precursors for compound, precursors in compound_prec_dict.items() if
                                          compound not in precursors}


            else:
                for compound in compounds:
                    if compound not in precursors:
                        compound_prec_dict.setdefault(compound, [])

                        for precursor in precursors:
                            try:
                                if precursor in graphs[map_id] and compound in graphs[map_id]:
                                    path = kegg_precursors_finder.find_path_from_source_to_target(compound, precursor)
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
    def generate_sdf(compound_prec_dict, compound_map_dict, graphs):
        writer = Chem.SDWriter('phenolics_and_terpenoids.sdf')

        precursor_ids = set(precursor for precursors in compound_prec_dict.values() for precursor in precursors)
        additional_precursor_ids = {'C00148', 'C00108', 'C00047', 'C00062', 'C00041', 'C00049', 'C01852', 'C00129',
                                    'C00135', 'C00187'}
        precursor_ids.update(additional_precursor_ids)

        for compound, precursors in compound_prec_dict.items():
            compound_structure = None

            for map_id, graph in graphs.items():
                if compound in graph.nodes:
                    compound_structure = graph.nodes[compound]["mol"]
                    break

            if compound_structure:
                mol = Chem.MolFromMolBlock(compound_structure)

                mol.SetProp("Compound_ID", compound)

                for precursor_id in precursor_ids:
                    flag = "1" if precursor_id in precursors else "0"
                    precursor_prop_name = precursor_id
                    mol.SetProp(precursor_prop_name, flag)

                mol.SetProp("Map_IDs", ";".join(compound_map_dict.get(compound, [])))

                writer.write(mol)

                if all(precursor_id == '0' for precursor_id in precursors):
                    print(f"Compound {compound} has zero precursors")

            else:
                result = KeggNetworkGenerator.get_compound_structure(compound.strip())
                if result:
                    mol = Chem.MolFromMolBlock(result)

                    mol.SetProp("Compound_ID", compound)

                    for precursor_id in precursor_ids:
                        flag = "1" if precursor_id in precursors else "0"
                        precursor_prop_name = precursor_id
                        mol.SetProp(precursor_prop_name, flag)

                    mol.SetProp("Map_IDs", ";".join(compound_map_dict.get(compound, [])))

                    writer.write(mol)

                    if all(precursor_id == '0' for precursor_id in precursors):
                        print(f"Compound {compound} has zero precursors")
                else:
                    print(f"Failed to retrieve structure for compound {compound}")

        writer.close()

