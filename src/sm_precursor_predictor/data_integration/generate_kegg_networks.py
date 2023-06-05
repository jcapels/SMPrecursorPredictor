import re

from sm_precursor_predictor.data_integration.kegg_api import KeggApi
import networkx as nx
from urllib.error import HTTPError
from Bio.KEGG import REST as kegg_api
from rdkit import Chem
from tqdm import tqdm

class KeggNetworkGenerator:

   
    @staticmethod
    def get_compound_structure(compound_id):
        try:
            result = kegg_api.kegg_get(compound_id, option='mol').read()
            return result
        except HTTPError:
            print(compound_id)
            return None


    @staticmethod
    def convert_to_smiles(compound_id):
   
        structure = KeggNetworkGenerator.get_compound_structure(compound_id)
        if structure is not None:
            mol = Chem.MolFromMolBlock(structure)

            if mol is not None:
                smiles = Chem.MolToSmiles(mol)
                return smiles
            else:
                return None

    @staticmethod
    def get_kegg_network(map_id, cofactor_list=None):
        """
        Get the KEGG network for the given map_id
        """
        to_ignore = [f"C{str(x).zfill(5)}" for x in range(1, 47)]

        if cofactor_list is None:
            cofactor_list = []
        df_reactions_in_map = KeggApi.get_links("reaction", f"path:{map_id}")
        reactions = df_reactions_in_map.iloc[:, 1]

        G = nx.DiGraph()
        mol_attr = {}
        progress_bar = tqdm(total=len(reactions), desc=map_id)
        for reaction in reactions:
            progress_bar.update(1)
            reaction = KeggApi.to_df(KeggApi.get(reaction))
            for reaction_component in reaction.iloc[:, 0]:

                if "ENTRY" in reaction_component:
                    for entry in reaction_component.split(" "):
                        if entry.startswith("R"):
                            reaction_id = entry.strip()
                            if reaction_id == "R00462":
                                print()
                            break

                if "EQUATION" in reaction_component:

                    equation = reaction_component.split("  ")
                    equation = " ".join(equation[1:]).strip()
                    # get the all the identifiers in the form RXXXXX
                    substrates, products = equation.split("<=>")

                    substrates_identifiers = [x.strip() for x in substrates.split(" ") if x.startswith("C")]
                    products_identifiers = [x.strip() for x in products.split(" ") if x.startswith("C")]
                    for substrate in substrates_identifiers:
                        if substrate not in cofactor_list and substrate not in to_ignore:
                            if substrate in mol_attr.keys():
                                G.add_edge(substrate, reaction_id)
                                G.add_edge(reaction_id, substrate)
                            else:
                                mol = KeggNetworkGenerator.get_compound_structure(substrate)
                                if mol is not None:
                                    mol_attr[substrate] = {"mol": mol}
                                    G.add_edge(substrate, reaction_id)
                                    G.add_edge(reaction_id, substrate)

                    for product in products_identifiers:
                        if product not in cofactor_list and product not in to_ignore:
                            
                            if product in mol_attr.keys():
                                G.add_edge(reaction_id, product)
                                G.add_edge(product, reaction_id)
                            else:
                                mol = KeggNetworkGenerator.get_compound_structure(product)
                                if mol is not None:
                                    mol_attr[product] = {"mol": mol}
                                    G.add_edge(product, reaction_id)
                                    G.add_edge(reaction_id, product)
        
        nx.set_node_attributes(G, mol_attr)

        return G

    @staticmethod
    def generate_graphs_for_map_ids(map_ids):
        """
        Generate and store graphs for the given map IDs as a pickle file.
        """
        graph_dict = {}

        for map_id in map_ids:
            graph = KeggNetworkGenerator.get_kegg_network(map_id)
            graph_dict[map_id] = graph

        return graph_dict