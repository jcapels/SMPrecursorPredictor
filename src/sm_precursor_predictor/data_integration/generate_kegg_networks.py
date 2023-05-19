import re

from sm_precursor_predictor.data_integration.kegg_api import KeggApi
import networkx as nx


class KeggNetworkGenerator:

    @staticmethod
    def get_kegg_network(map_id, cofactor_list=None):
        """
        Get the KEGG network for the given map_id
        """
        to_ignore = [f"C{str(x).zfill(5)}" for x in range(1, 51)]
        if cofactor_list is None:
            cofactor_list = []
        df_reactions_in_map = KeggApi.get_links("reaction", f"path:{map_id}")
        reactions = df_reactions_in_map.iloc[:, 1]

        G = nx.DiGraph()
        for reaction in reactions:
            reaction = KeggApi.to_df(KeggApi.get(reaction))
            for reaction_component in reaction.iloc[:, 0]:

                if "ENTRY" in reaction_component:
                    for entry in reaction_component.split(" "):
                        if entry.startswith("R"):
                            reaction_id = entry.strip()
                            if reaction_id == "R02177":
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
                            G.add_edge(substrate, reaction_id)
                            G.add_edge(reaction_id, substrate)

                    for product in products_identifiers:
                        if product not in cofactor_list and product not in to_ignore:
                            G.add_edge(reaction_id, product)
                            G.add_edge(product, reaction_id)

        return G
