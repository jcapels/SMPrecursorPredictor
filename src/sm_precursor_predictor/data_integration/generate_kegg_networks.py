from sqlalchemy.orm.exc import MultipleResultsFound

from sm_precursor_predictor.data_integration.equilibrator_caller import EquilibratorCaller
from sm_precursor_predictor.data_integration.kegg_api import KeggApi
import networkx as nx


class KeggNetworkGenerator:

    @staticmethod
    def get_kegg_network(map_id):
        """
        Get the KEGG network for the given map_id
        """

        # equilibrator = EquilibratorCaller()

        df_reactions_in_map = KeggApi.get_links("reaction", f"path:{map_id}")
        reactions = df_reactions_in_map.iloc[:, 1]

        G = nx.DiGraph()
        for reaction in reactions:
            reaction = KeggApi.to_df(KeggApi.get(reaction))

            for reaction_component in reaction.iloc[:, 0]:

                # if "ENTRY" in reaction_component:
                #     for entry in reaction_component.split(" "):
                #         if entry.startswith("R"):
                #             reaction_id = entry.strip()
                #             break

                if "RCLASS" in reaction_component:
                    reaction_class = reaction_component.split("  ")
                    # get only the ones with CXXXXX
                    reaction_classes = [x.strip() for x in reaction_class if x.startswith("C")]
                    for reaction_class in reaction_classes:
                        substrate, product = reaction_class.split("_")
                        G.add_edge(substrate, product)

                # if "EQUATION" in reaction_component:
                #
                #     equation = reaction_component.split("  ")
                #     equation = " ".join(equation[1:]).strip()
                #     # get the all the identifiers in the form RXXXXX
                #     try:
                #         g_energy, std = equilibrator.get_reaction_reversibility(equation)
                #         if g_energy < -20:
                #             direction = "left_to_right"
                #         elif g_energy > 20:
                #             direction = "right_to_left"
                #         else:
                #             direction = "reversible"
                #     except:
                #         direction = "reversible"
                #
                #     if direction == "left_to_right":
                #         substrates, products = equation.split("<=>")
                #     elif direction == "right_to_left":
                #         products, substrates = equation.split("<=>")
                #     else:
                #         substrates, products = equation.split("<=>")
                #
                #     substrates_identifiers = [x.strip() for x in substrates.split(" ") if x.startswith("C")]
                #     products_identifiers = [x.strip() for x in products.split(" ") if x.startswith("C")]
                #
                #     for substrate in substrates_identifiers:
                #         G.add_edge(substrate, reaction_id)
                #
                #     for product in products_identifiers:
                #         G.add_edge(reaction_id, product)
                #
                #     if direction == "reversible":
                #         for substrate in substrates_identifiers:
                #             G.add_edge(reaction_id, substrate)
                #
                #         for product in products_identifiers:
                #             G.add_edge(product, reaction_id)

        return G
