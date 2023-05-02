from equilibrator_api import ComponentContribution


class EquilibratorCaller:

    def __init__(self):
        self.cc = ComponentContribution()

    def get_reaction_reversibility(self, reaction):
        result = self.cc.parse_reaction_formula(reaction)
        dG0_prime = self.cc.standard_dg_prime(result)
        return dG0_prime.magnitude.nominal_value, dG0_prime.magnitude.std_dev

