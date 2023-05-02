from unittest import TestCase


class TestEquilibrator(TestCase):

    def test_equilibrator(self):
        from equilibrator_api import ComponentContribution
        cc = ComponentContribution()
        result = cc.parse_reaction_formula("C02485 + C00005 + C00080 <=> C11951 + C00006")
        dG0_prime = cc.standard_dg_prime(result)
        print(dG0_prime.magnitude.nominal_value), print(dG0_prime.magnitude.std_dev)