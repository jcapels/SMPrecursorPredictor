from unittest import TestCase


class TestEquilibrator(TestCase):

    def test_equilibrator(self):
        from equilibrator_api import ComponentContribution
        cc = ComponentContribution()
        result = cc.parse_reaction_formula("C03470 + C00001 <=> C00398 + C01852")
        dG0_prime = cc.standard_dg_prime(result)
        print(dG0_prime)