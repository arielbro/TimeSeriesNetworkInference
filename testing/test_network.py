from unittest import TestCase
from graphs import Network, FunctionTypeRestriction
import random
import sympy


# noinspection PyPep8Naming
class TestNetwork(TestCase):

    def test_cnet_export_and_import(self):
        n = random.randint(1, 20)
        for restriction in [FunctionTypeRestriction.NONE, FunctionTypeRestriction.SYMMETRIC_THRESHOLD,
                            FunctionTypeRestriction.SIMPLE_GATES]:
            G = Network.generate_random(n_vertices=n, indegree_bounds=[1, 20],
                                        function_type_restriction=restriction)
            G.export_to_cnet("./temp_test_network.cnet")
            G_tag = Network.parse_cnet("./temp_test_network.cnet")
            self.assertEqual(G, G_tag)

    def test_union(self):
        # check toy examples
        G1 = Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C")],
                     vertex_functions=[sympy.And]*3)
        G2 = Network(vertex_names=["A", "X"], edges=[("A", "X")], vertex_functions=[None, sympy.Nor])
        sum_true = Network(vertex_names=["a_A", "a_B", "b_A", "b_X", "a_C"],  # order changed on purpose
                           edges=[("a_A", "a_B"), ("b_A", "b_X"), ("a_B", "a_C")],
                           vertex_functions=[sympy.And, sympy.And, sympy.Nor, sympy.Nor, sympy.And])
        sum_false_1 = Network(vertex_names=["a_A", "a_B", "b_A", "b_X", "a_C"],
                              edges=[("a_A", "a_B"), ("b_A", "b_X"), ("a_B", "a_C")],
                              vertex_functions=[sympy.And, sympy.Nor, sympy.And, sympy.Nor, sympy.And])
        sum_false_2 = Network(vertex_names=["a_B", "b_A", "b_X", "a_C"],
                              edges=[("b_A", "b_X"), ("a_B", "a_C")],
                              vertex_functions=[sympy.And, sympy.Nor, sympy.Nor, sympy.And])
        G_sum = G1 + G2  # order matters, for now.
        self.assertEqual(G_sum, sum_true)
        self.assertNotEqual(G_sum, sum_false_1)
        self.assertNotEqual(G_sum, sum_false_2)
