import os
import shutil
from unittest import TestCase
from graphs import Network, FunctionTypeRestriction
import random
import sympy
import utility

class TestNetwork(TestCase):

    def test_cnet_export_and_import(self):
        for _ in range(10):
            n = random.randint(1, 20)
            for restriction in [FunctionTypeRestriction.NONE, FunctionTypeRestriction.SYMMETRIC_THRESHOLD]:
                # TODO: implement and test for FunctionTypeRestriction.SIMPLE_GATES
                if random.choice([False, True]):
                    indegree_bounds = [0, 4]
                    indegree_geometric_p = None
                else:
                    indegree_bounds = None
                    indegree_geometric_p = 0.8
                G = Network.generate_random(n_vertices=n, indegree_bounds=indegree_bounds,
                                            function_type_restriction=restriction,
                                            indegree_geometric_p=indegree_geometric_p)
                G.export_to_cnet("./temp_test_network.cnet")
                G_tag = Network.parse_cnet("./temp_test_network.cnet")
                self.assertEqual(G, G_tag)

        Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
                       "\\Attractors - for Ariel\\BNS_Dubrova_2011\\MAPK_large2.cnet")
        Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
                                      "\\Attractors - for Ariel\\BNS_Dubrova_2011\\arabidopsis.cnet.txt")
        Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
                       "\\Attractors - for Ariel\\BNS_Dubrova_2011\\EGFR_man_with_inputs_all_zero.cnet")
        Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
                       "\\Attractors - for Ariel\\BNS_Dubrova_2011\\MAPK_large.cnet")
        Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
                       "\\Attractors - for Ariel\\BNS_Dubrova_2011\\thelper.cnet.txt")
        Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
                       "\\Attractors - for Ariel\\BNS_Dubrova_2011\\tcr.cnet")

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

    def test_pow_mult(self):
        # manual with toy examples
        G = Network(vertex_names=["a", "b"], edges=[("a", "a"), ("a", "b"), ("b", "a")],
                    vertex_functions=[sympy.And, sympy.And])
        G_squared = G * G
        # 00 -> 00 -> 00
        self.assertTrue(utility.is_same_state(G_squared.next_state([False, False], return_as_string=False),
                                              [False, False]))
        # 01 -> 00 -> 00
        self.assertTrue(utility.is_same_state(G_squared.next_state([False, True], return_as_string=False),
                                               [False, False]))
        # 10 -> 01 -> 00
        self.assertTrue(utility.is_same_state(G_squared.next_state([True, False], return_as_string=False),
                                              [False, False]))
        # 11 -> 11 -> 11
        self.assertTrue(utility.is_same_state(G_squared.next_state([True, True], return_as_string=False),
                                              [True, True]))

        G = Network(vertex_names=["a", "b"], edges=[("a", "a"), ("a", "b"), ("b", "a")],
                    vertex_functions=[sympy.Nand, sympy.And])
        G_squared = G * G
        # 00 -> 10 -> 11
        self.assertTrue(utility.is_same_state(G_squared.next_state([False, False], return_as_string=False),
                                              [True, True]))
        # 01 -> 10 -> 11
        self.assertTrue(utility.is_same_state(G_squared.next_state([False, True], return_as_string=False),
                                              [True, True]))
        # 10 -> 11 -> 01
        self.assertTrue(utility.is_same_state(G_squared.next_state([True, False], return_as_string=False),
                                              [False, True]))
        # 11 -> 01 -> 10
        self.assertTrue(utility.is_same_state(G_squared.next_state([True, True]),
                                              [True, False]))

        # random testing
        test_graphs = [Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
                            "\\Attractors - for Ariel\\BNS_Dubrova_2011\\MAPK_large2.cnet")] +\
                            [Network.generate_random(n_vertices=random.randint(1, 10), indegree_bounds=(0, 5))
                                for _ in range(5)]
        for G in test_graphs:
            G_squared = G ** 2
            G_cubed = G ** 3
            for i in range(20):
                starting_state = [bool(random.randint(0, 1)) for _ in range(len(G.vertices))]
                true_double_step = G.next_state(G.next_state(starting_state, return_as_string=False),
                                                return_as_string=False)
                true_triple_step = G.next_state(G.next_state(G.next_state(starting_state, return_as_string=False),
                                                return_as_string=False), return_as_string=False)
                self.assertTrue(utility.is_same_state(G_squared.next_state(starting_state, return_as_string=False),
                                                      true_double_step))
                self.assertTrue(utility.is_same_state(G_cubed.next_state(starting_state, return_as_string=False),
                                                      true_triple_step))

    def test_randomize_incoming_edges(self):
        for i in range(50):
            n = random.randint(1, 10)
            G = Network.generate_random(n_vertices=n, indegree_bounds=[1, 3])
            G_tag = G.copy()
            G_tag.randomize_incoming_edges(include_self_loops=False)
            for v, v_tag in zip(G.vertices, G_tag.vertices):
                self.assertTrue(v.function == v_tag.function)
                self.assertTrue(len(v.predecessors()) == len(v_tag.predecessors()))

        G = Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("A", "A"), ("B", "A")],
                    vertex_functions=[None] * 3)
        for i in range(100):
            G_tag = G.copy()
            G_tag.randomize_incoming_edges(include_self_loops=False)
            self.assertTrue(G_tag.get_vertex("A").predecessors() == [G_tag.get_vertex("B"), G_tag.get_vertex("C")])
            self.assertTrue(G_tag.get_vertex("B").predecessors() in [[G_tag.get_vertex("A")], [G_tag.get_vertex("C")]])
            self.assertTrue(G_tag.get_vertex("C").predecessors() == [])

        had_self_loop = False
        for i in range(100):
            G_tag = G.copy()
            G_tag.randomize_incoming_edges(include_self_loops=True)
            self.assertTrue(G_tag.get_vertex("A").predecessors() in [[G_tag.get_vertex("A"), G_tag.get_vertex("B")],
                                                                     [G_tag.get_vertex("B"), G_tag.get_vertex("C")],
                                                                     [G_tag.get_vertex("A"), G_tag.get_vertex("C")]])
            self.assertTrue(G_tag.get_vertex("B").predecessors() in [[G_tag.get_vertex("A")],
                                                                     [G_tag.get_vertex("B")],
                                                                     [G_tag.get_vertex("C")]])
            self.assertTrue(G_tag.get_vertex("C").predecessors() == [])
            had_self_loop = (G.vertices[0] in G.vertices[0].predecessors()) or \
                            (G.vertices[1] in G.vertices[1].predecessors()) or had_self_loop
        self.assertTrue(had_self_loop)

    def test_randomize_edges_by_switching(self):
        for i in range(10):
            n = random.randint(1, 10)
            G = Network.generate_random(n_vertices=n, indegree_bounds=[1, 3])
            G_tag = G.copy()
            G_tag.randomize_edges_by_switching(include_self_loops=False)
            for v, v_tag in zip(G.vertices, G_tag.vertices):
                self.assertTrue(v.function == v_tag.function)
                self.assertTrue(len(v.predecessors()) == len(v_tag.predecessors()))

        G = Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("A", "A"), ("B", "A"), ("B", "C")],
                    vertex_functions=[None] * 3)
        for i in range(10):
            G_tag = G.copy()
            G_tag.randomize_edges_by_switching(include_self_loops=False)
            self.assertTrue(G_tag.get_vertex("A").predecessors() == [G_tag.get_vertex("A"), G_tag.get_vertex("B")])
            self.assertTrue(G_tag.get_vertex("B").predecessors() == [G_tag.get_vertex("A")])
            self.assertTrue(G_tag.get_vertex("C").predecessors() == [G_tag.get_vertex("B")])

    def test_randomize_functions(self):
        self.assertTrue(False)  # TODO: implement...

    def test_truth_table_import_export(self):

        for model_dir in os.listdir("../cellcollective_models"):
            print model_dir
            try:
                G = Network.parse_boolean_tables(os.path.join("../cellcollective_models", model_dir))
            except ValueError as e:
                if e.message.startswith("Model export from cellcollective failed"):
                    print "\nWarning: {}\n".format(e.message)
                    pass
                else:
                    raise e

        for i in range(50):
            n = random.randint(1, 10)
            G = Network.generate_random(n_vertices=n, indegree_bounds=[1, 3])
            if os.path.exists("temp_dir"):
                shutil.rmtree("temp_dir", ignore_errors=True)
            G.export_to_boolean_tables(".", "temp_dir")
            G_tag = Network.parse_boolean_tables("temp_dir")
            self.assertTrue(G == G_tag)
            # TODO: permute vertex order in the species list file and make sure the import still goes well.
            # TODO: load corresponding models from Dubrova's cnet and this format and make sure they are equal.
