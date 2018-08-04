import logic
from unittest import TestCase
import graphs
import sympy
from collections import namedtuple
import random
from attractors import find_num_attractors_onestage, vertex_impact_scores, find_num_steady_states, \
    find_attractors_dubrova, find_attractors_onestage_enumeration

ILPAttractorExperimentParameters = namedtuple("AttractorExperimentParameters", "G T P n_attractors")
VertexImpactExperimentParameters = namedtuple("VertexImpactExperimentParameters", "G T P impacts")
DubrovaExperimentParameters = namedtuple("DubrovaExperimentParameters", "G mutate n_attractors")


class TestAttractors(TestCase):
    def test_num_attractors_onestage(self):
        return True
        experiments = []

        """test on known toy models"""
        # 0, 1
        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=1, n_attractors=0))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=3, n_attractors=1))
        # 2, 3
        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[logic.SymmetricThresholdFunction(signs=[-1], threshold=1)])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=1, n_attractors=0))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=3, n_attractors=1))
        # 4, 5
        G = graphs.Network(vertex_names=["A"], edges=[],
                           vertex_functions=[None])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=3, n_attractors=2))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=3, n_attractors=2))
        # 6, 7
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[logic.SymmetricThresholdFunction(signs=[1], threshold=1),
                                             None])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=5, n_attractors=4))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=5, n_attractors=4))
        # 8, 9
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[logic.SymmetricThresholdFunction(signs=[-1], threshold=1),
                                             None])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=1, n_attractors=0))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=3, n_attractors=2))
        # 10, 11
        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.And])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=2, n_attractors=2))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=3, P=1, n_attractors=1))
        # 12, 13
        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[None])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=3, n_attractors=2))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=3, n_attractors=2))
        # 14, 15
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[None, None])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=5, n_attractors=4))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=6, n_attractors=4))
        # 16, 17
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[None, True])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=5, n_attractors=2))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=6, n_attractors=2))
        # 18, 19, 20
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[sympy.Nand, sympy.And])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=3, n_attractors=0))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=4, P=2, n_attractors=1))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=4, P=1, n_attractors=1))
        # 21, 22, 23
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=3, n_attractors=2))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=3, n_attractors=3))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=15, P=15, n_attractors=3))
        # 24, 25
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[lambda x: True, lambda x: False])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=4, P=2, n_attractors=1))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=2, n_attractors=1))
        # 26, 27
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[None, sympy.And])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=4, n_attractors=3))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=4, n_attractors=2))
        # 28, 29
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[None, lambda _: True])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=3, n_attractors=1))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=4, P=2, n_attractors=1))
        # 30, 31
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[None, None])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=6, n_attractors=3))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=6, n_attractors=2))
        # 32, 33, 34, 35, 36
        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "A"), ("C", "A")],
                           vertex_functions=[logic.SymmetricThresholdFunction.from_function(sympy.Nand, 2),
                                             logic.SymmetricThresholdFunction.from_function(sympy.Nand, 1),
                                             logic.SymmetricThresholdFunction.from_function(sympy.Nand, 0)])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=3, n_attractors=3))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=4, n_attractors=3))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=3, n_attractors=3))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=4, n_attractors=4))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=3, P=4, n_attractors=4))
        # 37
        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "A"), ("C", "A")],
                           vertex_functions=[logic.SymmetricThresholdFunction.from_function(sympy.Nand, 2),
                                             logic.SymmetricThresholdFunction.from_function(sympy.Nand, 1),
                                             None])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=3, P=3, n_attractors=3))
        # 38, 39, 40
        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand]*3)
        experiments.append(ILPAttractorExperimentParameters(G=G, T=6, P=2, n_attractors=2))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=10, P=10, n_attractors=2))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=5, P=10, n_attractors=1))
        # 41, 42
        # acyclic, should have 2**#input_nodes attractors of length 1
        G = graphs.Network(vertex_names=["v1", "v2", "v3", "v4", "v5", "v6"],
                           edges=[("v1", "v4"), ("v2", "v4"), ("v1", "v5"), ("v4", "v6")],
                           vertex_functions=[sympy.Nand]*6)
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=10, n_attractors=8))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=6, P=10, n_attractors=8))
        # 43, 44, 45
        G = graphs.Network(vertex_names=["A1", "B1", "B2", "C1", "C2"],
                           edges=[("A1", "A1"), ("B1", "B2"), ("B2", "B1"), ("C1", "C2"), ("C2", "C1")],
                           vertex_functions=[sympy.And]*5)
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=10, n_attractors=8))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=18, n_attractors=18))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=3, P=40, n_attractors=20))  # offsets!
        # 46, 47, 48
        # a failed random graph added as a constant test
        G = graphs.Network(
            vertex_names=['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16',
                          '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31',
                          '32', '33', '34'],
            edges=[('1', '2'), ('2', '16'), ('3', '17'), ('5', '15'), ('6', '29'), ('7', '28'), ('8', '22'),
                   ('9', '28'), ('10', '18'), ('11', '15'), ('12', '24'), ('13', '14'), ('15', '18'), ('16', '26'),
                   ('17', '27'), ('18', '20'), ('19', '23'), ('20', '27'), ('23', '26'), ('24', '29'), ('25', '33'),
                   ('26', '30'), ('27', '32'), ('28', '32'), ('30', '32'), ('31', '34'), ('32', '33'), ('33', '34')],
            vertex_functions=[None, None, sympy.Nand, None, None, None, None, None, None, None, None, None, None, None,
                              sympy.Or, sympy.Nand,
                              sympy.Nand, sympy.Nand, sympy.Nand, None, sympy.Xor, None, sympy.And, sympy.Nand,
                              sympy.Xor, None, sympy.And, sympy.Nand, sympy.And, sympy.Xor, sympy.Or, None, sympy.Or,
                              sympy.And, sympy.And])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=6, n_attractors=6))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=10, n_attractors=10))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=10, n_attractors=10))
        # 49, 50, 51
        G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
               "\\Attractors - for Ariel\\BNS_Dubrova_2011\\MAPK_large2.cnet")
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=15, n_attractors=12))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=15, n_attractors=14))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=3, P=15, n_attractors=14))

        # for _ in range(5):
        #     size = 35
        #     G = graphs.Network(vertex_names=list(range(size)),
        #                        edges=[(i, random.choice(list(range(i+1, size)))) for i in range(size)
        #                               if random.random() < 0.8 and i != size-1],
        #                        vertex_functions=[random.choice([sympy.And, sympy.Nand, sympy.Or, sympy.Xor])
        #                                          for _ in range(size)])
        #     input_nodes = 0
        #     for v in G.vertices:
        #         is_input = True
        #         for e in G.edges:
        #             if e[1] == v:
        #                 is_input = False
        #                 break
        #         if is_input:
        #             input_nodes += 1
        #     attractor_number = 2**input_nodes
        #     experiments.append(ExperimentParameters(G=G, T=1, P=3, n_attractors=min(3, attractor_number)))
        #     experiments.append(ExperimentParameters(G=G, T=2, P=10, n_attractors=min(10, attractor_number)))
        #     experiments.append(ExperimentParameters(G=G, T=10, P=3, n_attractors=min(3, attractor_number)))

        # TODO: figure out how disjoint long attractors work together (multiplying doesn't account for offsets)
        # """test on basic semi-random networks: create connectivity components of acyclis networks and simple cycles"""
        # n_random_experiment = 0
        # while n_random_experiment < 10:
        #     n_components = random.randint(1, 3)
        #     attractor_number = 1
        #     max_attractor_len = 0
        #     cur_graph = None
        #     for n_component in range(n_components):  # TODO: change to graph union method
        #         comp_size = random.randint(1, 5)
        #         V = [i for i in range(comp_size)]
        #         E = []
        #         comp_type =random.choice(["cycle", "acyclic"])
        #         if comp_type == "acyclic":
        #             for i in range(len(V) - 1): # create only forward facing edges
        #                 for j in range(i+1, len(V)):
        #                     if random.random() <= 0.8:
        #                         E.append((V[i], V[j]))
        #             component_graph = graphs.Network(vertex_names=V, edges=E)
        #             restriction_level = random.choice([graphs.FunctionTypeRestriction.NONE,
        #                                                graphs.FunctionTypeRestriction.SYMMETRIC_THRESHOLD,
        #                                                graphs.FunctionTypeRestriction.SIMPLE_GATES])
        #             component_graph.randomize_functions(function_type_restriction=restriction_level)
        #             input_nodes = 0
        #             for v in V:
        #                 is_input = True
        #                 for e in E:
        #                     if e[1] == v:
        #                         is_input = False
        #                         break
        #                 if is_input:
        #                     input_nodes += 1
        #             attractor_number *= 2**input_nodes
        #             max_attractor_len = max(max_attractor_len, 1)
        #         elif comp_type == "cycle":
        #             """currently supports only a cycle of identity function, using a group theory theorem from
        #             https://www.quora.com/How-many-unique-binary-matrices-are-there-up-to-rotations-translations-and-flips
        #             , can later add negation cycles"""
        #             for i in range(len(V)):
        #                 E.append((V[i], V[(i + 1) % len(V)]))
        #             component_graph = graphs.Network(vertex_names=V, edges=E, vertex_functions=[sympy.And]*len(V))
        #             attractor_number *= binary_necklaces(len(V))
        #             max_attractor_len = max(max_attractor_len, len(V))
        #         cur_graph = component_graph if cur_graph is None else cur_graph + component_graph
        #     if attractor_number * len(cur_graph.vertices) * max_attractor_len <= 250:
        #         experiments.append(ExperimentParameters(G=cur_graph, T=max_attractor_len,
        #                                                 P=attractor_number + 1,
        #                                                 n_attractors=attractor_number))
        #         n_random_experiment += 1

        print "number of experiments (with keys)={}".format(len(experiments))
        for i, experiment in enumerate(experiments):
            print "experiment #{}".format(i)
            print "n={}, T={}, P={}, expected_n_attractors={}".format(len(experiment.G.vertices),
                                                                   experiment.T, experiment.P, experiment.n_attractors)
            # continue
            use_sampling = bool(random.randint(0, 1))
            use_sampling_for_mip_start = bool(random.randint(0, 1))
            simplify = bool(random.randint(0, 1))
            key_slice_size = random.randint(1, 15)
            print "key_slice_size={}".format(key_slice_size)
            n_attractors = find_num_attractors_onestage(G=experiment.G, max_len=experiment.T, max_num=experiment.P,
                                                        use_sat=False, verbose=False,
                                                        sampling_bounds=(3, 3) if use_sampling else None,
                                                        use_sampling_for_mip_start=use_sampling_for_mip_start,
                                                        simplify_general_boolean=simplify,
                                                        key_slice_size=key_slice_size)
            try:
                self.assertEqual(n_attractors, experiment.n_attractors)
            except AssertionError as e:
                print e
                print experiment.G
                raise e
            except Exception as e:
                raise e

        # print "number of experiments (without keys)={}".format(len(experiments))
        # for i, experiment in enumerate(experiments):
        #     print "experiment #{}".format(i)
        #     print "n={}, T={}, P={}, expected_n_attractors={}".format(len(experiment.G.vertices),
        #                                                            experiment.T, experiment.P, experiment.n_attractors)
        #     # continue
        #     n_attractors = find_num_attractors_onestage(G=experiment.G, max_len=experiment.T, max_num=experiment.P,
        #                                                 use_sat=False, verbose=False,
        #                                                 use_state_keys=False, require_result=experiment.n_attractors)
        #     try:
        #         self.assertEqual(n_attractors, experiment.n_attractors)
        #     except AssertionError as e:
        #         print e
        #         print experiment.G
        #         raise e

    def test_vertex_impact_scores(self):
        # TODO: also test the resulting models (assure they have the correct number of attractors)
        experiments = []

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand])
        experiments.append(VertexImpactExperimentParameters(G=G, T=3, P=3, impacts=[2]))

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand])
        experiments.append(VertexImpactExperimentParameters(G=G, T=1, P=3, impacts=[2]))

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand])
        experiments.append(VertexImpactExperimentParameters(G=G, T=1, P=1, impacts=[1]))

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.And])
        experiments.append(VertexImpactExperimentParameters(G=G, T=3, P=3, impacts=[2]))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[sympy.And, None])
        experiments.append(VertexImpactExperimentParameters(G=G, T=3, P=5, impacts=[4, 4]))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand, None])
        experiments.append(VertexImpactExperimentParameters(G=G, T=3, P=5, impacts=[4, 2]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.Nand])
        experiments.append(VertexImpactExperimentParameters(G=G, T=7, P=5, impacts=[4, 4, 4]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.Nand])
        experiments.append(VertexImpactExperimentParameters(G=G, T=3, P=5, impacts=[4, 4, 4]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.Nand])
        experiments.append(VertexImpactExperimentParameters(G=G, T=2, P=5, impacts=[2, 2, 2]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.Nand])
        experiments.append(VertexImpactExperimentParameters(G=G, T=1, P=5, impacts=[2, 2, 2]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "A"), ("C", "A")],
                           vertex_functions=[sympy.And, sympy.And, None])
        experiments.append(VertexImpactExperimentParameters(G=G, T=2, P=5, impacts=[6, 4, 4]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "A"), ("C", "A")],
                           vertex_functions=[sympy.And, sympy.And, None])
        experiments.append(VertexImpactExperimentParameters(G=G, T=1, P=5, impacts=[4, 3, 3]))

        print "number of experiments (with keys)={}".format(len(experiments))
        for i, experiment in enumerate(experiments):
            print "experiment #{}".format(i)
            print "n={}, T={}, P={}, expected_impacts={}".format(len(experiment.G.vertices),
                                                                      experiment.T, experiment.P,
                                                                      experiment.impacts)
            # continue
            impacts, _ = vertex_impact_scores(G=experiment.G, attractor_length_threshold=experiment.T,
                                              attractor_num_threshold=experiment.P,
                                              model_type_restriction=graphs.FunctionTypeRestriction.NONE)
            try:
                self.assertEqual(impacts, experiment.impacts)
            except AssertionError as e:
                print e
                print experiment.G
                raise e

    def test_find_num_steady_states(self):
        """test on known toy models"""
        # 0, 1
        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand])
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=False), 0)
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=True), 0)

        G = graphs.Network(vertex_names=["A"], edges=[],
                           vertex_functions=[None])
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=False), 2)
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=True), 2)

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.And])
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=False), 2)
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=True), 2)

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[sympy.Nand, sympy.And])
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=False), 0)
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=True), 0)

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand])
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=False), 2)

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[lambda x: True, lambda x: False])
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=False), 1)
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=True), 1)

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand]*3)
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=False), 0)

        G = graphs.Network(vertex_names=["A", "B", "C", "D"], edges=[("A", "B"), ("B", "C"), ("C", "D"), ("D", "A")],
                           vertex_functions=[sympy.Nand]*4)
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=False), 2)

        # acyclic, should have 2**#input_nodes attractors of length 1
        G = graphs.Network(vertex_names=["v1", "v2", "v3", "v4", "v5", "v6"],
                           edges=[("v1", "v4"), ("v2", "v4"), ("v1", "v5"), ("v4", "v6")],
                           vertex_functions=[sympy.Nand]*6)
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=False), 8)

        G = graphs.Network(vertex_names=["A1", "B1", "B2", "C1", "C2"],
                           edges=[("A1", "A1"), ("B1", "B2"), ("B2", "B1"), ("C1", "C2"), ("C2", "C1")],
                           vertex_functions=[sympy.And]*5)
        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand]*3)
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=False), 0)

        G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
               "\\Attractors - for Ariel\\BNS_Dubrova_2011\\MAPK_large2.cnet")
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=False), 12)

    def find_attractors_dubrova(self):
        experiments = []

        """test on known toy models"""
        # 0, 1
        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand])
        experiments.append(DubrovaExperimentParameters(G=G, mutate=False, n_attractors=1))
        experiments.append(DubrovaExperimentParameters(G=G, mutate=True, n_attractors=1))
        # 2
        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[logic.SymmetricThresholdFunction(signs=[-1], threshold=1)])
        experiments.append(DubrovaExperimentParameters(G=G, mutate=False, n_attractors=1))
        # 3, 4
        G = graphs.Network(vertex_names=["A"], edges=[],
                           vertex_functions=[None])
        experiments.append(DubrovaExperimentParameters(G=G, mutate=False, n_attractors=1))
        experiments.append(DubrovaExperimentParameters(G=G, mutate=True, n_attractors=2))
        # 5, 6
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[sympy.Nand, sympy.And])
        experiments.append(DubrovaExperimentParameters(G=G, mutate=False, n_attractors=1))
        experiments.append(DubrovaExperimentParameters(G=G, mutate=True, n_attractors=1))
        # 7, 8
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[lambda x: True, lambda x: False])
        experiments.append(DubrovaExperimentParameters(G=G, mutate=False, n_attractors=1))
        experiments.append(DubrovaExperimentParameters(G=G, mutate=True, n_attractors=1))
        # 9, 10
        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "A"), ("C", "A")],
                           vertex_functions=[logic.SymmetricThresholdFunction.from_function(sympy.Nand, 2),
                                             logic.SymmetricThresholdFunction.from_function(sympy.Nand, 1),
                                             logic.SymmetricThresholdFunction.from_function(sympy.Nand, 0)])
        experiments.append(DubrovaExperimentParameters(G=G, mutate=False, n_attractors=3))
        experiments.append(DubrovaExperimentParameters(G=G, mutate=True, n_attractors=3))
        # 11, 12
        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "A"), ("C", "A")],
                           vertex_functions=[logic.SymmetricThresholdFunction.from_function(sympy.Nand, 2),
                                             logic.SymmetricThresholdFunction.from_function(sympy.Nand, 1),
                                             True])
        experiments.append(DubrovaExperimentParameters(G=G, mutate=False, n_attractors=1))
        experiments.append(DubrovaExperimentParameters(G=G, mutate=True, n_attractors=3))

        # 13
        G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
               "\\Attractors - for Ariel\\BNS_Dubrova_2011\\MAPK_large.cnet")
        experiments.append(DubrovaExperimentParameters(G=G, mutate=True, n_attractors=16))

        print "number of experiments (with keys)={}".format(len(experiments))
        for i, experiment in enumerate(experiments):
            print "experiment #{}".format(i)
            print "n={}, mutate={}, expected_n_attractors={}".format(len(experiment.G.vertices),
                                                                     experiment.mutate, experiment.n_attractors)
            # continue
            n_attractors = len(find_attractors_dubrova(G=experiment.G,
                                                   dubrova_dir_path="C:/Users/ariel/Downloads/Attractors - for Ariel/"
                                                                    "Attractors - for Ariel/BNS_Dubrova_2011",
                                                   mutate_input_nodes=experiment.mutate))
            try:
                self.assertEqual(n_attractors, experiment.n_attractors)
            except AssertionError as e:
                print e
                print experiment.G
                raise e
            except Exception as e:
                raise e

    def test_find_attractors_enumerate(self):
        return True
        experiments = []

        """test on known toy models"""
        # 0, 1
        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=None, n_attractors=0))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=None, n_attractors=1))
        # 2, 3
        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[logic.SymmetricThresholdFunction(signs=[-1], threshold=1)])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=None, n_attractors=0))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=None, n_attractors=1))
        # 4, 5
        G = graphs.Network(vertex_names=["A"], edges=[],
                           vertex_functions=[None])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=None, n_attractors=2))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=None, n_attractors=2))
        # 6, 7
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[logic.SymmetricThresholdFunction(signs=[-1], threshold=1),
                                             None])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=None, n_attractors=0))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=None, n_attractors=2))
        # 8, 9
        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.And])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=None, n_attractors=2))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=3, P=None, n_attractors=2))
        # 10, 11
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[sympy.Nand, sympy.And])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=None, n_attractors=0))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=4, P=None, n_attractors=1))
        # 12, 13, 14
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=None, n_attractors=2))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=None, n_attractors=3))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=15, P=None, n_attractors=3))
        # 15, 16
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[lambda x: True, lambda x: False])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=4, P=None, n_attractors=1))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=None, n_attractors=1))
        # 17, 18, 19
        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "A"), ("C", "A")],
                           vertex_functions=[logic.SymmetricThresholdFunction.from_function(sympy.Nand, 2),
                                             logic.SymmetricThresholdFunction.from_function(sympy.Nand, 1),
                                             logic.SymmetricThresholdFunction.from_function(sympy.Nand, 0)])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=None, n_attractors=3))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=None, n_attractors=4))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=3, P=None, n_attractors=4))
        # 20
        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "A"), ("C", "A")],
                           vertex_functions=[logic.SymmetricThresholdFunction.from_function(sympy.Nand, 2),
                                             logic.SymmetricThresholdFunction.from_function(sympy.Nand, 1),
                                             None])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=3, P=None, n_attractors=4))
        # 21, 22, 23
        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand]*3)
        experiments.append(ILPAttractorExperimentParameters(G=G, T=6, P=None, n_attractors=2))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=10, P=None, n_attractors=2))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=5, P=None, n_attractors=1))
        # 24, 25
        # acyclic, should have 2**#input_nodes attractors of length 1
        G = graphs.Network(vertex_names=["v1", "v2", "v3", "v4", "v5", "v6"],
                           edges=[("v1", "v4"), ("v2", "v4"), ("v1", "v5"), ("v4", "v6")],
                           vertex_functions=[sympy.Nand]*6)
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=None, n_attractors=8))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=6, P=None, n_attractors=8))
        # 26, 27
        G = graphs.Network(vertex_names=["A1", "B1", "B2", "C1", "C2"],
                           edges=[("A1", "A1"), ("B1", "B2"), ("B2", "B1"), ("C1", "C2"), ("C2", "C1")],
                           vertex_functions=[sympy.And]*5)
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=None, n_attractors=8))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=3, P=None, n_attractors=20))  # offsets!
        # 28, 29
        # a failed random graph added as a constant test
        G = graphs.Network(
            vertex_names=['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16',
                          '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31',
                          '32', '33', '34'],
            edges=[('1', '2'), ('2', '16'), ('3', '17'), ('5', '15'), ('6', '29'), ('7', '28'), ('8', '22'),
                   ('9', '28'), ('10', '18'), ('11', '15'), ('12', '24'), ('13', '14'), ('15', '18'), ('16', '26'),
                   ('17', '27'), ('18', '20'), ('19', '23'), ('20', '27'), ('23', '26'), ('24', '29'), ('25', '33'),
                   ('26', '30'), ('27', '32'), ('28', '32'), ('30', '32'), ('31', '34'), ('32', '33'), ('33', '34')],
            vertex_functions=[None, None, sympy.Nand, None, None, None, None, None, None, None, None, None, None, None,
                              sympy.Or, sympy.Nand,
                              sympy.Nand, sympy.Nand, sympy.Nand, None, sympy.Xor, None, sympy.And, sympy.Nand,
                              sympy.Xor, None, sympy.And, sympy.Nand, sympy.And, sympy.Xor, sympy.Or, None, sympy.Or,
                              sympy.And, sympy.And])
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=None, n_attractors=2**17))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=None, n_attractors=2**17))
        # 30, 31, 32
        G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
               "\\Attractors - for Ariel\\BNS_Dubrova_2011\\MAPK_large2.cnet")
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=None, n_attractors=12))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=None, n_attractors=14))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=3, P=None, n_attractors=14))

        print "number of experiments (with keys)={}".format(len(experiments))
        for i, experiment in enumerate(experiments):
            print "experiment #{}".format(i)
            print "n={}, T={}, expected_n_attractors={}".format(len(experiment.G.vertices),
                                                                   experiment.T, experiment.n_attractors)
            # continue
            simplify = bool(random.randint(0, 1))
            key_slice_size = random.randint(1, 15)
            print "key_slice_size={}".format(key_slice_size)
            n_attractors = len(find_attractors_onestage_enumeration(G=experiment.G, max_len=experiment.T,
                                                                verbose=False,
                                                                simplify_general_boolean=simplify,
                                                                key_slice_size=key_slice_size))
            try:
                self.assertEqual(n_attractors, experiment.n_attractors)
            except AssertionError as e:
                print e
                print experiment.G
                raise e
            except Exception as e:
                raise e


# TODO: add dubrova v.s. ILP testing again.
