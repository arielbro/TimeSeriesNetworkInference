import numpy as np
import logic
from unittest import TestCase
import graphs
import sympy
from collections import namedtuple
import random
from attractors import find_num_attractors_onestage, \
    vertex_model_impact_scores, stochastic_vertex_model_impact_scores, find_num_steady_states, \
    find_attractors_dubrova, find_attractors_onestage_enumeration, ImpactType, \
    vertex_state_impact_scores, stochastic_vertex_state_impact_scores, graph_model_impact_score, \
    graph_state_impact_score, stochastic_graph_model_impact_score, stochastic_graph_state_impact_score

import attractors

dubrova_path = "../" + attractors.dubrova_path

ILPAttractorExperimentParameters = namedtuple("AttractorExperimentParameters", "G T P n_attractors")
VertexModelImpactExperimentParameters = namedtuple("VertexModelImpactExperimentParameters", "G current_attractors T P "
                                                                                  "impact_types relative_basins "
                                                                                  "maximal_bits "
                                                                                  "impacts")
VertexStateImpactExperimentParameters = namedtuple("VertexStateImpactExperimentParameters", "G current_attractors "
                                                                                  "relative_basins "
                                                                                  "max_transient_len "
                                                                                  "impacts")
StochasticVertexModelImpactExperimentParameters = namedtuple(
    "StochasticVertexModelImpactExperimentParameters", "G current_attractors "
    "bits_of_change relative_basins impact_type impacts")

StochasticVertexStateImpactExperimentParameters = namedtuple(
    "StochasticVertexStateImpactExperimentParameters", "G impacts")

GraphModelImpactExperimentParameters = namedtuple("GraphModelImpactExperimentParameters", "G current_attractors T P "
                                                                                  "impact_types relative_basins "
                                                                                  "maximal_bits "
                                                                                  "impact")
GraphStateImpactExperimentParameters = namedtuple("GraphStateImpactExperimentParameters", "G current_attractors "
                                                                                  "relative_basins "
                                                                                  "max_transient_len maximal_bits "
                                                                                  "impact")
StochasticGraphModelImpactExperimentParameters = namedtuple(
    "StochasticGraphModelImpactExperimentParameters", "G current_attractors "
    "bits_of_change relative_basins impact_type impact")

StochasticGraphStateImpactExperimentParameters = namedtuple(
    "StochasticGraphStateImpactExperimentParameters", "G bits_of_change impact")

DubrovaExperimentParameters = namedtuple("DubrovaExperimentParameters", "G mutate n_attractors")


class TestAttractors(TestCase):
    def test_num_attractors_onestage(self):
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
        # G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
        #        "\\Attractors - for Ariel\\BNS_Dubrova_2011\\MAPK_large2.cnet")
        # experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=15, n_attractors=12))
        # experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=15, n_attractors=14))
        # experiments.append(ILPAttractorExperimentParameters(G=G, T=3, P=15, n_attractors=14))
        G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
               "\\Attractors - for Ariel\\BNS_Dubrova_2011\\tcr.cnet")
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=15, n_attractors=8))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=6, P=15, n_attractors=9))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=7, P=15, n_attractors=9))

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
        #     print "experiment #{}".format(i)h
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

    def test_vertex_degeneracy_scores(self):
        self.assertTrue(False)  # TODO: write...

    def test_graph_state_impact_scores(self):
        experiments = []

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #0
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=1,
                                                                 maximal_bits=1,
                                                                 impact=0))
        # experiment #1
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=0,
                                                                 maximal_bits=1,
                                                                 impact=0))
        # experiment #2
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=30,
                                                                 maximal_bits=1,
                                                                 impact=0))
        # experiment #3
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=30,
                                                                 maximal_bits=10,
                                                                 impact=0))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand, None])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #4
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=1,
                                                                 maximal_bits=1,
                                                                 impact=0))
        # experiment #5
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=0,
                                                                 maximal_bits=1,
                                                                 impact=0))
        # experiment #6
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=30,
                                                                 maximal_bits=1,
                                                                 impact=0))
        # experiment #7
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=30,
                                                                 maximal_bits=10,
                                                                 impact=0))

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #8
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=0,
                                                                 maximal_bits=1,
                                                                 impact=1))
        # experiment #9
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=5,
                                                                 maximal_bits=1,
                                                                 impact=1))
        # experiment #10
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=5,
                                                                maximal_bits=5,
                                                                impact=1))
        # experiment #11
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=[0.1, 0.9],
                                                                max_transient_len=5,
                                                                maximal_bits=5,
                                                                impact=1))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[sympy.And, None])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #12
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=[0.1, 0.9],
                                                                max_transient_len=5,
                                                                maximal_bits=5,
                                                                impact=1))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.Nand])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #13
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=0,
                                                                maximal_bits=1,
                                                                impact=1))
        # experiment #14
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=[0.1, 0.9],
                                                                max_transient_len=5,
                                                                maximal_bits=5,
                                                                impact=1))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #15
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=0,
                                                                 impacts=[1, 1, 1]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, lambda _: True])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #16
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=0,
                                                                maximal_bits=1,
                                                                impact=0))
        # experiment #17
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=0,
                                                                maximal_bits=3,
                                                                impact=0))
        # experiment #18
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=5,
                                                                maximal_bits=2,
                                                                impact=0))

        G = graphs.Network(vertex_names=["A", "B", "C", "D"], edges=[("B", "A"), ("C", "A"), ("D", "A"),
                                                                     ("A", "B"), ("C", "B"), ("D", "B"),
                                                                     ("A", "C"), ("B", "C"), ("D", "C"),
                                                                     ("A", "D"), ("B", "D"), ("C", "D")],
                           vertex_functions=[lambda a, b, c: a + b + c > 1 for _ in range(4)])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # 0000 and 1111 are stable points, and attract everything with hamming distance <= 1,
        # where 2 bits of change land right into another attractor.
        # Other three two-state attractors are unstable under one bit change, with transient length of 1,
        # Or they can be switched between eachother/stables with 2 (same as 0000/1111 ones, if needed)
        # bits of change.
        # experiment #19
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=0,
                                                                maximal_bits=1,
                                                                impact=0))
        # experiment #20
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=1,
                                                                maximal_bits=1,
                                                                impact=3 / 5.0))
        # experiment #21
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=5,
                                                                maximal_bits=1,
                                                                impact=3 / 5.0))
        # experiment #22
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=0,
                                                                maximal_bits=2,
                                                                impact=1))
        # experiment #23
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=3,
                                                                maximal_bits=2,
                                                                impact=1))

        relative_basins = [5 / float(16) if len(attractor) == 1 else 2 / float(16) for
                           attractor in current_attractors]
        # experiment #24
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=relative_basins,
                                                                max_transient_len=5,
                                                                maximal_bits=1,
                                                                impact=6 / 16.0))
        # experiment #25
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=relative_basins,
                                                                max_transient_len=0,
                                                                maximal_bits=2,
                                                                impact=1))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "C")],
                           vertex_functions=[None, sympy.And, sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #19
        # 000, 110 and 111 are the steady states. First is stable, other can change on
        # right vertex change, B with one step and C immediately.
        # experiment #26
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=0,
                                                                maximal_bits=1,
                                                                impact=2 / 3.0))
        # experiment #27
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=0,
                                                                maximal_bits=2,
                                                                impact=2 / 3.0))
        # experiment #28
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=5,
                                                                maximal_bits=5,
                                                                impact=2 / 3.0))

        relative_len_decider = lambda attractor: 0.5 if [
                                int(s) for s in attractor[0]] == [0, 0, 0] else 3 / float(8) if [
                                int(s) for s in attractor[0]] == [1, 1, 0] else 1 / float(8)
        relative_basins = [relative_len_decider(att) for att in current_attractors]
        # experiment #29
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=relative_basins,
                                                                max_transient_len=5,
                                                                maximal_bits=2,
                                                                impact=0.5))
        # experiment #30
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=relative_basins,
                                                                max_transient_len=0,
                                                                maximal_bits=1,
                                                                impact=0.5))

        G = graphs.Network(vertex_names=["A", "B", "C", "D"], edges=[("A", "B"), ("B", "C"), ("C", "D"),
                                                                     ("D", "D")],
                           vertex_functions=[None, sympy.And, sympy.And, sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # Now 0000 is stable, 1110 changes immediently on last vertex change, 1111 can change in 2, 1, or 0
        # steps on change of second, third or last vertex.
        # experiment #31
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=0,
                                                                maximal_bits=1,
                                                                impact=2 / 3.0))
        # experiment #31
        experiments.append(GraphStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                max_transient_len=3,
                                                                maximal_bits=3,
                                                                impact=2 / 3.0))

        print "number of experiments (with keys)={}".format(len(experiments))
        for i, experiment in enumerate(experiments):
            print "experiment #{}".format(i)
            print "attractors:"
            print experiment.current_attractors
            print "n={}, relative_basins={}, expected_impacts={}".\
                format(len(experiment.G.vertices), experiment.relative_basins, experiment.impact)
            impact = graph_state_impact_score(G=experiment.G, current_attractors=experiment.current_attractors,
                                                 max_transient_len=experiment.max_transient_len,
                                                 relative_attractor_basin_sizes=experiment.relative_basins,
                                                 key_slice_size=15, maximal_bits_of_change=experiment.maximal_bits)

            # (from vertex version) got numeric problems with test #16 regardless of key_slice
            impact = round(impact, 5)
            experiment_impact = round(experiment.impact, 5)
            print "expected impact:"
            print impact
            print "got impact:"
            print experiment_impact
            try:
                self.assertEqual(impact, experiment_impact)
            except AssertionError as e:
                print e
                print experiment.G
                raise e


    def test_vertex_state_impact_scores(self):
        # TODO: test stochastic kind
        experiments = []

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #0
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=1,
                                                                 impacts=[0]))
        # experiment #1
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=0,
                                                                 impacts=[0]))
        # experiment #2
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=30,
                                                                 impacts=[0]))

        # experiment #3
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=[1],
                                                                 max_transient_len=30,
                                                                 impacts=[0]))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand, None])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #4
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=30,
                                                                 impacts=[0, np.nan]))

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #5
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=30,
                                                                 impacts=[1]))
        # experiment #6
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=30,
                                                                 impacts=[1]))
        # experiment #7
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=30,
                                                                 impacts=[1]))
        # experiment #8
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=[0.1, 0.9],
                                                                 max_transient_len=1,
                                                                 impacts=[1]))
        # experiment #9
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=[0.1, 0.9],
                                                                 max_transient_len=0,
                                                                 impacts=[1]))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[sympy.And, None])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #10
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=0,
                                                                 impacts=[1, np.nan]))
        # experiment #11
        experiments.append(VertexStateImpactExperimentParameters(G=G,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=[0.1, 0.4, 0.4, 0.1],
                                                                 max_transient_len=0,
                                                                 impacts=[1, np.nan]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.Nand])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #12
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=0,
                                                                 impacts=[1] * 3))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #13
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=0,
                                                                 impacts=[1, 1, 1]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, lambda _: True])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #14
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=0,
                                                                 impacts=[0, 0, 0]))

        G = graphs.Network(vertex_names=["A", "B", "C", "D"], edges=[("B", "A"), ("C", "A"), ("D", "A"),
                                                                     ("A", "B"), ("C", "B"), ("D", "B"),
                                                                     ("A", "C"), ("B", "C"), ("D", "C"),
                                                                     ("A", "D"), ("B", "D"), ("C", "D")],
                           vertex_functions=[lambda a, b, c: a + b + c > 1 for _ in range(4)])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #15
        # 0000 and 1111 are stable points, and attract everything with hamming distance <= 1.
        # Other three two-state attractors are unstable under one bit change, with transient length of 1.
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=0,
                                                                 impacts=[0] * 4))
        # experiment #16
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=1,
                                                                 impacts=[3 / 5.0] * 4))
        # experiment #17
        relative_basins = [5 / float(16) if len(attractor) == 1 else 2 / float(16) for
                           attractor in current_attractors]
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=relative_basins,
                                                                 max_transient_len=1,
                                                                 impacts=[6 / 16.0, 6 / 16.0,
                                                                          6 / 16.0, 6 / 16.0]))
        # experiment #18
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=relative_basins,
                                                                 max_transient_len=2,
                                                                 impacts=[6 / 16.0, 6 / 16.0,
                                                                          6 / 16.0, 6 / 16.0]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "C")],
                           vertex_functions=[None, sympy.And, sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #19
        # 000, 110 and 111 are the steady states. First is stable, other can change on
        # right vertex change, B with one step and C immediately.
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=0,
                                                                 impacts=[np.nan, 0, 2 / 3.0]))
        # experiment #20
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=1,
                                                                 impacts=[np.nan, 1 / 3.0, 2/ 3.0]))
        relative_len_decider = lambda attractor: 0.5 if [
                                int(s) for s in attractor[0]] == [0, 0, 0] else 3 / float(8) if [
                                int(s) for s in attractor[0]] == [1, 1, 0] else 1 / float(8)
        relative_basins = [relative_len_decider(att) for att in current_attractors]
        # experiment #21
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=relative_basins,
                                                                 max_transient_len=1,
                                                                 impacts=[np.nan, 1 / 8.0, 0.5]))

        G = graphs.Network(vertex_names=["A", "B", "C", "D"], edges=[("A", "B"), ("B", "C"), ("C", "D"),
                                                                     ("D", "D")],
                           vertex_functions=[None, sympy.And, sympy.And, sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # Now 0000 is stable, 1110 changes immediently on last vertex change, 1111 can change in 2, 1, or 0
        # steps on change of second, third or last vertex.
        # experiment #22
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=0,
                                                                 impacts=[np.nan, 0, 0, 2 / float(3)]))
        # experiment #23
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=1,
                                                                 impacts=[np.nan, 0, 1 / float(3),
                                                                          2 / float(3)]))
        # experiment #24
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=2,
                                                                 impacts=[np.nan, 1 / float(3), 1 / float(3),
                                                                          2 / float(3)]))
        # experiment #25
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=3,
                                                                 impacts=[np.nan, 1 / float(3), 1 / float(3),
                                                                          2 / float(3)]))
        # experiment #26
        experiments.append(VertexStateImpactExperimentParameters(G=G, current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 max_transient_len=30,
                                                                 impacts=[np.nan, 1 / float(3), 1 / float(3),
                                                                          2 / float(3)]))

        print "number of experiments (with keys)={}".format(len(experiments))
        for i, experiment in enumerate(experiments):
            print "experiment #{}".format(i)
            print "attractors:"
            print experiment.current_attractors
            print "n={}, relative_basins={}, expected_impacts={}".\
                format(len(experiment.G.vertices), experiment.relative_basins, experiment.impacts)
            impacts = vertex_state_impact_scores(G=experiment.G, current_attractors=experiment.current_attractors,
                                                 max_transient_len=experiment.max_transient_len,
                                                 relative_attractor_basin_sizes=experiment.relative_basins,
                                                 key_slice_size=15)

            # got numeric problems with test #16 regardless of key_slice
            impacts = [round(x, 5) if not np.isnan(x) else x for x in impacts]
            experiment_impacts = [round(x, 5) if not np.isnan(x) else x for x in experiment.impacts]
            print "expected impacts:"
            print impacts
            print "got impacts:"
            print experiment_impacts
            try:
                self.assertEqual(impacts, experiment_impacts)
            except AssertionError as e:
                print e
                print experiment.G
                raise e

    def test_graph_model_impact_scores(self):
        # TODO: also test the resulting models (assure they have the correct number of attractors)
        experiments = []

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #0
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Invalidation,
                                                                maximal_bits=1,
                                                                current_attractors=current_attractors,
                                                                relative_basins=None,
                                                                impact=1))
        # experiment #1
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Both,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact=1))
        # experiment #2
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact=2))
        # experiment #3
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Both,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact=1.5))
        # experiment #4
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=3, P=1, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact=1))
        # experiment #5
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Both,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=[1],
                                                                 impact=1.5))

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #6
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=1, P=1, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact=0.5))
        # experiment #7
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=1, P=1, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=[0.1, 0.9],
                                                                 impact=0.9))
        # experiment #8
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=1, P=1, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact=1))
        # experiment #9
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Both,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=[0.1, 0.9],
                                                                 impact=0.75))
        # experiment #10
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact=0.5))
        # experiment #11
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=1, P=1, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact=0))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[sympy.And, None])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #12
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=1, P=1, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact=1))
        # experiment #13
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Both,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=[0.1, 0.4, 0.4, 0.1],
                                                                 impact=0.75))
        # experiment #14
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact=0.5))
        # experiment #15
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=3, P=1, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact=0.25))
        # experiment #16
        experiments.append(GraphModelImpactExperimentParameters(G=G, T=1, P=1, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact=0))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.Nand])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #17
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1] * 3))
        # experiment #18
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1] * 3))
        # experiment #19
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=6, P=5, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[2] * 3))
        # experiment #20
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=6, P=3, impact_types=ImpactType.Both,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=[0.1, 0.9],
                                                                 impacts=[1.25] * 3))
        # experiment #21
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=6, P=5, impact_types=ImpactType.Both,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=[0.1, 0.9],
                                                                 impacts=[1.5] * 3))
        # experiment #22
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=6, P=2, impact_types=ImpactType.Addition,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0.5] * 3))
        # experiment #23
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=1, P=1, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0.5] * 3))
        # experiment #24
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=1, P=5, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1] * 3))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #25
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0.75, 0.75, 0.75]))
        # experiment #26
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1, 1, 1]))
        # experiment #27
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0.5, 0.5, 0.5]))
        # experiment #28
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Both,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0.75, 0.75, 0.75]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, lambda _: True])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #29
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1, 1, 1]))
        # experiment #30
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Addition,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1, 1, 3]))
        # experiment #31
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1, 1, 3]))

        print "number of experiments (with keys)={}".format(len(experiments))
        for i, experiment in enumerate(experiments):
            print "experiment #{}".format(i)
            print "n={}, T={}, P={}, maximal_bits={}, relative_basins={}, expected_impacts={}".\
                format(len(experiment.G.vertices),
                       experiment.T, experiment.P, experiment.maximal_bits, experiment.relative_basins,
                       experiment.impacts)
            print experiment.current_attractors
            impacts = vertex_model_impact_scores(G=experiment.G, current_attractors=experiment.current_attractors,
                                                 max_len=experiment.T,
                                                 max_num=experiment.P,
                                                 impact_types=experiment.impact_types,
                                                 relative_attractor_basin_sizes=experiment.relative_basins,
                                                 maximal_bits_of_change=experiment.maximal_bits)
            try:
                self.assertEqual(impacts, experiment.impacts)
            except AssertionError as e:
                print e
                print experiment.G
                raise e


    def test_vertex_model_impact_scores(self):
        # TODO: also test the resulting models (assure they have the correct number of attractors)
        # TODO: test stochastic kind
        experiments = []

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #0
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1]))
        # experiment #1
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Both,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1]))
        # experiment #2
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[2]))
        # experiment #3
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Both,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1.5]))
        # experiment #4
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=3, P=1, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1]))
        # experiment #5
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Both,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=[1],
                                                                 impacts=[1.5]))

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #6
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=1, P=1, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0.5]))
        # experiment #7
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=1, P=1, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=[0.1, 0.9],
                                                                 impacts=[0.9]))
        # experiment #8
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=1, P=1, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1]))
        # experiment #9
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Both,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=[0.1, 0.9],
                                                                 impacts=[0.75]))
        # experiment #10
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0.5]))
        # experiment #11
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=1, P=1, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0]))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[sympy.And, None])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #12
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=1, P=1, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1, np.nan]))
        # experiment #13
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Both,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=[0.1, 0.4, 0.4, 0.1],
                                                                 impacts=[0.75, np.nan]))
        # experiment #14
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=3, P=3, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0.5, np.nan]))
        # experiment #15
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=3, P=1, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0.25, np.nan]))
        # experiment #16
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=1, P=1, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0, np.nan]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.Nand])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #17
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1] * 3))
        # experiment #18
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1] * 3))
        # experiment #19
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=6, P=5, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[2] * 3))
        # experiment #20
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=6, P=3, impact_types=ImpactType.Both,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=[0.1, 0.9],
                                                                 impacts=[1.25] * 3))
        # experiment #21
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=6, P=5, impact_types=ImpactType.Both,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=[0.1, 0.9],
                                                                 impacts=[1.5] * 3))
        # experiment #22
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=6, P=2, impact_types=ImpactType.Addition,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0.5] * 3))
        # experiment #23
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=1, P=1, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0.5] * 3))
        # experiment #24
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=1, P=5, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1] * 3))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #25
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0.75, 0.75, 0.75]))
        # experiment #26
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1, 1, 1]))
        # experiment #27
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0.5, 0.5, 0.5]))
        # experiment #28
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Both,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[0.75, 0.75, 0.75]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, lambda _: True])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #29
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Invalidation,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1, 1, 1]))
        # experiment #30
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Addition,
                                                                 maximal_bits=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1, 1, 3]))
        # experiment #31
        experiments.append(VertexModelImpactExperimentParameters(G=G, T=7, P=5, impact_types=ImpactType.Addition,
                                                                 maximal_bits=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impacts=[1, 1, 3]))

        print "number of experiments (with keys)={}".format(len(experiments))
        for i, experiment in enumerate(experiments):
            print "experiment #{}".format(i)
            print "n={}, T={}, P={}, maximal_bits={}, relative_basins={}, expected_impacts={}".\
                format(len(experiment.G.vertices),
                       experiment.T, experiment.P, experiment.maximal_bits, experiment.relative_basins,
                       experiment.impacts)
            print experiment.current_attractors
            impacts = vertex_model_impact_scores(G=experiment.G, current_attractors=experiment.current_attractors,
                                                 max_len=experiment.T,
                                                 max_num=experiment.P,
                                                 impact_types=experiment.impact_types,
                                                 relative_attractor_basin_sizes=experiment.relative_basins,
                                                 maximal_bits_of_change=experiment.maximal_bits)
            try:
                self.assertEqual(impacts, experiment.impacts)
            except AssertionError as e:
                print e
                print experiment.G
                raise e

    def test_stochastic_graph_state_impact_scores(self):
        experiments = []

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand])
        # experiment #0
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=1, impact=0))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand, None])
        # experiment #1
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=1, impact=0))

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.And])
        # experiment #2
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=1, impact=1))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[sympy.And, None])
        # experiment #3
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=1, impact=0))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.Nand])
        # experiment #4
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=1, impact=0.5))
        # experiment #5
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=2, impact=0.5))
        # experiment #6
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=3, impact=0))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.And])
        # experiment #7
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=1, impact=1))
        # experiment #8
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=2, impact=0.5))
        # experiment #9
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=3, impact=1))


        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, lambda _: True])
        # experiment #10
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=1, impact=0))
        # experiment #11
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=2, impact=0))
        # experiment #12
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=3, impact=0))

        G = graphs.Network(vertex_names=["A", "B", "C", "D"], edges=[("B", "A"), ("C", "A"), ("D", "A"),
                                                                     ("A", "B"), ("C", "B"), ("D", "B"),
                                                                     ("A", "C"), ("B", "C"), ("D", "C"),
                                                                     ("A", "D"), ("B", "D"), ("C", "D")],
                           vertex_functions=[lambda a, b, c: a + b + c > 1 for _ in range(4)])
        # experiment #12
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=1, impact=3 / 8.0))
        # experiment #13
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=2, impact=1))
        # experiment #14
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=3, impact=1))
        # experiment #15
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=4, impact=0))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "C")],
                           vertex_functions=[None, sympy.And, sympy.And])
        # 000, 110 and 111 are the steady states. First is stable, other can change on
        # right vertex change, B with one step and C immediately.
        # experiment #16
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=1, impact=5 / 16.0))
        # experiment #17
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=2, impact=1 / 16.0))

        G = graphs.Network(vertex_names=["A", "B", "C", "D"], edges=[("A", "B"), ("B", "C"), ("C", "D"),
                                                                     ("D", "D")],
                           vertex_functions=[None, sympy.And, sympy.And, sympy.And])
        # Now 0000 is stable, 1110 changes immediently on last vertex change, 1111 can change in 2, 1, or 0
        # steps on change of second, third or last vertex.
        # experiment #18
        experiments.append(StochasticGraphStateImpactExperimentParameters(G=G, bits_of_change=1,
                                                                          impact=0.20833333333))

        print "number of experiments (with keys)={}".format(len(experiments))
        for i, experiment in enumerate(experiments):
            print "experiment #{}".format(i)
            print "n={}, expected_impacts={}".\
                format(len(experiment.G.vertices), experiment.impacts)

            for iteration in range(10):
                n_iter = random.randint(700, 1400)
                estimated_impacts = stochastic_vertex_state_impact_scores(G=experiment.G, n_iter=n_iter)
                print "estimated_impacts={}".format(estimated_impacts)
                self.assertTrue(len(experiment.impacts) == len(estimated_impacts))
                for calculated_impact, estimated_impact in zip(experiment.impacts, estimated_impacts):
                    if np.isnan(calculated_impact):
                        self.assertTrue(np.isnan(estimated_impact))
                    else:
                        self.assertTrue(abs(estimated_impact - calculated_impact) < 0.1)


    def test_stochastic_vertex_state_impact_scores(self):
        experiments = []

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand])
        # experiment #0
        experiments.append(StochasticVertexStateImpactExperimentParameters(G=G, impacts=[0]))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand, None])
        # experiment #1
        experiments.append(StochasticVertexStateImpactExperimentParameters(G=G, impacts=[0, np.nan]))

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.And])
        # experiment #2
        experiments.append(StochasticVertexStateImpactExperimentParameters(G=G, impacts=[1]))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[sympy.And, None])
        # experiment #3
        experiments.append(StochasticVertexStateImpactExperimentParameters(G=G, impacts=[1, np.nan]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.Nand])
        # experiment #4
        experiments.append(StochasticVertexStateImpactExperimentParameters(G=G, impacts=[0.5] * 3))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.And])
        # experiment #5
        experiments.append(StochasticVertexStateImpactExperimentParameters(G=G, impacts=[1, 1, 1]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, lambda _: True])
        # experiment #6
        experiments.append(StochasticVertexStateImpactExperimentParameters(G=G, impacts=[0, 0, 0]))

        G = graphs.Network(vertex_names=["A", "B", "C", "D"], edges=[("B", "A"), ("C", "A"), ("D", "A"),
                                                                     ("A", "B"), ("C", "B"), ("D", "B"),
                                                                     ("A", "C"), ("B", "C"), ("D", "C"),
                                                                     ("A", "D"), ("B", "D"), ("C", "D")],
                           vertex_functions=[lambda a, b, c: a + b + c > 1 for _ in range(4)])
        # experiment #7
        experiments.append(StochasticVertexStateImpactExperimentParameters(G=G, impacts=[3 / 8.0] * 4))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "C")],
                           vertex_functions=[None, sympy.And, sympy.And])
        # experiment #8
        # 000, 110 and 111 are the steady states. First is stable, other can change on
        # right vertex change, B with one step and C immediately.
        experiments.append(StochasticVertexStateImpactExperimentParameters(G=G,
                                                                           impacts=[np.nan, 1/8.0, 0.5]))
        G = graphs.Network(vertex_names=["A", "B", "C", "D"], edges=[("A", "B"), ("B", "C"), ("C", "D"),
                                                                     ("D", "D")],
                           vertex_functions=[None, sympy.And, sympy.And, sympy.And])
        # Now 0000 is stable, 1110 changes immediently on last vertex change, 1111 can change in 2, 1, or 0
        # steps on change of second, third or last vertex.
        # experiment #9
        experiments.append(StochasticVertexStateImpactExperimentParameters(G=G,
                                                                           impacts=[np.nan,
                                                                                    1/16.0, 1/16.0,
                                                                                    0.5]))

        print "number of experiments (with keys)={}".format(len(experiments))
        for i, experiment in enumerate(experiments):
            print "experiment #{}".format(i)
            print "n={}, expected_impacts={}".\
                format(len(experiment.G.vertices), experiment.impacts)

            for iteration in range(10):
                n_iter = random.randint(700, 1400)
                estimated_impacts = stochastic_vertex_state_impact_scores(G=experiment.G, n_iter=n_iter)
                print "estimated_impacts={}".format(estimated_impacts)
                self.assertTrue(len(experiment.impacts) == len(estimated_impacts))
                for calculated_impact, estimated_impact in zip(experiment.impacts, estimated_impacts):
                    if np.isnan(calculated_impact):
                        self.assertTrue(np.isnan(estimated_impact))
                    else:
                        self.assertTrue(abs(estimated_impact - calculated_impact) < 0.1)

    def test_stochastic_vertex_model_impact_scores(self):
        # TODO: also test the resulting models (assure they have the correct number of attractors)
        experiments = []

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.Nand])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #0
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G,
                                                                 bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Invalidation,
                                                                 impacts=[1]))
        # experiment #1
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G,
                                                                 bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Invalidation,
                                                                 impacts=[1]))
        # experiment #2
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G,
                                                                 bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Addition,
                                                                 impacts=[1]))
        # experiment #3
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G,
                                                                 bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Addition,
                                                                 impacts=[2]))
        # experiment #4
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G,
                                                                 bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Both,
                                                                 impacts=[1.5]))

        G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
                           vertex_functions=[sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #5
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Invalidation,
                                                                 impacts=[0.5]))
        # experiment #6
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=[0.1, 0.9],
                                                                 impact_type=ImpactType.Invalidation,
                                                                 impacts=[0.5]))
        # experiment #7
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Invalidation,
                                                                 impacts=[1]))
        # experiment #8
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Addition,
                                                                 impacts=[0]))
        # experiment #9
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Addition,
                                                                 impacts=[0.5]))
        # experiment #10
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Both,
                                                                 impacts=[0.75]))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "A")],
                           vertex_functions=[sympy.And, None])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #11
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Invalidation,
                                                                 impacts=[0.5, np.nan]))
        # experiment #12
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Invalidation,
                                                                 impacts=[1, np.nan]))
        # experiment #13
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Addition,
                                                                 impacts=[0.5, np.nan]))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.Nand])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #14
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G,bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Invalidation,
                                                                 impacts=[1] * 3))
        # experiment #15
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G,bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Invalidation,
                                                                 impacts=[1] * 3))
        # experiment #16
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G,bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Addition,
                                                                 impacts=[0.5] * 3))
        # experiment #17
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G,bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Addition,
                                                                 impacts=[2] * 3))
        # experiment #18
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G,bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Both,
                                                                 impacts=[0.75] * 3))
        # experiment #19
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G,bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Both,
                                                                 impacts=[1.5] * 3))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #20
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Invalidation,
                                                                 impacts=[3 / 4.0] * 3))
        # experiment #21
        basin_sizes = [3 / 8.0 if len(att) > 1 else 1 / 8.0 for att in current_attractors]
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=basin_sizes,
                                                                 impact_type=ImpactType.Invalidation,
                                                                 impacts=[7 / 8.0] * 3))
        # experiment #22
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Invalidation,
                                                                 impacts=[1, 1, 1]))
        # experiment #23
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Addition,
                                                                 impacts=[0] * 3))
        # experiment #24
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Addition,
                                                                 impacts=[0.5, 0.5, 0.5]))
        # experiment #25
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Both,
                                                                 impacts=[7 / 16.0] * 3))
        # experiment #26
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Both,
                                                                 impacts=[0.75] * 3))

        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
                           vertex_functions=[sympy.Nand, sympy.Nand, lambda _: True])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #27
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Invalidation,
                                                                 impacts=[0.5] * 3))
        # experiment #28
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Invalidation,
                                                                 impacts=[1] * 3))
        # experiment #29
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Addition,
                                                                 impacts=[0.5, 0.5, 2.5]))
        # experiment #30
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=2,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Addition,
                                                                 impacts=[1, 1, 1]))
        # experiment #31
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=1,
                                                                 current_attractors=current_attractors,
                                                                 relative_basins=None,
                                                                 impact_type=ImpactType.Both,
                                                                 impacts=[0.5, 0.5, 1.5]))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A"), ("B", "B")],
                           vertex_functions=[sympy.And, sympy.And])
        current_attractors = find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=True)
        # experiment #32
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=1,
                                                                           current_attractors=current_attractors,
                                                                           relative_basins=None,
                                                                           impact_type=ImpactType.Invalidation,
                                                                           impacts=[0.5, 0.25]))
        # experiment #33
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=1,
                                                                           current_attractors=current_attractors,
                                                                           relative_basins=[0.1, 0.9],
                                                                           impact_type=ImpactType.Invalidation,
                                                                           impacts=[0.5, 0.25]))
        # experiment #34
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=2,
                                                                           current_attractors=current_attractors,
                                                                           relative_basins=None,
                                                                           impact_type=ImpactType.Invalidation,
                                                                           impacts=[1, 0.5]))
        # experiment #35
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=1,
                                                                           current_attractors=current_attractors,
                                                                           relative_basins=None,
                                                                           impact_type=ImpactType.Addition,
                                                                           impacts=[0.25, 1 / 8.0]))
        # experiment #36
        experiments.append(StochasticVertexModelImpactExperimentParameters(G=G, bits_of_change=2,
                                                                           current_attractors=current_attractors,
                                                                           relative_basins=None,
                                                                           impact_type=ImpactType.Addition,
                                                                           impacts=[0.5, 1 / 4.0]))

        print "number of experiments (with keys)={}".format(len(experiments))
        for i, experiment in enumerate(experiments):
            print "experiment #{}".format(i)
            print "n={}, bits_of_change={}, relative_basins={}, impact_type={}, expected_impacts={}".\
                format(len(experiment.G.vertices),
                       experiment.bits_of_change, experiment.relative_basins, experiment.impact_type,
                       experiment.impacts)
            print experiment.current_attractors

            for use_dubrova in [False, True]:
                n_iter = random.randint(400, 440)
                attractor_estimation_n_iter = random.randint(30, 35)

                estimated_impacts = stochastic_vertex_model_impact_scores(
                    G=experiment.G, current_attractors=experiment.current_attractors, n_iter=n_iter, use_dubrova=use_dubrova,
                    bits_of_change=experiment.bits_of_change,
                    relative_attractor_basin_sizes=experiment.relative_basins,
                    attractor_estimation_n_iter=attractor_estimation_n_iter,
                    impact_type=experiment.impact_type,
                    cur_dubrova_path=dubrova_path)

                self.assertTrue(len(experiment.impacts) == len(estimated_impacts))
                print "estimated_impacts={}".format(estimated_impacts)
                for calculated_impact, estimated_impact in zip(experiment.impacts, estimated_impacts):
                    if np.isnan(calculated_impact):
                        self.assertTrue(np.isnan(estimated_impact))
                    else:
                        self.assertTrue(abs(estimated_impact - calculated_impact) < 0.15)

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
               "\\Attractors - for Ariel\\BNS_Dubrova_2011\\tcr.cnet")
        self.assertEqual(find_num_steady_states(G, verbose=False, simplify_general_boolean=False), 8)

    def test_find_attractors_dubrova(self):
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
                                             True])
        experiments.append(DubrovaExperimentParameters(G=G, mutate=False, n_attractors=3))
        experiments.append(DubrovaExperimentParameters(G=G, mutate=True, n_attractors=3))
        # 11, 12
        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "A"), ("C", "A")],
                           vertex_functions=[logic.SymmetricThresholdFunction.from_function(sympy.Nand, 2),
                                             logic.SymmetricThresholdFunction.from_function(sympy.Nand, 1),
                                             False])
        experiments.append(DubrovaExperimentParameters(G=G, mutate=False, n_attractors=1))
        experiments.append(DubrovaExperimentParameters(G=G, mutate=True, n_attractors=1))

        # 13, 14
        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "A"), ("C", "A")],
                           vertex_functions=[logic.SymmetricThresholdFunction.from_function(sympy.Nand, 2),
                                             logic.SymmetricThresholdFunction.from_function(sympy.Nand, 1),
                                             None])
        experiments.append(DubrovaExperimentParameters(G=G, mutate=False, n_attractors=1))
        experiments.append(DubrovaExperimentParameters(G=G, mutate=True, n_attractors=4))

        # 15
        G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
               "\\Attractors - for Ariel\\BNS_Dubrova_2011\\tcr.cnet")
        # G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
        #        "\\Attractors - for Ariel\\BNS_Dubrova_2011\\MAPK_large.cnet")
        experiments.append(DubrovaExperimentParameters(G=G, mutate=False, n_attractors=9))

        print "number of experiments (with keys)={}".format(len(experiments))
        for i, experiment in enumerate(experiments):
            print "experiment #{}".format(i)
            print "n={}, mutate={}, expected_n_attractors={}".format(len(experiment.G.vertices),
                                                                     experiment.mutate, experiment.n_attractors)
            # continue
            attractors = find_attractors_dubrova(G=experiment.G,
                                                       dubrova_path="../bns_dubrova.exe",
                                                       mutate_input_nodes=experiment.mutate)
            n_attractors = len(attractors)
            try:
                self.assertEqual(n_attractors, experiment.n_attractors)
            except AssertionError as e:
                print e
                print experiment.G
                raise e
            except Exception as e:
                raise e

        print "testing state order in attractor"
        # TODO: expand? random graphs, compare ILP attractors with Dubrova's
        G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[sympy.And, sympy.Nand, True])
        desired_attractor = [[0, 0, 1], [0, 1, 1], [1, 1, 1], [1, 0, 1]]
        # repeat manually, (otherwise there's mutual dependence of tests).
        possible_attractors = [desired_attractor[shift:] + desired_attractor[:shift] for shift in range(4)]
        # print possible_attractors
        found_attractors = find_attractors_dubrova(G, dubrova_path="../bns_dubrova.exe", mutate_input_nodes=True)
        self.assertTrue(len(found_attractors) == 1)
        found_attractor = [[int(v) for v in state] for state in found_attractors[0]]
        # print found_attractor
        self.assertTrue(any(found_attractor == possible_attractors[i] for i in range(len(possible_attractors))))

        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[sympy.And, sympy.Nand])
        desired_attractor = [[0, 0], [0, 1], [1, 1], [1, 0]]
        # repeat manually, (otherwise there's mutual dependence of tests).
        possible_attractors = [desired_attractor[shift:] + desired_attractor[:shift] for shift in range(4)]
        # print possible_attractors
        found_attractors = find_attractors_dubrova(G, dubrova_path="../bns_dubrova.exe", mutate_input_nodes=True)
        self.assertTrue(len(found_attractors) == 1)
        found_attractor = [[int(v) for v in state] for state in found_attractors[0]]
        # print found_attractor
        self.assertTrue(any(found_attractor == possible_attractor for possible_attractor in possible_attractors))

    def test_find_attractors_enumerate(self):
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
        # 30, 31, 32, 33
        G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
               "\\Attractors - for Ariel\\BNS_Dubrova_2011\\tcr.cnet")
        experiments.append(ILPAttractorExperimentParameters(G=G, T=1, P=None, n_attractors=8))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=2, P=None, n_attractors=8))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=6, P=None, n_attractors=9))
        experiments.append(ILPAttractorExperimentParameters(G=G, T=8, P=None, n_attractors=9))

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
