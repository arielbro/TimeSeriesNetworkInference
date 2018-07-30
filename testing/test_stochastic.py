from unittest import TestCase
import stochastic
import graphs
import sympy
import utility
import random
import attractors


class TestStochastic(TestCase):

    def test_walk_to_attractor(self):
        # hard coded graphs with known attractors
        G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                           vertex_functions=[sympy.And, sympy.And])
        attractor = stochastic.walk_to_attractor(G, initial_state=(0, 0), max_walk=5, state_to_attractor_mapping=None)
        self.assertTrue(utility.is_same_attractor(attractor, [[0, 0]]))
        attractor = stochastic.walk_to_attractor(G, initial_state=(1, 0), max_walk=5, state_to_attractor_mapping=None)
        self.assertTrue(utility.is_same_attractor(attractor, [[1, 0], [0, 1]]))
        attractor = stochastic.walk_to_attractor(G, initial_state=(0, 1), max_walk=5, state_to_attractor_mapping=None)
        self.assertTrue(utility.is_same_attractor(attractor, [[1, 0], [0, 1]]))

        mapping = dict()
        attractor = stochastic.walk_to_attractor(G, initial_state=(0, 1), max_walk=5, state_to_attractor_mapping=mapping)
        self.assertTrue(utility.is_same_attractor(attractor, [[1, 0], [0, 1]]))
        unique_mapping_values = set(tuple(att) for att in mapping.values())
        self.assertTrue(utility.attractor_sets_equality(unique_mapping_values, {attractor}))
        attractor = stochastic.walk_to_attractor(G, initial_state=(0, 0), max_walk=5,
                                                 state_to_attractor_mapping=mapping)
        self.assertTrue(utility.is_same_attractor(attractor, [[0, 0]]))
        unique_mapping_values = set(tuple(att) for att in mapping.values())
        self.assertTrue(utility.attractor_sets_equality(unique_mapping_values, {[[0, 0]], [[0, 1], [1, 0]]}))

        # random graphs, assure what's found is an actual attractor.
        for test in range(50):
            n = random.randint(1, 20)
            G = graphs.Network.generate_random(n_vertices=n, indegree_bounds=[1, n - 1])
            attractor = stochastic.walk_to_attractor(G, stochastic.random_state(G))
            for t in range(len(attractor)):
                self.assertTrue(utility.is_same_state(attractor[(t + 1) % len(attractor)], G.next_state(attractor[t])))

    def test_estimate_attractors(self):
        for test in range(50):
            n = random.randint(1, 6)
            G = graphs.Network.generate_random(n_vertices=n, indegree_bounds=[1, 6])
        model_attractors = attractors.find_attractors_dubrova(G, attractors.dubrova_dir_path)
        estimated_attractors = stochastic.estimate_attractors(G, n_walks=1000, max_walk_len=1000, with_basins=False)
        print model_attractors
        print estimated_attractors
        self.assertTrue(utility.attractor_sets_equality(model_attractors, estimated_attractors))

