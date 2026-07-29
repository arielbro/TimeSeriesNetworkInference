import itertools
import random
from unittest import TestCase

import numpy as np

from attractor_learning.graphs import Network
from attractor_learning.logic import BooleanSymbolicFunc
from inference.benchmark_inference import (random_model_inference,
                                           exact_match_else_random_inference,
                                           linear_classifier_inference)

# The learned nodes and their true truth tables (inputs [a, b] with a as the most significant bit). "e" is
# ASYMMETRIC (a AND NOT b): its table changes under an input transposition, so these tests actually pin the
# first-predecessor-as-MSB bit convention - symmetric AND/OR alone would pass even with the bits reversed.
LEARNED_NODES = ["c", "d", "e"]


def _make_network():
    and_func = BooleanSymbolicFunc(input_names=["a", "b"], boolean_outputs=[False, False, False, True])   # a & b
    or_func = BooleanSymbolicFunc(input_names=["a", "b"], boolean_outputs=[False, True, True, True])       # a | b
    andnot_func = BooleanSymbolicFunc(input_names=["a", "b"], boolean_outputs=[False, False, True, False])  # a & ~b
    return Network(vertex_names=["a", "b", "c", "d", "e"],
                   edges=[("a", "c"), ("b", "c"), ("a", "d"), ("b", "d"), ("a", "e"), ("b", "e")],
                   vertex_functions=[None, None, and_func, or_func, andnot_func])


def _transition_matrices(network):
    """One clean 2-row (state, next_state) matrix per state, so every truth-table row of every node is
    observed exactly once as a genuine model transition."""
    return [np.array([list(state), list(network.next_state(state))], dtype=float)
            for state in itertools.product([0, 1], repeat=len(network))]


class TestBenchmarkInference(TestCase):
    def test_random_model_keeps_topology_and_input_nodes(self):
        random.seed(0)
        network = _make_network()
        inferred = random_model_inference([], network)
        self.assertEqual({(u.name, v.name) for u, v in network.edges},
                         {(u.name, v.name) for u, v in inferred.edges})
        # input nodes keep a None function; nodes with predecessors get a fair-coin truth table
        self.assertIsNone(inferred.get_vertex("a").function)
        self.assertIsNone(inferred.get_vertex("b").function)
        self.assertEqual(len(inferred.get_vertex("c").function.boolean_outputs), 4)

    def test_exact_match_recovers_observed_rows(self):
        random.seed(1)
        network = _make_network()
        inferred = exact_match_else_random_inference(_transition_matrices(network), network)
        # every row of each learned node was observed once, so consensus == the true output for every row
        for name in LEARNED_NODES:
            self.assertEqual(tuple(bool(o) for o in inferred.get_vertex(name).function.boolean_outputs),
                             network.get_vertex(name).function.boolean_outputs)

    def test_exact_match_leaves_unobserved_rows_alone(self):
        network = _make_network()
        # a single transition from (a=1, b=1); only truth-table row 3 is ever observed
        matrix = np.array([[1, 1, 0, 0, 0], list(network.next_state((1, 1, 0, 0, 0)))], dtype=float)

        random.seed(2)
        base = random_model_inference([matrix], network)
        base_outputs = {name: list(base.get_vertex(name).function.boolean_outputs) for name in LEARNED_NODES}

        random.seed(2)  # same seed -> same random base inside the inference, so unobserved rows must match it
        inferred = exact_match_else_random_inference([matrix], network)
        for name in LEARNED_NODES:
            idx = network.get_vertex(name).index
            outputs = list(inferred.get_vertex(name).function.boolean_outputs)
            self.assertEqual(outputs[3], bool(matrix[1, idx]))       # observed row overridden to the true output
            self.assertEqual(outputs[:3], base_outputs[name][:3])    # unobserved rows keep the random base

    def test_linear_classifier_recovers_separable_functions(self):
        random.seed(3)
        network = _make_network()
        inferred = linear_classifier_inference(_transition_matrices(network), network)
        # AND, OR and (a AND NOT b) are all linearly separable, so LR recovers them exactly
        for name in LEARNED_NODES:
            self.assertEqual(tuple(bool(o) for o in inferred.get_vertex(name).function.boolean_outputs),
                             network.get_vertex(name).function.boolean_outputs)

    def test_linear_classifier_constant_node(self):
        random.seed(4)
        network = _make_network()
        # inputs a,b never both 1, so c = a AND b is constant 0 in this data (single observed class)
        matrix = np.zeros((6, 5), dtype=float)
        matrix[:, 0] = [0, 1, 0, 1, 0, 1]
        matrix[:, 1] = [1, 0, 1, 0, 1, 0]
        inferred = linear_classifier_inference([matrix], network)
        self.assertEqual({bool(o) for o in inferred.get_vertex("c").function.boolean_outputs}, {False})

    def test_predictions_are_consistent_with_next_state(self):
        # end-to-end guard on row ordering: the recovered functions must reproduce the true transitions
        random.seed(5)
        network = _make_network()
        inferred = exact_match_else_random_inference(_transition_matrices(network), network)
        for state in itertools.product([0, 1], repeat=len(network)):
            expected = network.next_state(state)
            got = inferred.next_state(state)
            for name in LEARNED_NODES:
                idx = network.get_vertex(name).index
                self.assertEqual(got[idx], expected[idx])
