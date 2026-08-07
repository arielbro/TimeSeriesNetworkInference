import itertools
import random
from unittest import TestCase

import numpy as np

from attractor_learning.graphs import Network
from attractor_learning.logic import BooleanSymbolicFunc
from inference.benchmark_inference import (random_model_inference,
                                           exact_match_else_random_inference,
                                           linear_classifier_inference,
                                           reveal_inference, best_fit_inference,
                                           _select_reveal)

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


def _reveal_bestfit_network():
    # inputs a, b, c (indices 0,1,2); target g = a AND NOT b (asymmetric -> guards bit order); c is a
    # distractor input g does not depend on. Scaffold max in-degree = 2 (g's), so max_indegree=-1 gives k=2.
    g = BooleanSymbolicFunc(input_names=["a", "b"], boolean_outputs=[False, False, True, False])  # a & ~b
    return Network(vertex_names=["a", "b", "c", "g"], edges=[("a", "g"), ("b", "g")],
                   vertex_functions=[None, None, None, g])


def _all_transition_matrices(network):
    return [np.array([list(s), list(network.next_state(s))], dtype=float)
            for s in itertools.product([0, 1], repeat=len(network))]


class TestRevealBestFit(TestCase):
    def test_recover_regulators_and_positional_truth_table(self):
        for infer in (reveal_inference, best_fit_inference):
            random.seed(0)
            net = _reveal_bestfit_network()
            model = infer(_all_transition_matrices(net), net, max_indegree=-1)  # k = scaffold max in-degree = 2
            g = model.get_vertex("g")
            self.assertEqual([u.name for u in g.predecessors()], ["a", "b"], infer.__name__)
            # a AND NOT b over (a, b) with a as MSB -> (F,F,T,F); a reversed bit order would give (F,T,F,F)
            self.assertEqual(tuple(bool(o) for o in g.function.boolean_outputs),
                             (False, False, True, False), infer.__name__)

    def test_static_inputs_emitted_as_input_nodes_and_behaviour_matches(self):
        for infer in (reveal_inference, best_fit_inference):
            random.seed(1)
            net = _reveal_bestfit_network()
            model = infer(_all_transition_matrices(net), net, max_indegree=-1)
            for name in ("a", "b", "c"):  # source nodes hold their value -> emitted as input nodes
                self.assertEqual(len(model.get_vertex(name).predecessors()), 0, (infer.__name__, name))
                self.assertIsNone(model.get_vertex(name).function, (infer.__name__, name))
            # topology recovered exactly and dynamics reproduced
            self.assertEqual({(u.name, v.name) for u, v in model.edges}, {("a", "g"), ("b", "g")})
            for s in itertools.product([0, 1], repeat=len(net)):
                self.assertEqual(model.next_state(s), net.next_state(s), infer.__name__)

    def test_emit_static_off_gives_self_loops(self):
        random.seed(2)
        net = _reveal_bestfit_network()
        model = best_fit_inference(_all_transition_matrices(net), net, max_indegree=-1,
                                   emit_static_as_input=False)
        # with the toggle off, each source node holds its value via an identity self-loop
        self.assertIn(("a", "a"), {(u.name, v.name) for u, v in model.edges})

    def test_max_indegree_caps_regulators(self):
        random.seed(3)
        net = _reveal_bestfit_network()
        # g needs 2 regulators (a, b); a cap of 1 must force a single regulator
        model = best_fit_inference(_all_transition_matrices(net), net, max_indegree=1)
        self.assertEqual(len(model.get_vertex("g").predecessors()), 1)

    def test_reveal_is_parsimonious_on_noisy_data(self):
        # y = a AND NOT b with 10% flips; c,d,e,f are pure-noise distractors; budget k=5. REVEAL's MDL
        # relaxation must NOT grab the whole budget - it should keep the true small set {a, b}.
        rng = np.random.default_rng(0)
        n, N = 6, 800
        X = rng.integers(0, 2, size=(N, n)).astype(np.int8)
        y = ((X[:, 0] == 1) & (X[:, 1] == 0)).astype(np.int8)   # a AND NOT b
        flip = rng.random(N) < 0.10
        y = np.where(flip, 1 - y, y).astype(np.int8)
        reveal_set = _select_reveal(X, y, n, k=5)
        # A naive max-MI REVEAL would grab the whole budget (MI is monotone in added inputs); the MDL penalty
        # keeps it to exactly the true regulators. (Best-Fit's size here is noise/N-dependent - no assertion.)
        self.assertEqual(set(reveal_set), {0, 1})

    def test_scaffold_restriction_keeps_the_search_inside_the_given_topology(self):
        # g truly depends on a and b, but the scaffold only offers b and c. With added edges disallowed the
        # search may pick any subset of {b, c} and must never reach for a, however much better a would score.
        for infer in (reveal_inference, best_fit_inference):
            random.seed(6)
            net = _reveal_bestfit_network()  # true network: g = a AND NOT b
            scaffold = Network(vertex_names=["a", "b", "c", "g"], edges=[("b", "g"), ("c", "g")])
            model = infer(_all_transition_matrices(net), scaffold, max_indegree=-1,
                          emit_static_as_input=False, allow_additional_edges=False)
            self.assertLessEqual({(u.name, v.name) for u, v in model.edges},
                                 {("b", "g"), ("c", "g")}, infer.__name__)
            # a, b and c have no scaffold predecessors, so there is nothing to search: input nodes even with
            # emit_static_as_input off (which would otherwise have given them identity self-loops)
            for name in ("a", "b", "c"):
                self.assertEqual(len(model.get_vertex(name).predecessors()), 0, (infer.__name__, name))

    def test_timeout_degrades_to_single_regulators_but_still_returns_a_complete_model(self):
        for infer in (reveal_inference, best_fit_inference):
            net = _reveal_bestfit_network()   # g needs both a and b; k = 2
            data = _all_transition_matrices(net)

            random.seed(7)  # an already-spent budget: sizes above 1 are never swept
            spent = infer(data, net, max_indegree=-1, timeout_secs=0)
            g = spent.get_vertex("g")
            self.assertEqual(len(g.predecessors()), 1, infer.__name__)
            self.assertIsNotNone(g.function, infer.__name__)  # every gene still gets a function
            self.assertEqual(len(g.function.boolean_outputs), 2, infer.__name__)

            random.seed(7)  # a budget that cannot bind must leave the unbounded answer untouched
            timed = infer(data, net, max_indegree=-1, timeout_secs=60)
            random.seed(7)
            untimed = infer(data, net, max_indegree=-1)
            self.assertEqual({(u.name, v.name) for u, v in timed.edges},
                             {(u.name, v.name) for u, v in untimed.edges}, infer.__name__)

    def test_no_transition_data_returns_all_input_nodes(self):
        random.seed(4)
        net = _reveal_bestfit_network()
        model = reveal_inference([np.array([[0, 0, 0, 0]], dtype=float)], net, max_indegree=-1)  # single row
        self.assertEqual(model.edges, [])
        self.assertTrue(all(v.function is None for v in model.vertices))
