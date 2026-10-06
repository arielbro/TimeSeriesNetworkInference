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

    def test_linear_classifier_all_features_ignores_the_scaffold(self):
        # the scaffold offers c only the wrong input; with additional edges allowed every node is a feature,
        # so the lasso finds c's real regulators a, b anyway, and the source nodes a, b come out as inputs
        random.seed(6)
        network = _make_network()
        scaffold = Network(vertex_names=["a", "b", "c", "d", "e"],
                           edges=[("d", "c"), ("d", "e"), ("a", "d"), ("b", "d")])  # max in-degree 2
        inferred = linear_classifier_inference(_transition_matrices(network), scaffold,
                                               allow_additional_edges=True)
        for name in LEARNED_NODES:
            vertex = inferred.get_vertex(name)
            self.assertEqual([u.name for u in vertex.predecessors()], ["a", "b"], name)
            self.assertEqual(tuple(bool(o) for o in vertex.function.boolean_outputs),
                             network.get_vertex(name).function.boolean_outputs, name)
        for name in ("a", "b"):
            self.assertEqual(len(inferred.get_vertex(name).predecessors()), 0, name)
            self.assertIsNone(inferred.get_vertex(name).function, name)

    def test_linear_classifier_all_features_respects_max_indegree(self):
        random.seed(7)
        network = _make_network()
        inferred = linear_classifier_inference(_transition_matrices(network), network,
                                               allow_additional_edges=True, max_indegree=1)
        for name in LEARNED_NODES:
            vertex = inferred.get_vertex(name)
            self.assertEqual(len(vertex.predecessors()), 1, name)
            self.assertEqual(len(vertex.function.boolean_outputs), 2, name)

    def test_linear_classifier_scaffold_mode_is_the_default(self):
        # without the flag, only the scaffold's (wrong) input is available to c
        random.seed(8)
        network = _make_network()
        scaffold = Network(vertex_names=["a", "b", "c", "d", "e"], edges=[("d", "c"), ("a", "d"), ("b", "d")])
        inferred = linear_classifier_inference(_transition_matrices(network), scaffold)
        self.assertEqual([u.name for u in inferred.get_vertex("c").predecessors()], ["d"])

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


class TestBreadthFirstSubsetSearch(TestCase):
    """REVEAL and Best-Fit sweep regulator-set sizes breadth-first: every gene at size 1, then every gene at
    size 2, and so on. Depth-first per gene would spend a bounded budget on whichever genes come first and
    leave the rest at size 1, making the result depend on gene order."""

    def setUp(self):
        names = ["v%d" % i for i in range(4)]
        self.network = Network(vertex_names=names,
                               edges=[("v0", "v1"), ("v1", "v2"), ("v2", "v3"), ("v0", "v3")])
        for vertex in self.network.vertices:
            degree = len(vertex.predecessors())
            if degree:
                vertex.function = BooleanSymbolicFunc(
                    input_names=[u.name for u in vertex.predecessors()],
                    boolean_outputs=[bool(row % 2) for row in range(2 ** degree)])
        rng = random.Random(4)
        self.matrices = []
        for _ in range(6):
            state = [rng.randint(0, 1) for _ in names]
            rows = [list(state)]
            for _ in range(3):
                state = [int(bool(v)) for v in self.network.next_state(state)]
                rows.append(list(state))
            self.matrices.append(np.array(rows, dtype=float))

    def _scan_sizes(self, scan_name, inference, **kwargs):
        """The subset sizes the scans are called with, in order."""
        import inference.benchmark_inference as module
        sizes = []
        original = getattr(module, scan_name)

        def recording(X, y, candidates, size, best, deadline, self_column=None):
            sizes.append(size)
            return original(X, y, candidates, size, best, deadline, self_column)

        setattr(module, scan_name, recording)
        try:
            inference(self.matrices, self.network, timeout_secs=60, **kwargs)
        finally:
            setattr(module, scan_name, original)
        return sizes

    def test_reveal_sweeps_every_gene_at_each_size_before_the_next(self):
        sizes = self._scan_sizes("_scan_reveal", reveal_inference, max_indegree=3)
        self.assertEqual(sizes, sorted(sizes),
                         "sizes are not swept in nondecreasing order: {}".format(sizes))
        # 3 genes searched (v0 is a source node and is emitted as an input); each is offered the empty set
        # (size 0, an input node) before any regulator set, then sizes 1..3
        self.assertEqual(sizes, [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3])

    def test_best_fit_sweeps_every_gene_at_size_one_first(self):
        sizes = self._scan_sizes("_scan_best_fit", best_fit_inference, max_indegree=3)
        self.assertEqual(sizes, sorted(sizes),
                         "sizes are not swept in nondecreasing order: {}".format(sizes))
        self.assertEqual(sizes[:6], [0, 0, 0, 1, 1, 1],
                         "genes were not all swept at size 0, then all at size 1: {}".format(sizes))

    def test_an_expired_budget_leaves_every_gene_at_size_one(self):
        """With a budget already spent, no gene gets past size 1 - and every gene still gets size 1, so the
        model is complete rather than partly built."""
        for scan_name, inference in (("_scan_reveal", reveal_inference),
                                     ("_scan_best_fit", best_fit_inference)):
            with self.subTest(method=scan_name):
                import inference.benchmark_inference as module
                sizes = []
                original = getattr(module, scan_name)

                def recording(X, y, candidates, size, best, deadline, self_column=None, _o=original):
                    sizes.append(size)
                    return _o(X, y, candidates, size, best, deadline, self_column)

                setattr(module, scan_name, recording)
                try:
                    model = inference(self.matrices, self.network, max_indegree=3, timeout_secs=-1)
                finally:
                    setattr(module, scan_name, original)
                self.assertEqual(set(sizes), {0, 1}, "a gene was searched past size 1: {}".format(sizes))
                self.assertEqual(len(sizes), 6, "not every searched gene got its size-0 and size-1 sweeps")
                for vertex in model.vertices:
                    self.assertLessEqual(len(vertex.predecessors()), 1)


class TestInputNodeAsASearchOption(TestCase):
    """The empty regulator set is swept before any size-1 set, so a gene that (mostly) holds its value is
    emitted as an input node rather than given a regulator it does not need."""

    def _network_with_a_nearly_static_gene(self, exceptions):
        """v0 drives v1; v2 holds its value except in `exceptions` transitions, where it is flipped."""
        names = ["v0", "v1", "v2"]
        network = Network(vertex_names=names, edges=[("v0", "v1")])
        network.get_vertex("v1").function = BooleanSymbolicFunc(input_names=["v0"],
                                                               boolean_outputs=[False, True])
        rng = random.Random(12)
        matrices = []
        for _ in range(10):
            state = [rng.randint(0, 1) for _ in names]
            rows = [list(state)]
            for _ in range(3):
                state = [int(bool(v)) for v in network.next_state(state)]
                rows.append(list(state))
            matrices.append(np.array(rows, dtype=float))
        # flip v2's successor value in the first `exceptions` transitions, so it no longer holds its value
        # exactly and the static pre-check no longer fires for it
        flipped = 0
        for matrix in matrices:
            for t in range(1, matrix.shape[0]):
                if flipped < exceptions:
                    matrix[t, 2] = 1 - matrix[t, 2]
                    flipped += 1
        return network, matrices

    def test_a_nearly_static_gene_is_still_an_input_node(self):
        network, matrices = self._network_with_a_nearly_static_gene(exceptions=1)
        for inference in (reveal_inference, best_fit_inference):
            with self.subTest(method=inference.__name__):
                model = inference(matrices, network, max_indegree=2, timeout_secs=60)
                v2 = model.get_vertex("v2")
                self.assertEqual(len(v2.predecessors()), 0,
                                 "a gene that holds its value in all but one transition was given "
                                 "regulators: {}".format([u.name for u in v2.predecessors()]))
                self.assertIsNone(v2.function)

    def test_emit_static_as_input_off_removes_the_option(self):
        """Without the flag the raw methods' behaviour is restored: no empty set is offered, so the gene
        gets a regulator (in practice the self-loop the flag exists to avoid)."""
        network, matrices = self._network_with_a_nearly_static_gene(exceptions=1)
        for inference in (reveal_inference, best_fit_inference):
            with self.subTest(method=inference.__name__):
                model = inference(matrices, network, max_indegree=2, timeout_secs=60,
                                  emit_static_as_input=False)
                self.assertGreaterEqual(len(model.get_vertex("v2").predecessors()), 1)

    def test_a_genuinely_regulated_gene_still_gets_its_regulator(self):
        """The empty set winning on score, not by default: v1 is driven by v0 and must keep that edge."""
        network, matrices = self._network_with_a_nearly_static_gene(exceptions=0)
        for inference in (reveal_inference, best_fit_inference):
            with self.subTest(method=inference.__name__):
                model = inference(matrices, network, max_indegree=2, timeout_secs=60)
                self.assertEqual([u.name for u in model.get_vertex("v1").predecessors()], ["v0"])
