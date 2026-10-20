import os
import random
import itertools
import tempfile
from unittest import TestCase

from attractor_learning.graphs import Network
from attractor_learning.logic import SymmetricThresholdFunction, BooleanSymbolicFunc, SparseBooleanFunc


class TestNetworkSerialization(TestCase):

    def test_networkx_save_load_round_trip(self):
        # a high in-degree threshold node (would blow up cnet/truth-table export) plus a boolean node
        n_inputs = 25
        input_names = ["v{}".format(i) for i in range(n_inputs)]
        vertex_names = input_names + ["thr", "bs"]
        edges = [(name, "thr") for name in input_names] + [("v0", "bs"), ("v1", "bs")]
        signs = [random.choice([True, False]) for _ in range(n_inputs)]
        thr_func = SymmetricThresholdFunction(signs, threshold=13)
        bs_func = BooleanSymbolicFunc(input_names=["v0", "v1"], boolean_outputs=[False, True, True, False])
        functions = [None] * n_inputs + [thr_func, bs_func]
        G = Network(vertex_names=vertex_names, edges=edges, vertex_functions=functions)

        path = os.path.join(tempfile.gettempdir(), "test_rt_network.json")
        try:
            G.save(path)
            G_tag = Network.load(path)
        finally:
            if os.path.exists(path):
                os.remove(path)

        # compare behaviorally (== on networks would itself build truth tables) and by fields
        self.assertEqual({(u.name, v.name) for u, v in G.edges},
                         {(u.name, v.name) for u, v in G_tag.edges})
        loaded_thr = G_tag.get_vertex("thr").function
        self.assertIsInstance(loaded_thr, SymmetricThresholdFunction)
        self.assertEqual(loaded_thr.signs, thr_func.signs)
        self.assertEqual(loaded_thr.threshold, thr_func.threshold)
        for _ in range(2000):
            state = [random.randint(0, 1) for _ in range(len(vertex_names))]
            self.assertEqual(G.next_state(state), G_tag.next_state(state))

    def test_general_boolean_and_constant_input_round_trip(self):
        # BooleanSymbolicFunc built without input_names (as infer_known_topology_general does) gets generic
        # input_0/input_1 names that don't match the predecessors; plus a constant-True input node.
        or_func = BooleanSymbolicFunc(boolean_outputs=[False, True, True, True])  # 2-input OR, generic names
        G = Network(vertex_names=["a", "b", "c"], edges=[("a", "c"), ("b", "c")],
                    vertex_functions=[True, None, or_func])
        path = os.path.join(tempfile.gettempdir(), "test_general_network.json")
        try:
            G.save(path)
            G_tag = Network.load(path)
        finally:
            if os.path.exists(path):
                os.remove(path)
        self.assertEqual(bool(G_tag.get_vertex("a").function), True)  # constant input node preserved
        for state in itertools.product([0, 1], repeat=3):
            self.assertEqual(G.next_state(list(state)), G_tag.next_state(list(state)))

    def test_symmetric_threshold_compact_equality(self):
        # constant functions are independent of the signs, so different signs must compare/hash equal
        const_true_a = SymmetricThresholdFunction([True, False, True], 0)
        const_true_b = SymmetricThresholdFunction([False, False, False], 0)
        self.assertEqual(const_true_a, const_true_b)
        self.assertEqual(hash(const_true_a), hash(const_true_b))
        const_false_a = SymmetricThresholdFunction([True, True], 3)   # degree + 1
        const_false_b = SymmetricThresholdFunction([False, True], 3)
        self.assertEqual(const_false_a, const_false_b)
        self.assertEqual(hash(const_false_a), hash(const_false_b))
        self.assertNotEqual(const_true_a, const_false_a)
        # out-of-range threshold is rejected
        with self.assertRaises(AssertionError):
            SymmetricThresholdFunction([True, True], 4)


class TestSparseBooleanFunc(TestCase):
    """SparseBooleanFunc stores only the rows disagreeing with the majority output, so what needs covering
    is that it means the same function a full truth table would, and that it survives a file round trip."""

    def test_matches_dense_truth_table_on_every_row(self):
        n_inputs = 6
        input_names = ["v{}".format(i) for i in range(n_inputs)]
        rng = random.Random(17)
        # a strongly biased table, the case the sparse form exists for
        outputs = [rng.random() < 0.9 for _ in range(2 ** n_inputs)]
        dense = BooleanSymbolicFunc(input_names=input_names, boolean_outputs=outputs)
        sparse = SparseBooleanFunc(input_names=input_names, default_output=True,
                                   minority_rows=[row for row, out in enumerate(outputs) if not out])
        for row in itertools.product([False, True], repeat=n_inputs):
            self.assertEqual(bool(dense(*row)), bool(sparse(*row)), "disagreement on row {}".format(row))
        self.assertEqual(tuple(bool(o) for o in sparse.boolean_outputs), tuple(bool(o) for o in outputs))

    def test_network_round_trip_keeps_dynamics(self):
        # in-degree 12 is 4096 truth-table rows: fine sparsely, and the size the dense form is avoided at
        n_inputs = 12
        input_names = ["v{}".format(i) for i in range(n_inputs)]
        vertex_names = input_names + ["sparse"]
        edges = [(name, "sparse") for name in input_names]
        rng = random.Random(3)
        minority = {rng.randrange(2 ** n_inputs) for _ in range(40)}
        func = SparseBooleanFunc(input_names=input_names, minority_rows=minority, default_output=True)
        network = Network(vertex_names=vertex_names, edges=edges,
                          vertex_functions=[None] * n_inputs + [func])

        path = os.path.join(tempfile.gettempdir(), "test_sparse_network.json")
        try:
            network.save(path)
            reloaded = Network.load(path)
        finally:
            if os.path.exists(path):
                os.remove(path)
        self.assertIsInstance(reloaded.vertices[-1].function, SparseBooleanFunc)
        self.assertEqual(reloaded.vertices[-1].function, func)
        for _ in range(50):
            state = [rng.randint(0, 1) for _ in range(len(network))]
            self.assertEqual(list(reloaded.next_state(state)), list(network.next_state(state)))

    def test_canonical_form_makes_equality_and_hash_agree(self):
        a = SparseBooleanFunc(input_names=["x", "y", "z"], minority_rows=[1, 6], default_output=True)
        b = SparseBooleanFunc(input_names=["x", "y", "z"], minority_rows=[6, 1], default_output=True)
        c = SparseBooleanFunc(input_names=["x", "y", "z"], minority_rows=[1, 5], default_output=True)
        self.assertEqual(a, b)
        self.assertEqual(hash(a), hash(b))
        self.assertNotEqual(a, c)
        # a function with no disagreeing rows is the constant, and compares equal to it
        self.assertEqual(SparseBooleanFunc(input_names=["x"], minority_rows=[], default_output=True), True)

    def test_rejects_a_minority_that_is_not_the_minority(self):
        # the canonical form is what lets __eq__ skip the truth table, so the larger side is refused
        with self.assertRaises(ValueError):
            SparseBooleanFunc(input_names=["x", "y"], minority_rows=[0, 1, 2], default_output=True)
        with self.assertRaises(ValueError):
            SparseBooleanFunc(input_names=["x"], minority_rows=[5], default_output=True)
