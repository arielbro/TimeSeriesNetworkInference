"""Tests for infer_unknown_topology_symmetric, covering the candidate-set restriction.

The ILP searches a candidate input set per node. With allow_additional_edges the candidates are every
vertex, which is what the published (topology-free) behaviour needs; without it, inference may only keep or
drop scaffold edges, so the candidates are a node's scaffold predecessors and the model should be built at
that size rather than built over all n**2 pairs and constrained back down. These tests pin both the
behaviour (which must not change either way) and the size (which is the point of the restriction).
"""
import itertools
import random
import unittest

import numpy as np

from attractor_learning import graphs
from attractor_learning.logic import SymmetricThresholdFunction
from inference.binary_inference_ideas import infer_unknown_topology_symmetric, infer_symmetric_heuristic

TIMEOUT = 60


def threshold_network(vertex_names, edges, signs_and_thresholds):
    """A network whose every non-input node carries a symmetric threshold function.
    signs_and_thresholds maps a vertex name to (signs, threshold), signs ordered by predecessor index."""
    network = graphs.Network(vertex_names=vertex_names, edges=edges)
    for vertex in network.vertices:
        if vertex.name in signs_and_thresholds:
            signs, threshold = signs_and_thresholds[vertex.name]
            vertex.function = SymmetricThresholdFunction(signs, threshold)
    return network


def time_series(network, n_matrices, n_timepoints, seed):
    """Noise-free trajectories from random starting states, the shape inference expects."""
    rng = random.Random(seed)
    starts = [[rng.randint(0, 1) for _ in range(len(network))] for _ in range(n_matrices)]
    return _trajectories(network, starts, n_timepoints)


def exhaustive_time_series(network):
    """One two-row matrix per state of the network, i.e. every one-step transition there is. Inference fits
    the data, so a test that asks for the true model back has to hand it data that determines one - random
    starts do not, since a network's input nodes hold whatever value they started with."""
    starts = [list(state) for state in itertools.product([0, 1], repeat=len(network))]
    return _trajectories(network, starts, n_timepoints=2)


def _trajectories(network, starts, n_timepoints):
    matrices = []
    for state in starts:
        rows = [list(state)]
        for _ in range(n_timepoints - 1):
            state = [int(bool(v)) for v in network.next_state(state)]
            rows.append(list(state))
        matrices.append(np.array(rows, dtype=float))
    return matrices


def small_model():
    """A 5-node network: two inputs feeding an AND, an OR and a negation."""
    names = ["a", "b", "c", "d", "e"]
    edges = [("a", "c"), ("b", "c"), ("a", "d"), ("b", "d"), ("c", "e"), ("d", "e")]
    functions = {"c": ([True, True], 2),      # a AND b
                 "d": ([True, True], 1),      # a OR b
                 "e": ([True, False], 1)}     # c OR NOT d
    return threshold_network(names, edges, functions)


def edge_names(network):
    return {(u.name, v.name) for u, v in network.edges}


def infer(scaffold, matrices, **kwargs):
    kwargs.setdefault("timeout_secs", TIMEOUT)
    kwargs.setdefault("gurobi_threads", 1)
    return infer_unknown_topology_symmetric(matrices, scaffold, **kwargs)


class TestCandidateRestriction(unittest.TestCase):
    def setUp(self):
        self.true_network = small_model()
        self.matrices = time_series(self.true_network, n_matrices=4, n_timepoints=4, seed=7)

    def test_restricted_never_adds_an_edge_outside_the_scaffold(self):
        """allow_additional_edges=False: every inferred edge is a scaffold edge."""
        scaffold = graphs.Network(vertex_names=[v.name for v in self.true_network.vertices],
                                  edges=[(u.name, v.name) for u, v in self.true_network.edges])
        inferred = infer(scaffold, self.matrices, allow_additional_edges=False)
        self.assertTrue(edge_names(inferred) <= edge_names(scaffold),
                        "inferred edges outside the scaffold: {}".format(
                            edge_names(inferred) - edge_names(scaffold)))

    def test_restricted_can_drop_a_scaffold_edge(self):
        """The restriction bounds the candidates from above only - a superfluous scaffold edge can still be
        pruned, so this is not testing that the scaffold is copied verbatim."""
        scaffold_edges = [(u.name, v.name) for u, v in self.true_network.edges]
        scaffold_edges.append(("e", "c"))  # not in the true network
        scaffold = graphs.Network(vertex_names=[v.name for v in self.true_network.vertices],
                                  edges=scaffold_edges)
        inferred = infer(scaffold, self.matrices, allow_additional_edges=False)
        self.assertTrue(edge_names(inferred) <= set(scaffold_edges))
        self.assertIn(("a", "c"), edge_names(inferred))

    def test_unrestricted_may_use_a_non_scaffold_edge(self):
        """allow_additional_edges=True keeps every vertex a candidate, so an edge the scaffold omits can be
        recovered. Uses an empty scaffold, where nothing is available except added edges."""
        scaffold = graphs.Network(vertex_names=[v.name for v in self.true_network.vertices], edges=[])
        inferred = infer(scaffold, self.matrices, allow_additional_edges=True,
                         added_edges_relative_weight=-0.001)
        self.assertTrue(edge_names(inferred),
                        "unrestricted inference on an empty scaffold recovered no edges at all")

    def test_recovers_the_true_model_from_its_own_topology(self):
        """Given the true topology as the scaffold and noise-free data, the restricted search should
        reproduce the network's dynamics."""
        scaffold = graphs.Network(vertex_names=[v.name for v in self.true_network.vertices],
                                  edges=[(u.name, v.name) for u, v in self.true_network.edges])
        inferred = infer(scaffold, exhaustive_time_series(self.true_network),
                         allow_additional_edges=False)
        for state in itertools.product([0, 1], repeat=len(self.true_network)):
            self.assertEqual([int(bool(v)) for v in inferred.next_state(list(state))],
                             [int(bool(v)) for v in self.true_network.next_state(list(state))],
                             "dynamics differ at state {}".format(state))

    def test_node_without_scaffold_predecessors_holds_its_value(self):
        """A node the scaffold leaves isolated has no candidate inputs under the restriction, so it has no
        function and carries its value, the way an input node does."""
        scaffold = graphs.Network(vertex_names=["a", "b", "c", "d", "e"],
                                  edges=[("a", "c"), ("b", "c"), ("a", "d"), ("b", "d")])
        inferred = infer(scaffold, self.matrices, allow_additional_edges=False)
        isolated = inferred.get_vertex("e")
        self.assertEqual(len(isolated.predecessors()), 0)
        self.assertIsNone(isolated.function)
        state = [0, 1, 0, 1, 1]
        self.assertEqual(int(bool(inferred.next_state(state)[isolated.index])), state[isolated.index])

    def test_max_indegree_is_respected(self):
        scaffold = graphs.Network(vertex_names=[v.name for v in self.true_network.vertices],
                                  edges=[(u.name, v.name) for u, v in
                                         itertools.product(self.true_network.vertices, repeat=2)])
        inferred = infer(scaffold, self.matrices, allow_additional_edges=False, max_indegree=1)
        for vertex in inferred.vertices:
            self.assertLessEqual(len(vertex.predecessors()), 1, "node {} exceeded max_indegree".format(
                vertex.name))


class TestModelSize(unittest.TestCase):
    """The restriction exists to make the model smaller; correctness tests alone would not catch losing it."""

    def setUp(self):
        self.true_network = small_model()
        self.matrices = time_series(self.true_network, n_matrices=4, n_timepoints=4, seed=11)

    @staticmethod
    def _n_vars(scaffold, matrices, **kwargs):
        """Variable count of the model the inference builds, captured by intercepting the solve."""
        sizes = {}
        import gurobipy

        original = gurobipy.Model.optimize

        def record_then_solve(model, *args, **kwargs_):
            model.update()
            sizes['vars'] = model.NumVars
            sizes['constrs'] = model.NumConstrs
            return original(model, *args, **kwargs_)

        gurobipy.Model.optimize = record_then_solve
        try:
            infer(scaffold, matrices, **kwargs)
        finally:
            gurobipy.Model.optimize = original
        return sizes

    def test_restricted_model_scales_with_the_scaffold_not_with_n_squared(self):
        """A sparse scaffold on the same network must build a much smaller model than the unrestricted
        search, whose candidate set is every vertex."""
        names = [v.name for v in self.true_network.vertices]
        sparse = graphs.Network(vertex_names=names,
                                edges=[(u.name, v.name) for u, v in self.true_network.edges])
        complete = graphs.Network(vertex_names=names,
                                  edges=[(u.name, v.name) for u, v in
                                         itertools.product(self.true_network.vertices, repeat=2)])
        restricted = self._n_vars(sparse, self.matrices, allow_additional_edges=False)
        unrestricted = self._n_vars(sparse, self.matrices, allow_additional_edges=True)
        full_scaffold = self._n_vars(complete, self.matrices, allow_additional_edges=False)

        # 6 scaffold edges against 5*5=25 candidate pairs
        self.assertLess(restricted['vars'] * 2, unrestricted['vars'],
                        "restricted model ({} vars) is not appreciably smaller than the unrestricted one "
                        "({} vars)".format(restricted['vars'], unrestricted['vars']))
        self.assertLess(restricted['constrs'] * 2, unrestricted['constrs'])
        # restricting to a complete scaffold is the unrestricted problem, so the sizes should agree
        self.assertEqual(full_scaffold['vars'], unrestricted['vars'])

    @staticmethod
    def _largest_threshold_big_m(scaffold, matrices, **kwargs):
        """The largest coefficient appearing in the threshold-comparison constraints, i.e. the big-M the
        formulation chose. Not visible in any size counter, so it needs reading off the matrix."""
        import gurobipy

        original = gurobipy.Model.optimize
        found = {}

        def record_then_solve(model, *args, **kwargs_):
            model.update()
            coefficients = set()
            for constraint in model.getConstrs():
                if "threshold_function_path_constraint_>=" in constraint.ConstrName:
                    row = model.getRow(constraint)
                    coefficients.update(abs(row.getCoeff(i)) for i in range(row.size()))
            found['max'] = max(coefficients) if coefficients else None
            return original(model, *args, **kwargs_)

        gurobipy.Model.optimize = record_then_solve
        try:
            infer(scaffold, matrices, **kwargs)
        finally:
            gurobipy.Model.optimize = original
        return found['max']

    def test_max_indegree_tightens_the_big_m(self):
        """A node may use at most max_indegree of its candidates, so the threshold comparison is bounded by
        that rather than by the candidate count."""
        names = [v.name for v in self.true_network.vertices]
        complete = graphs.Network(vertex_names=names,
                                  edges=[(u.name, v.name) for u, v in
                                         itertools.product(self.true_network.vertices, repeat=2)])
        uncapped = self._largest_threshold_big_m(complete, self.matrices, allow_additional_edges=True)
        capped = self._largest_threshold_big_m(complete, self.matrices, allow_additional_edges=True,
                                               max_indegree=2)
        self.assertEqual(uncapped, len(names) + 1)
        self.assertEqual(capped, 3)

    def test_restriction_tightens_the_big_m(self):
        """Restricting the candidates to a sparse scaffold bounds a node's degree by its own candidate count,
        which tightens the same coefficient. max_indegree is set to n here so that the cap is not binding and
        the candidate count is what shows: left at -1 it defaults to the scaffold's max in-degree, which
        already tightens both paths to the same value."""
        names = [v.name for v in self.true_network.vertices]
        sparse = graphs.Network(vertex_names=names,
                                edges=[(u.name, v.name) for u, v in self.true_network.edges])
        restricted = self._largest_threshold_big_m(sparse, self.matrices, allow_additional_edges=False,
                                                   max_indegree=len(names))
        unrestricted = self._largest_threshold_big_m(sparse, self.matrices, allow_additional_edges=True,
                                                     max_indegree=len(names))
        self.assertEqual(unrestricted, len(names) + 1)
        self.assertLess(restricted, unrestricted)

    def test_max_indegree_does_not_change_the_variable_count(self):
        """max_indegree tightens bounds rather than removing candidates, so it must not silently shrink the
        model - if it ever does, the big-M tightening below is being applied to the wrong quantity."""
        names = [v.name for v in self.true_network.vertices]
        sparse = graphs.Network(vertex_names=names,
                                edges=[(u.name, v.name) for u, v in self.true_network.edges])
        without = self._n_vars(sparse, self.matrices, allow_additional_edges=False)
        with_cap = self._n_vars(sparse, self.matrices, allow_additional_edges=False, max_indegree=1)
        self.assertEqual(without['vars'], with_cap['vars'])


def ring3():
    """Repressilator a0 -| a1 -| a2 -| a0: a bijective state map, cycles of period 2 and 6 and no transients."""
    return threshold_network(["a0", "a1", "a2"], [("a2", "a0"), ("a0", "a1"), ("a1", "a2")],
                             {name: ([False], 1) for name in ["a0", "a1", "a2"]})


def copy_gate():
    """Input a (holds its value) copied into b: from (1, 0) one transient step to the fixed point (1, 1)."""
    return threshold_network(["a", "b"], [("a", "b")], {"b": ([True], 1)})


def fixed_function_vars(network):
    """The network's own signs and thresholds in add_path_to_model's function-variable shape, as constants,
    so the attractor components can be checked against a known model."""
    return [None if v.function is None else [[1 if s else -1 for s in v.function.signs], v.function.threshold]
            for v in network.vertices]


def run_states(network, state, n_steps):
    states = [list(state)]
    for _ in range(n_steps):
        states.append([int(bool(v)) for v in network.next_state(states[-1])])
    return states


class TestAttractorComponents(unittest.TestCase):
    """The attractor constraints on a known model and a fixed trajectory: satisfiable exactly when the
    trajectory has the claimed attractor property."""

    @staticmethod
    def _feasible(network, add_constraint):
        import gurobipy
        from attractor_learning.graphs import FunctionTypeRestriction
        model = gurobipy.Model()
        model.Params.OutputFlag = 0
        try:
            restrictions = [FunctionTypeRestriction.SYMMETRIC_THRESHOLD] * len(network)
            add_constraint(model, restrictions)
            model.optimize()
            return model.Status == gurobipy.GRB.OPTIMAL
        finally:
            model.dispose()

    def _returns(self, network, trajectory, max_path_len):
        from inference import ilp_components
        return self._feasible(network, lambda model, restrictions: ilp_components.add_attractor_return_path(
            network, model, trajectory, fixed_function_vars(network), max_path_len,
            function_type_restrictions=restrictions))

    def _repeats(self, network, trajectory):
        from inference import ilp_components
        return self._feasible(network, lambda model, _: ilp_components.add_attractor_repeat(model, trajectory))

    def test_return_path_needs_the_cycle_remainder(self):
        """Three states of ring3's 6-cycle: getting from the third back to the first takes 4 more steps."""
        network = ring3()
        trajectory = run_states(network, [0, 0, 1], 2)
        self.assertTrue(self._returns(network, trajectory, max_path_len=4))
        self.assertTrue(self._returns(network, trajectory, max_path_len=6))
        self.assertFalse(self._returns(network, trajectory, max_path_len=3))

    def test_return_path_rejects_a_transient(self):
        network = copy_gate()
        trajectory = run_states(network, [1, 0], 1)          # (1, 0) -> (1, 1), never back
        self.assertFalse(self._returns(network, trajectory, max_path_len=5))
        self.assertTrue(self._returns(network, [[1, 1], [1, 1]], max_path_len=1))

    def test_repeat_needs_a_closed_cycle(self):
        network = ring3()
        self.assertTrue(self._repeats(network, run_states(network, [0, 0, 1], 6)))    # the full 6-cycle, closed
        self.assertFalse(self._repeats(network, run_states(network, [0, 0, 1], 5)))   # one state short
        gate = copy_gate()
        self.assertTrue(self._repeats(gate, run_states(gate, [1, 0], 2)))     # transient, then the fixed point twice
        self.assertFalse(self._repeats(gate, run_states(gate, [1, 0], 1)))    # ends before repeating


class TestAttractorMethods(unittest.TestCase):
    def setUp(self):
        from inference.binary_inference_ideas import (infer_unknown_topology_symmetric_in_attractors,
                                                      infer_unknown_topology_symmetric_contains_attractors)
        self.methods = [infer_unknown_topology_symmetric_in_attractors,
                        infer_unknown_topology_symmetric_contains_attractors]
        self.network = ring3()
        self.scaffold = graphs.Network(vertex_names=[v.name for v in self.network.vertices],
                                       edges=[(u.name, v.name) for u, v in self.network.edges])

    def test_anchored_and_flipless_modes_are_rejected(self):
        matrices = [np.array(run_states(self.network, [0, 0, 1], 6), dtype=float)]
        for method in self.methods:
            for no_anchoring, allow_input_flips in [(False, True), (True, False), (False, False)]:
                with self.assertRaises(ValueError, msg="{} accepted no_anchoring={}, flips={}".format(
                        method.__name__, no_anchoring, allow_input_flips)):
                    method(matrices, self.scaffold, no_anchoring=no_anchoring,
                           allow_input_flips=allow_input_flips, timeout_secs=TIMEOUT, gurobi_threads=1)

    def test_attractor_data_recovers_the_model(self):
        """Clean data covering both of ring3's attractors in full: the true model satisfies either
        assumption, so it is still the best fit."""
        matrices = [np.array(run_states(self.network, start, 6), dtype=float)
                    for start in ([0, 0, 1], [0, 0, 0])]
        for method in self.methods:
            inferred = method(matrices, self.scaffold, no_anchoring=True, allow_input_flips=True,
                              flip_penalty=2.0, allow_additional_edges=False, timeout_secs=TIMEOUT,
                              gurobi_threads=1, attractor_max_path_len=6)
            for state in itertools.product([0, 1], repeat=3):
                self.assertEqual([int(bool(v)) for v in inferred.next_state(list(state))],
                                 [int(bool(v)) for v in self.network.next_state(list(state))],
                                 "{}: dynamics differ at {}".format(method.__name__, state))


class TestSymmetricHeuristic(unittest.TestCase):
    """infer_symmetric_heuristic: the per-node fit symmetric_topology warm-starts from, as a method."""

    def setUp(self):
        self.true_network = small_model()
        self.names = [v.name for v in self.true_network.vertices]
        self.true_edges = [(u.name, v.name) for u, v in self.true_network.edges]

    def test_recovers_the_true_model_from_its_own_topology(self):
        scaffold = graphs.Network(vertex_names=self.names, edges=self.true_edges)
        inferred = infer_symmetric_heuristic(exhaustive_time_series(self.true_network), scaffold,
                                             included_edges_relative_weight=-0.01)
        for state in itertools.product([0, 1], repeat=len(self.true_network)):
            self.assertEqual([int(bool(v)) for v in inferred.next_state(list(state))],
                             [int(bool(v)) for v in self.true_network.next_state(list(state))],
                             "dynamics differ at state {}".format(state))

    def test_restricted_never_adds_an_edge_outside_the_scaffold(self):
        scaffold = graphs.Network(vertex_names=self.names, edges=self.true_edges + [("e", "c")])
        inferred = infer_symmetric_heuristic(time_series(self.true_network, 4, 4, seed=7), scaffold)
        self.assertTrue(edge_names(inferred) <= edge_names(scaffold))

    def test_unrestricted_may_use_a_non_scaffold_edge(self):
        scaffold = graphs.Network(vertex_names=self.names, edges=[])
        inferred = infer_symmetric_heuristic(exhaustive_time_series(self.true_network), scaffold,
                                             allow_additional_edges=True, added_edges_relative_weight=-0.001)
        self.assertTrue(edge_names(inferred))

    def test_max_indegree_is_respected(self):
        scaffold = graphs.Network(vertex_names=self.names,
                                  edges=list(itertools.product(self.names, repeat=2)))
        inferred = infer_symmetric_heuristic(exhaustive_time_series(self.true_network), scaffold,
                                             max_indegree=1)
        for vertex in inferred.vertices:
            self.assertLessEqual(len(vertex.predecessors()), 1)

    def test_runner_resolves_the_method_name(self):
        from inference.run_inference_on_data import resolve_inference_method, inference_method_name
        self.assertIs(resolve_inference_method("symmetric_heuristic"), infer_symmetric_heuristic)
        self.assertEqual(inference_method_name(infer_symmetric_heuristic), "symmetric_heuristic")


if __name__ == "__main__":
    unittest.main()
