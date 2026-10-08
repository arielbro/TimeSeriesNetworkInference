import attractor_learning.graphs
import attractor_learning.ilp
import gurobipy
from inference import ilp_components
from attractor_learning.graphs import FunctionTypeRestriction, Network
from inference.ilp_components import ModelAdditionType, AttractorAssumption, get_value_of_gurobi_entity
from attractor_learning.logic import BooleanSymbolicFunc, SymmetricThresholdFunction
import numpy as np


# TODO: rethink the scattered way I do config, that made it worthwhile to have these
# TODO: "overload" functions for easy logging and usage.
def infer_known_topology_symmetric(*args, **kwargs):
    return infer_known_topology(*args, **kwargs,
                                function_type_restriction=FunctionTypeRestriction.SYMMETRIC_THRESHOLD)


def infer_known_topology_general(*args, **kwargs):
    return infer_known_topology(*args, **kwargs,
                                function_type_restriction=FunctionTypeRestriction.NONE)


def infer_known_topology(data_matrices, scaffold_network, function_type_restriction=None,
                         timeout_secs=None, log_file=None, allow_input_flips=False, flip_penalty=1.0,
                         no_anchoring=False, gurobi_threads=0, **kwargs):
    """
    Find a model with best fit to data_matrices, assuming that each node's inputs are defined by the scaffold network
    topology.
    The objective is the cell-wise agreement (fraction of (predicted state, node) cells the model explains correctly),
    optionally reduced by a penalty for treating observed input cells as noisy and flipping them.
    :param function_type_restriction:
    :param data_matrices:
    :param scaffold_network:
    :param timeout_secs: timelimit to pass to the solver.
    :param allow_input_flips: if true, the solver may flip bits of the (noisy) input data, paying flip_penalty per flip.
    :param flip_penalty: penalty per flipped input cell, on the same [0, 1] per-cell scale as the agreement score.
        flip_penalty < 1 denoises toward a self-consistent trajectory, > 1 only corrects errors that fix more than
        one cell, and = 1 (the default) is the break-even point at which the chosen flips may be under-determined.
    :param no_anchoring: if true, score a free-running trajectory (each step is fed the model's own previous
        prediction) instead of independent one-step transitions anchored to the data.
    :return:
    """
    # create function variables
    functions_variables = []
    model = gurobipy.Model()
    try:
        for i in range(len(scaffold_network)):
            degree = len(scaffold_network.vertices[i].predecessors())
            if degree == 0:
                functions_variables.append(None)
            elif function_type_restriction is None or function_type_restriction == FunctionTypeRestriction.NONE:
                functions_variables.append([model.addVar(lb=0, ub=1, vtype=gurobipy.GRB.INTEGER,
                                                     name="vertex_{}_row_{}_function_var".format(i, j))
                                            for j in range(2 ** degree)])
            elif function_type_restriction == FunctionTypeRestriction.SYMMETRIC_THRESHOLD:
                # use binary sign "pre-variables" so that the inputs _must_ be used (since we're assuming reference network
                # is gold truth)
                pre_signs = [model.addVar(lb=0, ub=1, vtype=gurobipy.GRB.INTEGER,
                                                     name="vertex_{}_input_{}_pre_sign_var".format(i, j))
                                for j in range(degree)]
                signs = [sign * 2 - 1 for sign in pre_signs]
                threshold = model.addVar(lb=0, ub=degree + 1, vtype=gurobipy.GRB.INTEGER,
                                         name="vertex_{}_threshold_var".format(i))
                functions_variables.append([signs, threshold])
            else:
                raise NotImplementedError()

        matrix_agreement_indicators, flip_cost_terms = ilp_components.add_matrices_as_model_paths(
           scaffold_network, model, data_matrices,
           function_vars=functions_variables,
           model_to_data_sample_rate_ratio=1,
           function_type_restrictions=[function_type_restriction] * len(scaffold_network),
           model_addition_type=ModelAdditionType.INDICATORS,
           per_cell_indicators=True, allow_input_flips=allow_input_flips, no_anchoring=no_anchoring)
        # cell-wise agreement normalized to [0, 1], so it shares a scale with the (per-cell) flip penalty
        n_cells = float(len(matrix_agreement_indicators))
        objective = gurobipy.quicksum(matrix_agreement_indicators) / n_cells
        if allow_input_flips and flip_cost_terms:
            objective = objective - flip_penalty * gurobipy.quicksum(flip_cost_terms) / n_cells

        if log_file is not None:
            model.Params.LogFile = log_file
        if timeout_secs is not None:
            model.Params.TimeLimit = timeout_secs
        if gurobi_threads:
            model.Params.Threads = gurobi_threads
        model.Params.MIPFocus = 1
        model.setObjective(objective, sense=gurobipy.GRB.MAXIMIZE)
        model.optimize()

        # set Boolean model found
        inferred_model = scaffold_network.copy()
        for i in range(len(inferred_model)):
            if len(inferred_model.vertices[i].predecessors()) == 0:
                func = None
            elif function_type_restriction is None or function_type_restriction == FunctionTypeRestriction.NONE:
                func = BooleanSymbolicFunc(boolean_outputs=[
                    get_value_of_gurobi_entity(functions_variables[i][j]) for j
                                                            in range(len(functions_variables[i]))])
            elif function_type_restriction == FunctionTypeRestriction.SYMMETRIC_THRESHOLD:
                signs, threshold = functions_variables[i]
                float_signs = [get_value_of_gurobi_entity(sign) for sign in signs]
                signs = [int(round(sign, 3)) for sign in float_signs]
                # gurobi can give "almost" integer values even for variables defined as integer type
                assert signs == [round(s, 3) for s in float_signs]
                threshold = get_value_of_gurobi_entity(threshold)
                func = SymmetricThresholdFunction(signs, threshold)
            else:
                raise NotImplementedError()
            inferred_model.vertices[i].function = func
    finally:
        model.dispose()
    return inferred_model


def _transition_rows(data_matrices):
    """Every one-step transition in the data as (X, Y) 0/1 arrays: X the state at t, Y the state at t+1."""
    xs, ys = [], []
    for matrix in data_matrices:
        matrix = np.asarray(matrix)
        if matrix.shape[0] >= 2:
            xs.append((matrix[:-1] != 0).astype(np.int8))
            ys.append((matrix[1:] != 0).astype(np.int8))
    return np.vstack(xs), np.vstack(ys)


def _fit_threshold_node(x, y, own, candidates, edge_costs, max_degree=None):
    """A symmetric threshold function for one node, fitted to anchored transitions without a solver.

    Signs start from each candidate input's correlation with the node's next value; the threshold is the
    best one for the signs in use, by enumeration; then single changes to one input's sign (+1, -1, or 0 to
    drop it) are taken while they improve agreed cells minus the edge costs (in cells). A node left with no
    input holds its value, as in the model. At most max_degree inputs are used, where given. Returns
    ({input index: sign}, threshold)."""
    columns = {i: x[:, i] for i in candidates}

    def value(signs):
        used = [i for i in candidates if signs[i] != 0]
        if max_degree is not None and len(used) > max_degree:
            return -np.inf, 1
        if not used:
            return int((own == y).sum()), 0
        agreeing = sum((columns[i] == 1) if signs[i] > 0 else (columns[i] == 0) for i in used).astype(int)
        best_score, best_t = -1, 1
        for t in range(1, len(used) + 1):
            score = int(((agreeing >= t) == y).sum())
            if score > best_score:
                best_score, best_t = score, t
        return best_score - sum(edge_costs[i] for i in used), best_t

    signs = {}
    for i in candidates:
        column = columns[i]
        if column.std() == 0 or y.std() == 0:
            signs[i] = 0
        else:
            signs[i] = 1 if np.corrcoef(column, y)[0, 1] >= 0 else -1
    if max_degree is not None:     # keep the strongest correlations within the cap
        strength = {i: abs(np.corrcoef(columns[i], y)[0, 1]) if signs[i] else 0.0 for i in candidates}
        for i in sorted(candidates, key=lambda i: -strength[i])[max_degree:]:
            signs[i] = 0
    current = value(signs)
    improved = True
    while improved:
        improved = False
        for i in candidates:
            for sign in (1, -1, 0):
                if sign == signs[i]:
                    continue
                trial = dict(signs)
                trial[i] = sign
                result = value(trial)
                if result[0] > current[0]:
                    signs, current, improved = trial, result, True
    return {i: s for i, s in signs.items() if s != 0}, current[1]


def _heuristic_warm_start_values(data_matrices, candidate_inputs, edge_costs, max_degree=None):
    """Per vertex a ({input vertex index: +1/-1 sign}, threshold) pair, fitted node by node on the anchored
    transitions (see _fit_threshold_node). Takes a fraction of a second, and on random NK data starts the
    free-run MIP far above what Gurobi's own heuristics find in that time. edge_costs[(i, j)] is edge i->j's
    cost in agreed cells."""
    X, Y = _transition_rows(data_matrices)
    warm_start = []
    for j, candidates in enumerate(candidate_inputs):
        if not candidates:
            warm_start.append(({}, 0))
            continue
        warm_start.append(_fit_threshold_node(X, Y[:, j], X[:, j], candidates,
                                              {i: edge_costs.get((i, j), 0.0) for i in candidates}, max_degree))
    return warm_start


def _candidate_inputs(scaffold_network, allow_additional_edges):
    """Candidate inputs per node, sorted by vertex index. With additional edges allowed this is every vertex -
    the topology-free search. Without them inference may only keep or drop scaffold edges, so the candidates
    are the node's scaffold predecessors."""
    if allow_additional_edges:
        return [list(range(len(scaffold_network))) for _ in range(len(scaffold_network))]
    return [sorted(u.index for u in v.predecessors()) for v in scaffold_network.vertices]


def _max_indegree_cap(scaffold_network, max_indegree):
    """Per-node in-degree cap: -1 means the scaffold's (global) max in-degree, computed per network. Clamped
    to >= 1 so a degenerate edgeless scaffold (whose max in-degree is 0) doesn't force degree <= 0, which
    would collapse the whole model to constants."""
    cap = max_indegree if max_indegree != -1 else scaffold_network.max_in_degree()
    return max(1, cap)


def _heuristic_edge_costs(data_matrices, scaffold_network, candidate_inputs, included_edges_relative_weight,
                          added_edges_relative_weight):
    """The symmetric_topology objective's per-edge terms, converted to agreed cells (its per-cell unit), so the
    heuristic weighs an edge against the cells it explains on the same scale the MIP does."""
    n_vertices = len(scaffold_network)
    n_cells = sum(max(np.asarray(m).shape[0] - 1, 0) for m in data_matrices) * n_vertices
    edge_norm = float(len(scaffold_network.edges)) or 1.0
    scaffold_edge_set = set(scaffold_network.edges)
    edge_costs = {}
    for j in range(n_vertices):
        for i in candidate_inputs[j]:
            in_scaffold = (scaffold_network.vertices[i], scaffold_network.vertices[j]) in scaffold_edge_set
            weight = included_edges_relative_weight if in_scaffold else added_edges_relative_weight
            edge_costs[(i, j)] = -weight * n_cells / edge_norm
    return edge_costs


def _network_from_signs_and_thresholds(vertex_names, signs_and_thresholds):
    """The Network a list of per-vertex ({input vertex index: +1/-1 sign}, threshold) pairs describes. A
    vertex with no input holds its value (no function)."""
    model = Network(vertex_names=vertex_names, edges=[], vertex_functions=[None] * len(vertex_names))
    for j, (signs_by_index, threshold) in enumerate(signs_and_thresholds):
        for i in sorted(signs_by_index):
            model.edges.append((model.vertices[i], model.vertices[j]))
    for vertex in model.vertices:
        vertex.precomputed_predecessors = None
    for j, (signs_by_index, threshold) in enumerate(signs_and_thresholds):
        if signs_by_index:
            model.vertices[j].function = SymmetricThresholdFunction(
                [signs_by_index[i] > 0 for i in sorted(signs_by_index)], threshold)
    return model


def infer_symmetric_heuristic(data_matrices, scaffold_network, allow_additional_edges=False,
                              included_edges_relative_weight=1, added_edges_relative_weight=-1, max_indegree=-1,
                              **kwargs):
    """A symmetric threshold model fitted node by node on the anchored transitions, without a solver: the
    heuristic symmetric_topology uses as its warm start (_heuristic_warm_start_values), returned as the model.
    Candidate inputs, the in-degree cap and the edge weights are as in infer_unknown_topology_symmetric. Each
    node maximizes the one-step cells it predicts on the observed (noisy) rows, minus its edges' costs; there
    are no input flips or free-run trajectories, so allow_input_flips, flip_penalty, no_anchoring and the
    timeout are ignored."""
    data_matrices = list(data_matrices)
    candidate_inputs = _candidate_inputs(scaffold_network, allow_additional_edges)
    edge_costs = _heuristic_edge_costs(data_matrices, scaffold_network, candidate_inputs,
                                       included_edges_relative_weight, added_edges_relative_weight)
    fitted = _heuristic_warm_start_values(data_matrices, candidate_inputs, edge_costs,
                                          _max_indegree_cap(scaffold_network, max_indegree))
    return _network_from_signs_and_thresholds([v.name for v in scaffold_network.vertices], fitted)


def _start_trajectories(trajectories, data_matrices, vertex_names, warm_start, free_run):
    """Start values for the data's (denoised) states: the observed rows, unflipped, and in free-run mode the
    states the warm-start model then predicts from the first of them - so the MIP start fixes the states as
    well as the functions, rather than leaving them for Gurobi to complete. trajectories as
    add_matrices_as_model_paths returns them."""
    model = _network_from_signs_and_thresholds(vertex_names, warm_start)
    matrices = [np.asarray(m) for m in data_matrices if np.asarray(m).shape[0] >= 2]
    for trajectory, matrix in zip(trajectories, matrices):
        observed = [[int(v != 0) for v in row] for row in matrix]
        # free-run: rows past the first are the model's own predictions; anchored: every row is a data row
        rows = [observed[0]] + [list(state) for state in model.next_states(observed[0], len(trajectory) - 1)[1:]]             if free_run else observed
        for state, values in zip(trajectory, rows):
            for var, value in zip(state, values):
                if isinstance(var, gurobipy.Var):
                    var.Start = int(value)


# Default bound on the path IN_ATTRACTORS lets the model take from a matrix's last state back to its first.
DEFAULT_ATTRACTOR_MAX_PATH_LEN = 10


def infer_unknown_topology_symmetric_in_attractors(data_matrices, scaffold_network, **kwargs):
    """infer_unknown_topology_symmetric, with every data matrix assumed to lie on an attractor: from its
    denoised last state the model returns to its denoised first state within attractor_max_path_len steps
    (AttractorAssumption.IN_ATTRACTORS). Free-run with input flips only."""
    return infer_unknown_topology_symmetric(data_matrices, scaffold_network,
                                            attractor_assumption=AttractorAssumption.IN_ATTRACTORS, **kwargs)


def infer_unknown_topology_symmetric_contains_attractors(data_matrices, scaffold_network, **kwargs):
    """infer_unknown_topology_symmetric, with every data matrix assumed to end in, and fully cover, an
    attractor: its denoised last state repeats an earlier one (AttractorAssumption.CONTAINS_ATTRACTORS).
    Free-run with input flips only."""
    return infer_unknown_topology_symmetric(data_matrices, scaffold_network,
                                            attractor_assumption=AttractorAssumption.CONTAINS_ATTRACTORS, **kwargs)


def infer_unknown_topology_symmetric(data_matrices, scaffold_network, allow_additional_edges=False,
                                   included_edges_relative_weight=1, added_edges_relative_weight=-1,
                                   timeout_secs=None, log_file=None, allow_input_flips=False, flip_penalty=1.0,
                                   no_anchoring=False, gurobi_threads=0, warm_start_heuristic=False, max_indegree=-1, attractor_assumption=None,
                                   attractor_max_path_len=DEFAULT_ATTRACTOR_MAX_PATH_LEN, **kwargs):
    """
    Find a symmetric threshold model with best fit to data_matrices and scaffold_network,
    by finding both the Boolean function and the incoming edges for each node.
    The objective is 1 * x + included_edges_relative_weight * y + added_edges_relative_weight * z - flip_penalty * w,
    where
    x is the proportion of cells in data_matrices (up to first rows) that are explained correctly by the model
    y is the proportion of scaffold_network edges that are included in the model (in [0, 1]).
    z is the number of edges not in scaffold_network that are included in the model, divided by
    the number of edges in scaffold_network (this normalizes added edges to the scaffold's scale, so z is NOT
    bounded by 1 - it is a deliberately stronger sparsity penalty than dividing by the number of possible edges).
    w is the proportion of (flippable) input cells the model chose to flip, on the same per-cell scale as x.
    :param allow_additional_edges: if false, doesn't allow any edges that weren't in the scaffold network originally.
    :param included_edges_relative_weight:
    :param added_edges_relative_weight:
    :param function_type_restriction:
    :param data_matrices:
    :param scaffold_network:
    :param timeout_secs: timelimit to pass to the solver.
    :param allow_input_flips: if true, the solver may flip bits of the (noisy) input data, paying flip_penalty per flip.
    :param flip_penalty: penalty per flipped input cell, on the same [0, 1] per-cell scale as x.
        flip_penalty < 1 denoises toward a self-consistent trajectory, > 1 only corrects errors that fix more than
        one cell, and = 1 (the default) is the break-even point at which the chosen flips may be under-determined.
    :param no_anchoring: if true, score a free-running trajectory (each step is fed the model's own previous
        prediction) instead of independent one-step transitions anchored to the data.
    :param warm_start_heuristic: if true, seed the MIP with symmetric
        threshold functions fitted node by node on the anchored transitions (_heuristic_warm_start_values),
        with the data rows unflipped and, in free-run mode, the states those functions predict - a complete
        start, found in under a second.
    :param attractor_assumption: None, or an ilp_components.AttractorAssumption imposed as a hard constraint
        on every data matrix's denoised trajectory. Requires no_anchoring and allow_input_flips. The
        symmetric_topology_in_attractors / _contains_attractors methods are this function with it set.
    :param attractor_max_path_len: for AttractorAssumption.IN_ATTRACTORS, the longest path the model may take
        from a matrix's last state back to its first.
    :return:
    """
    ilp_components.check_attractor_assumption_compatible(attractor_assumption, no_anchoring, allow_input_flips)
    data_matrices = list(data_matrices)  # iterated by the warm start and by the main solve below

    n_vertices = len(scaffold_network)
    # Candidate inputs per node (see _candidate_inputs). Building the model over just the scaffold's
    # predecessors when additional edges aren't allowed is what keeps it proportional to the scaffold rather
    # than to n**2 (the earlier formulation built every pair and constrained the disallowed ones to zero, which
    # costs the memory before presolve can remove them). Sorted by vertex index, because add_path_to_model
    # zips a node's sign variables against its predecessors in index order.
    candidate_inputs = _candidate_inputs(scaffold_network, allow_additional_edges)
    candidate_position = [{index: position for position, index in enumerate(candidates)}
                          for candidates in candidate_inputs]
    scaffold_edge_set = set(scaffold_network.edges)  # membership test inside the per-pair loops below

    max_indeg_cap = _max_indegree_cap(scaffold_network, max_indegree)

    warm_start = None
    if warm_start_heuristic:
        edge_costs = _heuristic_edge_costs(data_matrices, scaffold_network, candidate_inputs,
                                           included_edges_relative_weight, added_edges_relative_weight)
        warm_start = _heuristic_warm_start_values(data_matrices, candidate_inputs, edge_costs, max_indeg_cap)

    # create function variables
    functions_variables = []
    model = gurobipy.Model()
    try:
        for i in range(n_vertices):
            # use ternary signs, zero means the input isn't used
            signs = [model.addVar(lb=-1, ub=1, vtype=gurobipy.GRB.INTEGER,
                                  name="vertex_{}_input_{}_sign_var".format(i, j))
                                  for j in candidate_inputs[i]]
            # technically can have a constant False function with all nodes as input, which
            # will have a threshold of len(candidate_inputs[i]) + 1, but that function is
            # (better) representable with less inputs.
            threshold = model.addVar(lb=0, ub=max_indeg_cap + 1, vtype=gurobipy.GRB.INTEGER,
                                     name="vertex_{}_threshold_var".format(i))
            functions_variables.append([signs, threshold])
        model.update()

        # The edge indicators (|sign|) are built before the data constraints, because the per-node "unused"
        # indicator derived from them is needed while the transitions are added: a node that ends up using
        # no input holds its value there, the way an input node does in Network.next_state.
        included_edges_indicators = []
        added_edges_indicators = []
        edge_indicators = dict()
        unused_node_indicators = []
        for j in range(n_vertices):
            for i in candidate_inputs[j]:
                edge_indicators[i, j] = model.addVar(
                    lb=0, ub=1, vtype=gurobipy.GRB.INTEGER, name="edge_{}_{}_indicator_var".format(i, j))
        for j in range(n_vertices):
            # a node with no candidate inputs can only hold its value, which add_path_to_model does for it
            # directly (it has no predecessors there), so it needs no indicator
            unused_node_indicators.append(model.addVar(
                vtype=gurobipy.GRB.BINARY, name="vertex_{}_unused_indicator_var".format(j))
                if candidate_inputs[j] else None)
        model.update()  # repeat loop after update, so that there's one update
        for j in range(n_vertices):
            for i in candidate_inputs[j]:
                model.addGenConstrAbs(edge_indicators[i, j],
                                      functions_variables[j][0][candidate_position[j][i]],
                                      name="edge_{}_{}_abs_var".format(i, j))
                if (scaffold_network.vertices[i], scaffold_network.vertices[j]) in scaffold_edge_set:
                    included_edges_indicators.append(edge_indicators[i, j])
                else:
                    # only reachable with allow_additional_edges: otherwise the candidates are exactly the
                    # scaffold's edges, so a non-scaffold pair has no variable to constrain in the first place
                    added_edges_indicators.append(edge_indicators[i, j])
        for j in range(n_vertices):
            # unused_node_indicators[j] is 1 iff node j's in-degree is 0, and is handed to the transition
            # constraints as the third element of the node's function variables.
            functions_variables[j].append(unused_node_indicators[j])
            if not candidate_inputs[j]:
                continue
            degree = gurobipy.quicksum(edge_indicators[i, j] for i in candidate_inputs[j])
            model.addConstr(degree <= len(candidate_inputs[j]) * (1 - unused_node_indicators[j]),
                            name="node_{}_unused_indicator_constraint_<=".format(j))
            model.addConstr(degree >= 1 - unused_node_indicators[j],
                            name="node_{}_unused_indicator_constraint_>=".format(j))

        candidate_network = Network(
            vertex_names=[v.name for v in scaffold_network.vertices],
            edges=[(scaffold_network.vertices[i].name, v.name)
                   for j, v in enumerate(scaffold_network.vertices) for i in candidate_inputs[j]],
            vertex_functions=[None for v in scaffold_network.vertices])
        # add_path_to_model pairs a node's sign variables with its predecessors positionally, so the two
        # orderings have to agree; both are by vertex index, but assert it rather than rely on it.
        for j, vertex in enumerate(candidate_network.vertices):
            assert [u.index for u in vertex.predecessors()] == candidate_inputs[j], \
                "candidate order does not match predecessor order for vertex {}".format(vertex.name)
        function_type_restrictions = [FunctionTypeRestriction.SYMMETRIC_THRESHOLD] * n_vertices
        matrix_agreement_indicators, flip_cost_terms, trajectories = ilp_components.add_matrices_as_model_paths(
           candidate_network, model, data_matrices,
           function_vars=functions_variables,
           model_to_data_sample_rate_ratio=1,
           function_type_restrictions=function_type_restrictions,
           model_addition_type=ModelAdditionType.INDICATORS,
           per_cell_indicators=True, allow_input_flips=allow_input_flips, no_anchoring=no_anchoring,
           max_indegree=max_indeg_cap, return_trajectories=True)
        if attractor_assumption is not None:
            ilp_components.add_attractor_assumption(
                candidate_network, model, trajectories, attractor_assumption, functions_variables,
                attractor_max_path_len=attractor_max_path_len,
                function_type_restrictions=function_type_restrictions, max_indegree=max_indeg_cap)
        del candidate_network  # only needed to build the model; free it before the (heavy) solve
        n_cells = float(len(matrix_agreement_indicators))
        data_agreement = gurobipy.quicksum(matrix_agreement_indicators) / n_cells

        # we need to prevent a non-standard representation of a constant function that doesn't zero the sign
        # variables. With actual degree d, the threshold is limited to [1, d] when d is positive, and to 0
        # when d = 0 (a node using no input holds its value, so it needs no threshold at all). Over integers
        # the pair (threshold <= d, n * threshold >= d) gives that: at d = 0 both force 0, and at
        # d >= 1 the second forces threshold >= 1 while the first allows every value up to d - including
        # threshold == d, i.e. an AND over the node's inputs.
        for j in range(n_vertices):
            if not candidate_inputs[j]:
                continue
            degree = gurobipy.quicksum(edge_indicators[i, j] for i in candidate_inputs[j])
            threshold = functions_variables[j][1]
            model.addConstr(degree <= max_indeg_cap, name="node_{}_max_indegree_constraint".format(j))
            model.addConstr(threshold <= degree, name="node_{}_threshold_constraint_<=".format(j))
            # the multiplier only has to dominate the degree, which the node's own candidate count does
            model.addConstr(max_indeg_cap * threshold >= degree,
                            name="node_{}_threshold_constraint_>=".format(j))

        # both edge terms are normalized by the number of scaffold edges (see docstring); guard the rare
        # empty-scaffold case, where there are no included edges to normalize by and these terms are vacuous.
        n_scaffold_edges = float(len(included_edges_indicators))
        edge_norm = n_scaffold_edges if n_scaffold_edges > 0 else 1.0
        included_edges_agreement = included_edges_relative_weight * gurobipy.quicksum(included_edges_indicators) / edge_norm
        added_edges_agreement = added_edges_relative_weight * gurobipy.quicksum(added_edges_indicators) / edge_norm

        objective = data_agreement + included_edges_agreement + added_edges_agreement
        if allow_input_flips and flip_cost_terms:
            objective = objective - flip_penalty * gurobipy.quicksum(flip_cost_terms) / n_cells

        if warm_start is not None:
            # seed signs/thresholds and the trajectories' states from the warm start: a complete MIP start
            for i in range(n_vertices):
                signs_by_index, threshold_start = warm_start[i]
                signs_vars, threshold_var = functions_variables[i][0], functions_variables[i][1]
                # the warm start keys signs by vertex index; map through to this node's candidate positions
                for position, j in enumerate(candidate_inputs[i]):
                    signs_vars[position].Start = signs_by_index.get(j, 0)
                threshold_var.Start = min(threshold_start, len(candidate_inputs[i]))
            _start_trajectories(trajectories, data_matrices, [v.name for v in scaffold_network.vertices],
                                warm_start, free_run=no_anchoring)

        if log_file is not None:
            model.Params.LogFile = log_file
        if timeout_secs is not None:
            model.Params.TimeLimit = timeout_secs
        if gurobi_threads:
            model.Params.Threads = gurobi_threads
        model.Params.MIPFocus = 1

        model.setObjective(objective, sense=gurobipy.GRB.MAXIMIZE)
        model.optimize()

        # set Boolean model found
        inferred_model = Network(vertex_names=[v.name for v in scaffold_network.vertices], edges=[],
                                 vertex_functions=[None for v in scaffold_network.vertices])
        for i in range(len(inferred_model)):
            signs, threshold = functions_variables[i][0], functions_variables[i][1]
            float_signs = [get_value_of_gurobi_entity(sign) for sign in signs]
            signs = [int(round(sign, 3)) for sign in float_signs]
            # gurobi can give "almost" integer values even for variables defined as integer type
            assert signs == [round(s, 3) for s in float_signs]
            float_threshold = get_value_of_gurobi_entity(threshold)
            threshold = int(round(float_threshold, 3))
            assert threshold == round(float_threshold, 3)
            # assert threshold doesn't imply a constant function (that shouldn't be possible with the way we modelled this)

            assert threshold <= sum(abs(s) for s in signs)

            for position, j in enumerate(candidate_inputs[i]):
                if signs[position] != 0:
                    inferred_model.edges.append((inferred_model.vertices[j], inferred_model.vertices[i]))
            signs = [s for s in signs if s != 0]
            # threshold = max(-len(signs), min(threshold, len(signs) + 1))
            # a node that ended up using no input holds its value, which is what an absent function means in
            # Network.next_state - the same semantics the unused-node indicator gave it in the model
            inferred_model.vertices[i].function = (
                SymmetricThresholdFunction(signs, threshold) if signs else None)
    finally:
        model.dispose()
    return inferred_model
