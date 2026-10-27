import gurobipy
from inference import ilp_components
from attractor_learning.graphs import Network
from inference.ilp_components import AttractorAssumption
from attractor_learning.logic import SymmetricThresholdFunction
import numpy as np


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


# Share of the solve time limit spent in Gurobi's NoRel heuristic, which searches for solutions before the root
# LP is solved. On dense instances (an NK K = 15 network, the larger dynamic toy models with added edges) the root
# LP alone can take the whole limit, leaving only the warm start; on the dynamic toy models at 30 s, 25% gave a
# better solution in 10 of the 20 runs that timed out and a worse one in 3, at the cost of that time on runs that
# would have solved quickly.
NOREL_TIME_FRACTION = 0.25


def infer_unknown_topology_symmetric_norel(data_matrices, scaffold_network, **kwargs):
    """infer_unknown_topology_symmetric with the whole time limit given to Gurobi's NoRel heuristic
    (norel_time_fraction = 1) instead of NOREL_TIME_FRACTION of it, for instances whose root LP cannot be
    solved within the limit anyway; everything else, the warm start included, is unchanged."""
    return infer_unknown_topology_symmetric(data_matrices, scaffold_network, norel_time_fraction=1.0, **kwargs)


def infer_unknown_topology_symmetric(data_matrices, scaffold_network, allow_additional_edges=False,
                                   included_edges_relative_weight=1, added_edges_relative_weight=-1,
                                   timeout_secs=None, log_file=None, allow_input_flips=False, flip_penalty=1.0,
                                   no_anchoring=False, gurobi_threads=0, warm_start_heuristic=False, max_indegree=-1, attractor_assumption=None,
                                   attractor_max_path_len=DEFAULT_ATTRACTOR_MAX_PATH_LEN,
                                   norel_time_fraction=NOREL_TIME_FRACTION, **kwargs):
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

    Encoding: each candidate edge i -> j is two binaries, p (activating) and q (inhibiting) with
    p + q + u_j <= 1, where u_j is node j's unused indicator (1 iff it uses no input), so the sign is p - q and
    the edge indicator p + q. Node j's in-degree is capped by max_indeg_cap * (1 - u_j) and its threshold lies in
    [1 - u_j, degree]. Whether an input agrees with its edge's sign is p or q itself for a constant (data) input,
    and the convex hull of p*v + q*(1 - v) for a variable one; a predicted cell scored against a constant data
    cell enters the objective as the cell itself (or its complement). See ilp.add_path_to_model and
    ilp.add_state_agreement_terms. With a time limit, the first norel_time_fraction of it (NOREL_TIME_FRACTION by
    default) is spent in Gurobi's NoRel heuristic.
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
    # than to n**2. Sorted by vertex index.
    candidate_inputs = _candidate_inputs(scaffold_network, allow_additional_edges)
    scaffold_edge_set = set(scaffold_network.edges)  # membership test inside the per-pair loop below

    max_indeg_cap = _max_indegree_cap(scaffold_network, max_indegree)

    warm_start = None
    if warm_start_heuristic:
        edge_costs = _heuristic_edge_costs(data_matrices, scaffold_network, candidate_inputs,
                                           included_edges_relative_weight, added_edges_relative_weight)
        warm_start = _heuristic_warm_start_values(data_matrices, candidate_inputs, edge_costs, max_indeg_cap)

    model = gurobipy.Model()
    try:
        # per node: (positive edge binaries, negative edge binaries, threshold, unused indicator), the edge
        # binaries in candidate order
        functions = []
        included_edges, added_edges = [], []
        for j in range(n_vertices):
            candidates = candidate_inputs[j]
            positives = [model.addVar(vtype=gurobipy.GRB.BINARY, name="vertex_{}_input_{}_pos_var".format(j, i))
                         for i in candidates]
            negatives = [model.addVar(vtype=gurobipy.GRB.BINARY, name="vertex_{}_input_{}_neg_var".format(j, i))
                         for i in candidates]
            # threshold <= degree <= max_indeg_cap (constraints below); a node without candidates has no
            # function, so its threshold is fixed at 0
            threshold = model.addVar(lb=0, ub=min(len(candidates), max_indeg_cap), vtype=gurobipy.GRB.INTEGER,
                                     name="vertex_{}_threshold_var".format(j))
            # a node with no candidate inputs can only hold its value, which add_path_to_model does for it
            # directly, so it needs no indicator
            unused = model.addVar(vtype=gurobipy.GRB.BINARY, name="vertex_{}_unused_indicator_var".format(j)) \
                if candidates else None
            functions.append((positives, negatives, threshold, unused))
            if not candidates:
                continue
            for i, p, q in zip(candidates, positives, negatives):
                # one sign at most, and none at all on an unused node
                model.addConstr(p + q + unused <= 1, name="edge_{}_{}_unused_constraint".format(i, j))
                if (scaffold_network.vertices[i], scaffold_network.vertices[j]) in scaffold_edge_set:
                    included_edges.append(p + q)
                else:
                    # only reachable with allow_additional_edges: otherwise the candidates are exactly the
                    # scaffold's edges
                    added_edges.append(p + q)
            degree = gurobipy.quicksum(positives) + gurobipy.quicksum(negatives)
            # u_j is 1 iff node j's in-degree is 0
            model.addConstr(degree <= max_indeg_cap * (1 - unused), name="node_{}_max_indegree_constraint".format(j))
            model.addConstr(degree >= 1 - unused, name="node_{}_unused_indicator_constraint_>=".format(j))
            # we need to prevent a non-standard representation of a constant function that doesn't zero the
            # edges: the threshold lies in [1, degree] for a used node, and is 0 for an unused one (which holds
            # its value, so it needs no threshold at all). threshold == degree is an AND over the node's inputs.
            model.addConstr(threshold <= degree, name="node_{}_threshold_constraint_<=".format(j))
            model.addConstr(threshold >= 1 - unused, name="node_{}_threshold_constraint_>=".format(j))

        candidate_network = Network(
            vertex_names=[v.name for v in scaffold_network.vertices],
            edges=[(scaffold_network.vertices[i].name, v.name)
                   for j, v in enumerate(scaffold_network.vertices) for i in candidate_inputs[j]],
            vertex_functions=[None for v in scaffold_network.vertices])
        # add_path_to_model pairs a node's edge variables with its predecessors positionally, so the two
        # orderings have to agree; both are by vertex index, but assert it rather than rely on it.
        for j, vertex in enumerate(candidate_network.vertices):
            assert [u.index for u in vertex.predecessors()] == candidate_inputs[j], \
                "candidate order does not match predecessor order for vertex {}".format(vertex.name)
        agreement, flip_cost_terms, trajectories = ilp_components.add_matrices_as_model_paths(
            candidate_network, model, data_matrices, functions, max_indegree=max_indeg_cap,
            allow_input_flips=allow_input_flips, no_anchoring=no_anchoring)
        if attractor_assumption is not None:
            ilp_components.add_attractor_assumption(
                candidate_network, model, trajectories, attractor_assumption, functions,
                attractor_max_path_len=attractor_max_path_len, max_indegree=max_indeg_cap)
        del candidate_network  # only needed to build the model; free it before the (heavy) solve
        n_cells = float(len(agreement))
        data_agreement = gurobipy.quicksum(agreement) / n_cells

        # both edge terms are normalized by the number of scaffold edges (see docstring); guard the rare
        # empty-scaffold case, where there are no included edges to normalize by and these terms are vacuous.
        n_scaffold_edges = float(len(included_edges))
        edge_norm = n_scaffold_edges if n_scaffold_edges > 0 else 1.0
        included_edges_agreement = included_edges_relative_weight * gurobipy.quicksum(included_edges) / edge_norm
        added_edges_agreement = added_edges_relative_weight * gurobipy.quicksum(added_edges) / edge_norm

        objective = data_agreement + included_edges_agreement + added_edges_agreement
        if allow_input_flips and flip_cost_terms:
            objective = objective - flip_penalty * gurobipy.quicksum(flip_cost_terms) / n_cells

        if warm_start is not None:
            # seed the edges/thresholds and the trajectories' states from the warm start: a complete MIP start
            for j in range(n_vertices):
                signs_by_index, threshold_start = warm_start[j]
                positives, negatives, threshold, _ = functions[j]
                # the warm start keys signs by vertex index; map through to this node's candidate positions
                for i, p, q in zip(candidate_inputs[j], positives, negatives):
                    sign = signs_by_index.get(i, 0)
                    p.Start, q.Start = int(sign > 0), int(sign < 0)
                threshold.Start = min(threshold_start, len(candidate_inputs[j]))
            _start_trajectories(trajectories, data_matrices, [v.name for v in scaffold_network.vertices],
                                warm_start, free_run=no_anchoring)

        if log_file is not None:
            model.Params.LogFile = log_file
        if timeout_secs is not None:
            model.Params.TimeLimit = timeout_secs
            model.Params.NoRelHeurTime = norel_time_fraction * timeout_secs
        if gurobi_threads:
            model.Params.Threads = gurobi_threads
        # model.Params.MIPFocus = 1

        model.setObjective(objective, sense=gurobipy.GRB.MAXIMIZE)
        model.optimize()

        # set Boolean model found
        inferred_model = Network(vertex_names=[v.name for v in scaffold_network.vertices], edges=[],
                                 vertex_functions=[None for v in scaffold_network.vertices])
        for j in range(n_vertices):
            positives, negatives, threshold, _ = functions[j]
            signs = [int(round(p.X - q.X)) for p, q in zip(positives, negatives)]
            float_threshold = threshold.X
            threshold_value = int(round(float_threshold, 3))
            # gurobi can give "almost" integer values even for variables defined as integer type
            assert threshold_value == round(float_threshold, 3)
            assert threshold_value <= sum(abs(s) for s in signs)
            for i, sign in zip(candidate_inputs[j], signs):
                if sign != 0:
                    inferred_model.edges.append((inferred_model.vertices[i], inferred_model.vertices[j]))
            signs = [s for s in signs if s != 0]
            # a node that ended up using no input holds its value, which is what an absent function means in
            # Network.next_state - the same semantics the unused-node indicator gave it in the model
            inferred_model.vertices[j].function = SymmetricThresholdFunction(signs, threshold_value) if signs else None
    finally:
        model.dispose()
    return inferred_model
