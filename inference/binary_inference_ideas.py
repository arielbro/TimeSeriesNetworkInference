import attractor_learning.graphs
import attractor_learning.ilp
import gurobipy
from inference import ilp_components
from attractor_learning.graphs import FunctionTypeRestriction, Network
from inference.ilp_components import ModelAdditionType, get_value_of_gurobi_entity
from attractor_learning.logic import BooleanSymbolicFunc, SymmetricThresholdFunction
import itertools
import time


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


def _scaffold_warm_start_values(data_matrices, scaffold_network, timeout_secs, log_file,
                                allow_input_flips, flip_penalty, no_anchoring, gurobi_threads):
    """Solve the known-topology symmetric ILP (inputs fixed to the scaffold) and read back, per vertex, a
    ({input vertex index: +1/-1 sign}, threshold) pair to seed the unknown-topology model's variables.
    The sign of edge j->i is stored under index j (so it maps directly onto signs_vars[j] of vertex i).
    Thresholds are clamped into the unknown model's feasible [1, degree] range for nodes with inputs: the
    known model may return a constant function (threshold 0 or degree+1) while keeping nonzero signs, and
    feeding such a threshold back would make the whole seeded MIP start infeasible (Gurobi would discard it).
    Degree-0 (input) nodes get no signs and threshold 0."""
    scaffold_model = infer_known_topology_symmetric(
        data_matrices, scaffold_network, timeout_secs=timeout_secs, log_file=log_file,
        allow_input_flips=allow_input_flips, flip_penalty=flip_penalty,
        no_anchoring=no_anchoring, gurobi_threads=gurobi_threads)
    warm_start = []
    for vertex in scaffold_model.vertices:
        func = vertex.function
        if func is None:
            warm_start.append(({}, 0))
            continue
        predecessors = vertex.predecessors()
        # func.signs[k] aligns with predecessors()[k] (both ordered by vertex index in the known model)
        signs_by_index = {pred.index: (1 if sign else -1)
                          for pred, sign in zip(predecessors, func.signs)}
        degree = len(predecessors)
        threshold = min(max(func.threshold, 1), degree)
        warm_start.append((signs_by_index, threshold))
    return warm_start


def infer_unknown_topology_symmetric(data_matrices, scaffold_network, allow_additional_edges=False,
                                   included_edges_relative_weight=1, added_edges_relative_weight=-1,
                                   timeout_secs=None, log_file=None, allow_input_flips=False, flip_penalty=1.0,
                                   no_anchoring=False, gurobi_threads=0,
                                   warm_start_from_scaffold=False, warm_start_time_frac=0.2, **kwargs):
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
    :param warm_start_from_scaffold: if true, first solve the known-topology symmetric ILP (inputs fixed to the
        scaffold) and use the signs/thresholds it learns as a MIP start for this (unknown-topology) solve.
    :param warm_start_time_frac: fraction of timeout_secs given to the scaffold (warm-start) solve; this solve's
        actual elapsed time is then subtracted from timeout_secs to budget the main solve. Ignored when
        warm_start_from_scaffold is False or timeout_secs is None.
    :return:
    """
    warm_start = None
    main_timeout_secs = timeout_secs
    if warm_start_from_scaffold:
        data_matrices = list(data_matrices)  # iterated by both the scaffold solve and the main solve below
        scaffold_timeout = timeout_secs * warm_start_time_frac if timeout_secs is not None else None
        scaffold_start = time.time()
        warm_start = _scaffold_warm_start_values(
            data_matrices, scaffold_network, timeout_secs=scaffold_timeout, log_file=log_file,
            allow_input_flips=allow_input_flips, flip_penalty=flip_penalty,
            no_anchoring=no_anchoring, gurobi_threads=gurobi_threads)
        if timeout_secs is not None:
            # give the (harder) unknown-topology solve whatever of the budget the scaffold solve left unused
            main_timeout_secs = max(timeout_secs - (time.time() - scaffold_start), 0.0)

    # create function variables
    functions_variables = []
    model = gurobipy.Model()
    try:
        for i in range(len(scaffold_network)):
            # use ternary signs, zero means the input isn't used
            signs = [model.addVar(lb=-1, ub=1, vtype=gurobipy.GRB.INTEGER,
                                  name="vertex_{}_input_{}_sign_var".format(i, j))
                                  for j in range(len(scaffold_network))]
            # technically can have a constant False function with all nodes as input, which
            # will have a threshold of len(scaffold_network) + 1, but that function is
            # (better) representable with less inputs.
            threshold = model.addVar(lb=0, ub=len(scaffold_network), vtype=gurobipy.GRB.INTEGER,
                                     name="vertex_{}_threshold_var".format(i))
            functions_variables.append([signs, threshold])
        model.update()

        full_network = Network(vertex_names=[v.name for v in scaffold_network.vertices],
                               edges=[(u.name, v.name) for u, v in itertools.product(scaffold_network.vertices, repeat=2)],
                               vertex_functions=[None for v in scaffold_network.vertices])
        matrix_agreement_indicators, flip_cost_terms = ilp_components.add_matrices_as_model_paths(
           full_network, model, data_matrices,
           function_vars=functions_variables,
           model_to_data_sample_rate_ratio=1,
           function_type_restrictions=[FunctionTypeRestriction.SYMMETRIC_THRESHOLD] * len(scaffold_network),
           model_addition_type=ModelAdditionType.INDICATORS,
           per_cell_indicators=True, allow_input_flips=allow_input_flips, no_anchoring=no_anchoring)
        del full_network  # only needed to build the model; free it before the (heavy) solve
        n_cells = float(len(matrix_agreement_indicators))
        data_agreement = gurobipy.quicksum(matrix_agreement_indicators) / n_cells

        included_edges_indicators = []
        added_edges_indicators = []
        edge_indicators = dict()
        for i in range(len(scaffold_network)):
            for j in range(len(scaffold_network)):
                edge_indicators[i, j] = model.addVar(
                    lb=0, ub=1, vtype=gurobipy.GRB.INTEGER, name="edge_{}_{}_indicator_var".format(i, j))
        model.update()  # repeat loop after update, so that there's one update
        for i in range(len(scaffold_network)):
            for j in range(len(scaffold_network)):
                model.addGenConstrAbs(edge_indicators[i, j], functions_variables[j][0][i],
                                      name="edge_{}_{}_abs_var".format(i, j))
                if (scaffold_network.vertices[i], scaffold_network.vertices[j]) in scaffold_network.edges:
                    included_edges_indicators.append(edge_indicators[i, j])
                else:
                    if allow_additional_edges:
                        added_edges_indicators.append(edge_indicators[i, j])
                    else:
                        model.addConstr(edge_indicators[i, j] == 0)
                        added_edges_indicators.append(0)

        # we need to prevent a non-standard representation of a constant function that doesn't zero the sign variables.
        # With actual degree d, we want to limit the threshold to range [1, d] if d is positive
        # and [0, d + 1] if d=0 (to allow for the constant functions). Let n be the network size.
        # [d/n, d + (1 - d/n)], given an integer threshold, represents both cases, although it might result
        # in bad linear relaxations. An alternative (not implemented now) is to create some
        # indicator variable and conditional constraints.
        # BASICALLY - INDICATOR FOR WHETHER d=0 AND CONDITIONAL CONSTRAINTS.
        for j in range(len(scaffold_network)):
            degree = gurobipy.quicksum(edge_indicators[i, j] for i in range(len(scaffold_network)))
            threshold = functions_variables[j][1]
            model.addConstr(threshold <= degree + (1 - degree) / len(scaffold_network),
                            name="node_{}_threshold_constraint_<=".format(j))
            model.addConstr(len(scaffold_network) * threshold >= degree, name="node_{}_threshold_constraint_>=".format(j))

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
            # seed signs/threshold from the scaffold solve; Gurobi completes the rest of the partial MIP start
            for i in range(len(scaffold_network)):
                signs_by_index, threshold_start = warm_start[i]
                signs_vars, threshold_var = functions_variables[i]
                for j in range(len(scaffold_network)):
                    signs_vars[j].Start = signs_by_index.get(j, 0)
                threshold_var.Start = threshold_start

        if log_file is not None:
            model.Params.LogFile = log_file
        if main_timeout_secs is not None:
            model.Params.TimeLimit = main_timeout_secs
        if gurobi_threads:
            model.Params.Threads = gurobi_threads
        model.Params.MIPFocus = 1

        model.setObjective(objective, sense=gurobipy.GRB.MAXIMIZE)
        model.optimize()

        # set Boolean model found
        inferred_model = Network(vertex_names=[v.name for v in scaffold_network.vertices], edges=[],
                                 vertex_functions=[None for v in scaffold_network.vertices])
        for i in range(len(inferred_model)):
            signs, threshold = functions_variables[i]
            float_signs = [get_value_of_gurobi_entity(sign) for sign in signs]
            signs = [int(round(sign, 3)) for sign in float_signs]
            # gurobi can give "almost" integer values even for variables defined as integer type
            assert signs == [round(s, 3) for s in float_signs]
            float_threshold = get_value_of_gurobi_entity(threshold)
            threshold = int(round(float_threshold, 3))
            assert threshold == round(float_threshold, 3)
            # assert threshold doesn't imply a constant function (that shouldn't be possible with the way we modelled this)

            assert threshold <= sum(abs(s) for s in signs)

            for j in range(len(inferred_model.vertices)):
                if signs[j] != 0:
                    inferred_model.edges.append((inferred_model.vertices[j], inferred_model.vertices[i]))
            signs = [s for s in signs if s != 0]
            # threshold = max(-len(signs), min(threshold, len(signs) + 1))
            func = SymmetricThresholdFunction(signs, threshold)
            inferred_model.vertices[i].function = func
    finally:
        model.dispose()
    return inferred_model
