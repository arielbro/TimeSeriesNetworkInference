from attractor_learning import ilp
from enum import Enum
import logging
import gurobipy


class ModelAdditionType(Enum):
    CONSTRAINT = 1
    INDICATORS = 2


slice_size = 15

logger = logging.getLogger()
logger.info("slice_size={}".format(slice_size))


def _row_input_state(matrix, denoised_rows, t):
    """The state for data row t: its denoised variable list if that row is flippable, else the observed
    (constant) row."""
    if denoised_rows is not None and t in denoised_rows:
        return denoised_rows[t]
    return matrix[t]


class AttractorAssumption(Enum):
    """What the data is assumed to say about the model's attractors, beyond its transitions. Each is a hard
    constraint per data matrix on that matrix's denoised trajectory y_1, ..., y_T (see
    add_matrices_as_model_paths' return_trajectories):
      IN_ATTRACTORS        the matrix lies on an attractor: from y_T, the model returns to y_1 within
                           attractor_max_path_len steps (add_attractor_return_path).
      CONTAINS_ATTRACTORS  the matrix ends in, and fully covers, an attractor: y_T equals an earlier y_r
                           (add_attractor_repeat)."""
    IN_ATTRACTORS = 1
    CONTAINS_ATTRACTORS = 2


def check_attractor_assumption_compatible(attractor_assumption, no_anchoring, allow_input_flips):
    """The attractor assumptions are stated on the denoised trajectory the model itself produces, which only
    the free-run (no_anchoring) mode with input flips gives: y_1 is the denoised first row and y_2..y_T the
    model's own predictions from it. Anchored, the y's would be data rows the model's transitions are only
    softly scored against, so the assumption would constrain the data rather than the model; without flips,
    y_1 would be the observed (noisy) row, and a single noisy bit could make the assumption unsatisfiable."""
    if attractor_assumption is None:
        return
    if not no_anchoring:
        raise ValueError("attractor assumption {} requires free-run inference (no_anchoring=True); the anchored "
                         "mode is not supported".format(attractor_assumption.name))
    if not allow_input_flips:
        raise ValueError("attractor assumption {} requires allow_input_flips=True".format(attractor_assumption.name))


def _add_conditional_state_equality(model, indicator, first_state, second_state, name_prefix):
    """indicator = 1 -> first_state == second_state, cell by cell. Either state may hold variables or
    constants (a data row); all values are binary, so |a - b| <= 1 - indicator is exact."""
    for i, (first, second) in enumerate(zip(first_state, second_state)):
        first = first if isinstance(first, (gurobipy.Var, gurobipy.LinExpr)) else int(round(float(first)))
        second = second if isinstance(second, (gurobipy.Var, gurobipy.LinExpr)) else int(round(float(second)))
        model.addConstr(first - second <= 1 - indicator, name="{}_equal_<=_{}".format(name_prefix, i))
        model.addConstr(second - first <= 1 - indicator, name="{}_equal_>=_{}".format(name_prefix, i))


def add_attractor_return_path(graph, model, trajectory, function_vars, max_path_len,
                              function_type_restrictions=None, max_indegree=None, name_prefix=""):
    """IN_ATTRACTORS for one matrix: the model, run on from the trajectory's last state y_T, reaches its first
    state y_1 again after some l in [1, max_path_len] steps - so the whole trajectory lies on one cycle.

    The path z_1..z_L (L = max_path_len) is rolled out with add_path_to_model under the same function
    variables as the data, and one binary hit_l per step says z_l == y_1; at least one has to hold. The
    length is variable through that disjunction alone, so the attractor ILP's activity variables are not
    needed. Returns the hit indicators."""
    path = ilp.add_path_to_model(graph, model, path_len=max_path_len, first_state_vars=trajectory[-1],
                                 model_f_vars=function_vars, v_funcs_restrictions=function_type_restrictions,
                                 max_indegree=max_indegree, name_prefix="{}_return_path".format(name_prefix))
    hits = [model.addVar(vtype=gurobipy.GRB.BINARY, name="{}_returns_after_{}".format(name_prefix, l + 1))
            for l in range(max_path_len)]
    model.update()
    for l, (state, hit) in enumerate(zip(path, hits)):
        _add_conditional_state_equality(model, hit, state, trajectory[0],
                                        "{}_return_after_{}".format(name_prefix, l + 1))
    model.addConstr(gurobipy.quicksum(hits) >= 1, name="{}_returns_to_start".format(name_prefix))
    return hits


def add_attractor_repeat(model, trajectory, name_prefix=""):
    """CONTAINS_ATTRACTORS for one matrix: its last state y_T equals some earlier state y_r (r < T), so the
    trajectory has closed a cycle by its end - it ends in an attractor and covers all of it. One binary per
    candidate r, at least one of which has to hold. Returns those indicators."""
    hits = [model.addVar(vtype=gurobipy.GRB.BINARY, name="{}_repeats_row_{}".format(name_prefix, r))
            for r in range(len(trajectory) - 1)]
    model.update()
    for r, hit in enumerate(hits):
        _add_conditional_state_equality(model, hit, trajectory[-1], trajectory[r],
                                        "{}_repeat_row_{}".format(name_prefix, r))
    model.addConstr(gurobipy.quicksum(hits) >= 1, name="{}_ends_in_attractor".format(name_prefix))
    return hits


def add_attractor_assumption(graph, model, trajectories, attractor_assumption, function_vars,
                             attractor_max_path_len=None, function_type_restrictions=None, max_indegree=None):
    """Adds attractor_assumption's constraint for every trajectory (see AttractorAssumption)."""
    for index, trajectory in enumerate(trajectories):
        prefix = "matrix_{}_attractor".format(index)
        if attractor_assumption == AttractorAssumption.IN_ATTRACTORS:
            if attractor_max_path_len is None or attractor_max_path_len < 1:
                raise ValueError("IN_ATTRACTORS needs attractor_max_path_len >= 1, got {}".format(
                    attractor_max_path_len))
            add_attractor_return_path(graph, model, trajectory, function_vars, attractor_max_path_len,
                                      function_type_restrictions=function_type_restrictions,
                                      max_indegree=max_indegree, name_prefix=prefix)
        elif attractor_assumption == AttractorAssumption.CONTAINS_ATTRACTORS:
            add_attractor_repeat(model, trajectory, name_prefix=prefix)
        else:
            raise ValueError("Unrecognized attractor assumption {}".format(attractor_assumption))


def add_matrices_as_model_paths(graph, model, data_matrices, function_vars, model_to_data_sample_rate_ratio=1,
                                function_type_restrictions=None, model_addition_type=ModelAdditionType.CONSTRAINT,
                                per_cell_indicators=False, allow_input_flips=False, no_anchoring=False,
                                max_indegree=None, return_trajectories=False):
    """
    Adds the data matrices to the model as model transitions.

    With return_trajectories (INDICATORS mode only), a third element is returned: per matrix with at least two
    rows, the list of its states y_1..y_T as the model sees them. In free-run mode that is the (denoised, with
    flips) first row followed by the model's own predictions; anchored, it is the (denoised) data rows, with
    the last row the observed one, since anchored mode never feeds the last row to the model.

    In INDICATORS mode returns a tuple (agreement_indicators, flip_cost_terms):
      - agreement_indicators: binary indicators that are 1 iff the model correctly explains a transition.
        With per_cell_indicators, there is one indicator per (predicted state, node) cell (cell-wise scoring);
        otherwise one indicator per whole predicted state (legacy, row-wise scoring).
      - flip_cost_terms: when allow_input_flips, one term per flippable data cell, each equal to 1 iff the
        model chose to flip that observed bit (0 otherwise); empty otherwise. The caller is responsible for
        penalizing their sum in the objective.

    no_anchoring selects the transition structure:
      - anchored (default): each consecutive data pair (row t, row t+1) is an independent one-step transition,
        i.e. the model input at every step is the observed (or denoised) data state.
      - free-run: a single trajectory is rolled out from the first row, feeding the model its own predicted
        state at each subsequent step; predicted state t is compared against observed row t. Either way the
        cell count (and hence normalization) is (#rows - 1) * #nodes per matrix.

    allow_input_flips lets the solver treat the observed data as noisy: each cell that is used as a model
    input gets a "denoised" binary variable, with a flip cost when it differs from the observation. Flipping
    is only meaningful where the value feeds the dynamics, so denoising is restricted to model-input rows:
      - anchored: every row except the last (each interior denoised row is reused as both the source of its
        transition and the target of the previous one, so it cannot be flipped for free);
      - free-run: only the initial row, since later inputs are the model's own predictions and the remaining
        observed rows are pure comparison targets (flipping them would just trivially match the prediction).
    """
    if model_to_data_sample_rate_ratio != 1:
        raise NotImplementedError()
    if allow_input_flips and model_addition_type != ModelAdditionType.INDICATORS:
        raise NotImplementedError("input flips are only supported with indicator (soft) transitions")
    if no_anchoring and model_addition_type != ModelAdditionType.INDICATORS:
        raise NotImplementedError("free-run (no_anchoring) modeling is only supported with indicator (soft) "
                                  "transitions")

    if return_trajectories and model_addition_type != ModelAdditionType.INDICATORS:
        raise NotImplementedError("trajectories are only returned with indicator (soft) transitions")

    n = len(graph.vertices)
    indicators = []
    flip_cost_terms = []
    trajectories = []
    for matrix_index, matrix in enumerate(data_matrices):
        n_rows = len(matrix)
        if n_rows < 2:
            continue  # no transitions to add
        # rows that feed the model as inputs (and so may be denoised): in free-run only the initial row, in
        # anchored mode every row but the last.
        input_row_indices = [0] if no_anchoring else list(range(n_rows - 1))

        denoised_rows = None
        if allow_input_flips:
            denoised_rows = {t: [model.addVar(vtype=gurobipy.GRB.BINARY,
                                              name="denoised_matrix_{}_row_{}_node_{}".format(matrix_index, t, i))
                                 for i in range(n)]
                             for t in input_row_indices}
            model.update()
            for t in input_row_indices:
                for i in range(n):
                    observed = int(round(float(matrix[t][i])))
                    z = denoised_rows[t][i]
                    # flip cost is 1 iff the denoised value differs from the observed one
                    flip_cost_terms.append(z if observed == 0 else (1 - z))

        def add_agreement(predicted_state_vars, target_state, t):
            # note that it's a weak constraint indicator, i.e. indicator -> constraint
            indicator = ilp.add_state_equality_indicator(model, predicted_state_vars, target_state, force_equal=False,
                                  prefix="add_matrices_matrix_index_{}_row_{}".format(matrix_index, t),
                                  return_per_index=per_cell_indicators)
            if per_cell_indicators:
                indicators.extend(indicator)
            else:
                indicators.append(indicator)

        if model_addition_type == ModelAdditionType.CONSTRAINT:
            for t in range(n_rows - 1):
                ilp.add_path_to_model(graph, model, path_len=1,
                                      first_state_vars=_row_input_state(matrix, denoised_rows, t), model_f_vars=None,
                                      last_state_vars=matrix[t + 1], v_funcs_restrictions=function_type_restrictions, max_indegree=max_indegree,
                                      name_prefix="matrix_{}_time_{}".format(matrix_index, t))
        elif model_addition_type == ModelAdditionType.INDICATORS and no_anchoring:
            # single free-running trajectory: the model is fed its own predictions after the initial state
            predicted_states = ilp.add_path_to_model(graph, model, path_len=n_rows - 1,
                                  first_state_vars=_row_input_state(matrix, denoised_rows, 0), last_state_vars=None,
                                  v_funcs_restrictions=function_type_restrictions, max_indegree=max_indegree, model_f_vars=function_vars,
                                  name_prefix="matrix_{}_freerun".format(matrix_index))
            for t in range(1, n_rows):
                add_agreement(predicted_states[t - 1], matrix[t], t)
            trajectories.append([_row_input_state(matrix, denoised_rows, 0)] + list(predicted_states))
        elif model_addition_type == ModelAdditionType.INDICATORS:
            for t in range(n_rows - 1):
                next_step_vars = ilp.add_path_to_model(graph, model, path_len=1,
                                      first_state_vars=_row_input_state(matrix, denoised_rows, t),
                                      last_state_vars=None, v_funcs_restrictions=function_type_restrictions, max_indegree=max_indegree,
                                      model_f_vars=function_vars,
                                      name_prefix="matrix_{}_time_{}".format(matrix_index, t))[0]
                # the target is the denoised next row where it exists (i.e. the next row is itself a model
                # input), otherwise the observed constant row (the last row is only ever a target).
                add_agreement(next_step_vars, _row_input_state(matrix, denoised_rows, t + 1), t)
            trajectories.append([_row_input_state(matrix, denoised_rows, t) for t in range(n_rows)])
        else:
            raise ValueError("Unrecognized model addition type {}".format(model_addition_type))
    if model_addition_type == ModelAdditionType.INDICATORS:
        if return_trajectories:
            return indicators, flip_cost_terms, trajectories
        return indicators, flip_cost_terms


def add_scaffold_network_agreement_expression(graph, model, function_type_restrictions=None, missing_edge_cost=1,
                                              added_edge_cost=2):
    # TODO: add representation of unknown input nodes
    raise NotImplementedError()


def add_unknown_input_general_function_variables(graph, model, max_degrees):
    """
    Adds, for each node i,j, the variable IS_INPUT{i,j,k} determining whether i is the k'th input to j.
    We require that each node appears at most once as an input, so that \forall i,j.\sum_{k} IS_INPUT{i,j,k}\leq 1.
    We require that each input is defined uniquely, and so \forall k,j.\sum_{i} IS_INPUT{i,j,k}= 1.
    Finally, assume max_degrees is an array-like sequence of bounds for each node's degree (e.g. original degree in
    the graph + 1, or a constant maximal degree for all nodes). We will require \forall j.\sum_{i,k} \leq max_degrees[j]
    :param graph:
    :param model:
    :param max_degrees:
    :return: an array of length len(graph). For each node j, a multidimensional array of IS_INPUT{i,j,k}.
    """
    return NotImplementedError()


def add_unknown_input_threshold_function_variables(graph, model):
    """
    Adds
    :param graph:
    :param model:
    :param function_type_restrictions:
    :param v_funcs:
    :return:
    """
    # TODO: add documentation
    return NotImplementedError()

def get_value_of_gurobi_entity(v):
    """
    Gets the value of the entity in the current solution. Gurobi has a different interfece
    for querying variable and expression values, so need to separate them in code.
    :param entity:
    :return:
    """
    if isinstance(v, gurobipy.Var):
        return v.x
    elif isinstance(v, gurobipy.LinExpr):
        return v.getValue()
    else:
        raise NotImplementedError("Unknown gurobi entity type {}".format(type(v)))
