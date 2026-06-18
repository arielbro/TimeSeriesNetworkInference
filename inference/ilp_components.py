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


def add_matrices_as_model_paths(graph, model, data_matrices, function_vars, model_to_data_sample_rate_ratio=1,
                                function_type_restrictions=None, model_addition_type=ModelAdditionType.CONSTRAINT,
                                per_cell_indicators=False, allow_input_flips=False, no_anchoring=False):
    """
    Adds the data matrices to the model as model transitions.

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

    n = len(graph.vertices)
    indicators = []
    flip_cost_terms = []
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
                                      last_state_vars=matrix[t + 1], v_funcs_restrictions=function_type_restrictions,
                                      name_prefix="matrix_{}_time_{}".format(matrix_index, t))
        elif model_addition_type == ModelAdditionType.INDICATORS and no_anchoring:
            # single free-running trajectory: the model is fed its own predictions after the initial state
            predicted_states = ilp.add_path_to_model(graph, model, path_len=n_rows - 1,
                                  first_state_vars=_row_input_state(matrix, denoised_rows, 0), last_state_vars=None,
                                  v_funcs_restrictions=function_type_restrictions, model_f_vars=function_vars,
                                  name_prefix="matrix_{}_freerun".format(matrix_index))
            for t in range(1, n_rows):
                add_agreement(predicted_states[t - 1], matrix[t], t)
        elif model_addition_type == ModelAdditionType.INDICATORS:
            for t in range(n_rows - 1):
                next_step_vars = ilp.add_path_to_model(graph, model, path_len=1,
                                      first_state_vars=_row_input_state(matrix, denoised_rows, t),
                                      last_state_vars=None, v_funcs_restrictions=function_type_restrictions,
                                      model_f_vars=function_vars,
                                      name_prefix="matrix_{}_time_{}".format(matrix_index, t))[0]
                # the target is the denoised next row where it exists (i.e. the next row is itself a model
                # input), otherwise the observed constant row (the last row is only ever a target).
                add_agreement(next_step_vars, _row_input_state(matrix, denoised_rows, t + 1), t)
        else:
            raise ValueError("Unrecognized model addition type {}".format(model_addition_type))
    if model_addition_type == ModelAdditionType.INDICATORS:
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
