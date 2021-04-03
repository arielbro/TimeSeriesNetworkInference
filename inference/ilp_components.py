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


def add_matrices_as_model_paths(graph, model, data_matrices, function_vars, model_to_data_sample_rate_ratio=1,
                                function_type_restrictions=None, model_addition_type=ModelAdditionType.CONSTRAINT):
    if model_to_data_sample_rate_ratio != 1:
        raise NotImplementedError()

    indicators = []
    for matrix_index, matrix in enumerate(data_matrices):
        for t in range(len(matrix) - 1):
            if model_addition_type == ModelAdditionType.CONSTRAINT:
                ilp.add_path_to_model(graph, model, path_len=1, first_state_vars=matrix[t], model_f_vars=None,
                                      last_state_vars=matrix[t+1], v_funcs_restrictions=function_type_restrictions)
            elif model_addition_type == ModelAdditionType.INDICATORS:
                next_step_vars = ilp.add_path_to_model(graph, model, path_len=1, first_state_vars=matrix[t],
                                      last_state_vars=None, v_funcs_restrictions=function_type_restrictions,
                                      model_f_vars=function_vars)[0]
                # note that it's a weak constraint indicator, i.e. indicator -> constraint
                indicator = ilp.add_state_equality_indicator(model, next_step_vars, matrix[t + 1], force_equal=False,
                                      prefix="add_matrices_matrix_index_{}_row_{}".format(matrix_index, t))
                indicators.append(indicator)
            else:
                raise ValueError("Unrecognized model addition type {}".format(model_addition_type))
    if model_addition_type == ModelAdditionType.INDICATORS:
        return indicators


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
