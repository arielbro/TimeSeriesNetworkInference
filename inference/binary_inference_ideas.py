import attractor_learning.graphs
import attractor_learning.ilp
import gurobipy
from inference import ilp_components
from attractor_learning.graphs import FunctionTypeRestriction
from inference.ilp_components import ModelAdditionType, get_value_of_gurobi_entity
from attractor_learning.logic import BooleanSymbolicFunc, SymmetricThresholdFunction


# TODO: rethink the scattered way I do config, that made it worthwhile to have these
# TODO: "overload" functions for easy logging and usage.
def infer_known_topology_symmetric(*args, **kwargs):
    return infer_known_topology(*args, **kwargs,
                                function_type_restriction=FunctionTypeRestriction.SYMMETRIC_THRESHOLD)


def infer_known_topology_general(*args, **kwargs):
    return infer_known_topology(*args, **kwargs,
                                function_type_restriction=FunctionTypeRestriction.NONE)


def infer_known_topology(data_matrices, scaffold_network, function_type_restriction=None):
    """
    Find a model with best fit to data_matrices, assuming that each node's inputs are defined by the scaffold network
    topology.
    :param function_type_restriction:
    :param data_matrices:
    :param scaffold_network:
    :return:
    """
    # create function variables
    functions_variables = []
    model = gurobipy.Model()
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
            threshold = model.addVar(lb=-degree, ub=degree + 1, vtype=gurobipy.GRB.INTEGER,
                                     name="vertex_{}_threshold_var".format(i))
            functions_variables.append([signs, threshold])
        else:
            raise NotImplementedError()

    matrix_agreement_indicators = ilp_components.add_matrices_as_model_paths(
       scaffold_network, model, data_matrices,
       function_vars=functions_variables,
       model_to_data_sample_rate_ratio=1,
       function_type_restrictions=[function_type_restriction] * len(scaffold_network),
       model_addition_type=ModelAdditionType.INDICATORS)
    agreement = sum(matrix_agreement_indicators)

    model.setObjective(agreement, sense=gurobipy.GRB.MAXIMIZE)
    model.optimize()

    # set Boolean model found
    inferred_model = scaffold_network.copy()
    for i in range(len(inferred_model)):
        if len(inferred_model.vertices[i].predecessors()) == 0:
            inferred_model.vertices[i].function = None
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
    return inferred_model


