import attractor_learning.graphs
import attractor_learning.ilp
import gurobipy
from inference import ilp_components
from attractor_learning.graphs import FunctionTypeRestriction, Network
from inference.ilp_components import ModelAdditionType, get_value_of_gurobi_entity
from attractor_learning.logic import BooleanSymbolicFunc, SymmetricThresholdFunction
import itertools


# TODO: rethink the scattered way I do config, that made it worthwhile to have these
# TODO: "overload" functions for easy logging and usage.
def infer_known_topology_symmetric(*args, **kwargs):
    return infer_known_topology(*args, **kwargs,
                                function_type_restriction=FunctionTypeRestriction.SYMMETRIC_THRESHOLD)


def infer_known_topology_general(*args, **kwargs):
    return infer_known_topology(*args, **kwargs,
                                function_type_restriction=FunctionTypeRestriction.NONE)


def infer_known_topology(data_matrices, scaffold_network, function_type_restriction=None,
                         timeout_secs=None, log_file=None, **kwargs):
    """
    Find a model with best fit to data_matrices, assuming that each node's inputs are defined by the scaffold network
    topology.
    :param function_type_restriction:
    :param data_matrices:
    :param scaffold_network:
    :param timeout_secs: timelimit to pass to the solver.
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
            threshold = model.addVar(lb=0, ub=degree, vtype=gurobipy.GRB.INTEGER,
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

    if log_file is not None:
        model.Params.LogFile = log_file
    if timeout_secs is not None:
        model.Params.TimeLimit = timeout_secs
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


def infer_unknown_topology_symmetric(data_matrices, scaffold_network, allow_additional_edges=False,
                                   included_edges_relative_weight=1, added_edges_relative_weight=-1,
                                   timeout_secs=None, log_file=None, **kwargs):
    """
    Find a symmetric threshold model with best fit to data_matrices and scaffold_network,
    by finding both the Boolean function and the incoming edges for each node.
    The objective is 1 * x + included_edges_relative_weight * y + added_edges_relative_weight * z, where
    x is the proportion of cells in data_matrices (up to first rows) that are explained correctly by the model
    y is the proportion of scaffold_network edges that are included in the model.
    z is the number of edges not in scaffold_network that are included in the model, divided by
    the number of edges in scaffold_network.
    :param allow_additional_edges: if false, doesn't allow any edges that weren't in the scaffold network originally.
    :param included_edges_relative_weight:
    :param added_edges_relative_weight:
    :param function_type_restriction:
    :param data_matrices:
    :param scaffold_network:
    :param timeout_secs: timelimit to pass to the solver.
    :return:
    """
    # create function variables
    functions_variables = []
    model = gurobipy.Model()
    for i in range(len(scaffold_network)):
        # use ternary signs, zero means the input isn't used
        signs = [model.addVar(lb=-1, ub=1, vtype=gurobipy.GRB.INTEGER,
                              name="vertex_{}_input_{}_sign_var".format(i, j))
                              for j in range(len(scaffold_network))]
        threshold = model.addVar(lb=0, ub=len(scaffold_network), vtype=gurobipy.GRB.INTEGER,
                                 name="vertex_{}_threshold_var".format(i))
        functions_variables.append([signs, threshold])
    model.update()

    full_network = Network(vertex_names=[v.name for v in scaffold_network.vertices],
                           edges=[(u.name, v.name) for u, v in itertools.product(scaffold_network.vertices, repeat=2)],
                           vertex_functions=[None for v in scaffold_network.vertices])
    matrix_agreement_indicators = ilp_components.add_matrices_as_model_paths(
       full_network, model, data_matrices,
       function_vars=functions_variables,
       model_to_data_sample_rate_ratio=1,
       function_type_restrictions=[FunctionTypeRestriction.SYMMETRIC_THRESHOLD] * len(scaffold_network),
       model_addition_type=ModelAdditionType.INDICATORS)
    data_agreement = sum(matrix_agreement_indicators) / float(len(matrix_agreement_indicators))

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
    # that means that with actual degree d, we limit the threshold to range [1, d]
    for j in range(len(scaffold_network)):
        degree = sum(edge_indicators[i, j] for i in range(len(scaffold_network)))
        threshold = functions_variables[j][1]
        model.addConstr(threshold <= degree, name="node_{}_threshold_constraint_<=".format(j))
        # can be moved to lb in variable creation, but more readable here and gurobi surely can presolve this away.
        model.addConstr(threshold >= 1, name="node_{}_threshold_constraint_>=".format(j))

    included_edges_agreement = included_edges_relative_weight * sum(included_edges_indicators) / float(
                                          len(included_edges_indicators))
    added_edges_agreement = added_edges_relative_weight * sum(added_edges_indicators) / float(
                                          len(included_edges_indicators))

    if log_file is not None:
        model.Params.LogFile = log_file
    if timeout_secs is not None:
        model.Params.TimeLimit = timeout_secs
    model.setObjective(data_agreement + included_edges_agreement + added_edges_agreement, sense=gurobipy.GRB.MAXIMIZE)
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
        threshold = get_value_of_gurobi_entity(threshold)
        # assert threshold doesn't imply a constant function (that shouldn't be possible with the way we modelled this)
        assert threshold <= sum(abs(s) for s in signs)

        for j in range(len(inferred_model.vertices)):
            if signs[j] != 0:
                inferred_model.edges.append((inferred_model.vertices[j], inferred_model.vertices[i]))
        signs = [s for s in signs if s != 0]
        # threshold = max(-len(signs), min(threshold, len(signs) + 1))
        func = SymmetricThresholdFunction(signs, threshold)
        inferred_model.vertices[i].function = func
    return inferred_model
