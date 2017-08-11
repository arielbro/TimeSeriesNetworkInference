import sympy
import itertools


def pre_randomized_boolean_symbolic_func(boolean_outputs, *variables):
    formula = sympy.true
    for b_output, terms in zip(boolean_outputs, itertools.product(*[[~var, var] for var in variables])):
        if b_output: # TODO: Karnaugh maps?
            formula = formula | sympy.And(*terms)
    return formula


def pre_randomized_boolean_func(boolean_outputs):
    input_variables = [sympy.symbols("x" + str(i)) for i in range(len(boolean_outputs))]
    symbolic_func = pre_randomized_boolean_symbolic_func(boolean_outputs, *input_variables)
    return lambda *boolean_inputs: \
        symbolic_func.subs({input_variables[i]: boolean_inputs[i] for i in range(len(boolean_inputs))})


# TODO: create a pretty print option (e.g. a callable class creating partials with a homemade str()
def variable_bound_boolean_func(function_output_vars, *variables):
    formula = sympy.false
    for output_var, terms in zip(function_output_vars, itertools.product(*[[~var, var] for var in variables])):
        formula = formula | (sympy.And(*terms) & output_var)
    return formula


def formula_length(formula):
    # defined as the number of (non-unique) atoms in the formula
    if formula.is_Atom:
        return 1
    else:
        return sum(formula_length(arg) for arg in formula.args)