from utility import list_repr
import gurobipy
import sympy


def recursive_logic_to_var(formula, model, formulas_to_variables):
    if formula in formulas_to_variables:
        return formulas_to_variables[formula]
    assert not formula.is_Atom
    arg_vars = [recursive_logic_to_var(arg, model, formulas_to_variables) for arg in formula.args]
    if formula.func == sympy.And:
        return recursive_logic_to_var(~sympy.Or(*[~arg for arg in formula.args]), model, formulas_to_variables)
    elif formula.func == sympy.Or:
        name = "Or({})".format(list_repr(arg_vars))[:255] # gurobi len limit, no need to be unique though
        new_var = model.addVar(vtype=gurobipy.GRB.BINARY, name=name)
        formulas_to_variables[formula] = new_var
        model.update()
        model.addConstr(new_var <= sum(arg_vars), name=("constr: " + name + "<= sum")[:255]) #
        for var in arg_vars:
            model.addConstr(new_var >= var, name=("constr: " + name + ">= " +
                                                  var.VarName)[:255] if isinstance(var, gurobipy.Var)
                                                  else str(var)[:255])
        model.update()
        return new_var
    elif formula.func == sympy.Not:
        return 1 - arg_vars[0]
    elif formula.func == sympy.Implies:
        return recursive_logic_to_var(~formula.args[0] | formula.args[1], model, formulas_to_variables)
    elif formula.func == sympy.Equivalent:
        # return recursive_logic_to_var(formula.args[0] & formula.args[1] | ~formula.args[0] & ~formula.args[1],
        #                                model, formulas_to_variables)
        # save variables with a direct definition
        name = "not_Eq({})".format(list_repr(arg_vars))[:255] # gurobi len limit, no need to be unique though
        new_var = model.addVar(vtype=gurobipy.GRB.BINARY, name=name)
        formulas_to_variables[formula] = new_var
        model.update()
        a, b = arg_vars[0], arg_vars[1]
        model.addConstr(new_var >= a - b, name=("constr: " + name + ">= a - b")[:255])
        model.addConstr(new_var >= b - a, name=("constr: " + name + ">= b - 1")[:255])
        model.addConstr(new_var <= 2 -a - b, name=("constr: " + name + "<= 2 - a - b")[:255])
        model.addConstr(new_var <= a + b, name=("constr: " + name + "<= a + b")[:255])
        model.update()
        return 1 - new_var
    else:
        raise AssertionError("unexpected formula func encountered: {}".format(formula.func))


def logic_to_ilp(formula):
    symbols = formula.free_symbols
    model = gurobipy.Model()
    formulas_to_variables = {symbol: model.addVar(vtype=gurobipy.GRB.BINARY, name=symbol.name)
                             for symbol in symbols}
    model.update()
    top_formula_var = recursive_logic_to_var(formula, model, formulas_to_variables)  # create ILP vars/constrs
    model.addConstr(top_formula_var == 1, name="top level formula var")  # require top formula to be true
    return model
