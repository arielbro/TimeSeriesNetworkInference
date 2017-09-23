from unittest import TestCase
from ilp import create_state_keys_comparison_var, print_model_constraints
from utility import slice_int
import gurobipy
import random
import math


class TestIlp(TestCase):
    def test_state_key_comparison(self):
        for iteration in range(1000):
            a = random.randint(1, 2**5)
            b = random.randint(1, 2**5)
            if random.random() < 0.1:
                b = a
            model = gurobipy.Model()
            slice_size = random.randint(3, 3)
            domain_bits = int(math.floor(math.log(max(a, b), 2))) + 1
            include_equality = random.choice([False, True])
            variable_keys_tup = []
            for num in [a, b]:
                keys = slice_int(num, slice_size, domain_bits)
                variable_keys = [model.addVar(vtype=gurobipy.GRB.INTEGER) for key in keys]
                for var_key, key in zip(variable_keys, keys):
                    model.addConstr(var_key == key)
                variable_keys_tup.append(variable_keys)
            comparison_var = create_state_keys_comparison_var(model, variable_keys_tup[0],
                                                              variable_keys_tup[1],
                                                              include_equality=include_equality)
            model.setObjective(comparison_var, sense=gurobipy.GRB.MAXIMIZE)
            model.params.LogToConsole = 0
            print_model_constraints(model)
            model.optimize()
            print "a={}, b={}, domain_bits={}, slice_size={}, include_eq={}, res={}".format(a, b,
                domain_bits, slice_size, include_equality, model.ObjVal)
            print "\n"

            self.assertTrue(model.Status == gurobipy.GRB.OPTIMAL)
            if include_equality:
                self.assertEqual(a >= b, model.ObjVal == 1)
            else:
                self.assertEqual(a > b, model.ObjVal == 1)
