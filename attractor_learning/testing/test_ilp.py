from unittest import TestCase
from ilp import create_state_keys_comparison_var, print_model_constraints
from utility import slice_int
import gurobipy
import random
import math


class TestIlp(TestCase):
    def test_state_key_comparison(self):
        for iteration in range(1000):
            model = gurobipy.Model()
            slice_size = random.randint(1, 29)
            M = 2**slice_size
            a = random.randint(1, (M - 1)**3)
            b = random.randint(1, (M - 1)**3)
            if random.random() < 0.1:
                b = a
            domain_bits = int(math.floor(math.log(max(a, b), 2))) + 1
            include_equality = random.choice([False, True])
            variable_keys_tup = []
            for num in [a, b]:
                keys = slice_int(num, slice_size, domain_bits)
                variable_keys = [model.addVar(vtype=gurobipy.GRB.INTEGER) for _ in keys]
                for key in keys:
                    self.assertTrue(key <= M)
                    self.assertTrue(key >= 0)
                for var_key, key in zip(variable_keys, keys):
                    model.addConstr(var_key == key)
                variable_keys_tup.append(variable_keys)
            comparison_var = create_state_keys_comparison_var(model=model, first_state_keys=variable_keys_tup[0],
                                                              second_state_keys=variable_keys_tup[1],
                                                              include_equality=include_equality,
                                                              upper_bound=M,
                                                              name_prefix="test")
            model.setObjective(comparison_var, sense=gurobipy.GRB.MAXIMIZE)
            model.params.LogToConsole = 0
            print_model_constraints(model)
            model.optimize()
            print "a={}, b={}, domain_bits={}, slice_size={}, include_eq={}".format(a, b,
                domain_bits, slice_size, include_equality)
            print "res={}".format(model.ObjVal)
            print "\n"

            self.assertTrue(model.Status == gurobipy.GRB.OPTIMAL)
            if include_equality:
                self.assertEqual(a >= b, model.ObjVal == 1)
            else:
                self.assertEqual(a > b, model.ObjVal == 1)
