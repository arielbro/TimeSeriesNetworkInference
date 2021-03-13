from unittest import TestCase
import utility
import random
import graphs
import stochastic
import random
import sympy
import os
import attractors


class TestUtility(TestCase):

    def test_attractor_sets_equality(self):
        for test in range(50):
            n = random.randint(1, 6)
            G = graphs.Network.generate_random(n_vertices=n, indegree_bounds=[1, 5])
            first_attractors = list(
                stochastic.estimate_attractors(G, n_walks=500, max_walk_len=1000, with_basins=False))

            if len(first_attractors) == 0:
                continue

            second_attractors = first_attractors
            self.assertTrue(utility.attractor_sets_equality(first_attractors, second_attractors))
            second_attractors = list(s for s in first_attractors)
            random.shuffle(second_attractors)
            self.assertTrue(utility.attractor_sets_equality(first_attractors, second_attractors))

            # against dubrova
            second_attractors = attractors.find_attractors_dubrova(G, os.path.join("..", attractors.dubrova_path),
                                                                   mutate_input_nodes=True)
            self.assertTrue(utility.attractor_sets_equality(first_attractors, second_attractors))
            random.shuffle(second_attractors)
            self.assertTrue(utility.attractor_sets_equality(first_attractors, second_attractors))
            first_attractors = [utility.rotate(attractor, random.randint(0, 3)) for attractor in first_attractors]
            self.assertTrue(utility.attractor_sets_equality(first_attractors, second_attractors))
            first_attractors = [utility.rotate(attractor, random.randint(0, 3)) for attractor in first_attractors]
            self.assertFalse(utility.attractor_sets_equality(first_attractors, second_attractors[1:]))

            # no list conversion
            first_attractors = stochastic.estimate_attractors(G, n_walks=250, max_walk_len=1000, with_basins=False)
            self.assertTrue(utility.attractor_sets_equality(first_attractors, second_attractors))

            # from basins
            first_attractors_basin_pairs = list(
                stochastic.estimate_attractors(G, n_walks=250, max_walk_len=1000, with_basins=True))
            first_attractors = [p[0] for p in first_attractors_basin_pairs]
            self.assertTrue(utility.attractor_sets_equality(first_attractors, second_attractors))
            self.assertFalse(utility.attractor_sets_equality(first_attractors, second_attractors[:-1]))

            if len(first_attractors) == 1:
                continue

            second_attractors = second_attractors[:-1]
            self.assertFalse(utility.attractor_sets_equality(first_attractors, second_attractors))
            second_attractors = first_attractors[:-1]
            self.assertFalse(utility.attractor_sets_equality(first_attractors, second_attractors))
            first_attractors = first_attractors[1:]
            self.assertFalse(utility.attractor_sets_equality(first_attractors, second_attractors))

    def test_is_same_attractor(self):
        for test in range(50):
            n = random.randint(1, 20)
            G = graphs.Network.generate_random(n_vertices=n, indegree_bounds=[1, 5])
            attractor = stochastic.walk_to_attractor(G, stochastic.random_state(G))
            random_index = random.randint(0, len(attractor) - 1)
            shifted = utility.rotate(attractor, random_index)
            self.assertTrue(utility.is_same_attractor(attractor, shifted))

            perturbation_index = random.randint(0, n - 1)
            unaltered_state = attractor[random_index]
            perturbed_state = unaltered_state[:perturbation_index] + \
                              (1 - unaltered_state[perturbation_index],) + unaltered_state[perturbation_index + 1:]
            perturbed = attractor[:random_index] + (perturbed_state, ) + attractor[random_index + 1:]
            self.assertFalse(utility.is_same_attractor(attractor, perturbed))
            perturbed = utility.rotate(perturbed, random_index)
            self.assertFalse(utility.is_same_attractor(attractor, perturbed))

            format_changed = []
            for state in attractor:
                formatted_state = [random.choice([1, True, sympy.true]) if v_state else
                                   random.choice([0, False, sympy.false]) for v_state in state]
                if random.randint(0, 1):
                    formatted_state = tuple(formatted_state)
                format_changed.append(formatted_state)
            self.assertTrue(utility.is_same_attractor(attractor, format_changed))

    def test_is_same_state(self):
        for test in range(50):
            n = random.randint(1, 15)
            state = tuple(random.randint(0, 1) for _ in range(n))

            formatted_state = [random.choice([1, True, sympy.true]) if v_state else
                               random.choice([0, False, sympy.false]) for v_state in state]
            if random.randint(0, 1):
                formatted_state = tuple(formatted_state)
            self.assertTrue(utility.is_same_state(state, formatted_state))
            perturbation_index = random.randint(0, n - 1)
            perturbed_state = state[:perturbation_index
                              ] + (1 - state[perturbation_index],) + state[perturbation_index + 1:]
            self.assertFalse(utility.is_same_state(state, perturbed_state))

    def test_is_attractor_in_attractor_list(self):
        for test in range(50):
            n = random.randint(1, 6)
            G = graphs.Network.generate_random(n_vertices=n, indegree_bounds=[1, 5])
            first_attractors = list(
                stochastic.estimate_attractors(G, n_walks=500, max_walk_len=1000, with_basins=False))

            if len(first_attractors) == 0:
                continue

            second_attractors = first_attractors
            for attractor in first_attractors:
                self.assertTrue(utility.is_attractor_in_attractor_list(attractor, second_attractors))
            for attractor in second_attractors:
                self.assertTrue(utility.is_attractor_in_attractor_list(attractor, first_attractors))
            self.assertFalse(utility.is_attractor_in_attractor_list(first_attractors[0], second_attractors[1:]))

            second_attractors = list(s for s in first_attractors)
            random.shuffle(second_attractors)
            for attractor in first_attractors:
                self.assertTrue(utility.is_attractor_in_attractor_list(attractor, second_attractors))
            for attractor in second_attractors:
                self.assertTrue(utility.is_attractor_in_attractor_list(attractor, first_attractors))
            self.assertFalse(all(utility.is_attractor_in_attractor_list(a, first_attractors[1:])
                                 for a in second_attractors))

            # against dubrova
            second_attractors = attractors.find_attractors_dubrova(G, os.path.join("..", attractors.dubrova_path),
                                                                   mutate_input_nodes=True)
            for attractor in first_attractors:
                self.assertTrue(utility.is_attractor_in_attractor_list(attractor, second_attractors))
            for attractor in second_attractors:
                self.assertTrue(utility.is_attractor_in_attractor_list(attractor, first_attractors))
            self.assertFalse(all(utility.is_attractor_in_attractor_list(a, second_attractors[1:])
                                 for a in first_attractors))
            random.shuffle(second_attractors)
            self.assertFalse(all(utility.is_attractor_in_attractor_list(a, second_attractors[1:])
                                 for a in first_attractors))
            for attractor in first_attractors:
                self.assertTrue(utility.is_attractor_in_attractor_list(attractor, second_attractors))
            for attractor in second_attractors:
                self.assertTrue(utility.is_attractor_in_attractor_list(attractor, first_attractors))
            first_attractors = [utility.rotate(attractor, random.randint(0, 3)) for attractor in first_attractors]
            for attractor in first_attractors:
                self.assertTrue(utility.is_attractor_in_attractor_list(attractor, second_attractors))
            for attractor in second_attractors:
                self.assertTrue(utility.is_attractor_in_attractor_list(attractor, first_attractors))

            # no list conversion
            first_attractors = stochastic.estimate_attractors(G, n_walks=500, max_walk_len=1000, with_basins=False)
            for attractor in first_attractors:
                self.assertTrue(utility.is_attractor_in_attractor_list(attractor, second_attractors))
            for attractor in second_attractors:
                self.assertTrue(utility.is_attractor_in_attractor_list(attractor, first_attractors))
            self.assertFalse(all(utility.is_attractor_in_attractor_list(a, second_attractors[1:])
                                 for a in first_attractors))

            # from basins
            first_attractors_basin_pairs = list(
                stochastic.estimate_attractors(G, n_walks=250, max_walk_len=1000, with_basins=True))
            first_attractors = [p[0] for p in first_attractors_basin_pairs]
            for attractor in first_attractors:
                self.assertTrue(utility.is_attractor_in_attractor_list(attractor, second_attractors))
            for attractor in second_attractors:
                self.assertTrue(utility.is_attractor_in_attractor_list(attractor, first_attractors))
            self.assertFalse(all(utility.is_attractor_in_attractor_list(a, second_attractors[1:])
                                 for a in first_attractors))