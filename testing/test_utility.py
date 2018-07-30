from unittest import TestCase
import utility
import random
import graphs
import stochastic
import random
import sympy


class TestUtility(TestCase):

    def test_attractor_sets_equality(self):
        for test in range(50):
            n = random.randint(1, 20)
            G = graphs.Network.generate_random(n_vertices=n, indegree_bounds=[1,5])
            first_attractors = list(
                stochastic.estimate_attractors(G, n_walks=20, max_walk_len=1000, with_basins=False))

            second_attractors = first_attractors
            self.assertTrue(utility.attractor_sets_equality(first_attractors, second_attractors))
            second_attractors = first_attractors.copy()
            random.shuffle(second_attractors)
            self.assertTrue(utility.attractor_sets_equality(first_attractors, second_attractors))
            second_attractors = first_attractors[:-1]
            self.assertFalse(utility.attractor_sets_equality(first_attractors, second_attractors))
            first_attractors = first_attractors[1:]
            self.assertFalse(utility.attractor_sets_equality(first_attractors, second_attractors))

    def test_is_same_attractor(self):
        for test in range(50):
            n = random.randint(1, 20)
            G = graphs.Network.generate_random(n_vertices=n, indegree_bounds=[1,5])
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
