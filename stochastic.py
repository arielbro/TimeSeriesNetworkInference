import random
import sympy


def rotate(l, n):
    return l[n:] + l[:n]


def is_same_attractor(a1, a2):
    """
    :param a1: an attractor, represented as an ordered iterable of network states
    :param a2: ""
    :return: True iff the attractors have same states in same order, up to a shift.
    """
    assert type(a1) == type(a2)
    if len(a1) != len(a2):
        return False
    for shift in range(len(a1)):
        if a1 == rotate(a2, shift):
            return True
    return False


def next_network_state(G, current_state):
    """
    :param G: A network.
    :param current_state: a binary vector of length |G.vertices|
    :return: the next state after activation of G's boolean functions on current_state
    """
    bool_res = []
    for v, cur_value in zip(G.vertices, current_state):
        if len(v.predecessors()) != 0:
            bool_res.append(v.function(*[bool(current_state[u.index]) for u in v.predecessors()]))
        else:
            bool_res.append(cur_value)
    for b in bool_res:
        assert isinstance(b, bool) or b is sympy.true or b is sympy.false
    binary_res = tuple(0 if not b else 1 for b in bool_res)
    return binary_res


def estimate_attractors(G, n_walks, max_walk_len=None):
    """
    An attractor estimation algorithm based on sampling a starting state repeatedly and simulating the
    network until an attractor is reached.
    :param G:  A network.
    :param n_walks: The number of times the algorithm samples a starting state to walk from.
    :param max_walk_len: The max number of states the algorithm is willing to visit until an attractor if found
    :return: A list [(states, estimated_basin)] of attractors and their basins.
    """
    n = len(G.vertices)
    if max_walk_len is None:
        max_walk_len = 2**n
    attractors = list()
    # TODO: find how far we are from finding all (using sum(len(attractor)) for those found)
    state_to_attractor_mapping = dict()  # used to stop walks early as soon as you know who's basin it is.
    for _ in range(n_walks):
        # start walking until you hit a cycle or the attractor length bound.
        step_number = 1
        visited_states = []  # TODO: optimize with hash table while still keeping an ordered list alongside.
        initial_state = tuple(random.randint(0, 1) for _ in range(n))
        visited_states.append(initial_state)
        current_state = next_network_state(G, initial_state)
        while current_state not in visited_states and current_state not in state_to_attractor_mapping and\
                step_number < max_walk_len:
            visited_states.append(current_state)
            current_state = next_network_state(G, current_state)
            step_number += 1

        # stepped into the basin of one of the attractors (including inside the attractor)
        if current_state in state_to_attractor_mapping:
            for state in visited_states:
                state_to_attractor_mapping[state] = state_to_attractor_mapping[current_state]
            continue

        # either a new attractor, or an exhausting walk.
        cycle_start_index = -1 if step_number == max_walk_len else\
            min(i for i in range(len(visited_states)) if visited_states[i] == current_state)  # can optimize
        if cycle_start_index == -1 and next_network_state(G, visited_states[-1]) != visited_states[0]:
            # not an attractor!
            continue
        attractor = tuple(visited_states[cycle_start_index:])
        for existing_attractor in attractors:
            if is_same_attractor(attractor, existing_attractor):
                raise AssertionError("found an existing attractor without recognizing the basin first.")
        attractors.append(attractor)
        for state in visited_states:
            state_to_attractor_mapping[state] = attractor

    attractor_to_basin = dict()
    for state, attractor in state_to_attractor_mapping.items():
        attractor_to_basin[attractor] = attractor_to_basin.get(attractor, 0) + 1
    # print 'finished estimation'
    return attractor_to_basin.items()
