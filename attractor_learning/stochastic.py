import random
from . import utility
import numpy

def estimate_path_len_to_attractor(G, n_iter=1000):
    steps_nums = []
    for iteration in range(n_iter):
        if iteration and not iteration % 50:
            print("iteration #{}".format(iteration))
        step_number = 0
        initial_state = tuple(random.choice([False, True]) for _ in range(len(G.vertices)))
        visited_states_walk_count = {initial_state: 0}
        current_state = tuple(G.next_state(initial_state))
        while current_state not in visited_states_walk_count:
            step_number += 1
            visited_states_walk_count[current_state] = step_number
            current_state = G.next_state(current_state)
        steps_nums.append(visited_states_walk_count[current_state])
    print("mean number of steps to reach an attractor - {:.2f}. Max - {}.".format(numpy.mean(steps_nums),
                                                                                  max(steps_nums)))
    return steps_nums


def walk_to_attractor(G, initial_state, max_walk=None, state_to_attractor_mapping=None):
    """
    Simulates the network from an initial state until an attractor is reached, and returns the attractor.
    Optionally, a max_walk is supplied, and simulation terminates returning None if it is reached.
    If state_to_attractor_mapping is supplied, it is used to identify a basin's attractor without manually
    simulating a full loop on it. In this case, the mapping is also updated with the new states encountered.
    after which
    :param initial_state:
    :param max_walk:
    :param state_to_attractor_mapping:
    :return:
    """
    if state_to_attractor_mapping is None:
        state_to_attractor_mapping = dict()

    step_number = 0
    # states are used as dict keys here, and callers pass lists as well as tuples
    current_state = tuple(initial_state)
    # The list keeps the order, which the cycle needs; the index beside it answers "have I been here?" in
    # constant time. Scanning the list for that, as this did before, costs O(walk length ** 2) comparisons
    # of n-element tuples, which dominates the run time on a network whose transients are long - a chaotic
    # NK network, say, where a walk visits hundreds of states before it closes a cycle.
    visited_states = [current_state]
    visited_index = {current_state: 0}
    current_state = G.next_state(current_state)
    while current_state not in visited_index and current_state not in state_to_attractor_mapping and \
            (step_number < max_walk if max_walk is not None else True):
        visited_index[current_state] = len(visited_states)
        visited_states.append(current_state)
        current_state = G.next_state(current_state)
        step_number += 1

    if current_state in state_to_attractor_mapping:
        # stepped into the basin of a known attractor (including inside the attractor).
        attractor = state_to_attractor_mapping[current_state]
        for state in visited_states:
            state_to_attractor_mapping[state] = attractor
        return attractor
    elif (max_walk is not None) and step_number == max_walk:
        # An exhausting walk.
        return None
    elif current_state in visited_index:
        # A new attractor!
        attractor = tuple(visited_states[visited_index[current_state]:])
        # print("walked attractor: {}".format(attractor))
        # for state in attractor:
        #     print("source: {}. next: {}".format(state, G.next_state(state)))
        for state in visited_states:
            state_to_attractor_mapping[state] = attractor
        return attractor
    else:
        raise ValueError("Reached impossible case after network simulation")


def random_state(G):
    return tuple(random.randint(0, 1) for _ in range(len(G.vertices)))


def estimate_attractors(G, n_walks, max_walk_len=None, with_basins=True):
    """
    An attractor estimation algorithm based on sampling a starting state repeatedly and simulating the
    network until an attractor is reached.
    :param G:  A network.
    :param n_walks: The number of times the algorithm samples a starting state to walk from.
    :param max_walk_len: The max number of states the algorithm is willing to visit until an attractor if found
    :return: A list [(states, estimated_basin)] of attractors and their basins.
    """
    # TODO: write tests
    n = len(G.vertices)
    if max_walk_len is None:
        max_walk_len = 2**n
    # TODO: find how far we are from finding all (using sum(len(attractor)) for those found)
    # The shared mapping is what deduplicates: discovering a cycle records every one of its states, so a
    # later walk entering that cycle from any point returns the very same attractor object. A
    # rotation-invariant rescan of every attractor found so far used to run here too, costing
    # O(walks * attractors) comparisons, each trying every rotation of both cycles - and its result was
    # never read, since both returns below are built from the mapping instead.
    state_to_attractor_mapping = dict()  # used to stop walks early as soon as you know who's basin it is.
    for walk in range(n_walks):
        walk_to_attractor(G, random_state(G), max_walk=max_walk_len,
                          state_to_attractor_mapping=state_to_attractor_mapping)

    attractor_to_basin = dict()
    for state, attractor in state_to_attractor_mapping.items():
        attractor_to_basin[attractor] = attractor_to_basin.get(attractor, 0) + 1
    # print('finished estimation')
    if with_basins:
        return attractor_to_basin.items()
    return tuple(attractor_to_basin.keys())


def estimate_probability_to_reroll_attractor(G, n_walks, max_walk_len=None):
    """
    Using attractor estimation, estimate the probability of two uniform random choices of states to belong
    to same attractor. This is the sum of squares of probabilities to land in each attractor
    in the first place.
    :param G:
    :param n_walks:
    :param max_walk_len:
    :return:
    """
    attractor_to_basin_items = estimate_attractors(G, n_walks, max_walk_len, with_basins=True)
    total_basins = sum(basin_size for attractor, basin_size in attractor_to_basin_items)
    probabilities = [basin_size / float(total_basins) for attractor, basin_size in attractor_to_basin_items]
    return sum(p*p for p in probabilities)

