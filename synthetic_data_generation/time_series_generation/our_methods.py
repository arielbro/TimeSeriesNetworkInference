import enum
import random
from attractor_learning.stochastic import walk_to_attractor
import numpy as np

StateSampleType = enum.Enum("StateSampleType", "stable perturbed")
FrequencyHandling = enum.Enum("FrequencyHandling", "random floor")


def generate_one_experiment_data(model, **kwargs):
    """
    Given a model, generate one matrix of time-series data from the model. Data starts
    at some state (either a basin-weighted attractor state,
    or a perturbation of one, or a completely random state if not only_attractors),
    and advances according to the state space of the model.
    The method supports generating noisy data and sampling at different frequencies than model
    update.
    :param model:
    :param state_sample_type:
    :param timepoints_per_experiment:
    :param state_noise_chance:
    :param sample_to_model_freq_ratio:
    :param frequency_handling:
    :param frequency_noise_std:
    :param only_attractors:
    :return: an np.array matrix, corresponding to one experiment
    """

    n_nodes = len(model)
    data = np.zeros(shape=(kwargs['timepoints_per_experiment'], n_nodes))

    # sample starting state
    starting_point = [random.randint(0, 1) for _ in range(n_nodes)]
    if kwargs['only_attractors']:
        starting_attractor = walk_to_attractor(model, starting_point)
        # walk_to_attractor returns attractor states as tuples (next_state yields tuples); convert to a
        # list so the perturbation below can assign to an element.
        starting_point = list(random.choice(starting_attractor))

    if kwargs['state_sample_type'] == StateSampleType.perturbed:
        i = random.randint(0, n_nodes - 1)
        starting_point[i] = 1 - starting_point[i]
    elif kwargs['state_sample_type'] == StateSampleType.stable:
        pass
    else:
        raise ValueError("unkown StateSampleType {}".format(kwargs['state_sample_type']))

    if kwargs['frequency_noise_std'] != 0:
        raise NotImplementedError("frequency noise unimplemented")

    # fill in time-series
    model_states = [tuple(starting_point)]
    data[0, :] = model_states[0]
    q = kwargs['sample_to_model_freq_ratio']
    for t in range(1, kwargs['timepoints_per_experiment']):
        model_step = float(t) / q
        model_step_floor = int(model_step)
        for _ in range(model_step_floor + 1 - len(model_states)):
            model_states.append(model.next_state(model_states[-1]))
        state = np.zeros(shape=(n_nodes,))
        for i in range(n_nodes):
            if kwargs['frequency_handling'] == FrequencyHandling.floor:
                state[i] = model_states[model_step_floor][i]
            elif kwargs['frequency_handling'] == FrequencyHandling.random:
                if random.random() > model_step - model_step_floor:
                    state[i] = model_states[model_step_floor][i]
                else:
                    state[i] = model_states[model_step_floor + 1][i]
            else:
                raise ValueError("Unkown frequency handling mode {}".format(kwargs['frequency_handling']))
        data[t, :] = state

    noisy_data = data.copy()
    # add noise
    for t in range(kwargs['timepoints_per_experiment']):
        for i in range(n_nodes):
            if random.random() < kwargs['state_noise_chance']:
                noisy_data[t, i] = 1 - noisy_data[t, i]

    return data, noisy_data


def generate_experiments_data(model, **kwargs):
    for _ in range(kwargs['experiments_per_network']):
        yield generate_one_experiment_data(model, **kwargs)
