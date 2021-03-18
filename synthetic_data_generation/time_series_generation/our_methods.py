import enum
import random
from attractor_learning.stochastic import walk_to_attractor
import numpy as np

import logging

StateSampleType = enum.Enum("StateSampleType", "stable perturbed")
FrequencyHandling = enum.Enum("FrequencyHandling", "random floor")

timepoints_per_experiment = 10
state_sample_type = StateSampleType.stable
frequency_handling = FrequencyHandling.floor
sample_to_model_freq_ratio = 1.0
state_noise_chance = 0.0
frequency_noise_std = 0.0

have_logged = False

def generate_one_experiment_data(model, log=not have_logged):
    """
    Given a model, generate one matrix of time-series data from the model. Data starts
    at some state (either a basin-weighted attractor state, or a perturbation of one),
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
    :return: an np.array matrix, corresponding to one experiment
    """

    n_nodes = len(model)
    data = np.zeros(shape=(timepoints_per_experiment, n_nodes))

    # sample starting state
    starting_point = [random.randint(0, 1) for _ in range(n_nodes)]
    starting_attractor = walk_to_attractor(model, starting_point)

    starting_point = random.choice(starting_attractor)

    if state_sample_type == StateSampleType.perturbed:
        i = random.randint(0, n_nodes - 1)
        starting_point[i] = 1 - starting_point[i]
    elif state_sample_type == StateSampleType.stable:
        pass
    else:
        raise ValueError("unkown StateSampleType {}".format(state_sample_type))

    if frequency_noise_std != 0:
        raise NotImplementedError("frequency noise unimplemented")

    # fill in time-series
    model_states = [tuple(starting_point)]
    data[0, :] = model_states[0]
    q = sample_to_model_freq_ratio
    for t in range(1, timepoints_per_experiment):
        model_step = float(t) / q
        model_step_floor = int(model_step)
        for _ in range(model_step_floor + 1 - len(model_states)):
            model_states.append(model.next_state(model_states[-1]))
        state = np.zeros(shape=(n_nodes,))
        for i in range(n_nodes):
            if frequency_handling == FrequencyHandling.floor:
                state[i] = model_states[model_step_floor][i]
            elif frequency_handling == FrequencyHandling.random:
                if random.random() > model_step - model_step_floor:
                    state[i] = model_states[model_step_floor][i]
                else:
                    state[i] = model_states[model_step_floor + 1][i]
            else:
                raise ValueError("Unkown frequency handling mode {}".format(frequency_handling))
        data[t, :] = state

    # add noise
    for t in range(timepoints_per_experiment):
        for i in range(n_nodes):
            if random.random() < state_noise_chance:
                data[t, i] = 1 - data[t, i]

    return data


def generate_experiments_data(model, n_experiments):
    for _ in range(n_experiments):
        yield generate_one_experiment_data(model)


def log_params():
    logger = logging.getLogger(__name__)
    logger.info("timepoints_per_experiment={}".format(timepoints_per_experiment))
    logger.info("sample_to_model_freq_ratio={}".format(sample_to_model_freq_ratio))
    logger.info("state_noise_chance={}".format(state_noise_chance))
    logger.info("state_sample_type={}".format(state_sample_type))
    logger.info("frequency_handling={}".format(frequency_handling))