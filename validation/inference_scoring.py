import os
from attractor_learning import graphs
import numpy as np
import shutil
import itertools


def models_to_edge_vectors(reference_model, inference_model):
    assert {v.name for v in reference_model.vertices} == {v.name for v in inference_model.vertices}
    y_true, y_pred = np.zeros(shape=(len(reference_model) ** 2, )), \
                     np.zeros(shape=(len(reference_model) ** 2, ))
    for i, edge in enumerate(itertools.product([v.name for v in reference_model.vertices], repeat=2)):
        y_true[i] = 1 if edge in reference_model.edges else 0
        y_pred[i] = 1 if edge in inference_model.edges else 0
    return y_true, y_pred


def model_dirs_to_edge_vectors_list(reference_dir, inference_dir):
    reference_models = {f.name: graphs.Network.parse_cnet(os.path.join(f.path, "true_network")) for
                        f in os.scandir(reference_dir) if f.is_dir()}
    inference_models = {f.name: graphs.Network.parse_cnet(os.path.join(f.path, "inferred_network"))
                        for f in os.scandir(inference_dir) if f.is_dir()}
    assert (set(reference_models.keys()) == set(inference_models.keys()))

    ref_vectors = []
    inf_vectors = []
    for model_name, ref_model in reference_models.items():
        ref_model = reference_models[model_name]
        inf_model = inference_models[model_name]
        ref_vector, inf_vector = models_to_edge_vectors(ref_model, inf_model)
        ref_vectors.append(ref_vector)
        inf_vectors.append(inf_vector)
    return ref_vectors, inf_vectors


def model_dirs_to_timeseries_vectors(reference_dir, inference_dir):
    reference_model_names = {f.name for f in os.scandir(reference_dir) if f.is_dir()}
    inference_model_names = {f.name for f in os.scandir(inference_dir) if f.is_dir()}
    assert(reference_model_names == inference_model_names)

    for model_name in reference_model_names:
        ref_matrices = np.load(os.path.join(reference_dir, model_name, "matrices.npz"))
        inf_matrices = np.load(os.path.join(inference_dir, model_name, "matrices.npz"))

        assert(set(ref_matrices.keys()) == set(inf_matrices.keys()))

        ref_vecs = [mat[1:, ].flatten() for mat in ref_matrices]
        inf_vecs = [mat[1:, ].flatten() for mat in inf_matrices]
    return ref_vecs, inf_vecs


def model_dirs_to_boolean_function_vectors(reference_dir, inference_dir):
    raise NotImplementedError()