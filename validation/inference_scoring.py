import os
from attractor_learning import graphs
import numpy as np
import shutil
import itertools


def models_to_edge_vectors(reference_model, inference_model):
    assert {v.name for v in reference_model.vertices} == {v.name for v in inference_model.vertices}
    y_true, y_pred = np.zeros(shape=(len(reference_model) ** 2, )), \
                     np.zeros(shape=(len(reference_model) ** 2, ))
    for i, (u, v) in enumerate(itertools.product(reference_model.vertices, repeat=2)):
        y_true[i] = 1 if (u, v) in reference_model.edges else 0
        # different model, work by name
        edge_equiv = inference_model.get_vertex(u.name), inference_model.get_vertex(v.name)
        y_pred[i] = 1 if edge_equiv in inference_model.edges else 0
    assert i == len(y_true) - 1
    return y_true, y_pred


def model_dirs_to_edge_vectors_list(reference_dir, inference_dir):
    reference_models = {f.name: graphs.Network.parse_cnet(os.path.join(f.path, "true_network.cnet")) for
                        f in os.scandir(reference_dir) if f.is_dir()}
    inference_models = {f.name: graphs.Network.parse_cnet(os.path.join(f.path, "inferred_network.cnet"))
                        for f in os.scandir(inference_dir) if f.is_dir()}
    assert (set(reference_models.keys()) == set(inference_models.keys()))

    ref_vectors = []
    pred_vectors = []
    for model_name, ref_model in reference_models.items():
        ref_model = reference_models[model_name]
        pred_model = inference_models[model_name]
        ref_vector, pred_vector = models_to_edge_vectors(ref_model, pred_model)
        ref_vectors.append(ref_vector)
        pred_vectors.append(pred_vector)
    return {'ref': ref_vectors, 'pred': pred_vectors}


def model_dirs_to_timeseries_vectors(reference_dir, inference_dir):
    reference_model_names = {f.name for f in os.scandir(reference_dir) if f.is_dir()}
    inference_model_names = {f.name for f in os.scandir(inference_dir) if f.is_dir()}
    assert(reference_model_names == inference_model_names)

    ref_train_vecs = []
    ref_test_vecs = []
    pred_train_vecs = []
    pred_test_vecs = []
    for model_name in reference_model_names:
        ref_matrices = np.load(os.path.join(reference_dir, model_name, "matrices.npz"))

        pred_train_matrices = np.load(os.path.join(inference_dir, model_name, "train_matrices.npz"))
        pred_test_matrices = np.load(os.path.join(inference_dir, model_name, "test_matrices.npz"))

        ref_train_matrices = {i: mat for i, mat in ref_matrices.items() if i in pred_train_matrices}
        ref_test_matrices = {i: mat for i, mat in ref_matrices.items() if i in pred_test_matrices}

        assert(set(ref_matrices.keys()) ==
               (set(pred_train_matrices.keys()) | set(pred_test_matrices.keys())))

        ref_train_vecs.append(np.concatenate([mat[1:, ].flatten() for mat in ref_train_matrices.values()]))
        ref_test_vecs.append(np.concatenate([mat[1:, ].flatten() for mat in ref_test_matrices.values()]))
        pred_train_vecs.append(np.concatenate([mat[1:, ].flatten() for mat in pred_train_matrices.values()]))
        pred_test_vecs.append(np.concatenate([mat[1:, ].flatten() for mat in pred_test_matrices.values()]))
    return {'ref_train': ref_train_vecs, 'ref_test': ref_test_vecs,
            'pred_train': pred_train_vecs, 'pred_test': pred_test_vecs}


def model_dirs_to_boolean_function_vectors(reference_dir, inference_dir):
    raise NotImplementedError()


def model_dirs_to_time_taken_vector(inference_dir):
    timings = []
    for model_name in {f.name for f in os.scandir(inference_dir) if f.is_dir()}:
        timing = float(np.load(os.path.join(inference_dir, model_name, "inference_time.npy")))
        timings.append(timing)
    return timings


def aggregate_classification_metric(reference_vectors, prediction_vectors, metric):
    assert(len(reference_vectors) == len(prediction_vectors))
    has_constant_vecs = False
    for vec in reference_vectors:
        if len(np.unique(vec)) == 1:
            print("Warning: reference vector has only one value, {}".format(np.unique(vec)[0]))
            has_constant_vecs = True
    if has_constant_vecs:
        res_metrics = []
        for ref_vec, pred_vec in zip(reference_vectors, prediction_vectors):
            if len(np.unique(ref_vec)) > 1:
                res_metrics.append(metric(ref_vec, pred_vec))
        print(len(res_metrics))
        return res_metrics
    return [metric(y_ref, y_pred) for (y_ref, y_pred) in zip(reference_vectors, prediction_vectors)]
