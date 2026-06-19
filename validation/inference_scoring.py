import os
from attractor_learning import graphs
import numpy as np
from scipy import sparse
import shutil
from sklearn.metrics import accuracy_score
import itertools


def models_to_edge_vectors(reference_model, inference_model, use_sparse=True):
    assert {v.name for v in reference_model.vertices} == {v.name for v in inference_model.vertices}
    if use_sparse:
        y_true = sparse.lil_matrix((1, len(reference_model) ** 2))
        y_pred = sparse.lil_matrix((1, len(reference_model) ** 2))
        # use order in reference_model
        for vec, model in zip([y_true, y_pred], [reference_model, inference_model]):
            for edge in model.edges:
                u_index = reference_model.get_vertex(edge[0].name).index
                v_index = reference_model.get_vertex(edge[1].name).index
                index = u_index * len(reference_model) + v_index
                vec[0, index] = 1
    else:
        y_true, y_pred = np.zeros(shape=(len(reference_model) ** 2, )), \
                         np.zeros(shape=(len(reference_model) ** 2, ))
        for i, (u, v) in enumerate(itertools.product(reference_model.vertices, repeat=2)):
            y_true[i] = 1 if (u, v) in reference_model.edges else 0
            # different model, work by name
            edge_equiv = inference_model.get_vertex(u.name), inference_model.get_vertex(v.name)
            y_pred[i] = 1 if edge_equiv in inference_model.edges else 0
    return y_true, y_pred


def model_dirs_to_edge_vectors_list(reference_dir, inference_dir, use_sparse=True):
    model_names = {f.name for f in os.scandir(reference_dir) if f.is_dir()}
    assert model_names == {f.name for f in os.scandir(inference_dir) if f.is_dir()}

    for name in model_names:
        ref_model = graphs.Network.load(os.path.join(reference_dir, name, "true_network.json"))
        inferred_model = graphs.Network.load(os.path.join(inference_dir, name, "inferred_network.json"))
        ref_vector, pred_vector = models_to_edge_vectors(ref_model, inferred_model,
                                                         use_sparse=use_sparse)
        yield ref_vector, pred_vector


def model_dirs_to_network_sizes(reference_dir, inference_dir=None, use_sparse=True, with_edges=False):
    model_names = {f.name for f in os.scandir(reference_dir) if f.is_dir()}
    res = []

    for name in model_names:
        try:
            ref_model = graphs.Network.load(os.path.join(reference_dir, name, "true_network.json"))
        except FileNotFoundError as e:
            ref_model = graphs.Network.load(os.path.join(reference_dir, name, "inferred_network.json"))
        if with_edges:
            res.append(len(ref_model) + len(ref_model.edges))
        else:
            res.append(len(ref_model))
    return res


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
        assert((set(pred_train_matrices.keys()) & set(pred_test_matrices.keys())) == set())

        # iterate in a consistent way over train and test matrices
        train_keys = list(pred_train_matrices.keys())
        test_keys = list(pred_test_matrices.keys())
        ref_train_vecs.append(np.concatenate([ref_train_matrices[i][1:, ].flatten() for i in train_keys]))
        ref_test_vecs.append(np.concatenate([ref_test_matrices[i][1:, ].flatten() for i in test_keys]))
        pred_train_vecs.append(np.concatenate([pred_train_matrices[i][1:, ].flatten() for i in train_keys]))
        pred_test_vecs.append(np.concatenate([pred_test_matrices[i][1:, ].flatten() for i in test_keys]))
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


def aggregate_classification_metric(ref_inf_vector_iterator, metric):
    metric_fine_with_constant_vecs = True
    try:
        metric([1, 1], [1, 1])
    except ValueError as e:
        metric_fine_with_constant_vecs = False

    res_metrics = []
    i = 0
    for ref_vec, pred_vec in ref_inf_vector_iterator:
        if sparse.issparse(ref_vec):
            if ref_vec.getnnz() == 0:
                unique = [0]
            elif ref_vec.getnnz() == ref_vec.shape[1]:
                unique = np.unique(ref_vec.data[0])
            else:
                unique = [0] + list(np.unique(ref_vec.data[0]))
        else:
            unique = np.unique(ref_vec.data)
        if len(unique) == 1:
            print("Warning: reference vector has only one value, {}".format(unique[0]))
        if metric_fine_with_constant_vecs or (len(unique) > 1):
            if len(ref_vec.shape) == 2:
                # TODO: find way to get correct results from sklearn.metrics.accuracy_score with sparse one-row matrices
                res_metrics.append(metric(ref_vec.toarray()[0], pred_vec.toarray()[0]))
            else:
                res_metrics.append(metric(ref_vec, pred_vec))
        i += 1
        # if not i % 10:
        #     print(i)
    return res_metrics


def sparse_accuracy_score(x, y):
    """
    Returns the accuracy (hit rate) of equal length binary x and y, i.e. the fraction of
    positions where x_i == y_i. Agreement on zeros counts as well as agreement on ones.
    x and y can be dense (list, tuple, np.array) or sparse, in which case the computation
    exploits sparsity but still scores over the full length of the vectors.
    :param x:
    :param y:
    :return:
    """
    try:
        assert(len(x) == len(y))
    except TypeError:
        assert(x.shape == y.shape)

    if sparse.issparse(x) and sparse.issparse(y):
        n = float(max(x.shape))
        common = x.dot(y.T)
        assert(common.shape == (1, 1))
        n_common = common[0, 0]  # positions where both are 1 (binary vectors)
        nnz_x = x.getnnz()
        nnz_y = y.getnnz()
        # agreements = (both 1) + (both 0) = n_common + (n - nnz_x - nnz_y + n_common)
        n_agree = n - nnz_x - nnz_y + 2 * n_common
        return n_agree / n
    else:
        return accuracy_score(x, y)


def sparse_jaccard_score(x, y):
    """
    Returns the Jaccard index (intersection-over-union) of equal length binary x and y, i.e.
    |{i: x_i == y_i == 1}| / |{i: x_i == 1 or y_i == 1}|. True negatives (positions where both
    are 0) are ignored, making it suitable for sparse vectors such as adjacency matrices where
    accuracy would be dominated by the zeros.
    x and y can be dense (list, tuple, np.array) or sparse. An empty union (both all-zero)
    yields a score of 1.0.
    :param x:
    :param y:
    :return:
    """
    if sparse.issparse(x) and sparse.issparse(y):
        assert(x.shape == y.shape)
        intersection = x.dot(y.T)
        assert(intersection.shape == (1, 1))
        n_intersection = intersection[0, 0]  # positions where both are 1 (binary vectors)
        nnz_x = x.getnnz()
        nnz_y = y.getnnz()
    else:
        x = np.asarray(x)
        y = np.asarray(y)
        assert(x.shape == y.shape)
        n_intersection = int(np.sum((x != 0) & (y != 0)))
        nnz_x = int(np.count_nonzero(x))
        nnz_y = int(np.count_nonzero(y))

    n_union = nnz_x + nnz_y - n_intersection
    if n_union == 0:
        return 1.0
    return n_intersection / float(n_union)
