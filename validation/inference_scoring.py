import os
import random
from attractor_learning import graphs
from attractor_learning.logic import BooleanSymbolicFunc
import numpy as np
from scipy import sparse
import shutil
from sklearn.metrics import accuracy_score
import itertools


def prune_ignored_inputs(model):
    """Remove, in place, each edge feeding an input that a node's learned truth table does not actually depend
    on, so an edge score reflects the model's effective topology rather than the scaffold it was handed.

    Only nodes whose function is a BooleanSymbolicFunc are touched - i.e. the truth-table-based inference
    methods (general, random_model, exact_match_else_random, linear_classifier); threshold-function nodes
    (symmetric / symmetric_topology) and input nodes are left alone. An input is "ignored" when flipping it
    never changes the output for any assignment of the node's other inputs (probed by calling the function, so
    this is independent of the truth-table bit convention); the truth table is rebuilt over the remaining
    (relevant) inputs. A node whose function depends on no input at all (a learned constant) keeps a single,
    randomly chosen one of its inputs and stays a (constant) non-input node emitting that same value, with its
    other incoming edges removed. Every rewrite preserves the node's next_state output. Returns the number of
    edges removed.
    """
    predecessors_by_index = {v.index: list(v.predecessors()) for v in model.vertices}
    edges_to_remove = set()   # (predecessor_index, vertex_index) pairs
    replacement_functions = {}  # vertex_index -> new function (BooleanSymbolicFunc or None)

    for v in model.vertices:
        func = v.function
        predecessors = predecessors_by_index[v.index]
        degree = len(predecessors)
        if degree == 0 or not isinstance(func, BooleanSymbolicFunc):
            continue
        combinations = list(itertools.product([False, True], repeat=degree))
        outputs = {combo: bool(func(*combo)) for combo in combinations}

        # input k is relevant iff flipping it changes the output for some assignment of the other inputs
        relevant = []
        for k in range(degree):
            for combo in combinations:
                if not combo[k] and outputs[combo] != outputs[combo[:k] + (True,) + combo[k + 1:]]:
                    relevant.append(k)
                    break
        if len(relevant) == degree:
            continue  # depends on every input; nothing to prune

        if not relevant:
            # depends on no input (learned constant): keep one randomly chosen input and stay a constant,
            # non-input node emitting that same value; drop the node's other incoming edges.
            kept = random.choice(predecessors)
            for pred in predecessors:
                if pred.index != kept.index:
                    edges_to_remove.add((pred.index, v.index))
            constant_value = outputs[combinations[0]]  # identical for every assignment
            replacement_functions[v.index] = BooleanSymbolicFunc(
                input_names=[kept.name], boolean_outputs=[constant_value, constant_value])
        else:
            for k in range(degree):
                if k not in relevant:
                    edges_to_remove.add((predecessors[k].index, v.index))
            # ignored inputs are truly irrelevant, so fixing them (to False) while enumerating the relevant
            # ones reproduces the function over its support. reduced_names is in vertex-index order (a subset
            # of predecessors), matching the post-prune predecessor order the rebuilt function is called with.
            reduced_names = [predecessors[k].name for k in relevant]
            reduced_outputs = []
            for reduced_combo in itertools.product([False, True], repeat=len(relevant)):
                full = [False] * degree
                for position, k in enumerate(relevant):
                    full[k] = reduced_combo[position]
                reduced_outputs.append(outputs[tuple(full)])
            replacement_functions[v.index] = BooleanSymbolicFunc(input_names=reduced_names,
                                                                 boolean_outputs=reduced_outputs)

    if edges_to_remove:
        model.edges = [(a, b) for (a, b) in model.edges if (a.index, b.index) not in edges_to_remove]
        for vertex_index, new_func in replacement_functions.items():
            model.vertices[vertex_index].function = new_func
        for v in model.vertices:  # adjacency changed; drop caches so predecessors()/successors() recompute
            v.precomputed_predecessors = None
            v.precomputed_successors = None
    return len(edges_to_remove)


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


def varying_columns(matrix):
    """Boolean mask over the columns of one time-series matrix, True where the node is not constant down
    the matrix's rows - i.e. neither all zeros nor all ones, the two constant cases a binary matrix has.

    A node that holds the same value at every scored timepoint says nothing about whether the model got
    the dynamics right: any model that happens to hold it there scores it perfectly. In a sparse regime
    those nodes are most of the matrix, so they dominate the full-matrix accuracy; masking to the varying
    ones scores only the positions where the prediction had something to get wrong.
    """
    arr = np.asarray(matrix)
    if arr.shape[0] == 0:
        return np.zeros(arr.shape[1], dtype=bool)
    return np.any(arr != arr[0], axis=0)


def timeseries_score_vectors(ref_matrices, pred_matrices, keys=None, varying_only=False):
    """Flattened (y_true, y_pred) over the predicted part of a group of time-series matrices, the pair
    every time-series accuracy in this project is computed from.

    Scoring starts at row 1: row 0 is the state the model was seeded with, not something it predicted.
    With varying_only, each matrix additionally keeps only its varying columns (see varying_columns),
    taken from the GROUND TRUTH matrix and over the same rows that are scored - so the varying score
    always covers a subset of the positions the full score covers, and a node counts as varying only if
    it varies where the model is being graded. The mask is per matrix, as one trajectory's constant node
    is another's varying one.

    ref_matrices and pred_matrices are anything mapping a key to a (timepoints x nodes) array - a dict or
    an open .npz. keys defaults to the predicted matrices' own keys, which is what identifies the group
    (train/test) a set of predictions was made for.
    """
    keys = list(pred_matrices.keys()) if keys is None else list(keys)
    refs, preds = [], []
    for key in keys:
        ref = np.asarray(ref_matrices[key])[1:, ]
        pred = np.asarray(pred_matrices[key])[1:, ]
        assert ref.shape == pred.shape, \
            "reference and prediction matrices differ in shape for key {}: {} vs {}".format(
                key, ref.shape, pred.shape)
        if varying_only:
            mask = varying_columns(ref)
            ref, pred = ref[:, mask], pred[:, mask]
        refs.append(ref.flatten())
        preds.append(pred.flatten())
    if not refs:
        return np.empty(0), np.empty(0)
    return np.concatenate(refs), np.concatenate(preds)


def timeseries_accuracy_score(ref_matrices, pred_matrices, keys=None, varying_only=False):
    """Accuracy of the predicted time series against the reference one, over the full matrices or over
    the varying nodes alone (see timeseries_score_vectors).

    NaN when nothing is left to score: with varying_only, a group whose every node is constant in every
    trajectory has no positions the metric is defined over. That is a real (if uninformative) outcome,
    not a failure, so it is kept apart from the 0.0 the analysis scores an unfinished run as.
    """
    y_true, y_pred = timeseries_score_vectors(ref_matrices, pred_matrices, keys=keys,
                                              varying_only=varying_only)
    if y_true.size == 0:
        return float('nan')
    return sparse_accuracy_score(y_true, y_pred)


def timeseries_varying_score_filename(name, group):
    """Name (no extension) of the file holding one varying-node accuracy score, beside the full-matrix
    score it accompanies: timeseries_<name>_accuracy_score_<group>.npy gains a _varying counterpart.
    name is the comparison ('real', 'reference', 'real_start'), group is 'train' or 'test'."""
    return "timeseries_{}_varying_accuracy_score_{}".format(name, group)


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
