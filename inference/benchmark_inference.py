"""Baseline inference methods to benchmark the ILP-based inference against.

All three assume the scaffold network's topology as gold truth (they only decide each node's Boolean
function, keeping the scaffold's edges), and share the standard inference interface
``method(data_matrices, scaffold_network, **kwargs) -> Network``. ``data_matrices`` is an iterable of
time-series matrices, each of shape (timepoints, n_nodes) with column i holding vertex-index-i's values
(the same column<->index alignment the ILP path and scoring rely on); a transition is a consecutive pair
of rows within one matrix.

Row/bit convention: a node's inputs are its scaffold predecessors in vertex-index order, and truth-table
row j (index into ``BooleanSymbolicFunc.boolean_outputs``) decodes with the first predecessor as the most
significant bit - matching how ``BooleanSymbolicFunc`` builds its formula and how ``next_state`` calls it.
"""

import itertools
import random
import warnings
import numpy as np
from attractor_learning.logic import BooleanSymbolicFunc

# Inverse L1-regularization strength for the logistic-regression baseline (sklearn's C; smaller = stronger
# lasso). The task doesn't specify a strength, so we use sklearn's default.
DEFAULT_LASSO_C = 1.0


def _row_index(input_values):
    """Truth-table row index for predecessor values given in vertex-index order (first predecessor = MSB)."""
    row = 0
    for value in input_values:
        row = (row << 1) | int(bool(value))
    return row


def _fit_l1_logistic(X, y, C):
    """Fit an L1-regularized (lasso) logistic regression, spanning sklearn's API change around the `penalty`
    parameter. `penalty='l1'` is valid and correct on every version below 1.10 (liblinear has supported L1
    for years); 1.10 removes the parameter, so constructing with it raises TypeError and we fall back to its
    replacement, l1_ratio=1.0 (both keep the fast liblinear solver). The 1.8-1.9 deprecation/inconsistency
    warnings about `penalty` are silenced - the L1 behavior on those versions is verified correct."""
    from sklearn.linear_model import LogisticRegression
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*penalty.*")
        try:
            classifier = LogisticRegression(penalty='l1', solver='liblinear', C=C)  # sklearn < 1.10
        except TypeError:
            classifier = LogisticRegression(l1_ratio=1.0, solver='liblinear', C=C)  # sklearn >= 1.10
        classifier.fit(X, y)
    return classifier


def random_model_inference(data_matrices, scaffold_network, **kwargs):
    """A random Boolean model over the scaffold: each node's truth-table row output is an i.i.d. fair coin.

    The data is ignored; this is the "no information used" baseline (and the base the exact-match baseline
    starts from).
    """
    inferred_model = scaffold_network.copy()
    inferred_model.randomize_functions()  # NONE restriction: fair-coin output per truth-table row
    return inferred_model


def exact_match_else_random_inference(data_matrices, scaffold_network, **kwargs):
    """Start from a random model, then for every node override each truth-table row that some training
    transition lands on with the consensus (majority) observed output for that row (ties broken by a fair
    coin). Rows never observed in the training data keep their random value.
    """
    inferred_model = random_model_inference(data_matrices, scaffold_network, **kwargs)
    data_matrices = [np.asarray(matrix) for matrix in data_matrices]
    for vertex in inferred_model.vertices:
        predecessors = vertex.predecessors()
        if len(predecessors) == 0:
            continue  # input node: function stays None (value carried over by next_state)
        pred_indices = [u.index for u in predecessors]
        true_counts = {}
        total_counts = {}
        for matrix in data_matrices:
            for t in range(matrix.shape[0] - 1):
                row = _row_index(matrix[t, pred_indices])
                total_counts[row] = total_counts.get(row, 0) + 1
                true_counts[row] = true_counts.get(row, 0) + int(bool(matrix[t + 1, vertex.index]))
        boolean_outputs = list(vertex.function.boolean_outputs)  # random base, overridden below where seen
        for row, total in total_counts.items():
            n_true = true_counts[row]
            n_false = total - n_true
            if n_true > n_false:
                boolean_outputs[row] = True
            elif n_false > n_true:
                boolean_outputs[row] = False
            else:
                boolean_outputs[row] = random.choice([False, True])
        vertex.function = BooleanSymbolicFunc(input_names=[u.name for u in predecessors],
                                              boolean_outputs=boolean_outputs)
    return inferred_model


def linear_classifier_inference(data_matrices, scaffold_network, lasso_C=DEFAULT_LASSO_C, **kwargs):
    """For each node, fit a logistic regression (cross-entropy loss, L1/lasso regularization) predicting the
    node's next value from its scaffold predecessors' current values, then realize the learned classifier as
    a full truth table (predicting all 2**degree input rows).

    Degenerate cases keep the classifier's spirit: a node whose training outputs are all one class becomes
    the corresponding constant function; a node with no training transitions falls back to a fair-coin table.
    """
    inferred_model = scaffold_network.copy()
    data_matrices = [np.asarray(matrix) for matrix in data_matrices]
    for vertex in inferred_model.vertices:
        predecessors = vertex.predecessors()
        degree = len(predecessors)
        if degree == 0:
            vertex.function = None  # input node
            continue
        pred_indices = [u.index for u in predecessors]
        feature_blocks = [matrix[:-1, pred_indices] for matrix in data_matrices if matrix.shape[0] >= 2]
        target_blocks = [matrix[1:, vertex.index] for matrix in data_matrices if matrix.shape[0] >= 2]

        rows = list(itertools.product([False, True], repeat=degree))  # same row order as boolean_outputs
        row_features = np.array([[int(bit) for bit in row] for row in rows], dtype=int)
        if feature_blocks:
            X = np.vstack(feature_blocks).astype(int)
            y = np.concatenate(target_blocks).astype(int)
        else:
            X, y = np.empty((0, degree), dtype=int), np.empty((0,), dtype=int)

        classes = np.unique(y)
        if len(classes) == 0:  # no training transitions for this node
            outputs = [random.choice([False, True]) for _ in rows]
        elif len(classes) == 1:  # only one output ever observed: constant function
            outputs = [bool(classes[0])] * len(rows)
        else:
            classifier = _fit_l1_logistic(X, y, lasso_C)
            outputs = [bool(prediction) for prediction in classifier.predict(row_features)]
        vertex.function = BooleanSymbolicFunc(input_names=[u.name for u in predecessors],
                                              boolean_outputs=outputs)
    return inferred_model
