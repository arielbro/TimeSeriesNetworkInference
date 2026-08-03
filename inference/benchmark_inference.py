"""Baseline / published inference methods to benchmark the ILP-based inference against.

All share the standard inference interface ``method(data_matrices, scaffold_network, **kwargs) -> Network``.
``data_matrices`` is an iterable of time-series matrices, each of shape (timepoints, n_nodes) with column i
holding vertex-index-i's values (the same column<->index alignment the ILP path and scoring rely on); a
transition is a consecutive pair of rows within one matrix.

Two families live here:
  * Scaffold-topology baselines (random_model, exact_match_else_random, linear_classifier): keep the
    scaffold's edges and only decide each node's Boolean function.
  * Topology-free published methods (REVEAL, Best-Fit): infer each node's regulators from scratch by
    searching regulator sets of size up to an in-degree bound; the scaffold is used only to derive that
    bound when max_indegree == -1 (its max in-degree). See reveal_inference / best_fit_inference.

Row/bit convention: a node's inputs are ordered by vertex index, and truth-table row j (index into
``BooleanSymbolicFunc.boolean_outputs``) decodes with the first input as the most significant bit - matching
how ``BooleanSymbolicFunc`` builds its formula and how ``next_state`` calls it.
"""

import itertools
import random
import warnings
import numpy as np
from attractor_learning.graphs import Network
from attractor_learning.logic import BooleanSymbolicFunc

# Inverse L1-regularization strength for the logistic-regression baseline (sklearn's C; smaller = stronger
# lasso). The task doesn't specify a strength, so we use sklearn's default.
DEFAULT_LASSO_C = 1.0

# REVEAL/Best-Fit are topology-free: on the noise-holding data, an input/source node just carries its value
# (next == current every step), which any search reads as a perfect self-loop identity. When this flag is on,
# a gene whose next value equals its current value in every observed transition (a source node, or a
# globally-constant gene) is instead emitted as an input node (no regulators). This is data-driven (it never
# looks at the true topology) and keeps these two methods from scoring a false-positive self-loop on every
# source node - which would otherwise single them out in the edge (Jaccard) comparison, since the
# scaffold-topology baselines get input nodes right for free. Turn off to reproduce the raw methods' behavior
# (input nodes become self-loops). A genuinely-regulated node that never changes in-sample can be mislabeled.
EMIT_STATIC_GENES_AS_INPUTS = True


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


# ---------------------------------------------------------------------------------------------------------
# Topology-free published methods: REVEAL and Best-Fit. Both search, per gene, over candidate regulator sets
# of size up to an in-degree bound k (max_indegree, or the scaffold's max in-degree when -1), score each set,
# and read off the chosen set's truth table by majority vote per observed input pattern.
# ---------------------------------------------------------------------------------------------------------

def _transition_table(data_matrices):
    """Stack every one-step transition into (X, Y) int8 arrays of shape (n_transitions, n_nodes): X = states
    at t, Y = states at t+1, binarised (nonzero -> 1). Returns (None, None) if there are no transitions."""
    xs, ys = [], []
    for matrix in data_matrices:
        matrix = np.asarray(matrix)
        if matrix.shape[0] < 2:
            continue
        xs.append((matrix[:-1] != 0).astype(np.int8))
        ys.append((matrix[1:] != 0).astype(np.int8))
    if not xs:
        return None, None
    return np.vstack(xs), np.vstack(ys)


def _pattern_codes(X, subset):
    """Integer code in [0, 2**|subset|) for each transition's values on `subset` (given ascending, first =
    MSB) - the truth-table row that transition lands on."""
    columns = list(subset)
    weights = (1 << np.arange(len(columns) - 1, -1, -1)).astype(np.int64)
    return X[:, columns].astype(np.int64) @ weights


def _pattern_counts(codes, y, n_patterns):
    """(total, ones) per input pattern: transitions landing on each pattern, and how many output 1."""
    total = np.bincount(codes, minlength=n_patterns)
    ones = np.bincount(codes, weights=y.astype(float), minlength=n_patterns)
    return total, ones


def _conditional_entropy_bits(total, ones):
    """H(Y | S) in bits from per-pattern (total, ones) counts, weighted by pattern frequency."""
    n = total.sum()
    if n == 0:
        return 0.0
    mask = total > 0
    tp = total[mask].astype(float)
    p1 = ones[mask] / tp
    with np.errstate(divide='ignore', invalid='ignore'):
        h = np.where((p1 > 0) & (p1 < 1), -p1 * np.log2(p1) - (1 - p1) * np.log2(1 - p1), 0.0)
    return float(np.sum((tp / n) * h))


def _select_best_fit(X, y, n, k):
    """Best-Fit Extension: the regulator set (size 1..k) whose best Boolean function misclassifies the fewest
    observed transitions (the sum of the minority count over input patterns), ties broken toward fewer
    regulators. Note that with no complexity penalty this tends to use the full budget k on noisy data (adding
    regulators can only lower training error) - that is the method's known behavior, which is why k matters.
    Returns the chosen subset (ascending tuple), or None if k < 1."""
    y = y.astype(float)
    best = None  # (error, size, subset)
    for size in range(1, k + 1):
        n_patterns = 1 << size
        for subset in itertools.combinations(range(n), size):
            total, ones = _pattern_counts(_pattern_codes(X, subset), y, n_patterns)
            error = float(np.minimum(ones, total - ones).sum())
            if best is None or (error, size) < (best[0], best[1]):
                best = (error, size, subset)
        if best is not None and best[0] == 0.0:
            break  # perfect fit at this (smallest) size; a larger set can't beat zero error
    return best[2] if best else None


def _select_reveal(X, y, n, k):
    """REVEAL with an MDL/BIC relaxation for noisy data (ours): the regulator set (size 1..k) minimizing the
    description length  N*H(Y|S) + 0.5 * 2**|S| * log2(N)  (residual coding cost + model cost), ties broken
    toward fewer regulators. On noise-free data H(Y|S)=0 is achievable and this reduces to REVEAL's original
    criterion (the smallest fully-determining set). The complexity term counters the monotonicity of mutual
    information (H(Y|S) never increases as inputs are added), so pure-noise regulators are not accepted -
    without it, maximizing MI would always grab the largest allowed set. Returns the subset, or None if k<1."""
    y = y.astype(float)
    n_transitions = len(y)
    log2n = np.log2(n_transitions) if n_transitions > 1 else 0.0
    best = None  # (description_length, size, subset)
    for size in range(1, k + 1):
        n_patterns = 1 << size
        model_cost = 0.5 * n_patterns * log2n
        for subset in itertools.combinations(range(n), size):
            total, ones = _pattern_counts(_pattern_codes(X, subset), y, n_patterns)
            description_length = n_transitions * _conditional_entropy_bits(total, ones) + model_cost
            if best is None or (description_length, size) < (best[0], best[1]):
                best = (description_length, size, subset)
    return best[2] if best else None


def _truth_table_from_subset(X, y, subset):
    """boolean_outputs (length 2**|subset|) for regulators `subset` (ascending, first = MSB): majority output
    per observed input pattern, with ties and unobserved patterns filled by a fair coin (the relaxation that
    makes the learned function total)."""
    n_patterns = 1 << len(subset)
    total, ones = _pattern_counts(_pattern_codes(X, subset), y, n_patterns)
    outputs = []
    for row in range(n_patterns):
        if total[row] == 0:
            outputs.append(random.choice([False, True]))
        else:
            n_true = ones[row]
            n_false = total[row] - n_true
            if n_true > n_false:
                outputs.append(True)
            elif n_false > n_true:
                outputs.append(False)
            else:
                outputs.append(random.choice([False, True]))
    return outputs


def _infer_by_subset_search(data_matrices, scaffold_network, max_indegree, select_regulators,
                            emit_static_as_input):
    """Shared driver for the topology-free methods (REVEAL, Best-Fit). For each gene it (optionally) emits a
    static gene as an input node, else searches regulator sets of size up to the in-degree bound k with
    select_regulators and builds the gene's BooleanSymbolicFunc from the chosen set. k = max_indegree, or the
    scaffold's max in-degree when max_indegree == -1; k < 1 (e.g. an edgeless scaffold with -1) yields an
    all-input-node model."""
    data_matrices = list(data_matrices)
    names = [v.name for v in scaffold_network.vertices]
    n = len(names)
    X, Y = _transition_table(data_matrices)
    k = max_indegree if max_indegree != -1 else scaffold_network.max_in_degree()
    k = min(k, n)  # can't select more regulators than there are genes

    edges, functions = [], [None] * n
    if X is not None and k >= 1:
        for i in range(n):
            y = Y[:, i]
            if emit_static_as_input and np.all(y == X[:, i]):
                continue  # only ever holds its value (source node / constant) -> input node
            subset = select_regulators(X, y, n, k)
            if not subset:
                continue
            subset = tuple(sorted(subset))  # ascending -> first predecessor is MSB (BooleanSymbolicFunc order)
            functions[i] = BooleanSymbolicFunc(input_names=[names[j] for j in subset],
                                               boolean_outputs=_truth_table_from_subset(X, y, subset))
            edges.extend((names[j], names[i]) for j in subset)
    return Network(vertex_names=names, edges=edges, vertex_functions=functions)


def reveal_inference(data_matrices, scaffold_network, max_indegree=-1,
                     emit_static_as_input=EMIT_STATIC_GENES_AS_INPUTS, **kwargs):
    """REVEAL (Liang, Fuhrman & Somogyi 1998), topology-free, with an MDL relaxation for noisy data (see
    _select_reveal). max_indegree bounds each gene's regulator-set size (-1 = the scaffold's max in-degree)."""
    return _infer_by_subset_search(data_matrices, scaffold_network, max_indegree, _select_reveal,
                                   emit_static_as_input)


def best_fit_inference(data_matrices, scaffold_network, max_indegree=-1,
                       emit_static_as_input=EMIT_STATIC_GENES_AS_INPUTS, **kwargs):
    """Best-Fit Extension (Lähdesmäki, Shmulevich & Yli-Harja 2003), topology-free, with unobserved truth-table
    rows filled by a fair coin (the relaxation to a total function). max_indegree bounds each gene's
    regulator-set size (-1 = the scaffold's max in-degree)."""
    return _infer_by_subset_search(data_matrices, scaffold_network, max_indegree, _select_best_fit,
                                   emit_static_as_input)
