"""Baseline / published inference methods to benchmark the ILP-based inference against.

All share the standard inference interface ``method(data_matrices, scaffold_network, **kwargs) -> Network``.
``data_matrices`` is an iterable of time-series matrices, each of shape (timepoints, n_nodes) with column i
holding vertex-index-i's values (the same column<->index alignment the ILP path and scoring rely on); a
transition is a consecutive pair of rows within one matrix.

Two families live here:
  * Scaffold-topology baselines (random_model, exact_match_else_random, linear_classifier): keep the
    scaffold's edges and only decide each node's Boolean function.
  * Published topology-inferring methods (REVEAL, Best-Fit): infer each node's regulators by searching
    regulator sets of size up to an in-degree bound. With allow_additional_edges on (the published,
    topology-free behaviour) the scaffold is used only to derive that bound when max_indegree == -1 (its
    max in-degree); with it off, the search is confined to each node's scaffold predecessors, so these
    methods can prune scaffold edges but not invent new ones - the same restriction the ILP methods obey.
    See reveal_inference / best_fit_inference.

Row/bit convention: a node's inputs are ordered by vertex index, and truth-table row j (index into
``BooleanSymbolicFunc.boolean_outputs``) decodes with the first input as the most significant bit - matching
how ``BooleanSymbolicFunc`` builds its formula and how ``next_state`` calls it.
"""

import itertools
import random
import time
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


def all_constants_inference(data_matrices, scaffold_network, **kwargs):
    """The trivial "nothing ever changes" baseline: an edgeless model, so every node is an input node and
    holds its value (Network.next_state), and a predicted trajectory is the initial state repeated.

    Both the data and the scaffold's edges are ignored - only its vertex names are kept. This is the score
    a zero-edge inference result already gets, made explicit as a benchmark: on data sampled at or near an
    attractor most cells never change, so this baseline alone can reach a high time-series accuracy, and a
    method is only informative to the extent it beats it. Its edge set is empty, so it scores zero edge
    overlap by construction.
    """
    return Network(vertex_names=[v.name for v in scaffold_network.vertices], edges=[],
                   vertex_functions=[None for _ in scaffold_network.vertices])


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
# Topology-free published methods: REVEAL and Best-Fit. Both search over candidate regulator sets of size up
# to an in-degree bound k (max_indegree, or the scaffold's max in-degree when -1), score each set, and read
# off the chosen set's truth table by majority vote per observed input pattern. The search runs breadth-first
# over set size - every gene at size 1, then every gene at size 2 - so a run that spends its budget leaves
# every gene searched to the same depth (see _infer_by_subset_search).
#
# Two knobs bound that search, both taken from the standard inference kwargs:
#   * allow_additional_edges - when False the candidates for a gene are only its scaffold predecessors, so
#     the methods pick the best subset of the given topology rather than ranging over every gene;
#   * timeout_secs - a wall-clock budget for the whole call, after which the search stops enlarging
#     regulator sets and keeps what it has (see _infer_by_subset_search).
# ---------------------------------------------------------------------------------------------------------

# Subsets scanned between clock reads inside a size sweep. The sweep body is a few microseconds, so polling
# the clock every iteration would be a measurable share of the runtime; every 1024 keeps the overshoot past
# the deadline negligible without that cost.
_DEADLINE_CHECK_EVERY = 1024


def _out_of_time(deadline):
    """True once the wall-clock budget is spent. deadline is an absolute time.time() value, or None for no
    limit (in which case the search always runs to completion)."""
    return deadline is not None and time.time() >= deadline

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


def _scan_best_fit(X, y, cands, size, best, deadline, self_column=None):
    """Best-Fit Extension over one gene's regulator sets of exactly `size`: the set whose best Boolean
    function misclassifies the fewest observed transitions (the sum of the minority count over input
    patterns). Returns (best, done, timed_out), where best is the better of the incoming (error, size,
    subset) and anything found here - ties broken toward fewer regulators - and done says the gene needs no
    larger set, since a perfect fit cannot be beaten. Note that with no complexity penalty this tends to use
    the full size budget on noisy data; that is the method's known behavior, which is why the cap matters."""
    y = y.astype(float)
    if size == 0:
        # the empty set is the gene as an input node: it holds its value, so its error is the number of
        # transitions where the next value differs from the current one
        error = float(np.count_nonzero(y != self_column.astype(float)))
        if best is None or (error, size) < (best[0], best[1]):
            best = (error, size, ())
        return best, best[0] == 0.0, False
    n_patterns = 1 << size
    for scanned, subset in enumerate(itertools.combinations(cands, size)):
        if size > 1 and scanned % _DEADLINE_CHECK_EVERY == 0 and _out_of_time(deadline):
            return best, False, True
        total, ones = _pattern_counts(_pattern_codes(X, subset), y, n_patterns)
        error = float(np.minimum(ones, total - ones).sum())
        if best is None or (error, size) < (best[0], best[1]):
            best = (error, size, subset)
    return best, best is not None and best[0] == 0.0, False


def _scan_reveal(X, y, cands, size, best, deadline, self_column=None):
    """REVEAL with an MDL/BIC relaxation for noisy data (ours) over one gene's regulator sets of exactly
    `size`: the set minimizing the description length  N*H(Y|S) + 0.5 * 2**|S| * log2(N)  (residual coding
    cost + model cost), ties broken toward fewer regulators. On noise-free data H(Y|S)=0 is achievable and
    this reduces to REVEAL's original criterion (the smallest fully-determining set). The complexity term
    counters the monotonicity of mutual information (H(Y|S) never increases as inputs are added), so
    pure-noise regulators are not accepted - without it, maximizing MI would always grab the largest allowed
    set. Returns (best, done, timed_out) as _scan_best_fit does; REVEAL has no early exit, so done is False."""
    y = y.astype(float)
    n_transitions = len(y)
    log2n = np.log2(n_transitions) if n_transitions > 1 else 0.0
    if size == 0:
        # the empty set is the gene as an input node: its prediction is its own current value, so the
        # residual is H(Y | X_i) - the same residual as the identity self-loop {i}, but with no truth table
        # to encode and so no model cost, which is what makes an input node preferred over that self-loop
        total, ones = _pattern_counts(self_column.astype(np.int64), y, 2)
        description_length = n_transitions * _conditional_entropy_bits(total, ones)
        if best is None or (description_length, size) < (best[0], best[1]):
            best = (description_length, size, ())
        return best, False, False
    n_patterns = 1 << size
    model_cost = 0.5 * n_patterns * log2n
    for scanned, subset in enumerate(itertools.combinations(cands, size)):
        if size > 1 and scanned % _DEADLINE_CHECK_EVERY == 0 and _out_of_time(deadline):
            return best, False, True
        total, ones = _pattern_counts(_pattern_codes(X, subset), y, n_patterns)
        description_length = n_transitions * _conditional_entropy_bits(total, ones) + model_cost
        if best is None or (description_length, size) < (best[0], best[1]):
            best = (description_length, size, subset)
    return best, False, False


def _select_one_gene(scan_subsets, X, y, n, k, candidates=None, deadline=None, self_column=None):
    """One gene's search over sizes 1..k, returning the chosen subset (ascending tuple) or None if the
    effective cap is below 1. candidates restricts which gene indices may be chosen (default: all n); the
    size cap is lowered to len(candidates) when that is smaller than k. deadline is an absolute wall-clock
    time past which the search stops enlarging sets and returns the best found so far - size 1 is always
    swept in full, so the result is never None for a non-empty candidate list.

    _infer_by_subset_search does not use this: it drives the same scans itself so that it can sweep one size
    across every gene before moving to the next size. Kept for a caller that wants a single gene's answer."""
    cands = list(range(n)) if candidates is None else list(candidates)
    best = None
    # size 0 - the gene as an input node - is only an option when the caller supplies its own column to
    # score that against; without it the sweep starts at one regulator, as the raw methods do
    first_size = 0 if self_column is not None else 1
    for size in range(first_size, min(k, len(cands)) + 1):
        if size > 1 and _out_of_time(deadline):
            break  # keep the best set found so far rather than starting a larger, slower sweep
        best, done, timed_out = scan_subsets(X, y, cands, size, best, deadline, self_column)
        if timed_out or done:
            break
    return best[2] if best is not None else None


def _select_best_fit(X, y, n, k, candidates=None, deadline=None, self_column=None):
    """One gene's Best-Fit search over sizes 1..k, or 0..k when self_column is given (see _scan_best_fit,
    _select_one_gene)."""
    return _select_one_gene(_scan_best_fit, X, y, n, k, candidates, deadline, self_column)


def _select_reveal(X, y, n, k, candidates=None, deadline=None, self_column=None):
    """One gene's REVEAL search over sizes 1..k, or 0..k when self_column is given (see _scan_reveal,
    _select_one_gene)."""
    return _select_one_gene(_scan_reveal, X, y, n, k, candidates, deadline, self_column)


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


def _infer_by_subset_search(data_matrices, scaffold_network, max_indegree, scan_subsets,
                            emit_static_as_input, allow_additional_edges=True, timeout_secs=None):
    """Shared driver for the topology-free methods (REVEAL, Best-Fit). For each gene it (optionally) emits a
    static gene as an input node, else searches regulator sets of size up to the in-degree bound k with
    scan_subsets and builds the gene's BooleanSymbolicFunc from the chosen set. k = max_indegree, or the
    scaffold's max in-degree when max_indegree == -1; k < 1 (e.g. an edgeless scaffold with -1) yields an
    all-input-node model.

    The search is breadth-first over regulator-set size: every gene is swept at size 1, then every gene at
    size 2, and so on. Sweeping one gene through all its sizes before starting the next would let the genes
    that happen to come first spend the whole budget and leave the rest at size 1, which on a large network
    is the normal outcome rather than an edge case - the size-2 sweep alone is n choose 2 sets per gene. This
    way an overrun leaves every gene searched to the same depth, so the result does not depend on gene order.

    allow_additional_edges=False confines a gene's candidate regulators to its own scaffold predecessors, so
    the methods choose the best subset of the given topology and can only ever invent no new edge. A gene
    with no scaffold predecessors then has nothing to search and stays an input node. (With the flag on -
    the published, topology-free behaviour - every gene is a candidate for every other.)

    timeout_secs bounds the whole call rather than any single gene: the deadline is set once here and passed
    down, and the sweep stops enlarging regulator sets once it has passed, keeping the best sets found so
    far. Size 1 is always swept in full for every gene, so the returned model is always complete - a run that
    overruns degrades toward single-regulator functions instead of failing.
    """
    data_matrices = list(data_matrices)
    vertices = list(scaffold_network.vertices)
    names = [v.name for v in vertices]
    n = len(names)
    X, Y = _transition_table(data_matrices)
    k = max_indegree if max_indegree != -1 else scaffold_network.max_in_degree()
    k = min(k, n)  # can't select more regulators than there are genes
    deadline = None if timeout_secs is None else time.time() + timeout_secs

    edges, functions = [], [None] * n
    if X is not None and k >= 1:
        # the genes that get a search at all, in index order, each with its candidates and running best
        searched = {}
        for i in range(n):
            y = Y[:, i]
            if emit_static_as_input and np.all(y == X[:, i]):
                continue  # only ever holds its value (source node / constant) -> input node
            if allow_additional_edges:
                candidates = list(range(n))
            else:
                candidates = sorted(u.index for u in vertices[i].predecessors())
                if not candidates:
                    continue  # no scaffold parents and no new edges allowed -> input node
            searched[i] = {'y': y, 'candidates': candidates, 'best': None, 'done': False,
                           'self_column': X[:, i]}

        max_size = min(k, max((len(gene['candidates']) for gene in searched.values()), default=0))
        # size 0 - the gene as an input node, holding its value - is swept before any regulator set, so a
        # gene that mostly holds its value is not given a regulator it does not need. emit_static_as_input
        # off reproduces the raw methods, which have no such option (an input node becomes a self-loop).
        first_size = 0 if emit_static_as_input else 1
        for size in range(first_size, max_size + 1):
            if size > 1 and _out_of_time(deadline):
                break  # every gene is searched to the same depth, so stop the whole level rather than some
            timed_out = False
            for gene in searched.values():
                if gene['done'] or size > len(gene['candidates']):
                    continue
                gene['best'], gene['done'], timed_out = scan_subsets(
                    X, gene['y'], gene['candidates'], size, gene['best'], deadline, gene['self_column'])
                if timed_out:
                    break
            if timed_out:
                break

        for i, gene in searched.items():
            if gene['best'] is None:
                continue
            # ascending -> first predecessor is MSB (BooleanSymbolicFunc order)
            subset = tuple(sorted(gene['best'][2]))
            if not subset:
                continue  # the search chose the empty set: an input node, holding its value

            functions[i] = BooleanSymbolicFunc(input_names=[names[j] for j in subset],
                                               boolean_outputs=_truth_table_from_subset(X, gene['y'], subset))
            edges.extend((names[j], names[i]) for j in subset)
    return Network(vertex_names=names, edges=edges, vertex_functions=functions)


def reveal_inference(data_matrices, scaffold_network, max_indegree=-1,
                     emit_static_as_input=EMIT_STATIC_GENES_AS_INPUTS,
                     allow_additional_edges=True, timeout_secs=None, **kwargs):
    """REVEAL (Liang, Fuhrman & Somogyi 1998), topology-free, with an MDL relaxation for noisy data (see
    _select_reveal). max_indegree bounds each gene's regulator-set size (-1 = the scaffold's max in-degree);
    allow_additional_edges=False confines the search to each gene's scaffold predecessors, and timeout_secs
    caps the whole call (both described in _infer_by_subset_search). The defaults keep the published,
    unbounded behaviour; the pipeline passes both explicitly."""
    return _infer_by_subset_search(data_matrices, scaffold_network, max_indegree, _scan_reveal,
                                   emit_static_as_input, allow_additional_edges, timeout_secs)


def best_fit_inference(data_matrices, scaffold_network, max_indegree=-1,
                       emit_static_as_input=EMIT_STATIC_GENES_AS_INPUTS,
                       allow_additional_edges=True, timeout_secs=None, **kwargs):
    """Best-Fit Extension (Lähdesmäki, Shmulevich & Yli-Harja 2003), topology-free, with unobserved truth-table
    rows filled by a fair coin (the relaxation to a total function). max_indegree bounds each gene's
    regulator-set size (-1 = the scaffold's max in-degree); allow_additional_edges=False confines the search
    to each gene's scaffold predecessors, and timeout_secs caps the whole call (both described in
    _infer_by_subset_search). The defaults keep the published, unbounded behaviour; the pipeline passes both
    explicitly."""
    return _infer_by_subset_search(data_matrices, scaffold_network, max_indegree, _scan_best_fit,
                                   emit_static_as_input, allow_additional_edges, timeout_secs)
