"""START-BIT adapted to unordered single-cell snapshots: a heuristic solution of the snapshot-tree ILP
(snapshot_tree_ilp.py, joint_orientation_ilp.tex) without a solver, used both as a benchmark method and as a
complete MIP start. Run from the repository root.

The heuristic alternates between a tree step and a function step, from several starts:

1. Tree step. With the network F fixed and the latent states equal to the data, the cost of making a the
   predecessor of c is the number of genes on which c differs from a although F leaves that gene unchanged at
   a: the bits the ILP would have to flip. Placing c under the pseudo-root costs rho * n bits, the ILP's root
   penalty in bits. Candidate predecessors are the symmetrized k nearest neighbours of
   snapshot_tree_ilp.candidate_predecessors. The tree is a minimum-cost spanning arborescence rooted at the
   pseudo-root (Edmonds, arborescence.py). Ties are broken by a seeded jitter on the edge costs, far smaller than one bit.
2. Function step (weighted, the default). With the tree fixed, the ILP's terms separate by gene. For state a
   and gene j, let m_aj be the number of a's children that differ from a on j. Making j unstable at a costs
   gamma and explains all m_aj changes; keeping it stable costs m_aj flips. So F_j is fitted to one sample per
   cell a: target "unstable" (1 - X_aj) with weight m_aj - gamma if that is positive, else target "stable"
   (X_aj) with weight gamma - m_aj. An edge costs beta ln N bits, as in the ILP objective.
   The fit is START-BIT's greedy search (signs from weighted correlation, best threshold by enumeration, then
   single sign changes) with sample weights. Candidate inputs are the regulators and the node itself, at most d
   inputs, and a fit that uses a negative self-loop is refitted without the node's own input.
3. Repair, once, on the returned round only. Going down the tree from the pseudo-root, each violating bit of a
   child is set to its parent's value, which makes every transition constraint hold, for the ILP's MIP start and
   for an objective comparable with the ILP's. The pseudo-root's state is the per-gene majority of its children's
   detection calls.

Both steps minimize the same surrogate cost: the ILP objective with the latent states equal to the data and each
violation counted as one flip, so a lasting change costs one flip whatever the size of the subtree below it. Each
start begins with the tree step for F with every gene stable everywhere, where an edge costs the Hamming distance
between the two cells; it then alternates function and tree steps for up to L rounds (--iterations), stopping when
a round reproduces the previous tree. Each of the S starts (--starts) uses its own tie-breaking seed, and so its
own initial minimum Hamming tree. The round with the best surrogate over all starts is repaired and returned, and
the spread of the starts' best surrogates is reported.

With symmetric Hamming costs the arborescence is a minimum spanning forest of the neighbour graph, one component
per root, and any cell of a component can be its root at equal cost. The orientation of the initial tree is
therefore set by the tie-breaking alone; the first function step is what breaks the symmetry.

The earlier variants are kept as options. --init all-pairs replaces the initial Hamming tree by a first function
step fitted (unweighted) to every candidate pair in both directions, optionally only on the pairs where each gene
changes (--changes-only). --refit synchronous refits F to the tree's parent-to-child transitions unweighted,
counting every gene of every edge. Fitting all pairs is not neutral about stability: for gene j, a pair on which
j does not change is scored correct iff j is stable at the pair's source state, so the all-pairs fit is the
changes-only fit plus a stability bonus at each state weighted by the number of its neighbours that agree with it
on j, 26.5 on average on the LMPP branch (k = 20), i.e. an effective gamma of about 26 against the ILP's 0.5;
--changes-only is gamma = 0. Details and results: "Benchmark: START-BIT on candidate pairs" in
joint_orientation_ilp.tex.
"""

import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arborescence import minimum_arborescence  # noqa: E402
from inference.binary_inference_ideas import _fit_threshold_node, _network_from_signs_and_thresholds  # noqa: E402
from snapshot_tree_ilp import (PSEUDO_ROOT, candidate_predecessors, hamey_regulators, lmpp_cells,  # noqa: E402
                               load_detection_data)

# total tie-breaking jitter on a tree's cost stays below this many bits
JITTER_BITS = 1e-3


def smooth_majority(detected, size=3):
    """Each cell's detection calls replaced by the per-gene majority over the cell and its size - 1 nearest
    cells in Hamming distance (ties broken by index). Hamey et al. smooth with the mode of 3 consecutive cells in
    pseudotime; without an order, the nearest cells take their place. size must be odd; 1 leaves the data as
    is."""
    if size % 2 == 0:
        raise ValueError("size must be odd, got {}".format(size))
    if size == 1:
        return detected.copy()
    distances = (detected[:, None, :] != detected[None, :, :]).sum(axis=2)
    np.fill_diagonal(distances, -1)  # the cell itself comes first
    nearest = np.argsort(distances, axis=1, kind="stable")[:, :size]
    return (detected[nearest].sum(axis=1) * 2 > size).astype(detected.dtype)


def transition_pairs(predecessors):
    """Every ordered (a, c) with a in P(c), cells numbered from 1."""
    return [(a, c) for c, cells in sorted(predecessors.items()) for a in cells]


def tree_pairs(parent):
    """The tree's cell-to-cell edges as ordered (parent, child) pairs; edges from the pseudo-root, which has no
    data row, are left out."""
    return [(a, c) for c, a in sorted(parent.items()) if a != PSEUDO_ROOT]


def _candidate_inputs(genes, regulators, j):
    return sorted({genes.index(gene) for gene in regulators if gene in genes} | {j})


def fit_network(detected, genes, regulators, pairs, edge_penalty=1.0, max_indegree=None,
                forbid_self_repression=True, changes_only=False):
    """Unweighted START-BIT: per node a ({input index: +1/-1}, threshold) pair, fitted to the transitions in
    pairs, scoring agreed transition cells per node. An edge costs edge_penalty * ln N bits as in the ILP,
    scaled by len(pairs) / N transitions per cell into transition cells. With
    changes_only, node j is fitted only on the pairs where j changes, with the same edge cost."""
    n_genes = detected.shape[1]
    max_indegree = max_indegree or len(set(regulators) & set(genes))
    sources = np.array([a - 1 for a, _ in pairs])
    targets = np.array([c - 1 for _, c in pairs])
    x, y = detected[sources], detected[targets]
    n_cells = detected.shape[0]
    edge_cost = edge_penalty * np.log(n_cells) * len(pairs) / n_cells
    fitted = []
    for j in range(n_genes):
        candidates = _candidate_inputs(genes, regulators, j)
        rows = (x[:, j] != y[:, j]) if changes_only else np.ones(len(x), dtype=bool)
        if not rows.any():
            fitted.append(({}, 0))
            continue
        xj, yj = x[rows], y[rows, j]
        signs, threshold = _fit_threshold_node(xj, yj, xj[:, j], candidates, {i: edge_cost for i in candidates},
                                               max_indegree)
        if forbid_self_repression and signs.get(j, 0) < 0:
            candidates = [i for i in candidates if i != j]
            signs, threshold = _fit_threshold_node(xj, yj, xj[:, j], candidates,
                                                   {i: edge_cost for i in candidates}, max_indegree)
        fitted.append((signs, threshold))
    return fitted


def fit_node_weighted(x, target, weights, own, candidates, edge_cost, max_degree=None):
    """START-BIT's greedy search (_fit_threshold_node) with sample weights: maximizes the total weight of the
    samples whose target the function reproduces, minus edge_cost per input. A node with no input holds its
    value (own). Returns ({input index: sign}, threshold)."""
    columns = {i: x[:, i] for i in candidates}
    total = weights.sum()

    def value(signs):
        used = [i for i in candidates if signs[i] != 0]
        if max_degree is not None and len(used) > max_degree:
            return -np.inf, 1
        if not used:
            return float(weights[own == target].sum()), 0
        agreeing = sum((columns[i] == 1) if signs[i] > 0 else (columns[i] == 0) for i in used).astype(int)
        best_score, best_t = -1.0, 1
        for t in range(1, len(used) + 1):
            score = float(weights[(agreeing >= t) == target].sum())
            if score > best_score:
                best_score, best_t = score, t
        return best_score - edge_cost * len(used), best_t

    def weighted_correlation(column):
        if total == 0:
            return 0.0
        mx, my = (weights * column).sum() / total, (weights * target).sum() / total
        cov = (weights * (column - mx) * (target - my)).sum()
        sx = np.sqrt((weights * (column - mx) ** 2).sum())
        sy = np.sqrt((weights * (target - my) ** 2).sum())
        return 0.0 if sx == 0 or sy == 0 else cov / (sx * sy)

    correlation = {i: weighted_correlation(columns[i].astype(float)) for i in candidates}
    signs = {i: (0 if correlation[i] == 0 else (1 if correlation[i] > 0 else -1)) for i in candidates}
    if max_degree is not None:
        for i in sorted(candidates, key=lambda i: -abs(correlation[i]))[max_degree:]:
            signs[i] = 0
    current = value(signs)
    improved = True
    while improved:
        improved = False
        for i in candidates:
            for sign in (1, -1, 0):
                if sign == signs[i]:
                    continue
                trial = dict(signs)
                trial[i] = sign
                result = value(trial)
                if result[0] > current[0]:
                    signs, current, improved = trial, result, True
    return {i: s for i, s in signs.items() if s != 0}, current[1]


def _weighted_score(x, target, weights, own, signs, threshold, edge_cost):
    """fit_node_weighted's score of one function: weight of reproduced targets minus edge_cost per input."""
    if not signs:
        return float(weights[own == target].sum())
    agreeing = sum((x[:, i] == 1) if sign > 0 else (x[:, i] == 0) for i, sign in signs.items())
    return float(weights[(agreeing >= threshold) == target].sum()) - edge_cost * len(signs)


def fit_network_weighted(detected, genes, regulators, parent, edge_penalty=1.0, stability_penalty=0.5,
                         max_indegree=None, forbid_self_repression=True, previous=None):
    """The function step on the ILP's own terms, given the tree (see the module docstring, step 2). The greedy
    search is local, so each node keeps the best of its greedy fit, no inputs (holding its value) and, if
    given, its previous function; the step then never worsens the surrogate cost."""
    n_cells, n_genes = detected.shape
    max_indegree = max_indegree or len(set(regulators) & set(genes))
    changed_children = np.zeros(detected.shape)
    for a, c in tree_pairs(parent):
        changed_children[a - 1] += detected[c - 1] != detected[a - 1]
    edge_cost = edge_penalty * np.log(n_cells)
    fitted = []
    for j in range(n_genes):
        gain = changed_children[:, j] - stability_penalty
        target = np.where(gain > 0, 1 - detected[:, j], detected[:, j])
        weights = np.abs(gain)
        candidates = _candidate_inputs(genes, regulators, j)
        signs, threshold = fit_node_weighted(detected, target, weights, detected[:, j], candidates, edge_cost,
                                             max_indegree)
        if forbid_self_repression and signs.get(j, 0) < 0:
            candidates = [i for i in candidates if i != j]
            signs, threshold = fit_node_weighted(detected, target, weights, detected[:, j], candidates, edge_cost,
                                                 max_indegree)
        options = [(signs, threshold), ({}, 0)] + ([previous[j]] if previous is not None else [])
        scores = [_weighted_score(detected, target, weights, detected[:, j], option[0], option[1], edge_cost)
                  for option in options]
        fitted.append(options[int(np.argmax(scores))])
    return fitted


def apply_network(fitted, states):
    """F applied to every row of states (a 0/1 array), vectorized. A node with no input holds its value; F=None
    holds every value."""
    if fitted is None:
        return states.copy()
    outputs = np.empty_like(states)
    for j, (signs, threshold) in enumerate(fitted):
        if not signs:
            outputs[:, j] = states[:, j]
            continue
        agreeing = sum((states[:, i] == 1) if sign > 0 else (states[:, i] == 0) for i, sign in signs.items())
        outputs[:, j] = (agreeing >= threshold).astype(states.dtype)
    return outputs


def violations(parent_states, parent_outputs, child_states):
    """Per pair, the genes on which the child differs from its parent where the parent is stable."""
    return (child_states != parent_states) & (parent_outputs == parent_states)


def build_tree(detected, fitted, predecessors, root_penalty, seed=None):
    """{c: predecessor} (PSEUDO_ROOT for roots), a minimum-cost spanning arborescence of the candidate graph
    rooted at the pseudo-root (arborescence.minimum_arborescence), costs in flipped bits on the detection calls.
    fitted=None is F with every gene stable, for which an edge costs the Hamming distance. With seed, each edge
    cost gets a uniform jitter whose total over a tree stays below JITTER_BITS, which breaks ties at random
    without changing the optimum."""
    n_cells, n_genes = detected.shape
    outputs = apply_network(fitted, detected)
    src, dst, weight = [], [], []
    for c, cells in sorted(predecessors.items()):
        parents = np.array(cells) - 1
        costs = violations(detected[parents], outputs[parents], detected[c - 1][None, :]).sum(axis=1)
        src.extend([PSEUDO_ROOT] + list(cells))
        dst.extend([c] * (len(cells) + 1))
        weight.extend([root_penalty * n_genes] + [float(cost) for cost in costs])
    src, dst, weight = np.array(src), np.array(dst), np.array(weight)
    if seed is not None:
        # the same seed gives every edge the same jitter in every round, so ties are broken consistently and a
        # round can reproduce the previous tree
        weight = weight + np.random.default_rng(seed).uniform(0, JITTER_BITS / n_cells, size=len(weight))
    chosen = minimum_arborescence(n_cells + 1, src, dst, weight, root=PSEUDO_ROOT)
    return {int(dst[e]): int(src[e]) for e in chosen}


def repair_states(detected, fitted, parent):
    """(latent states with the pseudo-root as row 0, pseudotime per cell): walking down the tree, set each
    child's violating bits to its parent's (final) values, so every transition constraint holds."""
    n_cells, n_genes = detected.shape
    states = np.vstack([np.zeros(n_genes, dtype=detected.dtype), detected.copy()])
    roots = [c for c, a in parent.items() if a == PSEUDO_ROOT]
    states[PSEUDO_ROOT] = (detected[np.array(roots) - 1].mean(axis=0) >= 0.5).astype(detected.dtype)
    children = {}
    for c, a in parent.items():
        children.setdefault(a, []).append(c)
    pseudotime = {}
    frontier = [(PSEUDO_ROOT, 0)]
    while frontier:
        a, depth = frontier.pop()
        output = apply_network(fitted, states[a][None, :])[0]
        for c in children.get(a, []):
            bad = violations(states[a], output, states[c])
            states[c][bad] = states[a][bad]
            pseudotime[c] = depth + 1
            frontier.append((c, depth + 1))
    assert len(pseudotime) == n_cells, "the tree does not reach every cell"
    return states, pseudotime


def objective(detected, fitted, states, parent, edge_penalty, stability_penalty, root_penalty):
    """The ILP's objective for this solution, and its terms."""
    n_cells, n_genes = detected.shape
    bits = float(n_cells * n_genes)
    agreement = (states[1:] == detected).sum() / bits
    edges = 0 if fitted is None else sum(len(signs) for signs, _ in fitted)
    unstable = int((apply_network(fitted, states) != states).sum())
    roots = sum(a == PSEUDO_ROOT for a in parent.values())
    value = (agreement - edge_penalty * np.log(n_cells) / bits * edges - stability_penalty / bits * unstable
             - root_penalty / n_cells * roots)
    return value, {"agreement": agreement, "edges": edges, "unstable": unstable, "roots": roots,
                   "flipped_bits": int((states[1:] != detected).sum())}


def surrogate(detected, fitted, parent, edge_penalty, stability_penalty, root_penalty):
    """The cost both steps minimize, as an objective: the ILP objective with the latent states equal to the
    data, each tree-edge violation counted as one flip, and instability counted at the cells' data states."""
    n_cells, n_genes = detected.shape
    bits = float(n_cells * n_genes)
    outputs = apply_network(fitted, detected)
    pairs = tree_pairs(parent)
    flips = 0 if not pairs else int(violations(detected[[a - 1 for a, _ in pairs]],
                                               outputs[[a - 1 for a, _ in pairs]],
                                               detected[[c - 1 for _, c in pairs]]).sum())
    edges = 0 if fitted is None else sum(len(signs) for signs, _ in fitted)
    roots = sum(a == PSEUDO_ROOT for a in parent.values())
    return (1 - flips / bits - edge_penalty * np.log(n_cells) / bits * edges
            - stability_penalty / bits * int((outputs != detected).sum()) - root_penalty / n_cells * roots)


def run_start(detected, genes, regulators, predecessors, seed, iterations, edge_penalty, stability_penalty,
              root_penalty, max_indegree, init="hamming", refit="weighted", changes_only=False, log=None):
    """One start: the best round by the surrogate as (fitted, parent, surrogate, terms), and the per-round terms.
    Round 0 is the initial tree (init "hamming": F with every gene stable; init "all-pairs": the
    unweighted fit to all candidate pairs and its tree)."""
    if init == "hamming":
        fitted = None
    else:
        fitted = fit_network(detected, genes, regulators, transition_pairs(predecessors), edge_penalty,
                             max_indegree, changes_only=changes_only)
    parent, previous = build_tree(detected, fitted, predecessors, root_penalty, seed), None
    best, history = None, []
    for round_index in range(iterations + 1):
        if round_index > 0:
            if refit == "weighted":
                fitted = fit_network_weighted(detected, genes, regulators, parent, edge_penalty, stability_penalty,
                                              max_indegree, previous=fitted)
            else:
                fitted = fit_network(detected, genes, regulators, tree_pairs(parent), edge_penalty, max_indegree)
            previous, parent = parent, build_tree(detected, fitted, predecessors, root_penalty, seed)
        value = surrogate(detected, fitted, parent, edge_penalty, stability_penalty, root_penalty)
        terms = {"start": seed, "round": round_index, "surrogate": value,
                 "edges": 0 if fitted is None else sum(len(signs) for signs, _ in fitted),
                 "roots": sum(a == PSEUDO_ROOT for a in parent.values()),
                 "parents_changed": None if previous is None else sum(parent[c] != previous[c] for c in parent)}
        history.append(terms)
        if log is not None:
            log(terms)
        if best is None or value > best[2]:
            best = (fitted, parent, value, terms)
        if previous == parent:
            break
    return best, history


def solve(detected, genes, regulators, k=20, edge_penalty=1.0, stability_penalty=0.5, root_penalty=None,
          max_indegree=None, changes_only=False, iterations=10, starts=5, init="hamming", refit="weighted",
          seed=0, log=None):
    """The full heuristic over `starts` starts with seeds seed, seed + 1, ...: for the round with the best surrogate,
    repaired, (fitted network, {c: predecessor}, states with the pseudo-root as row 0, {c: pseudotime}, ILP objective,
    objective terms). terms["history"] holds every round of every start, and terms["starts"] each start's best
    surrogate with their mean, standard deviation, minimum and maximum. Defaults as in
    snapshot_tree_ilp.build_model."""
    root_penalty = 6.0 / detected.shape[1] if root_penalty is None else root_penalty
    predecessors = candidate_predecessors(detected, k)
    best, history, start_values = None, [], []
    for start in range(starts):
        result, start_history = run_start(detected, genes, regulators, predecessors, seed + start, iterations,
                                          edge_penalty, stability_penalty, root_penalty, max_indegree, init, refit,
                                          changes_only, log)
        history.extend(start_history)
        start_values.append(result[2])
        if best is None or result[2] > best[2]:
            best = result
    fitted, parent, _, round_terms = best
    states, pseudotime = repair_states(detected, fitted, parent)
    value, terms = objective(detected, fitted, states, parent, edge_penalty, stability_penalty, root_penalty)
    terms = dict(terms, objective=value, max_pseudotime=max(pseudotime.values()), start=round_terms["start"],
                 round=round_terms["round"], surrogate=round_terms["surrogate"])
    values = np.array(start_values)
    spread = {"best_per_start": [float(v) for v in values], "mean": float(values.mean()),
              "std": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
              "min": float(values.min()), "max": float(values.max())}
    return fitted, parent, states, pseudotime, value, dict(terms, history=history, starts=spread)


def set_mip_start(handles, fitted, parent, states, pseudotime):
    """Start values for snapshot_tree_ilp.build_model's variables: the functions, latent states, predecessors
    and pseudotimes. Gurobi completes the rest (function outputs, instability, auxiliaries), which they
    determine."""
    for j, (positives, negatives, threshold, unused) in enumerate(handles["functions"]):
        signs, value = fitted[j] if fitted is not None else ({}, 0)
        for i, p, q in zip(handles["candidate_inputs"][j], positives, negatives):
            p.Start, q.Start = int(signs.get(i, 0) > 0), int(signs.get(i, 0) < 0)
        threshold.Start, unused.Start = (value, 0) if signs else (0, 1)
    for row, variables in zip(states, handles["states"]):
        for var, value in zip(variables, row):
            var.Start = int(value)
    for (a, c), indicator in handles["pi"].items():
        indicator.Start = int(parent[c] == a)
    for c, tau in enumerate(handles["tau"]):
        if c != PSEUDO_ROOT:
            tau.Start = pseudotime[c]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cells", choices=["all", "lmpp"], default="lmpp")
    parser.add_argument("--max-cells", type=int, default=None)
    parser.add_argument("--k", type=int, default=20)
    parser.add_argument("--edge-penalty", type=float, default=1.0,
                        help="beta: an edge costs beta ln N changed data values")
    parser.add_argument("--stability-penalty", type=float, default=0.5)
    parser.add_argument("--root-penalty", type=float, default=None)
    parser.add_argument("--max-indegree", type=int, default=None)
    parser.add_argument("--iterations", type=int, default=10,
                        help="L, rounds of alternating function and tree steps per start")
    parser.add_argument("--starts", type=int, default=5, help="S, starts with different tie-breaking seeds")
    parser.add_argument("--seed", type=int, default=0, help="seed of the first start")
    parser.add_argument("--init", choices=["hamming", "all-pairs"], default="hamming")
    parser.add_argument("--refit", choices=["weighted", "synchronous"], default="weighted")
    parser.add_argument("--changes-only", action="store_true",
                        help="with --init all-pairs, fit each gene only on the pairs where it changes")
    parser.add_argument("--smooth", type=int, default=1,
                        help="majority-smooth each cell over itself and its nearest cells, this many in all (odd; "
                             "Hamey et al. use 3); the smoothed calls then replace the data throughout")
    parser.add_argument("--out", default=os.path.join("single_cell_inference", "results_heuristic"))
    options = parser.parse_args()

    genes, names, types, detected = load_detection_data()
    if options.cells == "lmpp":
        rows = [names.index(name) for name in lmpp_cells()]
        names, types, detected = [names[r] for r in rows], [types[r] for r in rows], detected[rows]
    if options.max_cells:
        names, types, detected = names[:options.max_cells], types[:options.max_cells], detected[:options.max_cells]

    if options.smooth > 1:
        raw = detected
        detected = smooth_majority(raw, options.smooth)
        nearest = lambda x: ((x[:, None, :] != x[None, :, :]).sum(axis=2) + np.eye(len(x), dtype=int) * 99).min(axis=1)
        print("smoothed over {} cells: {:.1%} of values changed; mean nearest-neighbour distance {:.2f} -> {:.2f}; "
              "distinct states {} -> {}".format(options.smooth, (detected != raw).mean(), nearest(raw).mean(),
                                                nearest(detected).mean(), len(np.unique(raw, axis=0)),
                                                len(np.unique(detected, axis=0))), flush=True)

    def log(terms):
        print({key: (round(float(v), 5) if key in ("objective", "agreement", "surrogate") else v)
               for key, v in terms.items()}, flush=True)

    started = time.time()
    fitted, parent, states, pseudotime, value, terms = solve(
        detected, genes, hamey_regulators(), options.k, options.edge_penalty, options.stability_penalty,
        options.root_penalty, options.max_indegree, options.changes_only, options.iterations, options.starts,
        options.init, options.refit, options.seed, log)
    network = _network_from_signs_and_thresholds(genes, fitted if fitted is not None else [({}, 0)] * len(genes))
    unstable = apply_network(fitted, states) != states
    os.makedirs(options.out, exist_ok=True)
    network.save(os.path.join(options.out, "network.json"))
    cells = {name: {"predecessor": None if parent[c] == PSEUDO_ROOT else names[parent[c] - 1],
                    "pseudotime": pseudotime[c], "state": [int(v) for v in states[c]],
                    "unstable": [g for g, u in zip(genes, unstable[c]) if u], "cell_type": cell_type}
             for c, (name, cell_type) in enumerate(zip(names, types), start=1)}
    terms = {key: v.item() if isinstance(v, np.generic) else v for key, v in terms.items()}
    with open(os.path.join(options.out, "cells.json"), "w") as cells_file:
        json.dump({"genes": genes, "pseudo_root_state": [int(v) for v in states[PSEUDO_ROOT]], "cells": cells,
                   "objective": value, "terms": terms, "options": vars(options)}, cells_file, indent=1, default=float)
    spread = terms["starts"]
    print("{} cells x {} genes in {:.0f}s: objective {:.5f} (start {}, round {}); surrogate over {} starts mean "
          "{:.5f}, sd {:.5f}, min {:.5f}, max {:.5f} -> {}".format(
              len(names), len(genes), time.time() - started, value, terms["start"], terms["round"],
              len(spread["best_per_start"]), spread["mean"], spread["std"], spread["min"], spread["max"], options.out))


if __name__ == "__main__":
    main()
