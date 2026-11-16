"""Loading and plotting behind paper_figures.ipynb.

Reads the per-network score files an inference experiment wrote under inferred_models/generated/<experiment>,
and draws three figures per experiment: the accuracy panels (hit-rate, edge-set similarity, edge-sign
agreement), the run-time panels (training time, fraction of networks reaching the time limit), and the fraction
of networks whose inference failed. Each panel is a bar per method, or a line per method along a swept parameter (K, network size).

The scoring rules follow post_generals_analysis.ipynb, so the numbers match its figures:
  - a network a method has no output for, while another method of the same condition does, counts as a failure:
    0 on every score, and the condition's longest successful run time as its time;
  - a varying-node score is NaN when nothing varies, and such a network is left out of that metric;
  - Constant Nodes has no edges, so it is left out of the edge metric;
  - a run reached the time limit when it never finished or took at least 95% of the longest successful run in
    the whole experiment;
  - a network failed when its inference never finished (no inference_time.npy).
Score files older runs did not write - the varying-node hit-rates and the edge-sign agreement - are computed
here by compute_missing_scores, as post_generals_analysis.ipynb's backfill and edge-sign cells do.
"""
import itertools
import json
import os
import re

import matplotlib
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns
from matplotlib import pyplot as plt

from attractor_learning import graphs
from attractor_learning.logic import SymmetricThresholdFunction
from validation import inference_scoring

# ---- methods -----------------------------------------------------------------------------------------------
# The methods drawn, in legend order. Colours and the font sizes below are the ones post_generals_analysis uses.
METHODS = ["all_constants", "random_model", "linear_classifier", "best_fit", "reveal",
           "symmetric_heuristic", "symmetric_topology"]

METHOD_DISPLAY_NAMES = {
    "all_constants": "Constant Nodes",
    "random_model": "Random Model on PKN",
    "linear_classifier": "Logistic Regression",
    "best_fit": "Best-Fit",
    "reveal": "REVEAL",
    "symmetric_heuristic": "START-BIT (heuristic)",
    "symmetric_topology": "STOP-BIT (ILP)",
}

METHOD_COLORS = {
    "all_constants": "#a3a3a3",
    "random_model": "#4d4d4d",
    "linear_classifier": "#56b4e9",
    "best_fit": "#009e73",
    "reveal": "#08306b",
    "symmetric_heuristic": "#ffb000",
    "symmetric_topology": "#b0207f",
}

# older runs named a variant directory by the inference function rather than the config name
_METHOD_ALIASES = {
    "infer_unknown_topology_symmetric": "symmetric_topology",
    "all_constants_inference": "all_constants",
    "random_model_inference": "random_model",
    "linear_classifier_inference": "linear_classifier",
    "reveal_inference": "reveal",
    "best_fit_inference": "best_fit",
}

EDGELESS_METHODS = ("all_constants",)
LEGEND_TITLE = "Method"

plt.rcParams.update({
    'figure.titlesize': 18,
    'axes.titlesize': 16,
    'axes.labelsize': 15,
    'xtick.labelsize': 13,
    'ytick.labelsize': 13,
    'legend.fontsize': 13,
    'legend.title_fontsize': 14,
})
VALUE_LABEL_FONTSIZE = 11

# ---- metrics -----------------------------------------------------------------------------------------------
# metric -> (panel title, value-axis label, logarithmic value axis)
METRICS = {
    'hit_rate': ("Hit-Rate on All Nodes", "Hit-rate on test set", False),
    'hit_rate_varying': ("Hit-Rate on Varying Nodes", "Hit-rate on test set", False),
    'edge': ("Edge Set Jaccard Similarity", "Edge sets Jaccard similarity", False),
    'sign': ("Edge Sign Agreement", "Fraction of shared edges with the true sign", False),
    'timing': ("Training Time", "training time (s)", True),
    'time_limit': ("Reaching the Time Limit", "fraction of networks reaching the time limit", False),
    'failed': ("Failed Inference", "fraction of networks with failed inference", False),
}
ACCURACY_PANELS = ['hit_rate', 'hit_rate_varying', 'edge', 'sign']
TIME_PANELS = ['timing', 'time_limit']
FAILED_PANELS = ['failed']

# the score each metric is read from; hit-rates are rolled out from the real first state of each test series
SCORE_FILES = {
    'hit_rate': 'timeseries_real_start_accuracy_score_test.npy',
    'hit_rate_varying': 'timeseries_real_start_varying_accuracy_score_test.npy',
    'edge': 'edge_accuracy_score.npy',
    'sign': 'edge_sign_agreement_score.npy',
}
TIME_FILE = "inference_time.npy"
TIME_LIMIT_FRACTION = 0.95
_MODEL_MARKER_FILES = (TIME_FILE, "inferred_network.json", "edge_accuracy_score.npy",
                       "timeseries_real_accuracy_score_train.npy")


# ---- loading -----------------------------------------------------------------------------------------------
def canonical_method(name):
    return _METHOD_ALIASES.get(name, name)


def _dir_name_tokens(name):
    """(params, plain tokens) of one directory name. A variant directory is named <comb_str>-<data folder>, its
    parts separated by commas, or by a dash that does not start a negative number."""
    tokens = [t for t in re.sub(r'-(?=\D)', ',', name).split(',') if t]
    params = dict(t.split('=', 1) for t in tokens if '=' in t)
    return params, [t for t in tokens if '=' not in t]


def _models(exp_dir):
    """(sub-experiment, params, network, model path) for every inferred model under exp_dir. Every key=value
    token on the path is a parameter; the plain names above the last such level (the data folder a run used)
    name the sub-experiment; the levels below it name the network."""
    for root, dirs, _ in os.walk(exp_dir):
        dirs.sort()
        if not any(os.path.exists(os.path.join(root, f)) for f in _MODEL_MARKER_FILES):
            continue
        dirs[:] = []
        parts = os.path.relpath(root, exp_dir).split(os.sep)
        with_params = [i for i, part in enumerate(parts[:-1]) if '=' in part]
        last = with_params[-1] if with_params else len(parts) - 2
        params, labels = {}, []
        for part in parts[:last + 1]:
            part_params, part_labels = _dir_name_tokens(part)
            params.update(part_params)
            labels.extend(part_labels)
        yield "/".join(labels), params, "/".join(parts[last + 1:]), root


def _load_score(path):
    return float(np.load(path)) if os.path.exists(path) else None


def load_experiment(exp_dir, sweep=None, methods=METHODS):
    """One row per (network, method, metric) of an experiment, as a DataFrame with columns
    sub, method, x (the swept parameter's value, NaN without one), network, metric, value, failed."""
    print("Loading", exp_dir)
    variants = {}       # (sub, method, params) -> {network: model path}
    for sub, params, network, path in _models(exp_dir):
        method = canonical_method(params.get('inference_method', ''))
        key = (sub, method, tuple(sorted((k, v) for k, v in params.items() if k != 'inference_method')))
        variants.setdefault(key, {})[network] = path
    if not any(method in methods for _, method, _ in variants):
        raise ValueError("no inferred models of {} under {}".format(", ".join(methods), exp_dir))
    # every method, drawn or not, counts towards which networks a condition should have, the failure run time
    # and the time limit; only the drawn ones get rows

    # parameters that differ between methods but never within one (a method-only option) don't separate
    # conditions; everything else does
    params_df = pd.DataFrame([dict(p) for _, _, p in variants]).assign(method=[m for _, m, _ in variants])
    method_specific = {c for c in params_df.columns if c != 'method' and params_df[c].nunique(dropna=False) > 1
                       and params_df.fillna("").groupby('method')[c].nunique().max() <= 1}

    def condition(key):
        sub, _, params = key
        return sub, tuple(p for p in params if p[0] not in method_specific)

    expected = {}
    for key, models in variants.items():
        expected.setdefault(condition(key), set()).update(models)
    produced = {(key[0], metric) for key, models in variants.items() for metric, f in SCORE_FILES.items()
                if any(os.path.exists(os.path.join(p, f)) for p in models.values())}
    times = {key: {n: _load_score(os.path.join(p, TIME_FILE)) for n, p in models.items()}
             for key, models in variants.items()}
    longest_by_sub = {}
    for (sub, _, _), by_network in times.items():
        finished = [t for t in by_network.values() if t is not None]
        if finished:
            longest_by_sub[sub] = max(longest_by_sub.get(sub, 0.0), max(finished))
    time_limit = TIME_LIMIT_FRACTION * max(longest_by_sub.values()) if longest_by_sub else None

    rows = []
    for key, models in sorted(variants.items()):
        sub, method, params = key
        if method not in methods:
            continue
        x = float(dict(params)[sweep]) if sweep else np.nan
        missing = sorted(expected[condition(key)] - set(models))
        if missing:
            print("  {} {}: {} network(s) with no output, counted as failed".format(
                sub, METHOD_DISPLAY_NAMES[method], len(missing)))

        def add(metric, network, value, failed=False):
            rows.append(dict(sub=sub, method=method, x=x, network=network, metric=metric, value=value,
                             failed=failed, variant=params))

        for metric, fname in SCORE_FILES.items():
            if metric in ('edge', 'sign') and method in EDGELESS_METHODS:
                continue
            if metric == 'sign' and not any(os.path.exists(os.path.join(p, fname)) for p in models.values()):
                continue        # a method without threshold functions has no sign score
            if (sub, metric) not in produced:
                continue
            scores = {n: _load_score(os.path.join(p, fname)) for n, p in models.items()}
            finite = [v for v in scores.values() if v is not None and not np.isnan(v)]
            if not finite and any(v is not None for v in scores.values()):
                continue        # undefined (NaN) on every network that got as far as scoring it
            for network, value in scores.items():
                if value is None:
                    add(metric, network, 0.0, failed=True)
                elif not np.isnan(value):
                    add(metric, network, value)
            for network in missing:
                add(metric, network, 0.0, failed=True)

        stand_in = longest_by_sub.get(sub)
        for network in list(models) + missing:
            t = times[key].get(network)
            if t is not None:
                add('timing', network, t)
            elif stand_in is not None:
                add('timing', network, stand_in, failed=True)
            if time_limit is not None:
                add('time_limit', network, float(t is None or t >= time_limit), failed=t is None)
            add('failed', network, float(t is None), failed=t is None)

    df = pd.DataFrame(rows)
    pooled = df.groupby(['sub', 'method', 'x', 'metric'], dropna=False)['variant'].nunique()
    if (pooled > 1).any():
        print("  pooling several parameter combinations per bar/point:",
              sorted({v for v in df['variant'].unique()}))
    return df.drop(columns='variant')


def varying_hit_rate_differs(df, threshold):
    """(whether some method's mean hit-rate on varying nodes differs from its hit-rate on all nodes by at least
    threshold, at some sub-experiment and sweep value; the largest such difference)."""
    means = (df[df['metric'].isin(['hit_rate', 'hit_rate_varying'])]
             .groupby(['sub', 'method', 'x', 'metric'], dropna=False)['value'].mean().unstack('metric'))
    if 'hit_rate_varying' not in means:
        return False, np.nan
    gap = float((means['hit_rate'] - means['hit_rate_varying']).abs().max())
    return gap >= threshold, gap


# ---- computing missing score files -------------------------------------------------------------------------
# Older runs lack the varying-node hit-rates and the edge-sign agreement. Both are computed from each model's
# saved output and its data folder's ground truth, and written into the model's own directory, which
# load_experiment reads. Whatever is already on disk is left alone, so a second run is a directory scan.
#
# Varying-node scores: the same comparison as the full-matrix score, restricted per matrix to the columns the
# ground truth moves over the scored timepoints (inference_scoring.varying_columns), or to the nodes a data
# folder's varying_nodes.json names. real_start predictions are not saved, so they are re-simulated from
# inferred_network.json. The full-matrix score is recomputed first and checked against the stored one; a
# mismatch means the data folder is not what the run used, and that comparison is reported and skipped.
#
# Edge-sign agreement, for models whose inferred functions are threshold functions (START-BIT, STOP-BIT): the
# fraction of edges both the inferred and the true network sign on which the signs agree, NaN when they share
# none. The true signs come from the data folder's true_signs.json, written from true_network.json when absent:
# a threshold function gives each input its sign; a truth table of up to TRUE_SIGN_MAX_TABLE_INPUTS inputs gives
# +1 to an input that never switches the output off, -1 to one that never switches it on; a constant function
# gives none. A topology-only ground truth uses scaffold_signs.json when the scaffold has the same edges.
# Inferred networks over MAX_INFERRED_NETWORK_BYTES (truth tables of several GB) are not loaded.

VARYING_COMPARISONS = {         # comparison -> (ground-truth matrix file, predictions re-simulated)
    'real': ("real_matrices.npz", False),
    'reference': ("noisy_matrices.npz", False),
    'real_start': ("real_matrices.npz", True),
}
TRUE_SIGNS_FILE = "true_signs.json"
SIGN_SCORE_FILE = SCORE_FILES['sign']
TRUE_SIGN_MAX_TABLE_INPUTS = 16
MAX_INFERRED_NETWORK_BYTES = 100 * 2 ** 20


def _source_dir_for_model(rel_parts, data_exp_dir):
    """The data directory of one inferred model's source network, or None. A model sits at
    <experiment>/[<levels>/]<comb_str>-<data folder>/<the network's path inside that data folder>, so the data
    path is the same path with one component replaced by the data folder's name. Every component is tried,
    shallowest first, matched against the data directories that exist beside it (comb_str can hold dashes)."""
    for i in range(len(rel_parts) - 1):
        base = os.path.join(data_exp_dir, *rel_parts[:i])
        if not os.path.isdir(base):
            break
        part = rel_parts[i]
        names = [d.name for d in os.scandir(base) if d.is_dir()]
        for name in sorted((n for n in names if part == n or part.endswith("-" + n)), key=len, reverse=True):
            candidate = os.path.join(base, name, *rel_parts[i + 1:])
            if os.path.isdir(candidate):
                return candidate
    return None


def _data_dir_resolver(exp_dir, data_exp_dir):
    """model path -> its network's data directory, or None; cached per parent directory."""
    cache = {}

    def resolve(model_path):
        rel_parts = os.path.relpath(model_path, exp_dir).split(os.sep)
        parent_key = tuple(rel_parts[:-1])
        if parent_key not in cache:
            resolved = _source_dir_for_model(rel_parts, data_exp_dir)
            cache[parent_key] = os.path.dirname(resolved) if resolved else None
        if cache[parent_key] is None:
            return None
        net_dir = os.path.join(cache[parent_key], rel_parts[-1])
        return net_dir if os.path.isdir(net_dir) else None

    return resolve


def _varying_scores_for_model(model_path, data_net_dir):
    """Writes the varying-node score files one model is missing. (files written, problem messages)."""
    varying_mask = None
    varying_file = os.path.join(data_net_dir, "varying_nodes.json")
    if os.path.exists(varying_file):
        with open(varying_file) as f:
            varying_names = set(json.load(f))
        scaffold = graphs.Network.load(os.path.join(data_net_dir, "scaffold_network.json"))
        varying_mask = np.array([v.name in varying_names for v in scaffold.vertices])

    targets = []
    for group in ('train', 'test'):
        if not os.path.exists(os.path.join(model_path, "{}_matrices.npz".format(group))):
            continue
        for name, (truth_file, resimulate) in VARYING_COMPARISONS.items():
            out_name = inference_scoring.timeseries_varying_score_filename(name, group)
            full_name = "timeseries_{}_accuracy_score_{}.npy".format(name, group)
            out_path = os.path.join(model_path, out_name + ".npy")
            if os.path.exists(out_path) and not (varying_mask is not None and np.isnan(float(np.load(out_path)))):
                continue
            if os.path.exists(os.path.join(model_path, full_name)):
                targets.append((name, group, truth_file, resimulate, out_name, full_name))
    if not targets:
        return 0, []

    written, problems = 0, []
    truths, predictions, network = {}, {}, []
    for name, group, truth_file, resimulate, out_name, full_name in targets:
        try:
            if truth_file not in truths:
                with np.load(os.path.join(data_net_dir, truth_file)) as f:
                    truths[truth_file] = {k: f[k] for k in f.files}
            truth = truths[truth_file]
            if (group, resimulate) not in predictions:
                with np.load(os.path.join(model_path, "{}_matrices.npz".format(group))) as f:
                    keys = list(f.files)
                    saved = None if resimulate else {k: f[k] for k in keys}
                if resimulate:
                    if not network:
                        network.append(graphs.Network.load(os.path.join(model_path, "inferred_network.json")))
                    saved = {k: np.asarray(network[0].next_states(truth[k][0, :], truth[k].shape[0] - 1))
                             for k in keys}
                predictions[(group, resimulate)] = (keys, saved)
            keys, predicted = predictions[(group, resimulate)]
            recomputed = inference_scoring.timeseries_accuracy_score(truth, predicted, keys=keys)
            stored = float(np.load(os.path.join(model_path, full_name)))
            if not np.isclose(recomputed, stored, rtol=0, atol=1e-9):
                problems.append("{}: recomputed {} {} accuracy {:.6f} != stored {:.6f}; data folder {} does not "
                                "match this run".format(model_path, name, group, recomputed, stored, data_net_dir))
                continue
            np.save(os.path.join(model_path, out_name),
                    inference_scoring.timeseries_accuracy_score(truth, predicted, keys=keys, varying_only=True,
                                                                varying_mask=varying_mask))
            written += 1
        except Exception as e:
            problems.append("{}: {} {} - {}: {}".format(model_path, name, group, type(e).__name__, e))
    return written, problems


def _function_input_signs(vertex, thresholds_only=False):
    """{predecessor name: 1 | -1} for one vertex's function (see the comment above)."""
    function, predecessors = vertex.function, vertex.predecessors()
    if not predecessors or not callable(function):
        return {}
    if isinstance(function, SymmetricThresholdFunction):
        if not 1 <= function.threshold <= len(function.signs):
            return {}
        return {u.name: 1 if sign else -1 for u, sign in zip(predecessors, function.signs)}
    if thresholds_only or len(predecessors) > TRUE_SIGN_MAX_TABLE_INPUTS:
        return {}
    table = np.array([bool(function(*row)) for row in itertools.product([False, True], repeat=len(predecessors))],
                     dtype=int).reshape((2,) * len(predecessors))
    signs = {}
    for i, u in enumerate(predecessors):
        change = np.take(table, 1, axis=i) - np.take(table, 0, axis=i)
        if change.any() and (change >= 0).all():
            signs[u.name] = 1
        elif change.any() and (change <= 0).all():
            signs[u.name] = -1
    return signs


def _network_signs(network, thresholds_only=False):
    return {(u, v.name): sign for v in network.vertices
            for u, sign in _function_input_signs(v, thresholds_only).items()}


def _true_signs(data_net_dir):
    """A data folder's true signs, writing true_signs.json when absent and derivable; None when not."""
    path = os.path.join(data_net_dir, TRUE_SIGNS_FILE)
    if not os.path.exists(path):
        true_network = graphs.Network.load(os.path.join(data_net_dir, "true_network.json"))
        if any(v.function is not None for v in true_network.vertices):
            signs = _network_signs(true_network)
        else:
            scaffold_signs_path = os.path.join(data_net_dir, "scaffold_signs.json")
            if not os.path.exists(scaffold_signs_path):
                return None
            scaffold = graphs.Network.load(os.path.join(data_net_dir, "scaffold_network.json"))
            if {(u.name, v.name) for u, v in true_network.edges} != {(u.name, v.name) for u, v in scaffold.edges}:
                return None
            with open(scaffold_signs_path) as f:
                signs = {(u, v): int(sign) for u, v, sign in json.load(f)}
        with open(path, "w") as f:
            json.dump([[u, v, sign] for (u, v), sign in sorted(signs.items())], f)
    with open(path) as f:
        return {(u, v): int(sign) for u, v, sign in json.load(f)}


def _sign_score_for_model(model_path, true_signs):
    """Writes one model's edge-sign score if it has none and infers threshold functions; whether it wrote one."""
    out_path = os.path.join(model_path, SIGN_SCORE_FILE)
    network_path = os.path.join(model_path, "inferred_network.json")
    if os.path.exists(out_path) or not os.path.exists(network_path):
        return False
    inferred = graphs.Network.load(network_path)
    if not any(isinstance(v.function, SymmetricThresholdFunction) for v in inferred.vertices):
        return False
    inferred_signs = _network_signs(inferred, thresholds_only=True)
    shared = [edge for edge in inferred_signs if edge in true_signs]
    np.save(out_path, np.mean([inferred_signs[e] == true_signs[e] for e in shared]) if shared else np.nan)
    return True


def compute_missing_scores(exp_dir, data_exp_dir, methods=METHODS, varying=True, signs=True,
                           max_problems_shown=5, progress_every=200):
    """Writes the varying-node and edge-sign score files the models of `methods` under exp_dir are missing,
    reading the ground truth from data_exp_dir (data/generated/<experiment>), and prints what it did."""
    print("Computing missing scores for", exp_dir)
    if not os.path.isdir(data_exp_dir):
        print("  no data folder {}; nothing computed".format(data_exp_dir))
        return
    resolve = _data_dir_resolver(exp_dir, data_exp_dir)
    models = [path for _, params, _, path in _models(exp_dir)
              if canonical_method(params.get('inference_method', '')) in methods]
    counts = dict(varying=0, sign=0, unsourced=0, unsigned=0, oversized=0)
    problems, true_signs = [], {}
    for checked, model_path in enumerate(models, start=1):
        data_net_dir = resolve(model_path)
        if data_net_dir is None:
            counts['unsourced'] += 1
            continue
        if varying:
            written, model_problems = _varying_scores_for_model(model_path, data_net_dir)
            counts['varying'] += written
            problems.extend(model_problems)
        if signs:
            network_path = os.path.join(model_path, "inferred_network.json")
            if os.path.exists(network_path) and os.path.getsize(network_path) > MAX_INFERRED_NETWORK_BYTES:
                counts['oversized'] += 1
            else:
                try:
                    if data_net_dir not in true_signs:
                        true_signs[data_net_dir] = _true_signs(data_net_dir)
                    if true_signs[data_net_dir] is None:
                        counts['unsigned'] += 1
                    else:
                        counts['sign'] += _sign_score_for_model(model_path, true_signs[data_net_dir])
                except Exception as e:
                    problems.append("{}: edge signs - {}: {}".format(model_path, type(e).__name__, e))
        if progress_every and checked % progress_every == 0 and checked < len(models):
            print("    {}/{} models checked".format(checked, len(models)))
    notes = [(counts['unsourced'], "with no source network in the data folder"),
             (counts['unsigned'], "whose data has no true signs"),
             (counts['oversized'], "with an inferred network over {} MB, not sign-scored".format(
                 MAX_INFERRED_NETWORK_BYTES // 2 ** 20)),
             (len(problems), "problem(s)")]
    print("  {} model(s): wrote {} varying-node and {} edge-sign score file(s)".format(
        len(models), counts['varying'], counts['sign'])
          + "".join("; {} {}".format(n, text) for n, text in notes if n))
    for message in problems[:max_problems_shown]:
        print("    " + message)
    if len(problems) > max_problems_shown:
        print("    ... and {} more".format(len(problems) - max_problems_shown))


# ---- drawing -----------------------------------------------------------------------------------------------
def _value_text(value):
    if not np.isfinite(value):
        return ""
    if abs(value) >= 1000:
        return "{:,.0f}".format(value)
    if value != 0 and abs(value) < 0.005:
        return "{:.1e}".format(value)
    return "{:.2f}".format(value)


def _add_headroom(ax, top, fraction=0.16):
    """Raise the top of the value axis so labels written above `top` stay inside the axes."""
    if top is None or not np.isfinite(top):
        return
    lo, hi = ax.get_ylim()
    if ax.get_yscale() == 'log':
        needed = top * (top / lo) ** fraction if 0 < lo < top else hi
    else:
        needed = top + fraction * (top - lo)
    if needed > hi:
        ax.set_ylim(lo, needed)


def _log_value_axis(ax):
    if not np.isfinite(ax.dataLim.minposy):
        return
    ax.set_yscale('log')
    lo, hi = ax.get_ylim()
    if hi / lo < 100:
        ax.grid(True, axis='y', which='minor', alpha=0.4, linewidth=0.5)
    ax.grid(True, axis='y', which='major')


def _label_bars(ax):
    """Each bar's value above it (above its error bar where it has one), turned on its side."""
    error_tops = {}
    for line in ax.lines:
        xs = np.asarray(line.get_xdata(), dtype=float)
        ys = np.asarray(line.get_ydata(), dtype=float)
        if len(xs) == 2 and np.isclose(xs[0], xs[1]) and np.isfinite(ys).any():
            key = round(float(xs[0]), 6)
            error_tops[key] = max(error_tops.get(key, -np.inf), float(np.nanmax(ys)))
    top = None
    for container in ax.containers:
        if not isinstance(container, matplotlib.container.BarContainer):
            continue
        for patch in container:
            height = patch.get_height()
            if not np.isfinite(height):
                continue
            x = patch.get_x() + patch.get_width() / 2
            y = max(height, error_tops.get(round(x, 6), height))
            ax.annotate(_value_text(height), (x, y), xytext=(0, 2), textcoords='offset points',
                        ha='center', va='bottom', rotation=90, fontsize=VALUE_LABEL_FONTSIZE)
            top = y if top is None else max(top, y)
    _add_headroom(ax, top)


def _label_line_ends(ax, points, lift_pt=5, spacing_pt=13):
    """Each (x, y) point's value above it; labels at one x are stacked so they don't overlap."""
    by_x = {}
    for x, y in points:
        if np.isfinite(y):
            by_x.setdefault(x, []).append(y)
    if not by_x:
        return
    tallest = max(len(ys) for ys in by_x.values())
    _add_headroom(ax, max(max(ys) for ys in by_x.values()), fraction=0.05 + 0.045 * tallest)
    to_points = 72.0 / ax.figure.dpi
    for x, ys in by_x.items():
        below = None
        for y in sorted(ys):
            point_pt = ax.transData.transform((x, y))[1] * to_points
            label_pt = point_pt + lift_pt if below is None else max(point_pt + lift_pt, below + spacing_pt)
            below = label_pt
            ax.annotate(_value_text(y), (x, y), xytext=(0, label_pt - point_pt), textcoords='offset points',
                        ha='center', va='bottom', fontsize=VALUE_LABEL_FONTSIZE)


def _draw_bars(ax, data, methods, log_y):
    """Mean per method with a 95% confidence interval, one group of bars per sub-experiment."""
    subs = sorted(data['sub'].unique())
    sns.barplot(ax=ax, data=data, x='sub', y='value', order=subs, hue='method', hue_order=methods,
                palette=[METHOD_COLORS[m] for m in methods])
    if len(subs) == 1:
        ax.set_xticks([])
    ax.set_xlabel("")
    if log_y:
        _log_value_axis(ax)
    _label_bars(ax)


def _draw_lines(ax, data, methods, log_y, x_label, log_x_ratio=20):
    """Mean per method along the swept parameter, with a 95% t-interval band."""
    if data['sub'].nunique() > 1:
        raise ValueError("line panels take one sub-experiment, found {}".format(sorted(data['sub'].unique())))
    xs = sorted(data['x'].unique())
    ends = []
    for method in methods:
        stats = data[data['method'] == method].groupby('x')['value'].agg(['mean', 'sem', 'count'])
        ci = (stats['sem'] * scipy.stats.t.ppf(0.975, (stats['count'] - 1).clip(lower=1))
              ).where(stats['count'] > 1, 0)
        ax.plot(stats.index, stats['mean'], marker='o', markersize=6, color=METHOD_COLORS[method], label=method)
        ax.fill_between(stats.index, stats['mean'] - ci, stats['mean'] + ci, color=METHOD_COLORS[method], alpha=0.2)
        ends.extend(zip(stats.index[[0, -1]], stats['mean'].iloc[[0, -1]]))
    if xs[0] > 0 and xs[-1] / xs[0] >= log_x_ratio:
        ax.set_xscale('log')
        ax.minorticks_off()
    ax.set_xticks(xs)
    ax.set_xticklabels(["{:g}".format(x) for x in xs])
    ax.set_xlabel(x_label)
    ax.grid(True, axis='y')
    if log_y:
        _log_value_axis(ax)
    _label_line_ends(ax, set(ends))


def panels_figure(df, panels, title, x_label=None, panel_width=6.5, height=5.5):
    """One figure with a panel per metric in `panels` that has data: bars, or lines along the swept parameter
    when the experiment has one (x_label names it). One legend for all panels, outside the top-right corner of
    the last one. None when no panel has data."""
    panels = [m for m in panels if (df['metric'] == m).any()]
    if not panels:
        return None
    swept = df['x'].notna().any()
    fig, axes = plt.subplots(1, len(panels), figsize=(panel_width * len(panels) + 2.5, height),
                             constrained_layout=True, squeeze=False)
    handles = {}
    for ax, metric in zip(axes[0], panels):
        data = df[df['metric'] == metric]
        methods = [m for m in METHODS if (data['method'] == m).any()]
        panel_title, ylabel, log_y = METRICS[metric]
        if swept:
            _draw_lines(ax, data, methods, log_y, x_label)
        else:
            _draw_bars(ax, data, methods, log_y)
        ax.set_title(panel_title)
        ax.set_ylabel(ylabel)
        ax.set_axisbelow(True)
        ax.grid(True, axis='y')
        for handle, label in zip(*ax.get_legend_handles_labels()):
            handles.setdefault(label, handle)
        if ax.get_legend() is not None:
            ax.get_legend().remove()
    labels = [m for m in METHODS if m in handles]
    axes[0][-1].legend([handles[m] for m in labels], [METHOD_DISPLAY_NAMES[m] for m in labels],
                       title=LEGEND_TITLE, loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)
    fig.suptitle(title)
    return fig


def experiment_figures(df, title, x_label=None, show_varying=None, varying_gap_threshold=0.1):
    """(accuracy figure, run-time figure, failed-inference figure) of one experiment. The hit-rate on varying nodes is drawn when
    show_varying says so, or, when it is None, when it differs from the hit-rate on all nodes by at least
    varying_gap_threshold for some method (see varying_hit_rate_differs)."""
    if show_varying is None:
        show_varying, gap = varying_hit_rate_differs(df, varying_gap_threshold)
        print("  largest all-nodes vs varying-nodes hit-rate gap: {:.3f} -> varying-nodes panel {}".format(
            gap, "drawn" if show_varying else "left out"))
    accuracy_panels = [m for m in ACCURACY_PANELS if show_varying or m != 'hit_rate_varying']
    return (panels_figure(df, accuracy_panels, title, x_label),
            panels_figure(df, TIME_PANELS, title, x_label),
            panels_figure(df, FAILED_PANELS, title, x_label))


def save_figure(fig, figures_dir, name, formats=("png", "pdf")):
    if fig is None:
        return
    os.makedirs(figures_dir, exist_ok=True)
    for fmt in formats:
        fig.savefig(os.path.join(figures_dir, "{}.{}".format(name, fmt)), dpi=200, bbox_inches="tight")
