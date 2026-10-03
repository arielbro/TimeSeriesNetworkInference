"""Build the critical NK graphs_dirs: Kauffman NK networks held at criticality across a sweep of K.

The companion to data/random_nk_regime_sweep_models, and the set that separates the two things confounded
there. That one holds the bias at 0.5 and varies K, so K moves the regime (2*p*(1-p)*K = K/2) and at the
same time moves the in-degree, the truth-table size, the size of REVEAL's and Best-Fit's subset search and
the ILP's per-node input count - a difference between its K = 1 and K = 3 sets is therefore not
attributable to the dynamics. Here the function parameter moves with K instead, so the regime stays put and
K sweeps the size of the inference problem alone.

Criticality fixes p only up to reflection: 2*p*(1-p)*K = 1 has two roots, p = (1 +/- sqrt(1 - 2/K)) / 2,
and since p(1-p) is invariant under p -> 1-p they are equally critical. In the two fitted families every
node tosses a coin for which of the two it gets. That is not a cosmetic choice. Drawn entirely at the upper
root, a high-K set has nearly every truth table full of ones and its trajectories collapse onto the
all-ones corner - the K = 15 set measured a mean state density of 0.974, with models pinned at 1.000, which
is a near-static signal that `all_constants` already fits perfectly. Mixing the two roots leaves the
sensitivity exactly where it was and puts the density back near 0.5, so the generated time series actually
carry information.

Three function families, selected by --function-type, each into its own directory:

  general    data/random_nk_critical_models
             An i.i.d. truth table per node, at whichever root of 2*p*(1-p)*K = 1 the node's coin chose.
             Outside the class the symmetric / symmetric_topology methods search, so this is the set that
             can show what model misspecification costs.

  threshold  data/random_nk_critical_threshold_models
             Signed symmetric threshold functions: a uniformly random sign per input, and the threshold
             whose average sensitivity is closest to 1. The coin picks between that threshold t and its
             mirror K + 1 - t, which have complementary biases and identical sensitivity. Inside the class
             the symmetric methods search, and stored as signs plus a threshold, so there is no truth
             table and no in-degree ceiling.

  threshold-uniform
             data/random_nk_critical_threshold_uniform_models
             The same, except that nothing is fitted: each node draws its own threshold from
             Uniform[1, K], which is how make_random_scale_free_models draws its functions, so this set
             and that one differ only in topology. Summing K*C(K-1,t-1)/2**(K-1) over t = 1..K gives
             K*2**(K-1)/2**(K-1), so the mean sensitivity over the draw is exactly 1 at every K and the
             ensemble is critical by that identity - unbiased where the fitted set is off by as much as
             0.30, at the cost of scattering around 1 instead of landing on one value. No polarity coin,
             which would be a no-op: a uniform threshold already samples the biases symmetrically about
             0.5, and mirroring t to K + 1 - t returns the same distribution.

Note for the threshold directories: the signs affect neither the bias nor the sensitivity (the count of
inputs agreeing with their sign is Binomial(K, 1/2) whatever the signs are), only which inputs activate and
which inhibit. The ILP inference method does not accept threshold functions at all - see
SymmetricThresholdFunction. And criticality is only approached, not reached: the achievable sensitivities
are discrete, so the sweep runs 1.00 at K = 2 but 1.25 at K = 5, 0.70 at K = 10 and 0.92 at K = 15, and its
distance from criticality is not monotonic in K. Each model's README reports the regime it actually has.

K = 50 (p = (5 + 2*sqrt(6)) / 10) is generated in neither directory. A node cannot draw 50 distinct
regulators from the other 49 nodes, and even allowing self-loops it would make the complete digraph, with no
topology left to infer. Its truth table is also out of reach: the sparse representation stores only the rows
disagreeing with the majority output, which is what makes K = 15 cheap (about 1100 rows of 32768), but at
K = 50 that minority is still about 1e13 rows. The threshold family has no such ceiling - K = 50 there is a
one-line edit to K_VALUES - but it is left out so that the two families sweep the same grid.

Layout is `K=<k>/size=<n>/<name>`, matching the regime sweep, so every set arrives in the analysis notebook
with the same K and size parameter columns even though size is single-valued here.

The grid is the module constants SIZES, K_VALUES and MODELS_PER_K_AND_SIZE, edited per sweep rather than
passed in, and --output-dir names the directory to write it to. A run deletes its output directory before
writing, so a sweep on a grid other than the constants' current one must be given a directory of its own -
without it, the set already sitting in the function type's own directory is replaced by one built on a grid
its README does not describe. The directories named above hold K = 2, 5, 10, 15 at size 50, ten models per
cell; `data/random_nk_critical_threshold_uniform_n15_models` holds the threshold-uniform family over
K = 1, 3, 5, 8, 10 at size 15, five per cell, which is the grid the constants currently carry.

K = 1 is reachable in the threshold families only. There the single available threshold has sensitivity
exactly 1, so a uniform-threshold ensemble is critical exactly rather than in expectation, while an i.i.d.
truth table cannot be critical at K = 1 at any bias (2p(1-p)K <= K/2) - critical_bias raises there, so
--function-type general has nothing to generate at that K.

Regenerate with `python synthetic_data_generation/make_random_nk_critical_models.py [--function-type ...]
[--output-dir ...]`.
"""
import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from attractor_learning.graphs import Network
from synthetic_data_generation.nk_ensemble import (attractor_summary, check_round_trip, critical_bias,
                                                   critical_threshold, exact_threshold_sensitivity,
                                                   high_polarity_count, mean_density, mirrored_threshold,
                                                   nk_network, nk_threshold_network, sensitivity,
                                                   threshold_bias, threshold_sensitivity)

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# The grid. Not arguments: a sweep is a set that gets written once and then read by every experiment
# downstream, so which cells it holds belongs in the file next to the reasons for them rather than in
# whatever command line last ran. Change these together with --output-dir, so the sweep that was here
# before is not overwritten by one built on different cells - the directories in OUTPUT_DIRS hold
# K = 2, 5, 10, 15 at size 50, ten per cell, and their READMEs describe that grid.
SIZES = [15]
# in-degrees to sweep; the function parameter for each is derived, not chosen - see critical_bias and
# critical_threshold. K = 50 was asked for and is not here, for the reasons in the module docstring.
# K = 1 works for the threshold families only, where it is exactly critical, and raises for general.
K_VALUES = [1, 3, 5, 8, 10]
MODELS_PER_K_AND_SIZE = 5

# Offsets so that the critical sweeps and the regime sweep draw independent networks where their grids
# overlap. K = 2, size = 50, p = 0.5 is the same ensemble in both the general critical set and the regime
# sweep, and without an offset the two would hold nearly identical models: fine for a paired comparison,
# but it would silently double-count if a figure ever pooled the sets.
SEED_OFFSET = 500000000
THRESHOLD_SEED_OFFSET = 600000000
UNIFORM_THRESHOLD_SEED_OFFSET = 700000000


OUTPUT_DIRS = {"general": "random_nk_critical_models",
               "threshold": "random_nk_critical_threshold_models",
               "threshold-uniform": "random_nk_critical_threshold_uniform_models"}
SEED_OFFSETS = {"general": SEED_OFFSET,
                "threshold": THRESHOLD_SEED_OFFSET,
                "threshold-uniform": UNIFORM_THRESHOLD_SEED_OFFSET}
NAME_PREFIXES = {"general": "nk_crit", "threshold": "nk_crit_thr",
                 "threshold-uniform": "nk_crit_thru"}


def output_dir(function_type, dir_name=None):
    """Where the set is written: the directory this function type owns, or `dir_name` under data/ when a
    run sweeps a grid other than the default one and so must not overwrite it."""
    return os.path.join(DATA_DIR, dir_name or OUTPUT_DIRS[function_type])


def critical_bias_note(k):
    """critical_bias(k) formatted for a README, or, below k = 2 where it raises rather than returning a
    number, why there is no such bias: 2p(1-p)k <= k/2 < 1 caps an i.i.d. table's sensitivity short of 1
    at every bias. Only the threshold families have a k = 1 cell for this row to appear in."""
    if k < 2:
        return "none, at any bias: 2p(1-p)K <= K/2 < 1 at K = {}".format(k)
    return "{:.6f}".format(critical_bias(k))


def threshold_sensitivity_spread(k):
    """(mean, standard deviation) of one node's sensitivity when its threshold is drawn from Uniform[1, k].
    The mean is exactly 1 at every k - that identity is what makes the threshold-uniform ensemble critical
    without anything being fitted - and the sd is what a network of n such nodes scatters by, over
    sqrt(n). Zero at k = 1, where the draw has a single value to return."""
    values = [threshold_sensitivity(k, t) for t in range(1, k + 1)]
    mean = sum(values) / len(values)
    var = sum((v - mean) ** 2 for v in values) / len(values)
    return mean, var ** 0.5


def regenerate_command(function_type, dir_name):
    """The command line that reproduces this run. The directory is named whenever it is not the one the
    function type owns, since that is the argument without which the run would land somewhere else."""
    command = ("python synthetic_data_generation/make_random_nk_critical_models.py "
               "--function-type {}".format(function_type))
    if dir_name is not None and dir_name != OUTPUT_DIRS[function_type]:
        command += " --output-dir {}".format(dir_name)
    return command


def model_parameters(function_type, k):
    """(bias, threshold, expected sensitivity) for one in-degree, all for the higher-bias side of the coin;
    the lower side is 1 - bias, and for the threshold family K + 1 - threshold, which has the same
    sensitivity. The threshold is None for the general family, which has a truth table instead."""
    if function_type == "general":
        bias = critical_bias(k)
        return bias, None, 2 * bias * (1 - bias) * k
    if function_type == "threshold-uniform":
        # the mean of k*C(k-1,t-1)/2**(k-1) over t = 1..k is exactly 1, at every k, so this ensemble is
        # critical without any threshold being fitted - see nk_threshold_network
        return None, None, 1.0
    threshold = critical_threshold(k)
    return threshold_bias(k, threshold), threshold, threshold_sensitivity(k, threshold)


def build_network(function_type, n, k, bias, threshold, seed):
    if function_type == "general":
        return nk_network(n, k, bias, seed, mixed_polarity=True)
    if function_type == "threshold-uniform":
        # threshold None draws it per node; a polarity coin would be a no-op on a uniform threshold
        return nk_threshold_network(n, k, None, seed)
    return nk_threshold_network(n, k, threshold, seed, mixed_polarity=True)


def write_model_readme(path, name, function_type, k, n, bias, threshold, expected, seed, lam, density,
                       high, exact, attractors):
    if function_type == "threshold-uniform":
        if k == 1:
            realized = ("this network's own value is exactly 1 as well, since Uniform[1, 1] leaves the "
                        "draw nothing to choose and every node is a single signed input")
        else:
            realized = ("this network's own value is {:.4f}, since {} draws scatter around the "
                        "mean".format(exact, n))
        description = ("Kauffman NK network of signed symmetric threshold functions: every node has exactly "
                       "{} inputs, chosen uniformly without replacement from the other nodes, a uniformly "
                       "random sign per input, and its own threshold drawn uniformly from [1, {}]. Averaged "
                       "over that draw a node's sensitivity is exactly 1 at every K, so the ensemble is "
                       "critical without any threshold being fitted; {}. The signs decide which inputs "
                       "activate and which inhibit; the bias and the sensitivity follow from the threshold "
                       "alone.".format(k, k, realized))
    elif function_type == "general":
        description = ("Kauffman NK network at criticality: every node has exactly {} inputs, chosen "
                       "uniformly without replacement from the other nodes, and a truth table whose {} rows "
                       "are independent coin flips. The bias is one of the two roots of 2p(1-p)K = 1, "
                       "{:.6f} or {:.6f}, chosen by a coin toss per node - both are equally critical, since "
                       "p(1-p) does not change under p -> 1-p.".format(k, 2 ** k, bias, 1 - bias))
    else:
        description = ("Kauffman NK network of signed symmetric threshold functions: every node has exactly "
                       "{} inputs, chosen uniformly without replacement from the other nodes, a uniformly "
                       "random sign per input, and threshold {} or {} by a coin toss per node. The two "
                       "thresholds have complementary biases ({:.6f} and {:.6f}) and identical sensitivity "
                       "({:.4f}). The signs decide which inputs activate and which inhibit; the bias and "
                       "the sensitivity follow from the threshold alone.".format(
                           k, threshold, mirrored_threshold(k, threshold), bias, 1 - bias, expected))
    lines = ["# {}".format(name), "", description, "",
             "| property | value |", "|---|---|",
             "| nodes | {} |".format(n),
             "| edges | {} |".format(k * n),
             "| K | {} |".format(k)]
    if function_type == "threshold-uniform":
        lines += ["| threshold, per node | uniform in [1, {}] |".format(k),
                  "| bias, per node | one of {} |".format(
                      ", ".join("{:.4f}".format(threshold_bias(k, t)) for t in range(1, k + 1))),
                  "| bias needed for criticality of an i.i.d. table | {} |".format(
                      critical_bias_note(k))]
    else:
        lines.append("| bias, per node | {:.6f} or {:.6f} |".format(bias, 1 - bias))
        if threshold is not None:
            lines += ["| threshold, per node | {} or {} |".format(threshold,
                                                                  mirrored_threshold(k, threshold)),
                      "| bias needed for criticality of an i.i.d. table | {} |".format(
                          critical_bias_note(k))]
    lines.append("| nodes on the higher-bias side | {} of {} |".format(high, n))
    if exact is not None:
        lines.append("| exact sensitivity of this network | {:.4f} |".format(exact))
    if function_type == "threshold-uniform":
        # the ensemble is critical in expectation, but this network is one draw of n thresholds - report
        # what it actually came out as rather than the mean it was drawn around
        if k == 1:
            regime = "critical exactly (t = 1 is the only threshold there is, and its sensitivity is 1)"
        else:
            regime = "critical in expectation (this network: {:.2f})".format(exact)
    elif function_type == "general":
        # 2p(1-p)K = 1 holds exactly at either root, so the claim needs no hedging
        regime = "critical (by construction)"
    elif expected < 0.95:
        # the threshold family's sensitivities are discrete, so most K cannot land on 1 - say where it did
        # land rather than calling every set critical because it was meant to be
        regime = "ordered (expected sensitivity {:.2f}, below 1)".format(expected)
    elif expected > 1.05:
        regime = "chaotic (expected sensitivity {:.2f}, above 1)".format(expected)
    else:
        regime = "critical (expected sensitivity {:.2f})".format(expected)
    lines += ["| regime | {} |".format(regime),
              "| sensitivity (Derrida slope at distance 1) | {:.2f} |".format(lam),
              "| expected sensitivity | {:.4f} |".format(expected),
              "| mean state density after burn-in | {:.3f} |".format(density),
              "| seed | {} |".format(seed)]
    if attractors is not None:
        lines += ["| attractors (exact) | {} |".format(attractors[0]),
                  "| longest period | {} |".format(attractors[1])]
    else:
        lines.append("| attractors | not enumerated (2**{} states) |".format(n))
    lines += ["", "Stored as `network.json`. `Network.parse_model_dir` reads it, so this folder works "
                  "anywhere a graphs_dir is expected.", ""]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def index_header(function_type, command):
    """The prose and parameter table at the top of the graphs_dir's README. `command` is what regenerates
    it, which depends on where it is being written and so cannot be spelled out here."""
    total = len(K_VALUES) * len(SIZES) * MODELS_PER_K_AND_SIZE
    ks = ", ".join(str(k) for k in K_VALUES)
    sizes = ", ".join(str(s) for s in SIZES)
    sizes_text = " or ".join(str(s) for s in SIZES)
    # How far the uniform draw scatters, per K, and where on this grid the fitted-threshold rule would
    # miss criticality by most. Both were sentences of fixed numbers written for K = 2, 5, 10, 15 at size
    # 50, and neither survives a change of the grid, so they are computed from it instead. The second is a
    # statement about the rule on these K, not about data/random_nk_critical_threshold_models, which holds
    # its own grid and whose numbers cannot be derived from this one.
    spreads = {k: threshold_sensitivity_spread(k)[1] for k in K_VALUES}
    narrowest, widest = min(spreads, key=spreads.get), max(spreads, key=spreads.get)
    worst_fitted_k = max(K_VALUES, key=lambda k: abs(threshold_sensitivity(k, critical_threshold(k)) - 1))
    worst_fitted = threshold_sensitivity(worst_fitted_k, critical_threshold(worst_fitted_k))
    coin = ("Criticality fixes p only up to reflection: 2p(1-p)K = 1 has two roots, "
            "p = (1 +/- sqrt(1 - 2/K)) / 2, equally critical because p(1-p) is unchanged by p -> 1-p. "
            "Every node tosses a coin for which of the two it gets. Drawing every node at the upper root "
            "instead makes a high-K set nearly all ones and collapses its trajectories onto the all-ones "
            "corner (the K = 15 set measured density 0.974 that way, some models pinned at 1.000); mixing "
            "the roots leaves the sensitivity exactly where it was and returns the density to about 0.5. "
            "The realized split is Binomial(n, 1/2), so individual models sit well off even - each model's "
            "README reports its own, and that is what explains its density.")
    if function_type == "threshold-uniform":
        lines = ["# Random NK models, critical sweep, threshold drawn per node", "",
                 "A `graphs_dir` of {} Kauffman NK networks: {} at each combination of K in {} and size in "
                 "{}. Every node has exactly K inputs drawn uniformly without replacement from the other "
                 "nodes, a uniformly random sign per input, and its own threshold drawn uniformly from "
                 "[1, K].".format(total, MODELS_PER_K_AND_SIZE, ks, sizes), "",
                 "This is how `data/random_scale_free_models` draws its functions, on a fixed in-degree "
                 "instead of a power-law one, so the two differ only in topology. Nothing is fitted here: "
                 "summing the sensitivity K*C(K-1,t-1)/2**(K-1) over t = 1..K gives K*2**(K-1)/2**(K-1), so "
                 "its mean over a uniform threshold is exactly 1 at every K, and the ensemble is critical "
                 "by that identity alone.", "",
                 "The trade against `data/random_nk_critical_threshold_models`, which gives a network "
                 "one fitted threshold, is bias for variance. Fitting misses at most K, and would miss on "
                 "this grid by as much as {:.2f} at K = {}, but identically so for every network; this "
                 "one centres on 1 and scatters instead, since a network is {} draws from a per-node "
                 "sensitivity whose standard deviation runs from {:.2f} at K = {} to {:.2f} at K = {}. "
                 "Each model's README carries its own exact value, and the table below the spread it is "
                 "drawn from.".format(
                     worst_fitted, worst_fitted_k, sizes_text, spreads[narrowest], narrowest,
                     spreads[widest], widest), "",
                 "A polarity coin would do nothing here: the achievable biases are symmetric about 0.5 and "
                 "a uniform threshold already samples them symmetrically, so mirroring t to K + 1 - t "
                 "returns the same distribution. The density comes out near 0.5 without it."]
        if 1 in K_VALUES:
            lines += ["",
                      "K = 1 is the one cell where this family is critical exactly rather than in "
                      "expectation: Uniform[1, 1] has nothing to draw, every node is a single signed "
                      "input of sensitivity 1, and the spread is zero. No i.i.d. truth table is critical "
                      "there at any bias, since 2p(1-p)K <= K/2, so `--function-type general` has no "
                      "K = 1 cell to put beside it."]
        header = "| K | thresholds | sensitivity per node, min .. max | mean | sd |"
        rule = "|---|---|---|---|---|"
        for n in SIZES:
            header += " sd of a {}-node net |".format(n)
            rule += "---|"
        lines += ["", header, rule]
        for k in K_VALUES:
            values = [threshold_sensitivity(k, t) for t in range(1, k + 1)]
            mean = sum(values) / len(values)
            row = "| {} | 1 .. {} | {:.4f} .. {:.4f} | {:.4f} | {:.3f} |".format(
                k, k, min(values), max(values), mean, spreads[k])
            for n in SIZES:
                row += " {:.3f} |".format(spreads[k] / n ** 0.5)
            lines.append(row)
        lines += ["",
                  "The signs affect neither the bias nor the sensitivity - the count of inputs agreeing "
                  "with their sign is Binomial(K, 1/2) however the signs fall. The ILP inference method "
                  "does not accept threshold functions."]
    elif function_type == "general":
        lines = ["# Random NK models, critical sweep", "",
                 "A `graphs_dir` of {} Kauffman NK networks: {} at each combination of K in {} and size in "
                 "{}. Every node has exactly K inputs drawn uniformly without replacement from the other "
                 "nodes, and a truth table whose 2**K rows are independent coin flips.".format(
                     total, MODELS_PER_K_AND_SIZE, ks, sizes), "",
                 "Unlike `data/random_nk_regime_sweep_models`, the bias is not held fixed: it is derived "
                 "from K as a root of the criticality condition 2p(1-p)K = 1. Every model here is "
                 "therefore at the same critical point, and K sweeps the size of the inference problem "
                 "without moving the dynamical regime - the comparison the regime sweep cannot make, since "
                 "there K sets both at once.", "", coin, "",
                 "| K | bias p | 1 - p | 2p(1-p)K | truth-table rows | stored as |",
                 "|---|---|---|---|---|---|"]
        for k in K_VALUES:
            bias = critical_bias(k)
            lines.append("| {} | {:.6f} | {:.6f} | {:.3f} | {} | {} |".format(
                k, bias, 1 - bias, 2 * bias * (1 - bias) * k, 2 ** k,
                "full table" if 2 ** k <= 1024 else "minority rows only"))
        lines += ["",
                  "At K = 2 the two roots coincide at 0.5, so the coin has nothing to choose there and the "
                  "split reported per model is only the tables' own fluctuation.", "",
                  "Above 1024 rows a node's function is stored as the rows that disagree with its majority "
                  "output rather than as a full table, and is sampled the same way - the count from the "
                  "binomial, then that many distinct rows - so the cost is the disagreeing rows (about 1100 "
                  "of 32768 at K = 15) rather than 2**K. That is what puts K = 15 within reach: dense, it "
                  "is minutes and megabytes per model; sparse, it is a twentieth of a second and 0.4 MB.",
                  "",
                  "K = 50 (p = (5 + 2*sqrt(6)) / 10) is absent by necessity. A node cannot draw 50 distinct "
                  "regulators out of the 49 other nodes, and the sparse encoding does not rescue its truth "
                  "table either: the minority is still about 1e13 rows. See the generating script.", "",
                  "`data/random_nk_critical_threshold_models` is the same sweep drawn from signed symmetric "
                  "threshold functions instead. Those lie inside the class the symmetric and "
                  "symmetric_topology methods search, where these lie outside it, so the pair is what "
                  "separates the difficulty of a critical network from the cost of searching the wrong "
                  "function class."]
    else:
        lines = ["# Random NK models, critical sweep, signed symmetric threshold functions", "",
                 "A `graphs_dir` of {} Kauffman NK networks: {} at each combination of K in {} and size in "
                 "{}. Every node has exactly K inputs drawn uniformly without replacement from the other "
                 "nodes, a uniformly random sign per input, and a threshold.".format(
                     total, MODELS_PER_K_AND_SIZE, ks, sizes), "",
                 "Unlike `data/random_nk_critical_models`, these functions are symmetric threshold rather "
                 "than general Boolean, so they lie inside the class the symmetric and symmetric_topology "
                 "inference methods search - the two sets together are what separates the difficulty of a "
                 "critical network from the cost of searching the wrong function class. They are stored as "
                 "signs plus a threshold, so a model is a few hundred bytes per node at any K, with no "
                 "truth table and no in-degree ceiling. The ILP method does not accept threshold "
                 "functions.", "",
                 "The threshold is the one whose average sensitivity, K*C(K-1,t-1)/2**(K-1), comes closest "
                 "to 1, which is criticality itself rather than the bias that stands in for it when the "
                 "truth table is i.i.d. Each node then tosses a coin between that threshold t and its "
                 "mirror K + 1 - t, which has the complementary bias and exactly the same sensitivity, so "
                 "the coin is the threshold family's version of choosing between the two critical "
                 "roots.", "",
                 "The signs affect neither the bias nor the sensitivity - the count of inputs agreeing with "
                 "their sign is Binomial(K, 1/2) however the signs fall - so both are properties of the "
                 "threshold, and the achievable biases are the discrete upper tails of Binomial(K, 1/2).",
                 "",
                 "That discreteness is the caveat to read before using these. The achievable sensitivities "
                 "are K*C(K-1,t-1)/2**(K-1) and nothing between, so only K = 2 lands on 1 exactly: the "
                 "sweep runs 1.00, 1.25, 0.70 and 0.92 at K = 2, 5, 10 and 15, and its distance from "
                 "criticality is not monotonic in K. Each model's README reports the regime it actually "
                 "has rather than the one it was aimed at. Fitting the bias to a root of 2p(1-p)K = 1 "
                 "instead does worse, since the families disagree at equal bias - at these thresholds a "
                 "threshold function carries about half to two thirds of the sensitivity an i.i.d. table "
                 "of the same bias would - and at K = 15 it would pick t = 4, sensitivity 0.33.",
                 "", coin, "",
                 "| K | threshold | mirror | bias | 1 - bias | i.i.d. critical bias | sensitivity |",
                 "|---|---|---|---|---|---|---|"]
        for k in K_VALUES:
            t = critical_threshold(k)
            lines.append("| {} | {} | {} | {:.6f} | {:.6f} | {:.6f} | {:.4f} |".format(
                k, t, mirrored_threshold(k, t), threshold_bias(k, t), 1 - threshold_bias(k, t),
                critical_bias(k), threshold_sensitivity(k, t)))
    lines += ["",
              "The sensitivity column below is the measured Derrida slope at distance 1 (mean nodes "
              "differing one step after a one-bit perturbation). It is a mean of 200 draws, and on the "
              "n = 50 sets the per-model standard deviation of that estimator ran from 0.06 to 0.11 - "
              "more at smaller n, since the slope is an average over the network's own nodes - so an "
              "individual value a couple of tenths from the expected one is sampling noise rather than a "
              "model that missed its regime; judge a cell by its {} models, not by one. The density "
              "column is the mean fraction of nodes holding the value 1 along a trajectory, after 50 "
              "burn-in steps from a random state are run and discarded; the `high` column is how many of "
              "the model's nodes came out on the higher-bias side of the coin, which is what a density "
              "away from 0.5 tracks.".format(MODELS_PER_K_AND_SIZE)]
    if function_type == "threshold-uniform":
        lines += ["",
                  "There is no coin in this family, so `high` counts the nodes whose drawn threshold puts "
                  "their bias above 0.5, which is those with 2t < K + 1. At K = 1 that is none of them: "
                  "t = 1 gives bias exactly 0.5, so the column reads 0 and the density still comes out "
                  "there."]
    lines += ["",
              "Grouped `K=<k>/size=<n>`, one directory level per parameter, matching the regime sweep so "
              "every set reaches the analysis notebook with the same parameter columns.", "",
              "Regenerate with `{}`. The grid - K, size and models per cell - lives in the generating "
              "script's constants, so regenerating this set means restoring those to the values in the "
              "table above.".format(command), ""]
    return lines


def main(function_type, dir_name=None):
    out_dir = output_dir(function_type, dir_name)
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    seed_offset = SEED_OFFSETS[function_type]
    prefix = NAME_PREFIXES[function_type]
    index = index_header(function_type, regenerate_command(function_type, dir_name))

    for k in K_VALUES:
        bias, threshold, expected = model_parameters(function_type, k)
        if function_type == "threshold-uniform":
            heading = "## K = {} (threshold uniform in [1, {}])".format(k, k)
        elif threshold is None:
            heading = "## K = {} (p = {:.6f} or {:.6f})".format(k, bias, 1 - bias)
        else:
            heading = "## K = {} (p = {:.6f} or {:.6f}, threshold = {} or {})".format(
                k, bias, 1 - bias, threshold, mirrored_threshold(k, threshold))
        index += [heading, "",
                  "| model | nodes | edges | high | sensitivity | density | attractors | longest period |",
                  "|---|---|---|---|---|---|---|---|"]
        for n in SIZES:
            measured, densities = [], []
            for model_index in range(MODELS_PER_K_AND_SIZE):
                seed = seed_offset + 1000000 * k + 1000 * n + model_index
                name = "{}_k{}_n{:04d}_{:02d}".format(prefix, k, n, model_index)
                network = build_network(function_type, n, k, bias, threshold, seed)
                lam = sensitivity(network, seed)
                density = mean_density(network, seed)
                high = high_polarity_count(network)
                exact = exact_threshold_sensitivity(network) if function_type != "general" else None
                attractors = attractor_summary(network)
                measured.append(lam)
                densities.append(density)

                model_dir = os.path.join(out_dir, "K={}".format(k), "size={}".format(n), name)
                os.makedirs(model_dir)
                network.save(os.path.join(model_dir, Network.MODEL_JSON_NAME))
                write_model_readme(os.path.join(model_dir, "README.md"), name, function_type, k, n, bias,
                                   threshold, expected, seed, lam, density, high, exact, attractors)

                # the folder is only useful if the pipeline's readers can load it back unchanged
                check_round_trip(network, model_dir, seed)

                index.append("| [{}]({}) | {} | {} | {}/{} | {:.2f} | {:.3f} | {} | {} |".format(
                    name, "K={}/size={}/{}/README.md".format(k, n, name), n, k * n, high, n, lam, density,
                    attractors[0] if attractors is not None else "not enumerated",
                    attractors[1] if attractors is not None else "-"))
            if bias is None:
                described = "t~U[1,{}]".format(k)
            else:
                described = "p={:.6f}{}".format(
                    bias, "" if threshold is None else " t={}".format(threshold))
            print("K={:>2} size={:>4} {}: {} models, mean sensitivity {:.2f} (expected {:.2f}), "
                  "mean density {:.3f}".format(
                      k, n, described, MODELS_PER_K_AND_SIZE, sum(measured) / len(measured), expected,
                      sum(densities) / len(densities)))
        index.append("")

    with open(os.path.join(out_dir, "README.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(index) + "\n")
    print("\nwrote {} models to {}".format(len(K_VALUES) * len(SIZES) * MODELS_PER_K_AND_SIZE, out_dir))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--function-type", choices=["general", "threshold", "threshold-uniform"],
                        default="general",
                        help="general Boolean truth tables (default), or signed symmetric threshold "
                             "functions - the class the symmetric inference methods search - with the "
                             "threshold either fitted per K and mirrored by a coin (threshold), or drawn "
                             "per node from Uniform[1, K] (threshold-uniform)")
    parser.add_argument("--output-dir", default=None, metavar="NAME",
                        help="directory name under data/ to write the set to, deleting it first if it "
                             "exists (default: the directory this function type owns). Pass it whenever "
                             "SIZES, K_VALUES or MODELS_PER_K_AND_SIZE differ from the grid the default "
                             "directory already holds, so that set is not replaced by one its README does "
                             "not describe")
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    main(arguments.function_type, arguments.output_dir)
