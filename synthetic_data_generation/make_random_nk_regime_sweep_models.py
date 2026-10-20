"""Build data/random_nk_regime_sweep_models: classic Kauffman NK random Boolean networks, in three regimes.

Every node has exactly K inputs drawn uniformly without replacement from the other nodes, and a truth table
whose 2**K rows are independent coin flips with P(1) = BIAS. That is the original NK ensemble, and it
differs from data/random_scale_free_models in both respects: the degree is fixed rather than heavy-tailed,
and the functions are general Boolean rather than symmetric threshold. The second difference matters for
what these are for - the scale-free models are drawn from exactly the class the symmetric /
symmetric_topology methods search, so they cannot show what model misspecification costs.

The regimes come from the annealed approximation: a network is critical when 2*p*(1-p)*K = 1, ordered below
and chaotic above. At the fixed bias p = 0.5 that reduces to K/2, so K = 1, 2, 3 give 0.5, 1.0 and 1.5 -
ordered, critical and chaotic. Equivalently, a one-bit perturbation is expected to flip K/2 nodes in the
next step, which is what `sensitivity` measures here and what the README reports per model, as a check that
the ensemble is in the regime it claims rather than an assumption.

Bias is a constant rather than a swept parameter, which confounds the regime with the in-degree - see BIAS
below. data/random_nk_critical_models is the set that separates them: it moves the bias with K so that the
whole sweep sits at criticality, leaving K as the only thing that changes.

Layout is `K=<k>/size=<n>/nk_k<k>_n<n>_<index>`, one directory level per parameter, because data generation
and inference both mirror a graphs_dir's grouping into their output: K and size then arrive in the analysis
notebook as swept numeric parameters that a figure can put on an axis, rather than as opaque grouping.

Regenerate with `python synthetic_data_generation/make_random_nk_regime_sweep_models.py`.
"""
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from attractor_learning.graphs import Network
from synthetic_data_generation.nk_ensemble import (attractor_summary, check_round_trip, nk_network,
                                                   sensitivity)

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
OUTPUT_DIR = os.path.join(DATA_DIR, "random_nk_regime_sweep_models")

SIZES = [5, 10, 50, 250]
# K -> the regime it puts the ensemble in at BIAS = 0.5 (2*p*(1-p)*K, i.e. K/2 here)
REGIMES = {1: "ordered", 2: "critical", 3: "chaotic"}
MODELS_PER_K_AND_SIZE = 10

# P(a truth-table row outputs 1). Held at 0.5 because that is where the K values above give exactly the
# three regimes, and because it keeps K the only thing varying between them.
#
# Worth sweeping if you want to separate the regime from the in-degree. As it stands the two are
# confounded: K also sets the truth-table size (2**K rows), the size of REVEAL's and Best-Fit's subset
# search, and the ILP's per-node input count, so a difference between the K = 1 and K = 3 sets is not
# attributable to the dynamics alone. Since criticality only depends on 2*p*(1-p)*K, a biased ensemble at
# fixed K reaches a different regime with the same in-degree - K = 3 is critical at p = 0.211 or 0.789, and
# K = 2 is ordered at any p outside (0.146, 0.854). A pair of sets sharing K and differing in regime is the
# comparison that isolates the dynamics; data/random_nk_critical_models is the other half of it. Note also
# that a strong bias produces mostly-constant nodes and so a large frozen core, which is the near-static
# regime where all_constants is already a strong baseline.
BIAS = 0.5


def write_model_readme(path, name, k, n, bias, seed, lam, attractors):
    lines = ["# {}".format(name), "",
             "Kauffman NK network: every node has exactly {} input{}, chosen uniformly without replacement "
             "from the other nodes, and a truth table whose {} rows are independent coin flips with "
             "P(1) = {}.".format(k, "" if k == 1 else "s", 2 ** k, bias), "",
             "| property | value |", "|---|---|",
             "| nodes | {} |".format(n),
             "| edges | {} |".format(k * n),
             "| K | {} |".format(k),
             "| bias | {} |".format(bias),
             "| regime | {} |".format(REGIMES[k]),
             "| sensitivity (Derrida slope at distance 1) | {:.2f} |".format(lam),
             "| expected sensitivity, 2p(1-p)K | {:.2f} |".format(2 * bias * (1 - bias) * k),
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


def main():
    if os.path.isdir(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR)

    index = ["# Random NK models, regime sweep", "",
             "A `graphs_dir` of {} Kauffman NK networks: {} at each combination of K in {} and size in {}. "
             "Every node has exactly K inputs drawn uniformly without replacement from the other nodes, and "
             "a truth table whose 2**K rows are independent coin flips with P(1) = {}.".format(
                 len(REGIMES) * len(SIZES) * MODELS_PER_K_AND_SIZE, MODELS_PER_K_AND_SIZE,
                 ", ".join(str(k) for k in sorted(REGIMES)), ", ".join(str(s) for s in SIZES), BIAS), "",
             "The regimes follow the annealed approximation, critical at 2*p*(1-p)*K = 1: at p = {} that is "
             "K/2, so K = 1, 2, 3 are ordered, critical and chaotic. The sensitivity column is the measured "
             "Derrida slope at distance 1 (mean nodes differing one step after a one-bit perturbation), "
             "which should sit near K/2.".format(BIAS), "",
             "Unlike `data/random_scale_free_models`, the functions here are general Boolean rather than "
             "symmetric threshold, so these models are outside the class the symmetric / "
             "symmetric_topology methods search.", "",
             "Grouped `K=<k>/size=<n>`, one directory level per parameter, so both reach the analysis "
             "notebook as swept numeric parameters.", "",
             "`data/random_nk_critical_models` is the companion set: it moves the bias with K so the whole "
             "sweep stays critical, which separates the regime from the in-degree that is confounded with "
             "it here.", "",
             "Regenerate with `python synthetic_data_generation/make_random_nk_regime_sweep_models.py`.", ""]

    for k in sorted(REGIMES):
        index += ["## K = {} ({})".format(k, REGIMES[k]), "",
                  "| model | nodes | edges | sensitivity | attractors | longest period |",
                  "|---|---|---|---|---|---|"]
        for n in SIZES:
            measured = []
            for model_index in range(MODELS_PER_K_AND_SIZE):
                seed = 1000000 * k + 1000 * n + model_index
                name = "nk_k{}_n{:04d}_{:02d}".format(k, n, model_index)
                network = nk_network(n, k, BIAS, seed)
                lam = sensitivity(network, seed)
                attractors = attractor_summary(network)
                measured.append(lam)

                model_dir = os.path.join(OUTPUT_DIR, "K={}".format(k), "size={}".format(n), name)
                os.makedirs(model_dir)
                network.save(os.path.join(model_dir, Network.MODEL_JSON_NAME))
                write_model_readme(os.path.join(model_dir, "README.md"), name, k, n, BIAS, seed,
                                   lam, attractors)

                # the folder is only useful if the pipeline's readers can load it back unchanged
                check_round_trip(network, model_dir, seed)

                index.append("| [{}]({}) | {} | {} | {:.2f} | {} | {} |".format(
                    name, "K={}/size={}/{}/README.md".format(k, n, name), n, k * n, lam,
                    attractors[0] if attractors is not None else "not enumerated",
                    attractors[1] if attractors is not None else "-"))
            print("K={} size={:>4}: {} models, mean sensitivity {:.2f} (expected {:.2f})".format(
                k, n, MODELS_PER_K_AND_SIZE, sum(measured) / len(measured), 2 * BIAS * (1 - BIAS) * k))
        index.append("")

    with open(os.path.join(OUTPUT_DIR, "README.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(index) + "\n")
    print("\nwrote {} models to {}".format(
        len(REGIMES) * len(SIZES) * MODELS_PER_K_AND_SIZE, OUTPUT_DIR))


if __name__ == "__main__":
    main()
