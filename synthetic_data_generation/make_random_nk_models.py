"""Build data/random_nk_models: classic Kauffman NK random Boolean networks, in three dynamical regimes.

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

Bias is a constant rather than a swept parameter, deliberately - see BIAS below for when to change that.

Layout is `K=<k>/size=<n>/nk_k<k>_n<n>_<index>`, one directory level per parameter, because data generation
and inference both mirror a graphs_dir's grouping into their output: K and size then arrive in the analysis
notebook as swept numeric parameters that a figure can put on an axis, rather than as opaque grouping.

Regenerate with `python synthetic_data_generation/make_random_nk_models.py`.
"""
import json
import os
import random
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from attractor_learning.graphs import Network
from attractor_learning.logic import BooleanSymbolicFunc

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
OUTPUT_DIR = os.path.join(DATA_DIR, "random_nk_models")

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
# comparison that isolates the dynamics. Note also that a strong bias produces mostly-constant nodes and so
# a large frozen core, which is the near-static regime where all_constants is already a strong baseline.
BIAS = 0.5

ENUMERATE_UP_TO = 10        # exact attractor counts only where 2**n is small enough to walk
SENSITIVITY_SAMPLES = 200   # one-bit perturbations averaged for the Derrida estimate


def nk_network(n, k, bias, seed):
    """One NK network: exactly k inputs per node, uniform without replacement from the other nodes (no
    self-loops, matching the scale-free generator), and an i.i.d. bias-`bias` truth table per node."""
    rng = random.Random(seed)
    names = ["v{}".format(i) for i in range(n)]
    edges, functions = [], []
    for target in range(n):
        regulators = rng.sample([i for i in range(n) if i != target], k)
        regulators.sort()                       # ascending, so the first input is the truth table's MSB
        edges.extend((names[source], names[target]) for source in regulators)
        functions.append(BooleanSymbolicFunc(
            input_names=[names[source] for source in regulators],
            boolean_outputs=[rng.random() < bias for _ in range(2 ** k)]))
    return Network(vertex_names=names, edges=edges, vertex_functions=functions)


def sensitivity(network, seed, samples=SENSITIVITY_SAMPLES):
    """Average number of nodes that differ one step after a single-bit perturbation - the Derrida slope at
    distance 1. Below 1 is ordered, 1 is critical, above 1 is chaotic; the annealed approximation puts it
    at 2*p*(1-p)*K."""
    rng = random.Random(seed)
    n = len(network)
    total = 0
    for _ in range(samples):
        state = [rng.randint(0, 1) for _ in range(n)]
        perturbed = list(state)
        flipped = rng.randrange(n)
        perturbed[flipped] = 1 - perturbed[flipped]
        after = [int(bool(v)) for v in network.next_state(state)]
        after_perturbed = [int(bool(v)) for v in network.next_state(perturbed)]
        total += sum(a != b for a, b in zip(after, after_perturbed))
    return total / float(samples)


def attractor_summary(network):
    """(number of attractors, longest period) by enumerating the whole state space, or None when too big."""
    n = len(network)
    if n > ENUMERATE_UP_TO:
        return None
    successor = {}
    for bits in range(2 ** n):
        state = tuple((bits >> (n - 1 - i)) & 1 for i in range(n))
        successor[state] = tuple(int(bool(v)) for v in network.next_state(list(state)))
    attractors, seen = [], {}
    for state in successor:
        path, current = [], state
        while current not in seen and current not in path:
            path.append(current)
            current = successor[current]
        if current in path:                     # closed a new cycle
            cycle = path[path.index(current):]
            attractors.append(cycle)
            for member in cycle:
                seen[member] = len(attractors) - 1
        for member in path:
            seen.setdefault(member, seen.get(current, None))
    return len(attractors), max(len(cycle) for cycle in attractors)


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

    index = ["# Random NK models", "",
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
             "Regenerate with `python synthetic_data_generation/make_random_nk_models.py`.", ""]

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
                assert Network.model_dir_size(model_dir) == n
                reloaded = Network.parse_model_dir(model_dir)
                check = random.Random(seed)
                for _ in range(20):
                    state = [check.randint(0, 1) for _ in range(n)]
                    assert [int(bool(v)) for v in reloaded.next_state(state)] == \
                           [int(bool(v)) for v in network.next_state(state)], \
                        "{}: json round trip changed dynamics".format(name)

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
