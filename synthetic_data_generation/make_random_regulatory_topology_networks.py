"""Build data/random_regulatory_topology_networks: 30 random models at each of the sizes 25, 50, 100, 200 and 400,
with degree distributions fitted to RegulonDB, usable as a `graphs_dir`.

Degree distributions, fitted at run time to data/regulondb/network.tsv (self-loops dropped, duplicate pairs
merged: 2,809 genes, 7,260 edges):
  In-degree:  P(0) is RegulonDB's share of genes with no regulator; for k >= 1, k - 1 is negative binomial
              (maximum likelihood; it ties with a power law with exponential cutoff by AIC, and a pure power
              law is far worse).
  Out-degree: given that a node regulates anything, P(k) for 1 <= k < OUT_XMIN is RegulonDB's empirical
              frequency among regulators, and k >= OUT_XMIN follows a discrete power law (maximum likelihood
              exponent, ~1.89) carrying RegulonDB's share of regulators in that range. OUT_XMIN = 8 is the
              cutoff the Clauset-Shalizi-Newman KS criterion picks on this network.
  Cap:        a node in an n-node model can have at most n - 1 in- or out-edges (no self-loops). Both
              distributions are cut to 0..n-1 and renormalized, so the legal degrees keep their relative
              probabilities and no extra mass lands on the cap.
  Regulators: the share of nodes with out-degree >= 1 is E[in-degree] / E[out-degree | regulator], both under
              the cap, so that incoming and outgoing edges balance in expectation. A gene's number of
              regulators does not depend on the model's size, but the cap shrinks a regulator's mean
              out-degree, so the share rises at small sizes (0.40 at n=25, 0.19 at n=200; at RegulonDB's own
              size it gives 0.098, against the observed 0.11). Keeping RegulonDB's 11% at n=25 would leave
              ~18 outgoing edges for ~65 incoming ones.

For each model:
  1. Draw n in-degrees and n out-degrees independently from the above, and redraw both lists together until
     their sums are equal, the sum is positive, and the pair is digraphical (has a simple realization without
     self-loops). Redrawing only one list does not converge reliably: the out-degree sum is heavy-tailed, and
     after a large hub the expected number of in-degree redraws reaches 1e11 - 1e46.
  2. Wire with networkx's directed_configuration_model (a uniformly random stub matching).
  3. Remove self-loops and duplicate edges by switching: pick a defective edge (a, b) (a self-loop, or a
     second copy of an edge) and a uniformly random other edge (c, d), and replace them by (a, d) and (c, b)
     when neither is a self-loop and the number of defects does not grow. Every switch keeps every node's
     in- and out-degree. Switches that only move a duplicate are needed: requiring every switch to remove a
     defect got stuck on a hub of out-degree 98 at n=100, which must reach every node with in-degree >= 1;
     the nodes it missed had in-degree 1, fed only by two other hubs that already reached all of its
     duplicated targets. If 1000 attempts per edge leave defects, the matching is redrawn (at most 20
     times, then the run stops with an error). Rejecting defective matchings instead does not work at any of
     these sizes (0 of 300 matchings were simple, for each of 40 sampled degree sequences per size), since
     hubs near the cap almost always get a duplicate. No further switches are made once the graph is simple.
  4. Every node with in-degree d >= 1 gets a uniformly random sign per regulator and a threshold drawn
     uniformly from [1, d], i.e. a non-constant symmetric threshold function (as in
     make_random_scale_free_models.py). A node with no regulators is an input node.

Model k of size n is drawn from seed 1000 * n + k. Stored as network.json (signs + threshold per node).
Regenerate from the repository root with
    python synthetic_data_generation/make_random_regulatory_topology_networks.py [--sizes 400 ...]
--sizes (default: all of SIZES) regenerates only those sizes' directories and leaves the others as they are;
the top-level README is then rebuilt from every size directory present.
"""
import argparse
import csv
import os
import random
import shutil
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy import optimize, stats
from scipy.special import zeta

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from attractor_learning.graphs import Network
from attractor_learning.logic import SymmetricThresholdFunction

REGULONDB_NETWORK = os.path.join(ROOT, "data", "regulondb", "network.tsv")
OUTPUT_DIR = os.path.join(ROOT, "data", "random_regulatory_topology_networks")
SIZES = [25, 50, 100, 200, 400]
MODELS_PER_SIZE = 30
SIZE_DIR = "size={}"        # key=value, so the size reaches the analysis as a swept parameter
OUT_XMIN = 8
MAX_SEQUENCE_DRAWS = 10 ** 6
MAX_SWITCH_ATTEMPTS_PER_EDGE = 1000
MAX_WIRINGS = 20


def regulondb_degrees(path=REGULONDB_NETWORK):
    """In- and out-degree of every RegulonDB gene, with self-loops dropped and duplicate pairs merged."""
    with open(path, newline="", encoding="utf-8") as f:
        pairs = {(row["regulator"], row["target"]) for row in csv.DictReader(f, delimiter="\t")
                 if row["regulator"] != row["target"]}
    genes = sorted({g for pair in pairs for g in pair})
    index = {g: i for i, g in enumerate(genes)}
    in_degrees, out_degrees = np.zeros(len(genes), int), np.zeros(len(genes), int)
    for regulator, target in pairs:
        out_degrees[index[regulator]] += 1
        in_degrees[index[target]] += 1
    return in_degrees, out_degrees, len(pairs)


class DegreeModel(object):
    """The fitted in- and out-degree distributions, and their capped versions for a given model size."""

    def __init__(self, in_degrees, out_degrees):
        self.p_in_zero = float(np.mean(in_degrees == 0))
        shifted = in_degrees[in_degrees > 0] - 1
        fit = optimize.minimize(
            lambda p: -stats.nbinom.logpmf(shifted, np.exp(p[0]), 1 / (1 + np.exp(-p[1]))).sum(),
            [0.0, 0.0], method="Nelder-Mead", options={"xatol": 1e-8, "fatol": 1e-8, "maxiter": 10000})
        self.nb_r, self.nb_p = float(np.exp(fit.x[0])), float(1 / (1 + np.exp(-fit.x[1])))

        regulators = out_degrees[out_degrees > 0]
        tail = regulators[regulators >= OUT_XMIN]
        self.out_alpha = float(optimize.minimize_scalar(
            lambda a: len(tail) * np.log(zeta(a, OUT_XMIN)) + a * np.log(tail).sum(),
            bounds=(1.01, 6), method="bounded").x)
        self.out_tail_share = float(len(tail)) / len(regulators)
        self.out_low_frequencies = {k: float(np.mean(regulators == k)) for k in range(1, OUT_XMIN)}
        self.observed_regulator_share = float(np.mean(out_degrees > 0))

    def in_pmf(self, cap):
        k = np.arange(cap + 1)
        weights = np.where(k == 0, self.p_in_zero,
                           (1 - self.p_in_zero) * stats.nbinom.pmf(np.maximum(k - 1, 0), self.nb_r, self.nb_p))
        return k, weights / weights.sum()

    def regulator_out_pmf(self, cap):
        k = np.arange(1, cap + 1)
        tail = self.out_tail_share * k.astype(float) ** -self.out_alpha / zeta(self.out_alpha, OUT_XMIN)
        weights = np.array([self.out_low_frequencies[j] if j < OUT_XMIN else tail[j - 1] for j in k])
        return k, weights / weights.sum()

    def regulator_share(self, n):
        k_in, p_in = self.in_pmf(n - 1)
        k_out, p_out = self.regulator_out_pmf(n - 1)
        share = float((k_in * p_in).sum() / (k_out * p_out).sum())
        if share > 1:
            raise ValueError("n={}: in-edges cannot be balanced even if every node regulates".format(n))
        return share


def draw_degree_sequences(model, n, rng):
    """In- and out-degree lists with equal positive sums that some simple digraph realizes."""
    k_in, p_in = model.in_pmf(n - 1)
    k_out, p_out = model.regulator_out_pmf(n - 1)
    share = model.regulator_share(n)
    for draws in range(1, MAX_SEQUENCE_DRAWS + 1):
        in_degrees = rng.choice(k_in, n, p=p_in)
        out_degrees = np.where(rng.random(n) < share, rng.choice(k_out, n, p=p_out), 0)
        if in_degrees.sum() == out_degrees.sum() > 0 and \
                nx.is_digraphical(in_degrees.tolist(), out_degrees.tolist()):
            return in_degrees.tolist(), out_degrees.tolist(), draws
    raise RuntimeError("n={}: no matching degree sequences in {} draws".format(n, MAX_SEQUENCE_DRAWS))


def repair_by_switching(edges, rng):
    """Remove self-loops and duplicate edges from a stub matching by degree-preserving switches.

    Returns (simple edge list, accepted switches), or None if the switch budget ran out."""
    edges = list(edges)
    counts = {}
    for e in edges:
        counts[e] = counts.get(e, 0) + 1

    def excess(e, count):
        """Defects carried by `count` copies of edge e: every copy of a self-loop, every copy beyond the first."""
        return count if e[0] == e[1] else max(count - 1, 0)

    def defective(i):
        return edges[i][0] == edges[i][1] or counts[edges[i]] > 1

    bad = [i for i in range(len(edges)) if defective(i)]
    switches, attempts = 0, 0
    while bad:
        attempts += 1
        if attempts > MAX_SWITCH_ATTEMPTS_PER_EDGE * len(edges):
            return None
        i = bad[rng.randrange(len(bad))]
        j = rng.randrange(len(edges))
        if j == i:
            continue
        (a, b), (c, d) = edges[i], edges[j]
        new1, new2 = (a, d), (c, b)
        if a == d or c == b or new1 == new2:
            continue
        # Accept a switch unless it adds defects. Switches that keep the count (moving a duplicate from one
        # node to another) are needed: a hub whose out-degree is close to the cap must reach nearly every
        # node, and when the nodes it misses are fed only by other hubs that already reach all of its
        # duplicated targets, no single switch removes a defect, but moving the duplicate to the other hub
        # lets that hub's spare targets absorb it.
        after = dict((e, counts.get(e, 0)) for e in (edges[i], edges[j], new1, new2))
        before = sum(excess(e, k) for e, k in after.items())
        after[edges[i]] -= 1
        after[edges[j]] -= 1
        after[new1] += 1
        after[new2] += 1
        if sum(excess(e, k) for e, k in after.items()) > before:
            continue
        for e, k in after.items():
            if k:
                counts[e] = k
            else:
                counts.pop(e, None)
        edges[i], edges[j] = new1, new2
        switches += 1
        bad = [k for k in range(len(edges)) if defective(k)]
    return edges, switches


def random_regulatory_network(model, n, seed):
    rng = random.Random(seed)
    np_rng = np.random.default_rng(seed)
    in_degrees, out_degrees, draws = draw_degree_sequences(model, n, np_rng)
    for wirings in range(1, MAX_WIRINGS + 1):
        multigraph = nx.directed_configuration_model(in_degrees, out_degrees, seed=rng.randrange(2 ** 31))
        repaired = repair_by_switching([(u, v) for u, v, _ in sorted(multigraph.edges(keys=True))], rng)
        if repaired is not None:
            edges, switches = repaired
            break
    else:
        raise RuntimeError("n={}, seed={}: switching left defects in {} wirings".format(n, seed, MAX_WIRINGS))
    names = ["v{}".format(i) for i in range(n)]
    regulator_count = [0] * n
    for _, v in edges:
        regulator_count[v] += 1
    functions = []
    for v in range(n):
        d = regulator_count[v]
        functions.append(None if d == 0 else SymmetricThresholdFunction(
            signs=[rng.choice([True, False]) for _ in range(d)], threshold=rng.randint(1, d)))
    network = Network(vertex_names=names, edges=[(names[u], names[v]) for u, v in edges],
                      vertex_functions=functions)
    actual_in = [len(v.predecessors()) for v in network.vertices]
    actual_out = [len(v.successors()) for v in network.vertices]
    assert actual_in == in_degrees and actual_out == out_degrees, "switching changed a degree"
    assert len(set(edges)) == len(edges) and all(u != v for u, v in edges)
    return network, {"sequence_draws": draws, "wirings": wirings, "switches": switches}


def stats_of(network):
    in_degrees = [len(v.predecessors()) for v in network.vertices]
    out_degrees = [len(v.successors()) for v in network.vertices]
    return {"n": len(network), "edges": len(network.edges),
            "inputs": sum(d == 0 for d in in_degrees), "regulators": sum(d > 0 for d in out_degrees),
            "max_in": max(in_degrees), "max_out": max(out_degrees),
            "in_degrees": in_degrees, "out_degrees": out_degrees}


def draw(info, path, title):
    fig, ax = plt.subplots(1, 1, figsize=(5.2, 4.0))
    for degrees, colour, label in ((info["in_degrees"], 'tab:blue', 'in-degree'),
                                   (info["out_degrees"], 'tab:orange', 'out-degree')):
        values, counts = np.unique(degrees, return_counts=True)
        keep = values > 0
        ax.scatter(values[keep], counts[keep], s=18, color=colour, label=label, alpha=0.8)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('degree')
    ax.set_ylabel('number of nodes')
    ax.set_title('degree distribution (zeros omitted)', fontsize=8)
    ax.legend(fontsize=7)
    ax.grid(True, which='both', alpha=0.3)
    fig.suptitle("{}\n{} nodes, {} edges, {} regulators, {} input nodes".format(
        title, info["n"], info["edges"], info["regulators"], info["inputs"]), fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def write_model_readme(path, name, info, seed, log):
    lines = ["# {}".format(name), "",
             "A random model with RegulonDB-fitted degree distributions (seed {}); see ../../README.md for "
             "how it was drawn.".format(seed), "",
             "| property | value |", "|---|---|",
             "| nodes | {} |".format(info["n"]),
             "| edges | {} |".format(info["edges"]),
             "| regulators (out-degree >= 1) | {} |".format(info["regulators"]),
             "| input nodes (in-degree 0) | {} |".format(info["inputs"]),
             "| max in-degree | {} |".format(info["max_in"]),
             "| max out-degree | {} |".format(info["max_out"]),
             "| degree-sequence draws | {} |".format(log["sequence_draws"]),
             "| switches to remove self-loops and duplicates | {} |".format(log["switches"]),
             "", "![degree distribution](network.png)", ""]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def write_index_readme(model, regulondb_info, rows):
    doc = __doc__.split("\n\n", 1)[1].split("Model k of size n")[0].rstrip()
    present = [n for n in SIZES if os.path.isdir(os.path.join(OUTPUT_DIR, SIZE_DIR.format(n)))]
    lines = ["# Random regulatory-topology networks", "",
             "A `graphs_dir` of {} random models, {} at each of the sizes {}, whose degree distributions are "
             "fitted to RegulonDB, with a uniformly random non-constant symmetric threshold function per "
             "non-input node.".format(len(rows), MODELS_PER_SIZE, ", ".join(str(s) for s in present)), "",
             "## Implementation", "", "```", doc, "```", "",
             "## Fitted parameters", "",
             "From `data/regulondb/network.tsv`: {} genes, {} edges after dropping self-loops and merging "
             "duplicates.".format(*regulondb_info), "",
             "| parameter | value |", "|---|---|",
             "| P(in-degree = 0) | {:.4f} |".format(model.p_in_zero),
             "| in-degree - 1 ~ NegBinomial(r, p) | r = {:.4f}, p = {:.4f} (mean {:.3f}) |".format(
                 model.nb_r, model.nb_p, model.nb_r * (1 - model.nb_p) / model.nb_p),
             "| regulators' out-degree 1..{} | empirical: {} |".format(
                 OUT_XMIN - 1, ", ".join("{:.3f}".format(model.out_low_frequencies[k])
                                         for k in range(1, OUT_XMIN))),
             "| regulators' out-degree >= {} | power law, exponent {:.4f}, share {:.4f} |".format(
                 OUT_XMIN, model.out_alpha, model.out_tail_share),
             "| observed regulator share in RegulonDB | {:.4f} |".format(model.observed_regulator_share)]
    lines += ["", "| size | regulator share used |", "|---|---|"]
    lines += ["| {} | {:.4f} |".format(n, model.regulator_share(n)) for n in present]
    lines += ["", "## Layout", "",
              "One directory per size (`{}`), one per model inside it, each holding `network.json` (signs + "
              "threshold per node, read by `Network.parse_model_dir`), a README and a degree plot.".format(
                  SIZE_DIR.format("<nodes>")), "",
              "Regenerate with `python synthetic_data_generation/make_random_regulatory_topology_networks.py`; "
              "`--sizes 400` (say) regenerates only those sizes.",
              "", "## Models", "",
              "| model | nodes | edges | regulators | input nodes | max in-degree | max out-degree |",
              "|---|---|---|---|---|---|---|"]
    lines += rows
    with open(os.path.join(OUTPUT_DIR, "README.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def index_row(name, size, info):
    return "| [{0}]({1}/{0}/README.md) | {2} | {3} | {4} | {5} | {6} | {7} |".format(
        name, SIZE_DIR.format(size), info["n"], info["edges"], info["regulators"], info["inputs"],
        info["max_in"], info["max_out"])


def generate_size(model, size):
    size_dir = os.path.join(OUTPUT_DIR, SIZE_DIR.format(size))
    if os.path.isdir(size_dir):
        shutil.rmtree(size_dir)
    per_size = []
    for k in range(MODELS_PER_SIZE):
        seed = 1000 * size + k
        name = "n{:04d}_{:02d}".format(size, k)
        network, log = random_regulatory_network(model, size, seed)
        info = stats_of(network)
        model_dir = os.path.join(size_dir, name)
        os.makedirs(model_dir)
        network.save(os.path.join(model_dir, Network.MODEL_JSON_NAME))
        draw(info, os.path.join(model_dir, "network.png"), "random regulatory-topology model " + name)
        write_model_readme(os.path.join(model_dir, "README.md"), name, info, seed, log)

        # the folder is only useful if the pipeline's readers can load it back unchanged
        assert Network.model_dir_size(model_dir) == size
        reloaded = Network.parse_model_dir(model_dir)
        check = random.Random(seed)
        for _ in range(20):
            state = [check.randint(0, 1) for _ in range(size)]
            assert [int(bool(v)) for v in reloaded.next_state(state)] == \
                   [int(bool(v)) for v in network.next_state(state)], \
                "{}: json round trip changed dynamics".format(name)
        per_size.append((info, log))
    print("size {:>4}: edges {}-{}, regulators {}-{}, max out {}-{}, sequence draws max {}, switches max {}"
          .format(size, min(i["edges"] for i, _ in per_size), max(i["edges"] for i, _ in per_size),
                  min(i["regulators"] for i, _ in per_size), max(i["regulators"] for i, _ in per_size),
                  min(i["max_out"] for i, _ in per_size), max(i["max_out"] for i, _ in per_size),
                  max(l["sequence_draws"] for _, l in per_size), max(l["switches"] for _, l in per_size)),
          flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--sizes", type=int, nargs="+", default=SIZES, choices=SIZES,
                        help="sizes to (re)generate; other sizes' directories are left untouched")
    args = parser.parse_args()
    in_degrees, out_degrees, edge_count = regulondb_degrees()
    model = DegreeModel(in_degrees, out_degrees)
    print("fit: P(in=0)={:.4f}  NB r={:.4f} p={:.4f}  out alpha={:.4f} tail share={:.4f}".format(
        model.p_in_zero, model.nb_r, model.nb_p, model.out_alpha, model.out_tail_share))
    if not os.path.isdir(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    for size in args.sizes:
        generate_size(model, size)

    # index of every size present, regenerated or not
    rows = []
    for size in SIZES:
        size_dir = os.path.join(OUTPUT_DIR, SIZE_DIR.format(size))
        if not os.path.isdir(size_dir):
            continue
        for name in sorted(os.listdir(size_dir)):
            rows.append(index_row(name, size, stats_of(Network.parse_model_dir(os.path.join(size_dir, name)))))
    write_index_readme(model, (len(in_degrees), edge_count), rows)
    print("\nregenerated sizes {}; README indexes {} models in {}".format(args.sizes, len(rows), OUTPUT_DIR))


if __name__ == "__main__":
    main()
