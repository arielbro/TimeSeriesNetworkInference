"""Build data/regulondb_strict_regulators_random_models: 50 random models on the RegulonDB strict regulator
subnetwork, usable as a `graphs_dir` for generate.py.

The source graph is the one data/regulondb/make_inference_data.py --subnetwork strict_regulators uses as its
true network: starting from data/regulondb/gse20305_network.tsv, every gene with no outgoing edge to another
remaining gene is removed, repeatedly, until none is left (self-loops do not count), and the edges between
the kept genes are kept (115 genes, 446 edges, 72 of them self-loops), in gse20305_binary.json's gene order.

For each model, with m the number of edges:
  Topology: 10m successful degree-preserving switches. A switch picks two edges (a, b) and (c, d) uniformly
            at random among the edges that are not self-loops and replaces them by (a, d) and (c, b), unless
            either of those is a self-loop or already an edge. Every node keeps its in- and out-degree, and
            the self-loops (autoregulation) are the source graph's in every model.
  Functions: every node with in-degree d >= 1 gets a uniformly random sign per regulator and a threshold
            drawn uniformly from [1, d], i.e. a symmetric threshold function that is not constant (as in
            make_random_scale_free_models.py). A node with no regulators is an input node.

Model k is drawn from random.Random(k). The edge list is kept as a list in a fixed order, not as a set, so
the result does not depend on Python's string hashing (Network.randomize_edges_by_switching samples from a
set, so it is not reproducible across runs, and it counts attempts rather than switches).

Stored as network.json (signs + threshold per node): the largest in-degree is 24, and a truth table would
need 2**24 rows. Regenerate from the repository root with
    python synthetic_data_generation/make_regulondb_random_models.py
"""
import json
import os
import random
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from attractor_learning.graphs import Network  # noqa: E402
from attractor_learning.logic import SymmetricThresholdFunction  # noqa: E402

DATA_DIR = os.path.join(ROOT, "data", "regulondb")
OUTPUT_DIR = os.path.join(ROOT, "data", "regulondb_strict_regulators_random_models")
N_MODELS = 50
SWITCHES_PER_EDGE = 10
MAX_ATTEMPTS_PER_SWITCH = 1000  # only guards against a graph where almost no switch is possible


def read_tsv(path):
    with open(path, encoding="utf8") as in_file:
        return [line.rstrip("\r\n").split("\t") for line in in_file if line.strip()][1:]


def regulator_graph():
    """(genes, edges) of the strict regulator subnetwork, as make_inference_data.py --subnetwork
    strict_regulators builds it."""
    with open(os.path.join(DATA_DIR, "gse20305_binary.json")) as meta_file:
        genes = json.load(meta_file)["genes"]
    edges = [(row[0], row[1]) for row in read_tsv(os.path.join(DATA_DIR, "gse20305_network.tsv"))]
    regulators = {u for u, v in edges if u != v}
    while True:
        pruned = {u for u, v in edges if u != v and u in regulators and v in regulators}
        if pruned == regulators:
            break
        regulators = pruned
    return ([gene for gene in genes if gene in regulators],
            [(u, v) for u, v in edges if u in regulators and v in regulators])


def switch_edges(edges, n_switches, rng):
    """edges after n_switches successful degree-preserving switches of non-self-loop edges that create no
    self-loop, so the self-loops are left as they are."""
    self_loops = [(u, v) for u, v in edges if u == v]
    edges = [(u, v) for u, v in edges if u != v]
    edge_set = set(edges) | set(self_loops)
    done, attempts = 0, 0
    while done < n_switches:
        attempts += 1
        if attempts > MAX_ATTEMPTS_PER_SWITCH * n_switches:
            raise RuntimeError("only {} of {} switches succeeded".format(done, n_switches))
        i, j = rng.sample(range(len(edges)), 2)
        (a, b), (c, d) = edges[i], edges[j]
        if a == d or c == b or (a, d) in edge_set or (c, b) in edge_set:
            continue  # would create a self-loop or duplicate an edge (this includes a == c and b == d)
        edge_set -= {(a, b), (c, d)}
        edge_set |= {(a, d), (c, b)}
        edges[i], edges[j] = (a, d), (c, b)
        done += 1
    return edges + self_loops, attempts


def random_model(genes, edges, seed):
    rng = random.Random(seed)
    switched, attempts = switch_edges(edges, SWITCHES_PER_EDGE * len(edges), rng)
    network = Network(vertex_names=genes, edges=switched)
    functions = []
    for vertex in network.vertices:
        degree = len(vertex.predecessors())
        functions.append(None if degree == 0 else SymmetricThresholdFunction(
            signs=[rng.choice([True, False]) for _ in range(degree)], threshold=rng.randint(1, degree)))
    return Network(vertex_names=genes, edges=switched, vertex_functions=functions), switched, attempts


def degrees(genes, edges):
    in_degree = {gene: 0 for gene in genes}
    out_degree = {gene: 0 for gene in genes}
    for u, v in edges:
        out_degree[u] += 1
        in_degree[v] += 1
    return in_degree, out_degree


def main():
    genes, edges = regulator_graph()
    source_degrees = degrees(genes, edges)
    if os.path.isdir(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR)
    rows = []
    for seed in range(N_MODELS):
        name = "model_{:02d}".format(seed)
        network, switched, attempts = random_model(genes, edges, seed)
        assert degrees(genes, switched) == source_degrees and len(set(switched)) == len(edges)
        assert {(u, v) for u, v in switched if u == v} == {(u, v) for u, v in edges if u == v}
        model_dir = os.path.join(OUTPUT_DIR, name)
        os.makedirs(model_dir)
        network.save(os.path.join(model_dir, Network.MODEL_JSON_NAME))
        # the folder is only useful if the pipeline's readers load it back unchanged
        reloaded = Network.parse_model_dir(model_dir)
        check_rng = random.Random(seed)
        for _ in range(20):
            state = [check_rng.randint(0, 1) for _ in range(len(genes))]
            assert [int(bool(x)) for x in reloaded.next_state(state)] == \
                   [int(bool(x)) for x in network.next_state(state)], name
        kept = len(set(switched) & set(edges))
        self_loops = sum(1 for u, v in switched if u == v)
        rows.append("| {} | {} | {} | {} | {} |".format(name, seed, kept, self_loops, attempts))
        print("{}: {} of {} source edges kept, {} self-loops, {} attempts".format(
            name, kept, len(edges), self_loops, attempts))
    in_degree, out_degree = source_degrees
    lines = ["# Random models on the RegulonDB strict regulator subnetwork", "",
             "A `graphs_dir` of {} random models, each a degree-preserving randomization of the RegulonDB "
             "strict regulator subnetwork ({} genes, {} edges, {} self-loops; input nodes {}, largest in-degree {}, "
             "largest out-degree {}) with random symmetric threshold functions. Written by "
             "`synthetic_data_generation/make_regulondb_random_models.py`, whose docstring gives the "
             "procedure.".format(N_MODELS, len(genes), len(edges), sum(1 for u, v in edges if u == v),
                                 sum(1 for d in in_degree.values() if d == 0), max(in_degree.values()),
                                 max(out_degree.values())), "",
             "- Topology: {} x {} = {} successful edge switches per model; every gene keeps its in- and "
             "out-degree. Self-loops are never switched or created, so every model keeps the source's "
             "autoregulation.".format(
                 SWITCHES_PER_EDGE, len(edges), SWITCHES_PER_EDGE * len(edges)),
             "- Functions: a uniformly random sign per regulator and a threshold uniform in [1, in-degree].",
             "- Model k uses seed k.", "",
             "| model | seed | source edges kept | self-loops | switch attempts |", "|---|---|---|---|---|"] + rows
    with open(os.path.join(OUTPUT_DIR, "README.md"), "w", encoding="utf-8") as readme:
        readme.write("\n".join(lines) + "\n")
    print("wrote {} models to {}".format(N_MODELS, OUTPUT_DIR))


if __name__ == "__main__":
    main()
