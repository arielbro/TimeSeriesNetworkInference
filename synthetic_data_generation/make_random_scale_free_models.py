"""Build data/random_scale_free_models: ten completely random models at each of the sizes 5, 10, 50, 250
and 1000, usable as a `graphs_dir` (the role data/cellcollective_models plays for the real runs).

Topology: networkx's directed scale-free generator (a preferential-attachment process, so in- and
out-degrees are power-law and a few nodes are hubs), with self-loops dropped and parallel edges collapsed.
Functions: for every node with at least one regulator, a uniformly random sign per regulator and a
uniformly random threshold in [1, in-degree] - i.e. a uniformly random symmetric threshold function over
that node's inputs, excluding the two constant functions (threshold 0 or in-degree + 1), which the
unknown-topology ILP cannot represent anyway. A node with no regulators is an input node.

These are written as `network.json`, not as truth tables: a scale-free graph has hub nodes whose in-degree
runs into the dozens, and a cellcollective truth table needs 2**in-degree rows. The JSON form stores each
symmetric threshold function as its signs and threshold, so the file stays small - the readers used by data
generation and inference accept either format (Network.parse_model_dir / model_dir_size).
"""
import os
import random
import shutil
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from attractor_learning.graphs import Network
from attractor_learning.logic import SymmetricThresholdFunction

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data",
                          "random_scale_free_models")
SIZES = [5, 10, 50, 250, 1000]
MODELS_PER_SIZE = 10
ENUMERATE_UP_TO = 10        # full state space only where 2**n is small enough to walk
DRAW_GRAPH_UP_TO = 12       # above this a node-link picture says nothing; draw the degree profile instead


def random_scale_free_network(n, seed):
    """A scale-free digraph with a uniformly random symmetric threshold function per non-input node."""
    rng = random.Random(seed)
    multigraph = nx.scale_free_graph(n, seed=seed)
    graph = nx.DiGraph((u, v) for u, v in multigraph.edges() if u != v)   # simple, no self-loops
    graph.add_nodes_from(range(n))
    names = ["v{}".format(i) for i in range(n)]
    edges = [(names[u], names[v]) for u, v in graph.edges()]
    regulators = {}
    for source, target in edges:
        regulators.setdefault(target, []).append(source)
    functions = []
    for name in names:
        degree = len(regulators.get(name, []))
        if degree == 0:
            functions.append(None)                      # input node: holds its value
        else:
            functions.append(SymmetricThresholdFunction(
                signs=[rng.choice([True, False]) for _ in range(degree)],
                threshold=rng.randint(1, degree)))
    return Network(vertex_names=names, edges=edges, vertex_functions=functions)


def attractor_summary(network):
    """(number of attractors, sorted periods) by walking the whole state space; None when too large."""
    n = len(network)
    if n > ENUMERATE_UP_TO:
        return None
    transitions = {}
    for bits in range(2 ** n):
        state = tuple((bits >> (n - 1 - i)) & 1 for i in range(n))
        transitions[state] = tuple(int(bool(v)) for v in network.next_state(list(state)))
    cycle_id, cycles = {}, []
    for start in transitions:
        path, state = [], start
        while state not in cycle_id and state not in path:
            path.append(state)
            state = transitions[state]
        if state in path:
            cycles.append(path[path.index(state):])
            for member in cycles[-1]:
                cycle_id[member] = len(cycles) - 1
        target = cycle_id[state]
        for member in path:
            cycle_id.setdefault(member, target)
    return len(cycles), sorted(len(c) for c in cycles)


def stats(network):
    in_degrees = [len(v.predecessors()) for v in network.vertices]
    out_degrees = [len(v.successors()) for v in network.vertices]
    return {
        "n": len(network), "edges": len(network.edges),
        "inputs": sum(1 for d in in_degrees if d == 0),
        "max_in": max(in_degrees), "max_out": max(out_degrees),
        "mean_in": float(np.mean(in_degrees)),
        "in_degrees": in_degrees, "out_degrees": out_degrees,
    }


def draw(network, info, path, title):
    n = info["n"]
    if n <= DRAW_GRAPH_UP_TO:
        fig, axes = plt.subplots(1, 2, figsize=(9, 4.0), gridspec_kw={'width_ratios': [1.1, 1]})
        graph = nx.DiGraph()
        graph.add_nodes_from(v.name for v in network.vertices)
        signs = {}
        for v in network.vertices:
            if v.function is None:
                continue
            for u, sign in zip(v.predecessors(), v.function.signs):
                graph.add_edge(u.name, v.name)
                signs[(u.name, v.name)] = 1 if sign else -1
        positions = nx.circular_layout(graph)
        labels = {v.name: "{}\n".format(v.name) +
                          (r"$\theta$={}".format(v.function.threshold) if v.function is not None else "input")
                  for v in network.vertices}
        nx.draw_networkx_nodes(graph, positions, ax=axes[0], node_color='white', edgecolors='black',
                               node_size=1000)
        nx.draw_networkx_labels(graph, positions, labels=labels, ax=axes[0], font_size=7)
        nx.draw_networkx_edges(graph, positions, ax=axes[0],
                               edgelist=[e for e, s in signs.items() if s > 0], edge_color='tab:green',
                               arrowstyle='-|>', arrowsize=12, node_size=1000,
                               connectionstyle='arc3,rad=0.12')
        nx.draw_networkx_edges(graph, positions, ax=axes[0],
                               edgelist=[e for e, s in signs.items() if s < 0], edge_color='tab:red',
                               arrowstyle='-[', arrowsize=8, node_size=1000,
                               connectionstyle='arc3,rad=0.12')
        axes[0].margins(0.2)
        axes[0].axis('off')
        degree_ax = axes[1]
    else:
        fig, degree_ax = plt.subplots(1, 1, figsize=(5.2, 4.0))

    for degrees, colour, label in ((info["in_degrees"], 'tab:blue', 'in-degree'),
                                   (info["out_degrees"], 'tab:orange', 'out-degree')):
        values, counts = np.unique(degrees, return_counts=True)
        keep = values > 0
        degree_ax.scatter(values[keep], counts[keep], s=18, color=colour, label=label, alpha=0.8)
    degree_ax.set_xscale('log')
    degree_ax.set_yscale('log')
    degree_ax.set_xlabel('degree')
    degree_ax.set_ylabel('number of nodes')
    degree_ax.set_title('degree distribution (zeros omitted)', fontsize=8)
    degree_ax.legend(fontsize=7)
    degree_ax.grid(True, which='both', alpha=0.3)
    fig.suptitle("{}\n{} nodes, {} edges, {} input nodes, max in-degree {}".format(
        title, info["n"], info["edges"], info["inputs"], info["max_in"]), fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def write_readme(path, name, network, info, seed, attractors):
    lines = ["# {}".format(name), "",
             "A random scale-free model: topology from networkx's directed scale-free generator "
             "(seed {}), with a uniformly random symmetric threshold function per non-input node - random "
             "signs, and a threshold drawn uniformly from [1, in-degree].".format(seed), "",
             "| property | value |", "|---|---|",
             "| nodes | {} |".format(info["n"]),
             "| edges | {} |".format(info["edges"]),
             "| input nodes (in-degree 0) | {} |".format(info["inputs"]),
             "| mean in-degree | {:.2f} |".format(info["mean_in"]),
             "| max in-degree | {} |".format(info["max_in"]),
             "| max out-degree | {} |".format(info["max_out"])]
    if attractors is not None:
        count, periods = attractors
        lines.append("| attractors (full enumeration) | {} |".format(count))
        lines.append("| attractor periods | {} |".format(", ".join(str(p) for p in sorted(set(periods)))))
    else:
        lines.append("| attractors | not enumerated (2**{} states) |".format(info["n"]))
    lines += ["",
              "Stored as `network.json` (signs + threshold per node) rather than truth tables, which would "
              "need 2**{} rows for this model's largest node.".format(info["max_in"]), "",
              "![topology and degree distribution](network.png)", ""]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def main():
    if os.path.isdir(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR)
    index = ["# Random scale-free models", "",
             "A `graphs_dir` of {} completely random models: {} at each of the sizes {}. Topology from "
             "networkx's directed scale-free generator; every non-input node gets a uniformly random "
             "symmetric threshold function (random signs, threshold uniform in [1, in-degree]).".format(
                 len(SIZES) * MODELS_PER_SIZE, MODELS_PER_SIZE,
                 ", ".join(str(s) for s in SIZES)), "",
             "Stored as `network.json` (each function as signs + threshold), not as cellcollective truth "
             "tables: hub nodes here have in-degrees in the dozens, and a truth table needs 2**in-degree "
             "rows. `Network.parse_model_dir` / `Network.model_dir_size` read either format, so this "
             "folder works anywhere a graphs_dir is expected.", "",
             "Regenerate with `python synthetic_data_generation/make_random_scale_free_models.py`.", "",
             "| model | nodes | edges | input nodes | mean in-degree | max in-degree | attractors |",
             "|---|---|---|---|---|---|---|"]
    for size in SIZES:
        per_size = []
        for index_in_size in range(MODELS_PER_SIZE):
            seed = 1000 * size + index_in_size
            name = "n{:04d}_{:02d}".format(size, index_in_size)
            network = random_scale_free_network(size, seed)
            info = stats(network)
            attractors = attractor_summary(network)
            model_dir = os.path.join(OUTPUT_DIR, name)
            os.makedirs(model_dir)
            network.save(os.path.join(model_dir, Network.MODEL_JSON_NAME))
            draw(network, info, os.path.join(model_dir, "network.png"),
                 "random scale-free model {}".format(name))
            write_readme(os.path.join(model_dir, "README.md"), name, network, info, seed, attractors)

            # the folder is only useful if the pipeline's readers can load it back unchanged
            assert Network.model_dir_size(model_dir) == size
            reloaded = Network.parse_model_dir(model_dir)
            rng = random.Random(seed)
            for _ in range(20):
                state = [rng.randint(0, 1) for _ in range(size)]
                assert [int(bool(v)) for v in reloaded.next_state(state)] == \
                       [int(bool(v)) for v in network.next_state(state)], \
                    "{}: json round trip changed dynamics".format(name)
            index.append("| [{}]({}/README.md) | {} | {} | {} | {:.2f} | {} | {} |".format(
                name, name, info["n"], info["edges"], info["inputs"], info["mean_in"], info["max_in"],
                attractors[0] if attractors is not None else "not enumerated"))
            per_size.append(info)
        print("size {:>4}: {} models, edges {}-{}, max in-degree {}-{}".format(
            size, MODELS_PER_SIZE,
            min(i["edges"] for i in per_size), max(i["edges"] for i in per_size),
            min(i["max_in"] for i in per_size), max(i["max_in"] for i in per_size)))
    with open(os.path.join(OUTPUT_DIR, "README.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(index) + "\n")
    print("\nwrote {} models to {}".format(len(SIZES) * MODELS_PER_SIZE, OUTPUT_DIR))


if __name__ == "__main__":
    main()
