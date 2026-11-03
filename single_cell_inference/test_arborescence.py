"""Checks arborescence.minimum_arborescence against networkx.minimum_spanning_arborescence on random graphs:
the same edges and cost when the optimum is unique (continuous costs, or integer costs with a small jitter),
and the same cost when ties are left in. Run from the repository root:
python single_cell_inference/test_arborescence.py"""

import os
import sys

import networkx as nx
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arborescence import minimum_arborescence  # noqa: E402


def random_graph(rng, n_nodes, n_extra, integer, jitter):
    """A root (node 0) with an edge to every node, plus n_extra random edges; costs in [0, 10)."""
    src = [0] * (n_nodes - 1)
    dst = list(range(1, n_nodes))
    for _ in range(n_extra):
        u, v = rng.integers(1, n_nodes, size=2)
        if u != v:
            src.append(int(u))
            dst.append(int(v))
    weight = rng.integers(0, 10, size=len(src)).astype(float) if integer else rng.uniform(0, 10, size=len(src))
    weight[:n_nodes - 1] += 10  # root edges cost more, as in the snapshot tree
    if jitter:
        weight += rng.uniform(0, 1e-6, size=len(src))
    return np.array(src), np.array(dst), weight


def networkx_tree(src, dst, weight):
    graph = nx.DiGraph()
    for u, v, w in zip(src, dst, weight):
        if not graph.has_edge(u, v) or graph[u][v]["weight"] > w:
            graph.add_edge(int(u), int(v), weight=float(w))
    tree = nx.minimum_spanning_arborescence(graph, attr="weight")
    return {(u, v) for u, v in tree.edges()}, sum(d["weight"] for _, _, d in tree.edges(data=True))


def main():
    rng = np.random.default_rng(0)
    counts = {"same edges": 0, "same cost only": 0}
    for trial in range(300):
        n_nodes = int(rng.integers(3, 80))
        integer, jitter = trial % 3 != 0, trial % 3 == 1
        src, dst, weight = random_graph(rng, n_nodes, int(rng.integers(0, 6 * n_nodes)), integer, jitter)
        # duplicate (u, v) edges: keep the cheapest, as networkx's DiGraph does
        cheapest = {}
        for e, (u, v) in enumerate(zip(src, dst)):
            if (u, v) not in cheapest or weight[e] < weight[cheapest[u, v]]:
                cheapest[u, v] = e
        index = np.array(sorted(cheapest.values()))
        src, dst, weight = src[index], dst[index], weight[index]
        chosen = minimum_arborescence(n_nodes, src, dst, weight, root=0)
        ours = {(int(src[e]), int(dst[e])) for e in chosen}
        assert len(ours) == n_nodes - 1 and {v for _, v in ours} == set(range(1, n_nodes))
        theirs, their_cost = networkx_tree(src, dst, weight)
        our_cost = float(weight[chosen].sum())
        assert abs(our_cost - their_cost) < 1e-6, (trial, our_cost, their_cost)
        unique = not integer or jitter
        if unique:
            assert ours == theirs, trial
        counts["same edges" if ours == theirs else "same cost only"] += 1
    print("300 random graphs: equal cost in all;", counts,
          "(edge sets may differ only with integer costs and no jitter, i.e. tied optima)")


if __name__ == "__main__":
    main()
