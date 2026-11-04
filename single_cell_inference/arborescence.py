"""Minimum-cost spanning arborescence (Chu-Liu/Edmonds) on NumPy edge arrays.

networkx.minimum_spanning_arborescence stores a full copy of the working graph and branching at every cycle
contraction, so its memory grows as (number of contractions) x (number of edges): about 2 GB on the 531-cell
LMPP candidate graph. This version contracts every cycle of a level at once and runs the levels in a loop (the
2,167-cell graph needs more than 1,000 levels, too deep for recursion). For the expansion it keeps, per level,
only arrays over the nodes: which contracted node each original node belongs to, each node's cheapest incoming
edge, and whether the node lies on a cycle. Memory is O(E + levels x V). Both implementations are exact, so the
minimum cost is the same; with a unique optimum (e.g. after a random jitter on the costs) so is the arborescence.
"""

import numpy as np


def _best_incoming(n_nodes, dst, weight):
    """Per node, the position (in dst) of its cheapest incoming edge, -1 for a node with none."""
    order = np.lexsort((weight, dst))
    first = np.ones(len(order), dtype=bool)
    first[1:] = dst[order][1:] != dst[order][:-1]
    best = np.full(n_nodes, -1, dtype=np.int64)
    best[dst[order][first]] = order[first]
    return best


def _cycles(n_nodes, parent):
    """The cycles of the functional graph node -> parent[node] (-1 for none), as a cycle label per node (-1
    outside every cycle) and the number of cycles."""
    label = np.full(n_nodes, -1, dtype=np.int64)
    state = np.zeros(n_nodes, dtype=np.int8)  # 0 unvisited, 1 on the current path, 2 done
    n_cycles = 0
    for start in range(n_nodes):
        if state[start]:
            continue
        path, v = [], start
        while v != -1 and state[v] == 0:
            state[v] = 1
            path.append(v)
            v = parent[v]
        if v != -1 and state[v] == 1:  # closed a cycle at v
            u = v
            while True:
                label[u] = n_cycles
                u = parent[u]
                if u == v:
                    break
            n_cycles += 1
        for u in path:
            state[u] = 2
    return label, n_cycles


def minimum_arborescence(n_nodes, src, dst, weight, root=0):
    """A minimum-cost spanning arborescence of the directed graph on nodes 0..n_nodes-1 with edges
    src[e] -> dst[e] of cost weight[e], rooted at root: the array of chosen edge indices, one per non-root
    node, in increasing order. Every node must be reachable from the root. Self-loops are ignored."""
    src, dst = np.asarray(src, dtype=np.int64), np.asarray(dst, dtype=np.int64)
    alive = np.flatnonzero((src != dst) & (dst != root))  # original indices of the edges still in play
    reduced = np.asarray(weight, dtype=float)[alive].copy()  # their current (reduced) costs
    node_of = np.arange(n_nodes, dtype=np.int64)  # current node of each original node
    current_nodes, current_root = n_nodes, root
    levels = []  # per contraction: (node_of before it, best original edge per node, cycle label per node)

    while True:
        src_now, dst_now = node_of[src[alive]], node_of[dst[alive]]
        best = _best_incoming(current_nodes, dst_now, reduced)
        best[current_root] = -1
        missing = [v for v in range(current_nodes) if v != current_root and best[v] < 0]
        if missing:
            raise ValueError("{} nodes cannot be reached from the root".format(len(missing)))
        parent = np.where(best >= 0, src_now[np.maximum(best, 0)], -1)
        cycle, n_cycles = _cycles(current_nodes, parent)
        best_original = np.where(best >= 0, alive[np.maximum(best, 0)], -1)
        if n_cycles == 0:
            chosen = set(int(e) for v, e in enumerate(best_original) if v != current_root)
            break
        levels.append((node_of.copy(), best_original, cycle))

        # contract every cycle into one node; an edge into a cycle node v is charged relative to v's cycle edge,
        # which it would replace
        in_cycle = cycle >= 0
        entering = in_cycle[dst_now]
        reduced[entering] -= reduced[best[dst_now[entering]]]
        new_id = np.empty(current_nodes, dtype=np.int64)
        singles = np.flatnonzero(~in_cycle)
        new_id[singles] = n_cycles + np.arange(len(singles))
        new_id[in_cycle] = cycle[in_cycle]
        node_of = new_id[node_of]
        current_root = int(new_id[current_root])
        current_nodes = n_cycles + len(singles)
        keep = node_of[src[alive]] != node_of[dst[alive]]
        alive, reduced = alive[keep], reduced[keep]

    # expand, last contraction first: each cycle keeps its edges except the one into the node that the chosen
    # edge entering the cycle reaches
    for level_node_of, best_original, cycle in reversed(levels):
        result = set(int(best_original[v]) for v in np.flatnonzero(cycle >= 0))
        for e in chosen:
            v = level_node_of[dst[e]]
            if cycle[v] >= 0:
                result.discard(int(best_original[v]))
            result.add(e)
        chosen = result
    return np.array(sorted(chosen), dtype=np.int64)
