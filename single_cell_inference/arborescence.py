"""Minimum-cost spanning arborescence (Chu-Liu/Edmonds) on NumPy edge arrays.

networkx.minimum_spanning_arborescence stores a full copy of the working graph and branching at every cycle
contraction, so its memory grows as (number of contractions) x (number of edges): about 2 GB on the 531-cell
LMPP candidate graph. This version contracts every cycle of a level at once, so there are few levels, and keeps
only integer and float arrays per level. Both are exact, so the minimum cost is the same; with a unique optimum
(e.g. after a random jitter on the costs) so is the arborescence.
"""

import numpy as np


def _best_incoming(n_nodes, dst, weight, root):
    """Per node, the index of its cheapest incoming edge (-1 for the root or a node with none)."""
    order = np.lexsort((weight, dst))
    first = np.ones(len(order), dtype=bool)
    first[1:] = dst[order][1:] != dst[order][:-1]
    best = np.full(n_nodes, -1, dtype=np.int64)
    best[dst[order][first]] = order[first]
    best[root] = -1
    return best


def _cycles(n_nodes, parent, root):
    """The cycles of the functional graph node -> parent[node] (parent[root] = -1), as a component label per node
    (-1 outside every cycle) and the number of cycles."""
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


def _solve(n_nodes, src, dst, weight, root):
    """Indices of the edges of a minimum arborescence rooted at root, one per non-root node."""
    best = _best_incoming(n_nodes, dst, weight, root)
    missing = [v for v in range(n_nodes) if v != root and best[v] < 0]
    if missing:
        raise ValueError("nodes {} cannot be reached from the root".format(missing[:5]))
    parent = np.where(best >= 0, src[np.maximum(best, 0)], -1)
    parent[root] = -1
    cycle, n_cycles = _cycles(n_nodes, parent, root)
    if n_cycles == 0:
        return best[np.arange(n_nodes) != root]

    # contract: every cycle becomes one node, every other node keeps its own
    in_cycle = cycle >= 0
    new_id = np.empty(n_nodes, dtype=np.int64)
    singles = np.flatnonzero(~in_cycle)
    new_id[singles] = n_cycles + np.arange(len(singles))
    new_id[in_cycle] = cycle[in_cycle]
    keep = new_id[src] != new_id[dst]
    kept = np.flatnonzero(keep)
    entering_cycle = in_cycle[dst[kept]]
    new_weight = weight[kept].copy()
    # an edge into a cycle node v is charged relative to v's cycle edge, which it would replace
    new_weight[entering_cycle] -= weight[best[dst[kept][entering_cycle]]]
    chosen = kept[_solve(n_cycles + len(singles), new_id[src[kept]], new_id[dst[kept]], new_weight,
                         new_id[root])]

    # expand: each cycle keeps its edges except the one into the node its chosen entering edge reaches
    result = set(int(best[v]) for v in np.flatnonzero(in_cycle))
    for e in chosen:
        v = dst[e]
        if in_cycle[v]:
            result.discard(int(best[v]))
        result.add(int(e))
    return np.array(sorted(result), dtype=np.int64)


def minimum_arborescence(n_nodes, src, dst, weight, root=0):
    """A minimum-cost spanning arborescence of the directed graph on nodes 0..n_nodes-1 with edges
    src[e] -> dst[e] of cost weight[e], rooted at root: the array of chosen edge indices, one per non-root
    node. Every node must be reachable from the root. Self-loops are ignored."""
    src, dst = np.asarray(src, dtype=np.int64), np.asarray(dst, dtype=np.int64)
    weight = np.asarray(weight, dtype=float)
    usable = np.flatnonzero((src != dst) & (dst != root))
    return usable[_solve(n_nodes, src[usable], dst[usable], weight[usable], root)]
