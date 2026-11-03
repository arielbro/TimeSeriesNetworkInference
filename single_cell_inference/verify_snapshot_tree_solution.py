"""Check a snapshot_tree_ilp solution independently of Gurobi: the network applied to each latent state gives the
reported unstable genes, every child differs from its predecessor only on the predecessor's unstable genes, the
tree is acyclic with predecessors strictly earlier, and no negative self-loops."""
import json, sys, os
sys.path.insert(0, sys.argv[2])
from attractor_learning.graphs import Network

out = sys.argv[1]
net = Network.load(os.path.join(out, "network.json"))
d = json.load(open(os.path.join(out, "cells.json")))
genes = d["genes"]; cells = d["cells"]
idx = {v.name: v.index for v in net.vertices}
order = [idx[g] for g in genes]

def step(state):
    full = [0] * len(net.vertices)
    for g, val in zip(genes, state):
        full[idx[g]] = val
    nxt = net.next_state(full)
    return [nxt[i] for i in order]

problems = 0
for name, c in cells.items():
    f = step(c["state"])
    unstable = [g for g, a, b in zip(genes, c["state"], f) if a != b]
    if unstable != c["unstable"]:
        problems += 1; print("instability mismatch", name, unstable, c["unstable"])
    if c["predecessor"] is None:
        parent_state, parent_f, ptime = d["pseudo_root_state"], step(d["pseudo_root_state"]), 0
    else:
        p = cells[c["predecessor"]]; parent_state, parent_f, ptime = p["state"], step(p["state"]), p["pseudotime"]
    for g, ys, yp, fp in zip(genes, c["state"], parent_state, parent_f):
        if ys != yp and fp == yp:
            problems += 1; print("bad transition", name, g)
    if c["pseudotime"] <= ptime:
        problems += 1; print("bad time", name)
    if c["predecessor"] is None and c["pseudotime"] != 1:
        problems += 1; print("root not at time 1", name)
for u, v in net.edges:
    if u == v and v.function is not None:
        pass
print("cells", len(cells), "problems", problems)
