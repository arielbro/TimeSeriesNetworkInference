"""The joint binarization, cell ordering and signed threshold network ILP of joint_orientation_ilp.tex, on the
Hamey et al. haematopoietic stem and progenitor qRT-PCR data in single_cell_inference/data. Run from the
repository root.

Indexing follows the note: state 0 is the pseudo-root and states 1..N are the cells, so Y[c] is cell c's latent
state and Y[0] the pseudo-root's. The subprograms are one function each:

    add_threshold_functions   theta, p, q, u of the manuscript's Sections 3.3.1-3.3.2 (the same encoding as
                              inference.binary_inference_ideas.infer_unknown_topology_symmetric), plus the
                              no-negative-self-loop constraint
    add_states                Y, and the agreement terms C against the detection calls
    add_function_outputs      FUNCTIONS: F(Y[c]) for every state, through attractor_learning.ilp.add_path_to_model
    add_instability           INSTABILITY: U[c][j] = F(Y[c])_j xor Y[c][j]
    add_tree                  TREE: predecessor indicators pi and pseudotimes tau
    add_transitions           TRANSITIONS: a child differs from its predecessor only on genes unstable there
    build_model               all of the above and the objective

Candidate predecessors are the symmetrized k nearest neighbours of each cell in Hamming distance on the
detection calls, plus the pseudo-root. Candidate inputs of node j are the regulators and j itself.
"""

import argparse
import json
import os
import sys
import time

import gurobipy
import numpy as np
from gurobipy import GRB

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "data"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from attractor_learning import ilp  # noqa: E402
from attractor_learning.graphs import Network  # noqa: E402
from attractor_learning.logic import SymmetricThresholdFunction  # noqa: E402
from convert_hamey import EXPRESSION_FILE, ORDER_FILE  # noqa: E402
from normalize_hamey_qpcr import DATA_DIR, PROCESSED_FILE, RAW_FILE, UNDETECTED_CT, read_csv  # noqa: E402

PSEUDO_ROOT = 0
# normalizers, and genes Hamey et al. removed for technical reasons (SI)
DROPPED_GENES = ("Ubc", "Polr2a", "Eif2b1", "Egfl7", "Gfi1", "Sfpi1")


# ---------------------------------------------------------------------------------------------- data

def load_detection_data(data_dir=DATA_DIR, dropped_genes=DROPPED_GENES):
    """(genes, cell names, cell types, N x n 0/1 array): the cells of the processed file (those that passed
    quality control), each gene 1 iff detected (raw Ct below UNDETECTED_CT), Hamey et al.'s binarization."""
    processed = read_csv(os.path.join(data_dir, PROCESSED_FILE))
    raw = read_csv(os.path.join(data_dir, RAW_FILE))
    names = [row[1] for row in processed[1:]]
    types = [row[0] for row in processed[1:]]
    genes = [gene for gene in processed[0][2:] if gene not in dropped_genes]
    raw_rows = {row[0]: row[2:] for row in raw[1:]}
    raw_columns = [raw[0][2:].index(gene) for gene in genes]
    detected = np.array([[float(raw_rows[name][j]) < UNDETECTED_CT for j in raw_columns] for name in names])
    return genes, names, types, detected.astype(np.int8)


def hamey_regulators(source_dir=DATA_DIR):
    """The 31 genes Hamey et al. modelled as transcription factors: the columns of their LMPP binary file."""
    with open(os.path.join(source_dir, EXPRESSION_FILE)) as expression_file:
        return expression_file.readline().split()


def lmpp_cells(source_dir=DATA_DIR):
    """The cells of Hamey et al.'s HSC -> LMPP branch, in their pseudotime order."""
    with open(os.path.join(source_dir, ORDER_FILE)) as order_file:
        return [line.strip() for line in order_file if line.strip()]


def candidate_predecessors(detected, k):
    """P(c) for cells c = 1..N: the symmetrized k-nearest-neighbour sets in Hamming distance, ties broken by
    index, as {c: sorted list of cells}. a is in P(c) iff c is in P(a)."""
    n_cells = len(detected)
    distances = (detected[:, None, :] != detected[None, :, :]).sum(axis=2)
    np.fill_diagonal(distances, detected.shape[1] + 1)
    nearest = np.argsort(distances, axis=1, kind="stable")[:, :min(k, n_cells - 1)]
    neighbours = {c: set() for c in range(1, n_cells + 1)}
    for row, columns in enumerate(nearest):
        for column in columns:
            neighbours[row + 1].add(int(column) + 1)
            neighbours[int(column) + 1].add(row + 1)
    return {c: sorted(cells) for c, cells in neighbours.items()}


# ---------------------------------------------------------------------------------------------- subprograms

def add_threshold_functions(model, candidate_inputs, max_indegree, forbid_self_repression=True):
    """Per node j: (positives, negatives, threshold, unused), in the shape attractor_learning.ilp.add_path_to_model
    takes, with the edge binaries in candidate_inputs[j]'s order. Constraints as in the manuscript's Sections
    3.3.1-3.3.2 and infer_unknown_topology_symmetric: p + q + u <= 1 per edge, 1 - u <= degree <= d (1 - u),
    1 - u <= theta <= degree. Returns (functions, edge terms p + q)."""
    functions, edges = [], []
    for j, candidates in enumerate(candidate_inputs):
        positives = [model.addVar(vtype=GRB.BINARY, name="p_{}_{}".format(i, j)) for i in candidates]
        negatives = [model.addVar(vtype=GRB.BINARY, name="q_{}_{}".format(i, j)) for i in candidates]
        threshold = model.addVar(lb=0, ub=min(len(candidates), max_indegree), vtype=GRB.INTEGER,
                                 name="theta_{}".format(j))
        unused = model.addVar(vtype=GRB.BINARY, name="u_{}".format(j))
        for i, p, q in zip(candidates, positives, negatives):
            model.addConstr(p + q + unused <= 1, name="sign_{}_{}".format(i, j))
            if forbid_self_repression and i == j:
                # f_j = 1 - x_j makes j unstable at every state, explaining every change of j
                model.addConstr(q == 0, name="no_self_repression_{}".format(j))
            edges.append(p + q)
        degree = gurobipy.quicksum(positives) + gurobipy.quicksum(negatives)
        model.addConstr(degree <= max_indegree * (1 - unused), name="max_indegree_{}".format(j))
        model.addConstr(degree >= 1 - unused, name="unused_{}".format(j))
        model.addConstr(threshold <= degree, name="threshold_le_degree_{}".format(j))
        model.addConstr(threshold >= 1 - unused, name="threshold_ge_1_{}".format(j))
        functions.append((positives, negatives, threshold, unused))
    return functions, edges


def add_states(model, detected):
    """(Y, agreement terms): Y[c] for c = 0..N, binaries, and per cell and gene the agreement C of Y with the
    detection call, which for a constant target is Y itself or 1 - Y (ilp.add_state_agreement_terms). The
    pseudo-root has no data row and so no agreement terms."""
    n_cells, n_genes = detected.shape
    states = [[model.addVar(vtype=GRB.BINARY, name="Y_{}_{}".format(c, j)) for j in range(n_genes)]
              for c in range(n_cells + 1)]
    agreement = []
    for c in range(1, n_cells + 1):
        agreement.extend(ilp.add_state_agreement_terms(model, states[c], detected[c - 1], prefix="C_{}".format(c)))
    return states, agreement


def add_function_outputs(model, graph, states, functions, max_indegree):
    """FUNCTIONS: F(Y[c]) for every state c, as new binaries."""
    return [ilp.add_path_to_model(graph, model, path_len=1, first_state_vars=state, model_f_vars=functions,
                                  max_indegree=max_indegree, name_prefix="F_{}".format(c))[0]
            for c, state in enumerate(states)]


def add_instability(model, states, outputs):
    """INSTABILITY: U[c][j] = 1 iff F(Y[c])_j != Y[c][j], exactly (both directions), since U also relaxes the
    transition constraints."""
    unstable = []
    for c, (state, output) in enumerate(zip(states, outputs)):
        row = []
        for j, (y, f) in enumerate(zip(state, output)):
            u = model.addVar(vtype=GRB.BINARY, name="U_{}_{}".format(c, j))
            model.addConstr(u >= f - y, name="U_{}_{}_ge1".format(c, j))
            model.addConstr(u >= y - f, name="U_{}_{}_ge2".format(c, j))
            model.addConstr(u <= f + y, name="U_{}_{}_le1".format(c, j))
            model.addConstr(u <= 2 - f - y, name="U_{}_{}_le2".format(c, j))
            row.append(u)
        unstable.append(row)
    return unstable


def add_tree(model, predecessors, max_depth):
    """TREE: pi[(a, c)] for a in {0} + P(c), exactly one per cell; tau[c] in [1, max_depth] with
    tau_c >= tau_a + 1 - H (1 - pi_ac) for cell predecessors and tau_c <= 1 + (H - 1)(1 - pi_0c); tau_0 = 0.
    Returns (pi, tau), tau indexed by cell with tau[0] the constant 0."""
    n_cells = len(predecessors)
    tau = [0] + [model.addVar(lb=1, ub=max_depth, vtype=GRB.INTEGER, name="tau_{}".format(c))
                 for c in range(1, n_cells + 1)]
    pi = {}
    for c in range(1, n_cells + 1):
        for a in [PSEUDO_ROOT] + predecessors[c]:
            pi[a, c] = model.addVar(vtype=GRB.BINARY, name="pi_{}_{}".format(a, c))
        model.addConstr(gurobipy.quicksum(pi[a, c] for a in [PSEUDO_ROOT] + predecessors[c]) == 1,
                        name="one_predecessor_{}".format(c))
        for a in predecessors[c]:
            model.addConstr(tau[c] >= tau[a] + 1 - max_depth * (1 - pi[a, c]), name="order_{}_{}".format(a, c))
        model.addConstr(tau[c] <= 1 + (max_depth - 1) * (1 - pi[PSEUDO_ROOT, c]), name="root_time_{}".format(c))
    return pi, tau


def add_transitions(model, pi, states, unstable):
    """TRANSITIONS: if pi[(a, c)] then |Y[c][j] - Y[a][j]| <= U[a][j] for every gene j."""
    for (a, c), indicator in pi.items():
        for j in range(len(states[c])):
            model.addConstr(states[c][j] - states[a][j] <= unstable[a][j] + 1 - indicator,
                            name="transition_{}_{}_{}_up".format(a, c, j))
            model.addConstr(states[a][j] - states[c][j] <= unstable[a][j] + 1 - indicator,
                            name="transition_{}_{}_{}_down".format(a, c, j))


def build_model(detected, genes, regulators, k=20, max_indegree=None, max_depth=None, edge_penalty=1.0,
                stability_penalty=0.5, root_penalty=None, forbid_self_repression=True):
    """The full program of joint_orientation_ilp.tex: maximize
        sum C / (N n) - (beta ln N / (N n)) sum (p + q) - (gamma / (N n)) sum U - (rho / N) sum pi_0c
    so an edge costs beta ln N changed data values (a BIC-like growth with the number of cells), with
    beta = edge_penalty, gamma = stability_penalty, rho = root_penalty (default 6 / n: a cell is made a
    root only if attaching it costs more than about six flipped bits). max_indegree defaults to the number of
    regulators, max_depth (H) to N. Returns (model, handles), handles a dict of the variables by name."""
    n_cells, n_genes = detected.shape
    regulator_indices = [genes.index(gene) for gene in regulators if gene in genes]
    candidate_inputs = [sorted(set(regulator_indices) | {j}) for j in range(n_genes)]
    max_indegree = max_indegree or len(regulator_indices)
    max_depth = max_depth or n_cells
    root_penalty = 6.0 / n_genes if root_penalty is None else root_penalty
    predecessors = candidate_predecessors(detected, k)
    # the candidate network, for add_path_to_model, which pairs a node's edge variables with its predecessors
    graph = Network(vertex_names=list(genes),
                    edges=[(genes[i], genes[j]) for j in range(n_genes) for i in candidate_inputs[j]],
                    vertex_functions=[None] * n_genes)
    for j, vertex in enumerate(graph.vertices):
        assert [u.index for u in vertex.predecessors()] == candidate_inputs[j]

    model = gurobipy.Model("snapshot_tree_ilp")
    functions, edges = add_threshold_functions(model, candidate_inputs, max_indegree, forbid_self_repression)
    states, agreement = add_states(model, detected)
    outputs = add_function_outputs(model, graph, states, functions, max_indegree)
    unstable = add_instability(model, states, outputs)
    pi, tau = add_tree(model, predecessors, max_depth)
    add_transitions(model, pi, states, unstable)

    bits = float(n_cells * n_genes)
    objective = (gurobipy.quicksum(agreement) / bits
                 - edge_penalty * np.log(n_cells) / bits * gurobipy.quicksum(edges)
                 - stability_penalty / bits * gurobipy.quicksum(u for row in unstable for u in row)
                 - root_penalty / n_cells * gurobipy.quicksum(pi[PSEUDO_ROOT, c] for c in range(1, n_cells + 1)))
    model.setObjective(objective, GRB.MAXIMIZE)
    model.update()
    handles = {"genes": list(genes), "candidate_inputs": candidate_inputs, "predecessors": predecessors,
               "functions": functions, "states": states, "outputs": outputs, "unstable": unstable,
               "pi": pi, "tau": tau}
    return model, handles


def start_from_collapse(handles, detected):
    """A feasible MIP start: every cell under the pseudo-root with one shared latent state (the per-gene
    majority of the detection calls), and every node unused (holding its value), so nothing is unstable."""
    majority = (detected.mean(axis=0) >= 0.5).astype(int)
    for state in handles["states"]:
        for y, value in zip(state, majority):
            y.Start = int(value)
    for positives, negatives, threshold, unused in handles["functions"]:
        for var in positives + negatives:
            var.Start = 0
        threshold.Start, unused.Start = 0, 1
    for (a, c), indicator in handles["pi"].items():
        indicator.Start = int(a == PSEUDO_ROOT)
    for tau in handles["tau"][1:]:
        tau.Start = 1


# ---------------------------------------------------------------------------------------------- solution

def extract_solution(handles, names):
    """The solved program as plain data: the network (Network with SymmetricThresholdFunctions; a node that
    uses no input has no function, i.e. holds its value), and per cell its predecessor (None for the
    pseudo-root), pseudotime, latent state and unstable genes."""
    genes = handles["genes"]
    network = Network(vertex_names=genes, edges=[], vertex_functions=[None] * len(genes))
    for j, (positives, negatives, threshold, _) in enumerate(handles["functions"]):
        signs = [int(round(p.X - q.X)) for p, q in zip(positives, negatives)]
        for i, sign in zip(handles["candidate_inputs"][j], signs):
            if sign:
                network.edges.append((network.vertices[i], network.vertices[j]))
        signs = [sign for sign in signs if sign]
        network.vertices[j].function = SymmetricThresholdFunction(signs, int(round(threshold.X))) if signs else None
    cells = {}
    for (a, c), indicator in handles["pi"].items():
        if indicator.X > 0.5:
            state = [int(round(y.X)) for y in handles["states"][c]]
            cells[names[c - 1]] = {
                "predecessor": None if a == PSEUDO_ROOT else names[a - 1],
                "pseudotime": int(round(handles["tau"][c].X)),
                "state": state,
                "unstable": [genes[j] for j, u in enumerate(handles["unstable"][c]) if u.X > 0.5]}
    root_state = [int(round(y.X)) for y in handles["states"][PSEUDO_ROOT]]
    return network, cells, root_state


def read_config(path):
    """{option: value} from a config file in the format of configs_for_experiments/inference: one
    "name = value" per line, value a Python literal, # comments. Names are the command-line options with
    underscores (e.g. stability_penalty)."""
    import ast
    values = {}
    with open(path) as config_file:
        for line in config_file:
            line = line.split("#", 1)[0].strip()
            if line:
                name, value = (part.strip() for part in line.split("=", 1))
                values[name] = ast.literal_eval(value)
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=None,
                        help="config file whose values replace the defaults below; options given on the command "
                             "line override it")
    parser.add_argument("--data-dir", default=DATA_DIR,
                        help="folder holding hspc_qpcr_processed.csv, hspc_qpcr_raw.csv, binary_expression_LMPP.txt "
                             "(for the regulator names) and LMPP_trajectory_order.txt (for --cells lmpp)")
    parser.add_argument("--cells", choices=["all", "lmpp"], default="lmpp",
                        help="all 2,167 cells, or Hamey et al.'s 531-cell HSC -> LMPP branch")
    parser.add_argument("--max-cells", type=int, default=None, help="keep only the first this many cells (testing)")
    parser.add_argument("--k", type=int, default=20, help="nearest neighbours per cell before symmetrizing")
    parser.add_argument("--edge-penalty", type=float, default=1.0,
                        help="beta: an edge costs beta ln N changed data values")
    parser.add_argument("--stability-penalty", type=float, default=0.5, help="gamma, per unstable gene and state")
    parser.add_argument("--root-penalty", type=float, default=None, help="rho, per root (default 6 / n)")
    parser.add_argument("--max-indegree", type=int, default=None, help="d (default: all regulators)")
    parser.add_argument("--timeout", type=float, default=None, help="Gurobi time limit in seconds")
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--build-only", action="store_true", help="build the model and print its size")
    parser.add_argument("--heuristic-starts", type=int, default=5,
                        help="starts of the heuristic warm start (snapshot_tree_heuristic.solve)")
    parser.add_argument("--warm-start", choices=["heuristic", "collapse"], default="collapse",
                        help="MIP start: the START-BIT tree heuristic (snapshot_tree_heuristic.py), or every cell "
                             "collapsed into one state under the pseudo-root")
    parser.add_argument("--out", default=os.path.join("single_cell_inference", "results"),
                        help="output directory for the solution")
    options = parser.parse_args()
    if options.config:
        config = read_config(options.config)
        unknown = set(config) - set(vars(options))
        if unknown:
            raise ValueError("{}: unknown options {}".format(options.config, sorted(unknown)))
        parser.set_defaults(**config)
        options = parser.parse_args()

    genes, names, types, detected = load_detection_data(options.data_dir)
    if options.cells == "lmpp":
        rows = [names.index(name) for name in lmpp_cells(options.data_dir)]
        names, types, detected = [names[r] for r in rows], [types[r] for r in rows], detected[rows]
    if options.max_cells:
        names, types, detected = names[:options.max_cells], types[:options.max_cells], detected[:options.max_cells]
    regulators = hamey_regulators(options.data_dir)

    started = time.time()
    model, handles = build_model(detected, genes, regulators, k=options.k, max_indegree=options.max_indegree,
                                 edge_penalty=options.edge_penalty, stability_penalty=options.stability_penalty,
                                 root_penalty=options.root_penalty)
    print("{} cells x {} genes ({} regulators), k={} (mean |P(c)| {:.1f}): {} variables, {} constraints, "
          "built in {:.0f}s".format(len(names), len(genes), len(set(regulators) & set(genes)), options.k,
                                    np.mean([len(p) for p in handles["predecessors"].values()]),
                                    model.NumVars, model.NumConstrs, time.time() - started))
    if options.build_only:
        return
    if options.warm_start == "heuristic":
        import snapshot_tree_heuristic  # imports this module, so not at the top
        fitted, parent, states, pseudotime, value, terms = snapshot_tree_heuristic.solve(
            detected, genes, regulators, options.k, options.edge_penalty, options.stability_penalty,
            options.root_penalty, options.max_indegree, starts=options.heuristic_starts)
        snapshot_tree_heuristic.set_mip_start(handles, fitted, parent, states, pseudotime)
        print("heuristic start: objective {:.5f}, {}".format(value, terms))
    else:
        start_from_collapse(handles, detected)
    if options.timeout:
        model.Params.TimeLimit = options.timeout
    if options.threads:
        model.Params.Threads = options.threads
    model.optimize()
    if model.SolCount == 0:
        raise RuntimeError("no solution found (status {})".format(model.Status))

    network, cells, root_state = extract_solution(handles, names)
    os.makedirs(options.out, exist_ok=True)
    network.save(os.path.join(options.out, "network.json"))
    with open(os.path.join(options.out, "cells.json"), "w") as cells_file:
        json.dump({"genes": genes, "pseudo_root_state": root_state,
                   "cells": {name: dict(cells[name], cell_type=cell_type) for name, cell_type in zip(names, types)},
                   "objective": model.ObjVal, "bound": model.ObjBound, "options": vars(options)},
                  cells_file, indent=1)
    roots = sum(cell["predecessor"] is None for cell in cells.values())
    flips = int(sum((np.array(cells[name]["state"]) != row).sum() for name, row in zip(names, detected)))
    print("objective {:.5f} (bound {:.5f}); {} edges, {} roots, {} flipped bits, max pseudotime {} -> {}".format(
        model.ObjVal, model.ObjBound, len(network.edges), roots, flips,
        max(cell["pseudotime"] for cell in cells.values()), options.out))


if __name__ == "__main__":
    main()
