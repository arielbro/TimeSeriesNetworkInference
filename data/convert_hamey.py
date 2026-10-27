"""Turn the Hamey et al. 2017 LMPP pseudotime data into this repository's generated-data format: one
two-row matrix per input-output pair, with their partial-correlation network as the scaffold.

Source: Hamey et al. 2017, "Reconstructing blood stem cell regulatory network models from single-cell
molecular profiles", PNAS 114:5822-5829, doi:10.1073/pnas.1610609114 (main text and SI in documents/).
The three input files are copied unchanged from github.com/fionahamey/Pseudotime-network-inference into
data/hamey_pseudotime_inference:

    binary_expression_LMPP.txt   31 TFs x 531 cells of the HSC -> LMPP branch, already binary
    LMPP_trajectory_order.txt    the cells in Wanderlust pseudotime order, one name per line
    startingNetworkParCor.txt    "from.gene relation to.gene", relation 1 (activation) or -1 (repression):
                                 the top 100 partial correlations between TFs plus a self-activation for
                                 every gene - the candidate edges their rule search is restricted to

Input-output pairs follow the paper's Network Construction section (and booleanRules.py): for ordered
cells c_0..c_{n-1} and step k, pair t (k + 1 <= t <= n - 1) has input mode(c_{t-k-1}, c_{t-k}, c_{t-k+1})
per gene and output c_t. Consecutive pairs are not chained (pair t's input is not pair t-1's output), so
each pair is its own two-row time series rather than a row of one long one.

Pairs are held out by 5-fold cross-validation, each fold its own data directory and so its own
experiment in the inference output: <out>/contiguous_fold<k>/lmpp holds out the k-th of five consecutive
blocks of the pseudotime-sorted pairs, <out>/random_fold<k>/lmpp the k-th fold of a fixed random
permutation (seed 0). Neighbouring pairs share input cells, so random folds put near-copies of the test
inputs in training; contiguous folds leave that only at the block edges. Every network directory holds
the files run_inference_on_data.py reads:

    source_model.txt, conversion_notes.json
    scaffold_network.json   the partial-correlation network, unsigned (the scaffold format has no signs).
                            Self-loops included
    scaffold_signs.json     [[regulator, target, 1 | -1], ...], the network's signs, read by the hamey
                            inference method
    true_network.json       a copy of the scaffold. There is no ground truth: it is here because
                            run_inference_on_data.py finds networks by this file, and the edge scores it
                            yields only measure overlap with the scaffold
    real_matrices.npz       "t0004", "t0005", ... = the pair whose output is cell t (2 x 31, binary).
                            Keys sort in pseudotime order, so a contiguous split can be taken from them
    noisy_matrices.npz      a copy of real_matrices (the measurements have no hidden clean version)
    test_keys.json          the fold's held-out matrix keys, which run_inference_on_data.py uses in
                            place of its random train_size split

Note that booleanRules.py reads the order file with pandas' default header=0, so Hamey et al.'s own run
drops the first cell (HSC_150_043) and its pair indices are shifted by one; this script keeps every cell.
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from attractor_learning.graphs import Network  # noqa: E402

SOURCE_DIR = os.path.join("data", "hamey_pseudotime_inference")
EXPRESSION_FILE = "binary_expression_LMPP.txt"
ORDER_FILE = "LMPP_trajectory_order.txt"
NETWORK_FILE = "startingNetworkParCor.txt"
DEFAULT_OUT = os.path.join("data", "generated", "hamey")
NETWORK_NAME = "lmpp"
N_FOLDS = 5
RANDOM_FOLD_SEED = 0  # the random folds of the local LMPP comparisons
DEFAULT_STEP = 3  # the paper's k; the repository README's example uses 5


def read_table(path):
    """Tab-separated rows of a text file, header included, blank lines skipped."""
    with open(path) as table_file:
        return [line.rstrip("\r\n").split("\t") for line in table_file if line.strip()]


def read_expression(path):
    """(gene names, {cell name: 0/1 array over the genes}). The header row lists the genes only; every
    other row is a cell name followed by one value per gene."""
    rows = read_table(path)
    genes = rows[0]
    cells = {}
    for row in rows[1:]:
        if len(row) != len(genes) + 1:
            raise ValueError("{}: row {} has {} values for {} genes".format(path, row[0], len(row) - 1,
                                                                          len(genes)))
        values = np.array([int(value) for value in row[1:]])
        if not set(values) <= {0, 1}:
            raise ValueError("{}: row {} is not binary".format(path, row[0]))
        cells[row[0]] = values
    return genes, cells


def read_signed_network(path, genes):
    """[(regulator, target, sign)] in file order."""
    rows = read_table(path)
    if rows[0] != ["from.gene", "relation", "to.gene"]:
        raise ValueError("{}: unexpected header {}".format(path, rows[0]))
    edges = []
    for source, relation, target in rows[1:]:
        if source not in genes or target not in genes:
            raise ValueError("{}: edge {} -> {} names a gene not in the expression data".format(
                path, source, target))
        if relation not in ("1", "-1"):
            raise ValueError("{}: edge {} -> {} has relation {}".format(path, source, target, relation))
        edges.append((source, target, int(relation)))
    return edges


def input_output_pairs(ordered, step):
    """{t: (input, output)} for every t with a full window, input the per-gene majority of cells
    t-step-1, t-step, t-step+1 (the mode of three binary values) and output cell t."""
    pairs = {}
    for t in range(step + 1, len(ordered)):
        window = np.vstack(ordered[t - step - 1:t - step + 2])
        pairs[t] = ((window.sum(axis=0) >= 2).astype(int), ordered[t])
    return pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", default=SOURCE_DIR, help="folder holding the three Hamey et al. files")
    parser.add_argument("--out", default=DEFAULT_OUT, help="output data parent directory, one data directory per fold")
    parser.add_argument("--step", type=int, default=DEFAULT_STEP, help="pseudotime step k")
    options = parser.parse_args()

    genes, cells = read_expression(os.path.join(options.source, EXPRESSION_FILE))
    with open(os.path.join(options.source, ORDER_FILE)) as order_file:
        order = [line.strip() for line in order_file if line.strip()]
    if sorted(order) != sorted(cells):
        raise ValueError("the order file and the expression file list different cells")
    ordered = [cells[name] for name in order]
    signed_edges = read_signed_network(os.path.join(options.source, NETWORK_FILE), genes)

    pairs = input_output_pairs(ordered, options.step)
    matrices = {"t{:04d}".format(t): np.vstack(pair) for t, pair in pairs.items()}
    unchanged = sum(int((matrix[0] == matrix[1]).sum()) for matrix in matrices.values())
    total = len(matrices) * len(genes)

    scaffold = Network(vertex_names=genes, edges=[(source, target) for source, target, _ in signed_edges])
    keys = sorted(matrices)
    folds = {"contiguous": np.array_split(np.arange(len(keys)), N_FOLDS),
             "random": np.array_split(np.random.default_rng(RANDOM_FOLD_SEED).permutation(len(keys)), N_FOLDS)}

    notes = {
        "source_data": options.source.replace(os.sep, "/"),
        "step": options.step,
        "pairs": "input = per-gene mode of cells t-{0}-1, t-{0}, t-{0}+1; output = cell t".format(options.step),
        "genes": genes,
        "cell_order": order,
        "matrix_keys": {key: {"output_cell": order[t],
                              "input_cells": order[t - options.step - 1:t - options.step + 2]}
                        for key, t in zip(matrices, pairs)},
        "scaffold": "{} edges ({} self-loops), from {}".format(
            len(signed_edges), sum(source == target for source, target, _ in signed_edges), NETWORK_FILE),
        "scaffold_signs": [[source, target, sign] for source, target, sign in signed_edges],
        "true_network": "copy of the scaffold; there is no ground truth for this data",
        "output_equals_input": [unchanged, total],
    }
    lines = ["Hamey et al. 2017 HSC -> LMPP pseudotime branch, binary qPCR",
             "source: {} (github.com/fionahamey/Pseudotime-network-inference)".format(notes["source_data"]),
             "{} genes, {} ordered cells, {} input-output pairs with step {}".format(
                 len(genes), len(order), len(matrices), options.step),
             "scaffold: {}".format(notes["scaffold"]),
             "true_network.json is a copy of the scaffold (no ground truth)",
             "output equals input on {}/{} gene values (persistence baseline)".format(unchanged, total)]
    for split, parts in folds.items():
        for k, test in enumerate(parts):
            fold = "{} fold {} of {}".format(split, k, N_FOLDS)
            network_dir = os.path.join(options.out, "{}_fold{}".format(split, k), NETWORK_NAME)
            os.makedirs(network_dir, exist_ok=True)
            scaffold.save(os.path.join(network_dir, "scaffold_network.json"))
            scaffold.save(os.path.join(network_dir, "true_network.json"))
            with open(os.path.join(network_dir, "scaffold_signs.json"), "w") as signs_file:
                json.dump([[source, target, sign] for source, target, sign in signed_edges], signs_file)
            with open(os.path.join(network_dir, "test_keys.json"), "w") as test_keys_file:
                json.dump(sorted(keys[i] for i in test), test_keys_file)
            np.savez(os.path.join(network_dir, "real_matrices"), **matrices)
            np.savez(os.path.join(network_dir, "noisy_matrices"), **matrices)
            with open(os.path.join(network_dir, "conversion_notes.json"), "w") as notes_file:
                json.dump(dict(notes, held_out=fold), notes_file, indent=1)
            with open(os.path.join(network_dir, "source_model.txt"), "w") as source_file:
                source_file.write("\n".join(lines + ["held out: " + fold]) + "\n")
    print("\n".join(lines))
    print("{} folds each of contiguous and random splits under {}".format(N_FOLDS, options.out))


if __name__ == "__main__":
    main()
