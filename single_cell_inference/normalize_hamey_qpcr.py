"""Run from the repository root.

Continuous [0, 1] expression for the Hamey et al. 2017 single-cell qPCR, and its agreement with the
published LMPP binarisation.

Source: blood.stemcells.cam.ac.uk/single_cell_qpcr.html, the primary bone marrow HSPC files (2,167 cells
that passed quality control, 48 genes), downloaded into single_cell_inference/data if missing:

    hspc_qpcr_processed.csv   "Cell_type,Name,<genes>": delta Ct, each gene's Ct minus a per-cell housekeeping
                              offset (log2 cycles, lower = more expression). Undetected wells hold one fixed
                              value per gene
    hspc_qpcr_raw.csv         "Name,Cell_type,<genes>": raw Ct, 40 = undetected; also holds control wells
                              (names with NTC, 10X, 20X) and cells that failed quality control

Scaling, per gene over all 2,167 cells: a detected cell (raw Ct < 40) gets
(max dCt - dCt) / (max dCt - min dCt) over the gene's detected cells, so the most expressing cell is 1 and the
least expressing detected cell is 0; an undetected cell gets 0. A gene detected in no cell or at one level only
is set to 0 or 1 by detection alone.

Output: single_cell_inference/data/continuous_expression_hspc.txt, tab separated in the layout of binary_expression_LMPP.txt (header row
of genes, then a cell name and one value per gene), every cell and gene of the processed file. The LMPP files
used for the agreement check are read from the same folder (--source).
"""

import argparse
import os
import sys
import urllib.request

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data"))
from convert_hamey import EXPRESSION_FILE, ORDER_FILE, read_expression  # noqa: E402

URL = "http://blood.stemcells.cam.ac.uk/qpcr_data/"
PROCESSED_FILE = "hspc_qpcr_processed.csv"
RAW_FILE = "hspc_qpcr_raw.csv"
OUT_FILE = "continuous_expression_hspc.txt"
DATA_DIR = os.path.join("single_cell_inference", "data")
UNDETECTED_CT = 40
THRESHOLD = 0.5


def read_csv(path):
    """Rows of a comma-separated file whose lines may end in a bare carriage return (as these do)."""
    with open(path, newline="") as csv_file:
        text = csv_file.read().replace("\r\n", "\n").replace("\r", "\n")
    return [line.split(",") for line in text.split("\n") if line.strip()]


def fetch(name, folder):
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, name)
    if not os.path.exists(path):
        print("downloading", URL + name)
        urllib.request.urlretrieve(URL + name, path)
    return path


def scale(processed_path, raw_path):
    """(genes, cell names, cell types, cells x genes array in [0, 1])."""
    processed = read_csv(processed_path)
    if processed[0][:2] != ["Cell_type", "Name"]:
        raise ValueError("{}: unexpected header {}".format(processed_path, processed[0][:2]))
    genes = processed[0][2:]
    types = [row[0] for row in processed[1:]]
    names = [row[1] for row in processed[1:]]
    delta = np.array([[float(value) for value in row[2:]] for row in processed[1:]])

    raw = read_csv(raw_path)
    if raw[0][:2] != ["Name", "Cell_type"]:
        raise ValueError("{}: unexpected header {}".format(raw_path, raw[0][:2]))
    raw_rows = {row[0]: row[2:] for row in raw[1:]}
    raw_columns = [raw[0][2:].index(gene) for gene in genes]
    missing = [name for name in names if name not in raw_rows]
    if missing:
        raise ValueError("{} processed cells are not in {}, e.g. {}".format(len(missing), raw_path, missing[:3]))
    ct = np.array([[float(raw_rows[name][j]) for j in raw_columns] for name in names])
    detected = ct < UNDETECTED_CT

    expression = np.zeros_like(delta)
    for j in range(len(genes)):
        values = delta[detected[:, j], j]
        if values.size == 0:
            continue
        low, high = values.min(), values.max()
        expression[detected[:, j], j] = (high - values) / (high - low) if high > low else 1.0
    return genes, names, types, expression


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", default=DATA_DIR, help="folder holding the qPCR files and the output")
    parser.add_argument("--source", default=DATA_DIR, help="folder holding the Hamey et al. LMPP files")
    options = parser.parse_args()

    genes, names, types, expression = scale(fetch(PROCESSED_FILE, options.data), fetch(RAW_FILE, options.data))
    out_path = os.path.join(options.data, OUT_FILE)
    with open(out_path, "w") as out_file:
        out_file.write("\t".join(genes) + "\n")
        for name, row in zip(names, expression):
            out_file.write(name + "\t" + "\t".join("{:.4f}".format(value) for value in row) + "\n")
    print("{} cells x {} genes ({}) -> {}".format(len(names), len(genes), ", ".join(sorted(set(types))), out_path))

    binary_genes, binary_cells = read_expression(os.path.join(options.source, EXPRESSION_FILE))
    with open(os.path.join(options.source, ORDER_FILE)) as order_file:
        order = [line.strip() for line in order_file if line.strip()]
    rows = [names.index(name) for name in order]
    columns = [genes.index(gene) for gene in binary_genes]
    published = np.array([binary_cells[name] for name in order])
    continuous = expression[np.ix_(rows, columns)]
    thresholded = (continuous >= THRESHOLD).astype(int)
    detected = (continuous > 0).astype(int)
    print("LMPP branch, {} cells x {} genes, agreement with {}:".format(len(order), len(binary_genes), EXPRESSION_FILE))
    print("  threshold {}: {:.1%}   (detection, value > 0: {:.1%})".format(
        THRESHOLD, (thresholded == published).mean(), (detected == published).mean()))
    print("  per gene at {} (gene, agreement, published on-fraction, thresholded on-fraction):".format(THRESHOLD))
    for j in np.argsort((thresholded == published).mean(axis=0)):
        print("    {:9s} {:.2f}  {:.2f}  {:.2f}".format(binary_genes[j], (thresholded[:, j] == published[:, j]).mean(),
                                                    published[:, j].mean(), thresholded[:, j].mean()))


if __name__ == "__main__":
    main()
