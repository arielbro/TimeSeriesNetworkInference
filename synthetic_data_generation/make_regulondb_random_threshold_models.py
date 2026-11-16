"""Build 30 random symmetric threshold models on a fixed RegulonDB topology with RegulonDB's signs, usable as a
`graphs_dir`. --network picks the topology:

    full                data/regulondb/signs.tsv (the same gene pairs as network.tsv, duplicates already
                        merged by download.py): every gene that appears as a regulator or target, and every
                        edge. Written to data/regulondb_random_threshold_models.
    strict_regulators   the strict regulator subnetwork that data/regulondb/make_inference_data.py
                        --subnetwork strict_regulators uses as its true network: from gse20305_signs.tsv,
                        repeatedly remove every gene with no outgoing edge to another remaining gene, and keep
                        the edges between the remaining genes, in gse20305_binary.json's gene order. Written
                        to data/regulondb_strict_regulators_random_threshold_models.

Self-loops (autoregulation) are kept. The topology is the same in every model.

For each model, every gene with in-degree d >= 1 gets a symmetric threshold function over its regulators:
  Signs:     RegulonDB's sign where it is 1 (activation) or -1 (repression); where it is "dual" or "unknown",
             a uniformly random sign, drawn independently in each model.
  Threshold: uniform in [1, d], so the function is never constant.
A gene with no regulators is an input node.

Model k uses seed k. Stored as network.json (signs + threshold per node). Regenerate from the repository
root with
    python synthetic_data_generation/make_regulondb_random_threshold_models.py --network full
    python synthetic_data_generation/make_regulondb_random_threshold_models.py --network strict_regulators
"""
import argparse
import importlib.util
import json
import os
import random
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from attractor_learning.graphs import Network  # noqa: E402
from attractor_learning.logic import SymmetricThresholdFunction  # noqa: E402

REGULONDB_DIR = os.path.join(ROOT, "data", "regulondb")
NETWORKS = {
    "full": dict(signs="signs.tsv", output="regulondb_random_threshold_models",
                 title="the full RegulonDB network"),
    "strict_regulators": dict(signs="gse20305_signs.tsv",
                              output="regulondb_strict_regulators_random_threshold_models",
                              title="the RegulonDB strict regulator subnetwork"),
}
N_MODELS = 30


def read_signed_edges(path):
    with open(path, encoding="utf8") as in_file:
        rows = [line.rstrip("\r\n").split("\t") for line in in_file if line.strip()]
    assert rows[0][:3] == ["regulator", "target", "sign"], rows[0]
    signs = {}
    for regulator, target, sign in (row[:3] for row in rows[1:]):
        assert (regulator, target) not in signs, "duplicate pair {}".format((regulator, target))
        assert sign in ("1", "-1", "dual", "unknown"), sign
        signs[(regulator, target)] = sign
    return signs


def random_model(genes, edges, signs, seed):
    rng = random.Random(seed)
    topology = Network(vertex_names=genes, edges=edges)
    functions, drawn = [], 0
    for vertex in topology.vertices:
        regulators = vertex.predecessors()       # the order Network evaluates a function's inputs in
        if not regulators:
            functions.append(None)
            continue
        vertex_signs = []
        for regulator in regulators:
            sign = signs[(regulator.name, vertex.name)]
            if sign in ("1", "-1"):
                vertex_signs.append(sign == "1")
            else:
                vertex_signs.append(rng.choice([True, False]))
                drawn += 1
        functions.append(SymmetricThresholdFunction(signs=vertex_signs,
                                                    threshold=rng.randint(1, len(regulators))))
    return Network(vertex_names=genes, edges=edges, vertex_functions=functions), drawn


def strict_regulator_genes(edges):
    """make_inference_data.kept_genes(edges, "strict_regulators"), so both use the same gene set."""
    spec = importlib.util.spec_from_file_location(
        "make_inference_data", os.path.join(REGULONDB_DIR, "make_inference_data.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.kept_genes(edges, "strict_regulators")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--network", choices=sorted(NETWORKS), required=True)
    options = parser.parse_args()
    config = NETWORKS[options.network]
    output_dir = os.path.join(ROOT, "data", config["output"])
    signs_path = os.path.join(REGULONDB_DIR, config["signs"])
    signs = read_signed_edges(signs_path)
    if options.network == "full":
        edges = sorted(signs)
        genes = sorted({gene for edge in edges for gene in edge})
    else:
        kept = strict_regulator_genes(sorted(signs))
        with open(os.path.join(REGULONDB_DIR, "gse20305_binary.json")) as meta_file:
            genes = [gene for gene in json.load(meta_file)["genes"] if gene in kept]
        signs = {(u, v): sign for (u, v), sign in signs.items() if u in kept and v in kept}
        edges = sorted(signs)
    in_degree = {gene: 0 for gene in genes}
    out_degree = {gene: 0 for gene in genes}
    for u, v in edges:
        out_degree[u] += 1
        in_degree[v] += 1
    sign_counts = {s: sum(1 for x in signs.values() if x == s) for s in ("1", "-1", "dual", "unknown")}
    print("{} genes, {} edges ({} self-loops); signs {}".format(
        len(genes), len(edges), sum(u == v for u, v in edges), sign_counts), flush=True)

    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    rows = []
    for seed in range(N_MODELS):
        name = "model_{:02d}".format(seed)
        network, drawn = random_model(genes, edges, signs, seed)
        model_dir = os.path.join(output_dir, name)
        os.makedirs(model_dir)
        network.save(os.path.join(model_dir, Network.MODEL_JSON_NAME))

        # the folder is only useful if the pipeline's readers load it back unchanged, with the fixed signs
        reloaded = Network.parse_model_dir(model_dir)
        assert len(reloaded) == len(genes) and len(reloaded.edges) == len(edges)
        for vertex in reloaded.vertices:
            if vertex.function is None:
                continue
            for regulator, sign in zip(vertex.predecessors(), vertex.function.signs):
                fixed = signs[(regulator.name, vertex.name)]
                if fixed in ("1", "-1"):
                    assert sign == (fixed == "1"), (name, regulator.name, vertex.name)
        check_rng = random.Random(seed)
        for _ in range(20):
            state = [check_rng.randint(0, 1) for _ in range(len(genes))]
            assert [int(bool(x)) for x in reloaded.next_state(state)] == \
                   [int(bool(x)) for x in network.next_state(state)], name
        rows.append("| {} | {} | {} |".format(name, seed, drawn))
        print("{}: {} dual/unknown signs drawn".format(name, drawn), flush=True)

    lines = ["# Random threshold models on {}".format(config["title"]), "",
             "A `graphs_dir` of {} models sharing one RegulonDB topology and its signs, differing only in "
             "their thresholds (and in the signs RegulonDB leaves open). Written by "
             "`synthetic_data_generation/make_regulondb_random_threshold_models.py --network {}`.".format(
                 N_MODELS, options.network), "",
             "- Source: `data/regulondb/{}`{}: {} genes, {} edges, {} of them self-loops (kept); {} input "
             "nodes, largest in-degree {}, largest out-degree {}.".format(
                 config["signs"], "" if options.network == "full" else ", restricted to the strict regulators "
                 "as `data/regulondb/make_inference_data.py --subnetwork strict_regulators` does",
                 len(genes), len(edges), sum(u == v for u, v in edges),
                 sum(1 for d in in_degree.values() if d == 0), max(in_degree.values()),
                 max(out_degree.values())),
             "- Signs: RegulonDB's for the {} activating and {} repressing edges; a uniformly random sign, "
             "redrawn in each model, for the {} dual and {} unknown ones.".format(
                 sign_counts["1"], sign_counts["-1"], sign_counts["dual"], sign_counts["unknown"]),
             "- Thresholds: uniform in [1, in-degree] for every gene with a regulator, so no function is "
             "constant. Genes with no regulator are input nodes.",
             "- Model k uses seed k. Each `model_XX/network.json` stores signs + threshold per node.", "",
             "| model | seed | dual/unknown signs drawn |", "|---|---|---|"] + rows
    with open(os.path.join(output_dir, "README.md"), "w", encoding="utf-8") as readme:
        readme.write("\n".join(lines) + "\n")
    print("wrote {} models to {}".format(N_MODELS, output_dir))


if __name__ == "__main__":
    main()
