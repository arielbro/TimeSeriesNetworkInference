import os
from attractor_learning import graphs
import numpy as np
import shutil
import logging
import time
import datetime
import enum
import traceback
from concurrent.futures import ProcessPoolExecutor
from inference import dummy_inference, binary_inference_ideas, benchmark_inference
from inference.binary_inference_ideas import infer_known_topology_symmetric, infer_known_topology_general, infer_unknown_topology_symmetric
from sklearn.model_selection import train_test_split
from validation import inference_scoring
import time
import re
import itertools
import functools
import json
import sys
import random
import subprocess
import configargparse


def process_network(network_name, network_path, output_parent_dir, kwargs):
    inference_method = kwargs['inference_method']
    network_out_dir = os.path.join(output_parent_dir, network_name)
    os.makedirs(network_out_dir, exist_ok=True)
    network_log_file = os.path.join(network_out_dir, "inference_log.txt")
    network_started_iso = datetime.datetime.now().isoformat(timespec='seconds')
    network_wall_t0 = time.time()
    status = 'ok'

    try:
        scaffold_network = graphs.Network.load(os.path.join(network_path, "scaffold_network.json"))
        true_network = graphs.Network.load(os.path.join(network_path, "true_network.json"))
        with np.load(os.path.join(network_path, "noisy_matrices.npz")) as reference_matrices, \
                np.load(os.path.join(network_path, "real_matrices.npz")) as hidden_matrices:
            if kwargs['train_size'] != 1:
                train_keys, test_keys = train_test_split(list(reference_matrices.keys()),
                                                         train_size=kwargs['train_size'])
                reference_train = {i: reference_matrices[i] for i in train_keys}
                real_train = {i: hidden_matrices[i] for i in train_keys}
                reference_test = {i: reference_matrices[i] for i in test_keys}
                real_test = {i: hidden_matrices[i] for i in test_keys}
            else:
                reference_test = dict()
                reference_train = {i: reference_matrices[i] for i in reference_matrices.keys()}
                real_test = dict()
                real_train = {i: hidden_matrices[i] for i in hidden_matrices.keys()}

        start = time.time()
        inferred_model = inference_method(reference_train.values(), scaffold_network,
                                          timeout_secs=kwargs['model_inference_timeout_secs'],
                                          log_file=network_log_file,
                                          allow_additional_edges=kwargs['allow_additional_edges'],
                                          included_edges_relative_weight=kwargs['included_edges_relative_weight'],
                                          added_edges_relative_weight=kwargs['added_edges_relative_weight'],
                                          allow_input_flips=kwargs['allow_input_flips'],
                                          flip_penalty=kwargs['flip_penalty'],
                                          no_anchoring=kwargs['no_anchoring'],
                                          warm_start_from_scaffold=kwargs['warm_start_from_scaffold'],
                                          warm_start_time_frac=kwargs['warm_start_time_frac'],
                                          max_indegree=kwargs.get('max_indegree', -1),
                                          gurobi_threads=int(os.environ.get('SLURM_CPUS_PER_TASK') or 0))
        time_taken = time.time() - start
        del scaffold_network

        # Truth-table inference (general / random_model / exact_match_else_random / linear_classifier) keeps
        # every scaffold edge even when a node's learned truth table ignores an input. Prune those functionally
        # irrelevant inputs so the saved model, its time-series predictions, and the edge score all reflect the
        # model's effective topology (a node depending on no input keeps one randomly chosen input and stays
        # constant). Threshold-function models are untouched. Every rewrite preserves next_state output, so the
        # time-series predictions below are unchanged; done before the save so inferred_network.json is the
        # pruned graph. Kept out of inference_time above (post-processing, not solving).
        inference_scoring.prune_ignored_inputs(inferred_model)
        inferred_model.save(os.path.join(network_out_dir, "inferred_network.json"))

        for group, reference_data, real_data in zip(['train', 'test'], [reference_train, reference_test], [real_train, real_test]):
            inferred_matrices = {i: np.asarray(inferred_model.next_states(data_matrix[0, :], data_matrix.shape[0] - 1))
                                 for i, data_matrix in reference_data.items()}
            assert len(inferred_matrices) == len(reference_data)
            assert all(pred_mat.shape[0] == reference_data[i].shape[0] for i, pred_mat in inferred_matrices.items())
            np.savez(os.path.join(network_out_dir, "{}_matrices".format(group)), **inferred_matrices)

            # compare against non-noisy data
            pred_vec = np.concatenate([inferred_matrices[i][1:, ].flatten() for i in inferred_matrices.keys()])
            ref_vec = np.concatenate([real_data[i][1:, ].flatten() for i in inferred_matrices.keys()])
            timeseries_score = inference_scoring.sparse_accuracy_score(ref_vec, pred_vec)
            np.save(os.path.join(network_out_dir, "timeseries_real_accuracy_score_{}".format(group)),
                    timeseries_score)
            del ref_vec
            # compare against given reference data
            ref_vec = np.concatenate([reference_data[i][1:, ].flatten() for i in inferred_matrices.keys()])
            timeseries_score = inference_scoring.sparse_accuracy_score(ref_vec, pred_vec)
            np.save(os.path.join(network_out_dir, "timeseries_reference_accuracy_score_{}".format(group)),
                    timeseries_score)
            del ref_vec, pred_vec, inferred_matrices

        del reference_train, real_train, reference_test, real_test

        np.save(os.path.join(network_out_dir, "inference_time"), time_taken)

        y_true, y_pred = inference_scoring.models_to_edge_vectors(true_network, inferred_model, use_sparse=True)
        edge_score = inference_scoring.sparse_jaccard_score(y_true, y_pred)
        np.save(os.path.join(network_out_dir, "edge_accuracy_score"), edge_score)
        del true_network, inferred_model, y_true, y_pred
    except Exception:
        status = 'error'
        with open(network_log_file, 'a', encoding='utf-8') as f:
            f.write("\nERROR processing network={} pid={}\n".format(network_name, os.getpid()))
            f.write(traceback.format_exc())

    return (network_name, network_log_file, network_started_iso,
            time.time() - network_wall_t0, status)


def append_network_log_to_master(master_log_path, network_name, network_log_path,
                                  started_iso, elapsed_s, status):
    with open(master_log_path, 'a', encoding='utf-8') as out:
        out.write("\n>>> begin network={name} started={started} elapsed={elapsed:.2f}s status={status}\n".format(
            name=network_name, started=started_iso, elapsed=elapsed_s, status=status))
        if os.path.exists(network_log_path) and os.path.getsize(network_log_path) > 0:
            with open(network_log_path, 'r', encoding='utf-8') as src:
                shutil.copyfileobj(src, out)
            out.write("\n")
        else:
            out.write("(no per-network log produced)\n")
        out.write("<<< end network={name}\n".format(name=network_name))


def expected_output_filenames(kwargs):
    """The files process_network writes for one network on a successful run. The test-group files only
    exist when there is a test split (train_size != 1)."""
    names = ["inferred_network.json", "inference_time.npy", "edge_accuracy_score.npy",
             "train_matrices.npz", "timeseries_real_accuracy_score_train.npy",
             "timeseries_reference_accuracy_score_train.npy"]
    if kwargs.get('train_size') != 1:
        names += ["test_matrices.npz", "timeseries_real_accuracy_score_test.npy",
                  "timeseries_reference_accuracy_score_test.npy"]
    return names


def network_output_is_complete(network_out_dir, kwargs):
    """True iff network_out_dir already holds every (non-empty) output file inference would write, so the
    network can be skipped on a re-run."""
    return all(os.path.exists(os.path.join(network_out_dir, f)) and
               os.path.getsize(os.path.join(network_out_dir, f)) > 0
               for f in expected_output_filenames(kwargs))


def count_timeseries_matrices(network_path):
    """Number of time-series matrices for a network - the samples process_network splits into train/test
    (it splits on the keys of noisy_matrices.npz, line ~41). Returns None if the file can't be read, leaving
    that network's fate to the runtime error path rather than pre-judging it."""
    try:
        with np.load(os.path.join(network_path, "noisy_matrices.npz")) as data:
            return len(data.files)
    except Exception:
        return None


def split_precheck_reason(n_matrices, train_size):
    """None if process_network's train/test split would succeed for this many time-series matrices and this
    train_size, else the ValueError message explaining why it can't be split. Mirrors runtime exactly by
    invoking the same train_test_split (on a dummy range) rather than reimplementing sklearn's size math, so
    the precheck can't drift from process_network across sklearn versions. train_size==1 does no split."""
    if train_size == 1:
        return None
    try:
        train_test_split(list(range(n_matrices)), train_size=train_size)
        return None
    except ValueError as e:
        return str(e)


# argparse names that control HOW inference is run rather than inference parameters; excluded from the
# per-problem kwargs stored in a manifest.
_CONTROL_ARGS = {"config", "emit_manifest", "run_manifest", "task_id", "chunk_size", "n_processes",
                 "summarize_errors", "array_job_id"}


def parse_bool_option(value):
    """argparse/configargparse `type` for a boolean given as a value rather than a bare flag (`--x False`,
    or `x = False` in a config file). `type=bool` cannot be used for this: it just calls bool() on the
    string, so every non-empty value - "False" and "0" included - comes out True."""
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in ('true', 'yes', 'y', '1'):
        return True
    if normalized in ('false', 'no', 'n', '0'):
        return False
    raise ValueError("expected a boolean (true/false), got {!r}".format(value))


INFERENCE_METHODS = {
    "random": dummy_inference.dummy_inference_method,
    "general": infer_known_topology_general,
    "symmetric": infer_known_topology_symmetric,
    "symmetric_topology": infer_unknown_topology_symmetric,
    "random_model": benchmark_inference.random_model_inference,
    "exact_match_else_random": benchmark_inference.exact_match_else_random_inference,
    "linear_classifier": benchmark_inference.linear_classifier_inference,
    "reveal": benchmark_inference.reveal_inference,
    "best_fit": benchmark_inference.best_fit_inference,
}


def resolve_inference_method(name):
    if name not in INFERENCE_METHODS:
        raise ValueError("Unknown inference method {} (known: {})".format(
            name, ", ".join(sorted(INFERENCE_METHODS))))
    return INFERENCE_METHODS[name]


def inference_method_name(method):
    """The config name of a resolved inference method. Directory names are written with this rather than the
    function's __name__, so a folder reads inference_method=best_fit - the same spelling a config uses -
    instead of best_fit_inference."""
    for name, func in INFERENCE_METHODS.items():
        if func is method:
            return name
    return getattr(method, "__name__", str(method))


# Abbreviations for parameter names in output directory names. The full argparse names make a folder name
# unwieldy (and push toward path-length limits) as soon as a few parameters vary. 'inference_method' is
# deliberately NOT abbreviated: the analysis notebook keys off that exact name when picking a heatmap's
# categorical axis, so shortening it there would silently drop those figures.
COMB_STR_SHORTHANDS = {
    "allow_additional_edges": "addedges",
    "included_edges_relative_weight": "inclw",
    "added_edges_relative_weight": "addw",
    "model_inference_timeout_secs": "timeout",
    "warm_start_from_scaffold": "warmstart",
    "warm_start_time_frac": "warmfrac",
    "allow_input_flips": "inpflips",
    "flip_penalty": "flippen",
    "max_indegree": "maxindeg",
    "train_size": "trainsize",
}


def build_comb_str(options_combination):
    """Directory-name fragment for one inference-parameter combination (inference_method given as the
    function). Parameter names are abbreviated via COMB_STR_SHORTHANDS and the method is written with its
    short config name, so the folder stays readable when several parameters vary."""
    rendered = {}
    for key, value in options_combination.items():
        if key == 'data_parent_dir':
            continue
        if key == 'inference_method':
            value = inference_method_name(value)
        elif isinstance(value, enum.Enum):
            value = value.name
        rendered[COMB_STR_SHORTHANDS.get(key, key)] = value
    comb_str = str(rendered)
    return comb_str.translate(str.maketrans('', '', "'{}")).replace(": ", "=").replace(", ", ",")


def enumerate_problems(options):
    """All single-graph inference problems for a config: the cross product of inference-parameter combinations
    (the varying/appended options), subexperiments (data_dir) and graphs (network_X) in the data folder. Each
    problem is a JSON-serializable dict with the inference kwargs (inference_method kept as its name string),
    the network path/name, and the output_parent_dir - the same layout the whole-grid runner produces.
    Problems whose output is already fully populated are omitted, so the manifest (and the array sized from it)
    covers only missing work; a partially-written output counts as missing and is re-run from scratch.
    Problems that would fail process_network's train/test split (too few time-series matrices for the given
    train_size) are also kept out of the manifest and returned separately as prevalidation errors - this just
    moves detection that already happened at runtime earlier, so no SLURM task is wasted on a certain failure;
    summarize_errors still reports them. Returns (problems, prevalidation_errors)."""
    all_options = {k: v for (k, v) in options._get_kwargs()}
    constant_options = {k: v for k, v in all_options.items() if not isinstance(v, list)}
    variable_options = {k: v for k, v in all_options.items() if isinstance(v, list)}

    problems = []
    prevalidation_errors = []
    matrix_count_cache = {}  # net_path -> matrix count; the count is a property of the data, reused across combos
    for combo_values in itertools.product(*variable_options.values()):
        options_combination = dict(zip(variable_options.keys(), combo_values))
        resolved = dict(options_combination)
        resolved['inference_method'] = resolve_inference_method(resolved['inference_method'])
        comb_str = build_comb_str(resolved)

        kwargs = dict(options_combination)
        kwargs.update(constant_options)
        manifest_kwargs = {k: v for k, v in kwargs.items() if k not in _CONTROL_ARGS}

        data_parent_dir = kwargs['data_parent_dir']
        data_dir_partial_path = os.path.join(*os.path.normpath(data_parent_dir).split(os.sep)[-2:])
        for data_dir in os.scandir(data_parent_dir):
            if not data_dir.is_dir():
                continue
            output_parent_dir = os.path.join("inferred_models", data_dir_partial_path,
                                             "{}-{}".format(comb_str, data_dir.name))
            networks = [(f.name, f.path) for f in os.scandir(data_dir.path) if f.is_dir()]
            if not networks:
                raise ValueError("Did not find network directories in data_dir {}".format(data_dir.name))
            for net_name, net_path in networks:
                if network_output_is_complete(os.path.join(output_parent_dir, net_name), manifest_kwargs):
                    continue  # already fully generated; leave it out of the manifest
                # pre-check the train/test split so a network with too few time-series matrices is excluded
                # from the array (no wasted task) and recorded for the error summary instead. If the matrix
                # count can't be read, fall through to the manifest and let the runtime error path handle it.
                if net_path not in matrix_count_cache:
                    matrix_count_cache[net_path] = count_timeseries_matrices(net_path)
                n_matrices = matrix_count_cache[net_path]
                if n_matrices is not None:
                    reason = split_precheck_reason(n_matrices, manifest_kwargs.get('train_size'))
                    if reason is not None:
                        prevalidation_errors.append({"kwargs": manifest_kwargs,
                                                     "network_name": net_name,
                                                     "network_path": net_path,
                                                     "output_parent_dir": output_parent_dir,
                                                     "n_timeseries_matrices": n_matrices,
                                                     "reason": "train/test split precheck failed: {}".format(reason)})
                        continue
                problems.append({"kwargs": manifest_kwargs,
                                 "network_name": net_name,
                                 "network_path": net_path,
                                 "output_parent_dir": output_parent_dir})
    return problems, prevalidation_errors


def run_single_problem(problem):
    """Run inference for one manifest problem, skipping it if its output is already complete. Returns the
    process_network metadata tuple, or None if skipped."""
    kwargs = dict(problem["kwargs"])
    kwargs["inference_method"] = resolve_inference_method(kwargs["inference_method"])
    output_parent_dir = problem["output_parent_dir"]
    network_name = problem["network_name"]
    if network_output_is_complete(os.path.join(output_parent_dir, network_name), kwargs):
        print("Skipping {}/{}: output already fully populated".format(output_parent_dir, network_name))
        return None
    os.makedirs(output_parent_dir, exist_ok=True)
    return process_network(network_name, problem["network_path"], output_parent_dir, kwargs)


def error_parts_dir(manifest_path):
    """Directory where each array task drops a record of the problems it saw error, for summarize_errors to
    aggregate once the whole array has finished."""
    return manifest_path + ".errparts"


def prevalidation_path(manifest_path):
    """Sidecar file beside the manifest listing problems an emit-time precheck excluded from the array (e.g.
    too few time-series matrices to train/test split). Kept out of the manifest so no SLURM task is created
    for them, and out of error_parts_dir (which submit_inference_array.sh wipes); summarize_errors folds these
    into the error summary. Rewritten on every emit, so a prior submission's records never leak into this one."""
    return manifest_path + ".prevalidation.jsonl"


def summary_state_path(experiment_dir, summary_filename):
    """Durable per-experiment JSONL store backing the human-readable summary, keyed by
    (output_parent_dir, network_name). It accumulates prevalidation/error/incomplete records across runs so a
    later run updates the summary instead of replacing it, and marks an entry resolved (rather than dropping it)
    once its output is complete. Co-located with the summary in the experiment dir, so it persists across
    submissions (unlike the manifest/errparts/sidecar, which live beside the transient manifest)."""
    return os.path.join(experiment_dir, os.path.splitext(summary_filename)[0] + ".state.jsonl")


def run_manifest_chunk(manifest_path, task_id, chunk_size):
    """Run the contiguous slice [task_id*chunk_size, (task_id+1)*chunk_size) of the manifest (one SLURM array
    task). Appends the process_network metadata of each errored problem to a per-task file under
    error_parts_dir() as soon as it happens, so an end-of-array finalize job can collect them and a later kill
    in this chunk doesn't lose errors recorded earlier. Exits non-zero if any problem in the chunk errored, so
    SLURM marks the task failed - a resubmit skips the problems that did complete."""
    with open(manifest_path) as f:
        lines = [line for line in f if line.strip()]
    start = task_id * chunk_size
    end = min(start + chunk_size, len(lines))
    if start >= len(lines):
        print("task_id {} starts past the {} problems in {}; nothing to do".format(task_id, len(lines), manifest_path))
        return
    parts_dir = error_parts_dir(manifest_path)
    os.makedirs(parts_dir, exist_ok=True)
    part_path = os.path.join(parts_dir, "task_{}.jsonl".format(task_id))
    open(part_path, 'w', encoding='utf-8').close()  # truncate any stale file from a prior run of this task id
    had_error = False
    for i in range(start, end):
        problem = json.loads(lines[i])
        print("=== manifest problem {} ({} of {} in this task) ===".format(i, i - start + 1, end - start))
        meta = run_single_problem(problem)
        if meta is not None and meta[4] == 'error':
            had_error = True
            with open(part_path, 'a', encoding='utf-8') as out:  # append now, before any later problem can be killed
                out.write(json.dumps({"manifest_index": i,
                                      "output_parent_dir": problem["output_parent_dir"],
                                      "meta": meta}) + "\n")
    if had_error:
        sys.exit(1)


def query_slurm_task_states(array_job_id):
    """Map array-task index -> "STATE (exit CODE)" via sacct, so the summary can explain why a task was killed
    (e.g. TIMEOUT, OUT_OF_MEMORY, CANCELLED). Returns {} when no job id is given or sacct is unavailable/fails
    (e.g. not running on a SLURM node), so callers degrade gracefully."""
    if not array_job_id:
        return {}
    try:
        completed = subprocess.run(
            ["sacct", "-j", str(array_job_id), "--parsable2", "--noheader",
             "--format=JobID,State,ExitCode"],
            capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return {}
    if completed.returncode != 0:
        return {}
    states = {}
    for line in completed.stdout.splitlines():
        fields = line.split("|")
        if len(fields) < 3:
            continue
        job_id, state, exit_code = fields[0], fields[1], fields[2]
        if "." in job_id:  # skip job substeps like <A>_<a>.batch / .extern
            continue
        # main array-task line is "<A>_<a>"; sacct may collapse pending/cancelled ones into "<A>_[lo-hi]"
        m = re.match(r'^\d+_\[?(\d+)(?:-(\d+))?(?:%\d+)?\]?$', job_id)
        if not m:
            continue
        lo = int(m.group(1))
        hi = int(m.group(2)) if m.group(2) else lo
        for task_index in range(lo, hi + 1):
            states[task_index] = "{} (exit {})".format(state, exit_code)
    return states


def summarize_errors(manifest_path, summary_filename, chunk_size, array_job_id=None):
    """Update an errors summary named summary_filename in each experiment directory the manifest targets
    (inferred_models/<experiment> = the parent of a problem's output_parent_dir). Rather than replacing the
    summary with only this run's findings, it merges into a durable per-experiment store (summary_state_path,
    keyed by (output_parent_dir, network_name)) so results accumulate across submissions. Each entry is one of:
      1. errored - a network that raised a caught exception (recorded by the array tasks), rendered with its
         full per-network log/traceback,
      2. incomplete - a manifest problem still incomplete yet never recorded as an error, i.e. whose array task
         was killed/timed out (or died without a Python error). sacct is queried (via array_job_id, mapping
         manifest index -> task_id = index // chunk_size) for that task's terminal state, or
      3. pre-excluded - a problem an emit-time precheck kept out of the array (prevalidation sidecar), e.g. too
         few time-series matrices to train/test split.
    This run's findings update the matching entries (a network re-run to a different fate is reclassified);
    entries the run didn't touch are retained. Every entry (touched or not) is reconciled against on-disk
    output: one whose output is now complete is marked resolved and kept for history rather than dropped.
    Meant to run once after the whole array has finished (via a SLURM afterany job), and also inline by the
    submit script when there is no array to depend on."""
    chunk_size = max(1, int(chunk_size or 1))
    now_iso = datetime.datetime.now().isoformat(timespec='seconds')

    # full manifest (index -> problem) for kwargs lookup and completeness checks
    problems = []
    with open(manifest_path) as f:
        for index, line in enumerate(f):
            if line.strip():
                problems.append((index, json.loads(line)))
    problem_by_index = {index: problem for index, problem in problems}

    # caught-exception records the array tasks left behind this run
    error_records = []
    recorded_indices = set()
    parts_dir = error_parts_dir(manifest_path)
    if os.path.isdir(parts_dir):
        for entry in os.scandir(parts_dir):
            if not entry.name.endswith(".jsonl"):
                continue
            with open(entry.path, encoding='utf-8') as f:
                for line in f:
                    if line.strip():
                        rec = json.loads(line)
                        error_records.append(rec)
                        recorded_indices.add(rec["manifest_index"])

    # emit-time prevalidation exclusions for this run (never submitted, so absent from the manifest)
    prevalidation_records = []
    pv_path = prevalidation_path(manifest_path)
    if os.path.isfile(pv_path):
        with open(pv_path, encoding='utf-8') as f:
            for line in f:
                if line.strip():
                    prevalidation_records.append(json.loads(line))

    slurm_states = query_slurm_task_states(array_job_id)

    def key_of(output_parent_dir, network_name):
        return (os.path.normpath(output_parent_dir), network_name)

    # this run's findings, keyed by (output_parent_dir, network_name). Categories are disjoint by construction:
    # pre-excluded problems aren't in the manifest, and a recorded error's index is excluded from the incomplete
    # scan below - so a key lands in exactly one bucket.
    current = {}
    for rec in prevalidation_records:
        current[key_of(rec["output_parent_dir"], rec["network_name"])] = {
            "output_parent_dir": rec["output_parent_dir"], "network_name": rec["network_name"],
            "kwargs": rec.get("kwargs", {}), "category": "prevalidation",
            "reason": rec.get("reason", ""), "n_timeseries_matrices": rec.get("n_timeseries_matrices")}
    for rec in error_records:
        net = rec["meta"][0]
        problem = problem_by_index.get(rec["manifest_index"], {})
        current[key_of(rec["output_parent_dir"], net)] = {
            "output_parent_dir": rec["output_parent_dir"], "network_name": net,
            "kwargs": problem.get("kwargs", {}), "category": "error", "meta": rec["meta"]}
    for index, problem in problems:
        if index in recorded_indices:
            continue
        network_out_dir = os.path.join(problem["output_parent_dir"], problem["network_name"])
        if network_output_is_complete(network_out_dir, problem["kwargs"]):
            continue
        task_id = index // chunk_size
        current[key_of(problem["output_parent_dir"], problem["network_name"])] = {
            "output_parent_dir": problem["output_parent_dir"], "network_name": problem["network_name"],
            "kwargs": problem["kwargs"], "category": "incomplete", "task_id": task_id,
            "array_task": "{}_{}".format(array_job_id, task_id) if array_job_id else None,
            "slurm": slurm_states.get(task_id, "unknown (sacct unavailable)")}

    # experiments to (re)write: those this run touched, whether via the manifest or a current finding. An
    # experiment whose entries are ALL untouched this run isn't revisited - the run that last touched it already
    # reconciled it (a problem is in the manifest of the run that completes it, so its resolution is recorded then).
    experiment_dirs = {os.path.dirname(problem["output_parent_dir"]) for _, problem in problems}
    experiment_dirs |= {os.path.dirname(e["output_parent_dir"]) for e in current.values()}

    current_by_exp = {}
    for k, entry in current.items():
        current_by_exp.setdefault(os.path.dirname(entry["output_parent_dir"]), {})[k] = entry

    for experiment_dir in sorted(experiment_dirs):
        state_path = summary_state_path(experiment_dir, summary_filename)
        entries = {}  # key -> stored entry, carried over from prior runs
        if os.path.isfile(state_path):
            with open(state_path, encoding='utf-8') as f:
                for line in f:
                    if line.strip():
                        e = json.loads(line)
                        entries[key_of(e["output_parent_dir"], e["network_name"])] = e

        # overlay this run's findings (reclassify a touched entry, keep its first_seen)
        for k, entry in current_by_exp.get(experiment_dir, {}).items():
            prior = entries.get(k)
            merged = dict(entry)
            merged["first_seen"] = prior.get("first_seen", now_iso) if prior else now_iso
            merged["last_updated"] = now_iso
            entries[k] = merged

        # reconcile every entry against on-disk output: complete now -> resolved (kept, not dropped). A
        # carried-over entry uses its stored kwargs; if train_size (not part of comb_str) is edited in the
        # config mid-experiment, an untouched entry's completeness verdict could lag until the run re-touches it.
        for e in entries.values():
            complete = network_output_is_complete(
                os.path.join(e["output_parent_dir"], e["network_name"]), e.get("kwargs", {}))
            if complete and not e.get("resolved"):
                e["resolved"] = True
                e["resolved_at"] = now_iso
            elif not complete and e.get("resolved"):  # regressed since last resolved
                e["resolved"] = False
                e["resolved_at"] = None
            else:
                e.setdefault("resolved", False)
                e.setdefault("resolved_at", None)

        os.makedirs(experiment_dir, exist_ok=True)
        with open(state_path, 'w', encoding='utf-8') as out:
            for e in entries.values():
                out.write(json.dumps(e) + "\n")

        def in_cat(cat):
            return sorted((e for e in entries.values() if not e.get("resolved") and e["category"] == cat),
                          key=lambda e: (e["output_parent_dir"], e["network_name"]))
        errors, incompletes, prevalids = in_cat("error"), in_cat("incomplete"), in_cat("prevalidation")
        resolved = sorted((e for e in entries.values() if e.get("resolved")),
                          key=lambda e: (e["output_parent_dir"], e["network_name"]))

        summary_path = os.path.join(experiment_dir, summary_filename)
        with open(summary_path, 'w', encoding='utf-8') as out:
            out.write("Error summary for {}: {} errored, {} incomplete, {} pre-excluded, {} resolved "
                      "(updated {})\n".format(
                          manifest_path, len(errors), len(incompletes), len(prevalids), len(resolved), now_iso))
            out.write("\n==== resolved (output now complete; kept for history) ====\n")
            for e in resolved:
                out.write("- {dir}/{net}  [FIXED {when}]  (was {cat})\n".format(
                    dir=e["output_parent_dir"], net=e["network_name"], when=e.get("resolved_at"),
                    cat=e["category"]))
            if not resolved:
                out.write("(none)\n")
            out.write("\n==== pre-excluded problems (failed an emit-time precheck; never submitted) ====\n")
            for e in prevalids:
                out.write("- {dir}/{net}  ({n} time-series matrices)  {reason}\n".format(
                    dir=e["output_parent_dir"], net=e["network_name"],
                    n=e.get("n_timeseries_matrices"), reason=e.get("reason", "")))
            if not prevalids:
                out.write("(none)\n")
            out.write("\n==== incomplete but unrecorded problems (no caught error; likely killed/timed out) ====\n")
            for e in incompletes:
                out.write("- {dir}/{net}  task={task}{arr}  slurm={slurm}\n".format(
                    dir=e["output_parent_dir"], net=e["network_name"], task=e.get("task_id"),
                    arr="" if e.get("array_task") is None else " ({})".format(e.get("array_task")),
                    slurm=e.get("slurm")))
            if not incompletes:
                out.write("(none)\n")
            out.write("\n==== errored networks (caught exceptions; full per-network log below) ====\n")
            if not errors:
                out.write("(none)\n")
        for e in errors:
            name, network_log_path, started_iso, elapsed_s, status = e["meta"]
            append_network_log_to_master(summary_path, name, network_log_path, started_iso, elapsed_s, status)
        print("Wrote {} ({} errored, {} incomplete, {} pre-excluded, {} resolved)".format(
            summary_path, len(errors), len(incompletes), len(prevalids), len(resolved)))


def main():
    p = configargparse.ArgParser(default_config_files=['./config.txt'])
    p.add_argument('--config', is_config_file=True, help='config file path to override defaults')
    p.add_argument('--inference_method', required=False, type=str, action='append')
    p.add_argument('--data_parent_dir', required=False, type=str, action='append')
    p.add_argument('--train_size', required=False, type=float)
    p.add_argument('--model_inference_timeout_secs', required=False, type=float)
    # appendable (like inference_method) so a config can grid-search it - `allow_additional_edges = [False,
    # True]` runs both - and so it always lands in the varying options, hence in the output directory name.
    # default is applied after parsing, since argparse appends onto a list default instead of replacing it.
    p.add_argument('--allow_additional_edges', required=False, default=None, type=parse_bool_option,
                   action='append')
    p.add_argument('--included_edges_relative_weight', required=False, type=float, action='append')
    p.add_argument('--added_edges_relative_weight', required=False, type=float, action='append')
    p.add_argument('--allow_input_flips', required=False, default=False, type=bool)
    p.add_argument('--flip_penalty', required=False, default=1.0, type=float)
    p.add_argument('--no_anchoring', required=False, default=False, type=bool)
    p.add_argument('--warm_start_from_scaffold', required=False, default=False, type=bool)
    p.add_argument('--warm_start_time_frac', required=False, default=0.2, type=float)
    # max in-degree bound; only used by reveal / best_fit (regulator-set size cap) and symmetric_topology
    # (per-node MILP constraint). -1 = use the scaffold network's max in-degree, computed per network.
    p.add_argument('--max_indegree', required=False, default=-1, type=int)
    p.add_argument('--n_processes', required=False, type=int, default=1)
    # SLURM-array modes (each exits after running; otherwise a normal whole-grid in-process run happens):
    p.add_argument('--emit_manifest', required=False, type=str, default=None,
                   help='enumerate all single-graph problems for the config, write them (one JSON per line) '
                        'to this path, print the count, and exit')
    p.add_argument('--run_manifest', required=False, type=str, default=None,
                   help='run a chunk of problems from a manifest produced by --emit_manifest')
    p.add_argument('--task_id', required=False, type=int, default=0,
                   help='SLURM array task id; selects which chunk to run with --run_manifest')
    p.add_argument('--chunk_size', required=False, type=int, default=1,
                   help='number of manifest problems each array task runs')
    p.add_argument('--summarize_errors', required=False, type=str, default=None,
                   help='aggregate the per-task error records for the manifest given by --run_manifest into an '
                        'errors-only summary file with this name, written into each targeted experiment dir '
                        '(inferred_models/<experiment>); also lists incomplete-but-unrecorded problems, then '
                        'exit (run as an end-of-array finalize job)')
    p.add_argument('--array_job_id', required=False, type=str, default=None,
                   help='SLURM job id of the inference array, used by --summarize_errors to query sacct for the '
                        'terminal state of tasks that left problems incomplete')
    options = p.parse_args()
    if options.allow_additional_edges is None:
        options.allow_additional_edges = [False]  # kept a list so it stays a (single-valued) varying option

    if options.summarize_errors is not None:
        if options.run_manifest is None:
            p.error("--summarize_errors requires --run_manifest <manifest> to locate the error records")
        summarize_errors(options.run_manifest, options.summarize_errors, options.chunk_size, options.array_job_id)
        return
    if options.emit_manifest is not None:
        problems, prevalidation_errors = enumerate_problems(options)
        random.shuffle(problems)  # spread cheap/expensive problems across chunks for balanced array tasks
        with open(options.emit_manifest, 'w') as f:
            for problem in problems:
                f.write(json.dumps(problem) + "\n")
        # always (re)write the sidecar, even when empty, so a prior submission's prevalidation records don't
        # leak into this run's error summary (the manifest is likewise rewritten above)
        with open(prevalidation_path(options.emit_manifest), 'w', encoding='utf-8') as f:
            for rec in prevalidation_errors:
                f.write(json.dumps(rec) + "\n")
        print(len(problems))
        return
    if options.run_manifest is not None:
        run_manifest_chunk(options.run_manifest, options.task_id, options.chunk_size)
        return

    constant_options = {k: v for (k, v) in options._get_kwargs() if not isinstance(v, list)}
    variable_options = {k: v for (k, v) in options._get_kwargs() if isinstance(v, list)}
    options_combinations = (dict(zip(variable_options, x)) for x in itertools.product(*variable_options.values()))

    # parameter folders this run is responsible for, and the output bases they live under, so we can flag
    # pre-existing folders left over from a different config setup at the end of the run.
    expected_output_parent_dirs = set()
    output_bases = set()

    # run over different combinations of options as specified in the config
    for options_combination in options_combinations:
        print("Current parameters: ", options_combination)

        options_combination['inference_method'] = resolve_inference_method(options_combination['inference_method'])

        # represent the argument combination as a string to use in the output directory name
        comb_str = build_comb_str(options_combination)

        # kwargs = options_combination | constant_options (works on python>=3.9)
        kwargs = options_combination.copy()
        kwargs.update(constant_options)

        for data_dir in os.scandir(kwargs['data_parent_dir']):
            print("Processing data in {}".format(data_dir.path))
            if not os.path.isdir(data_dir):
                continue
            data_dir_partial_path = os.path.join(*os.path.normpath(kwargs['data_parent_dir']).split(os.sep)[-2:])
            output_base = os.path.join("inferred_models", data_dir_partial_path)
            output_parent_dir = os.path.join(output_base, "{}-{}".format(comb_str, data_dir.name))
            # don't wipe existing output: we pick up where a previous run left off (per-network skip below)
            os.makedirs(output_parent_dir, exist_ok=True)
            output_bases.add(os.path.normpath(output_base))
            expected_output_parent_dirs.add(os.path.normpath(output_parent_dir))

            # concatenate log for model generation with the inference log
            # with open(os.path.join(output_parent_dir, "log.txt"), 'wb') as new_log_file:
            #     with open(os.path.join(data_dir.path, "log.txt"), 'rb') as old_log_file:
            #         shutil.copyfileobj(old_log_file, new_log_file)

            logger = logging.getLogger()
            for h in list(logger.handlers):
                h.close()  # close the previous data_dir's file handler so its descriptor isn't leaked
                logger.removeHandler(h)
            logging.basicConfig(filename=os.path.join(output_parent_dir, "log.txt"),
                                filemode='a',
                                format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                                datefmt='%H:%M:%S',
                                level=logging.DEBUG)
            for var, val in kwargs.items():
                logger.info("{}={}".format(var, val))

            paths = [(f.name, f.path) for f in os.scandir(data_dir.path) if f.is_dir()]
            if len(paths) == 0:
                raise ValueError("Did not find network directories in data_dir {}".format(data_dir.name))

            # pick-up: skip networks whose output folder is already fully populated from a previous run
            paths_to_run = []
            for name, path in paths:
                if network_output_is_complete(os.path.join(output_parent_dir, name), kwargs):
                    msg = "Skipping network {}: output already fully populated".format(name)
                    print(msg)
                    logger.info(msg)
                else:
                    paths_to_run.append((name, path))

            n_processes = max(1, int(kwargs.get('n_processes', 1) or 1))
            worker = functools.partial(process_network,
                                       output_parent_dir=output_parent_dir,
                                       kwargs=kwargs)
            if not paths_to_run:
                network_log_meta = []
            elif n_processes > 1:
                with ProcessPoolExecutor(max_workers=n_processes) as ex:
                    network_log_meta = list(ex.map(worker,
                                                   [name for name, _ in paths_to_run],
                                                   [path for _, path in paths_to_run]))
            else:
                network_log_meta = [worker(name, path) for name, path in paths_to_run]

            master_log_path = os.path.join(output_parent_dir, "log.txt")
            for h in logger.handlers:
                h.flush()
            for meta in network_log_meta:
                append_network_log_to_master(master_log_path, *meta)

    # after everything is saved: any parameter folder under a touched base that this run did not produce
    # is a leftover from a different config setup (since this run only creates its own folders). Flag it so
    # the user knows old outputs weren't cleared before running a different config.
    stale_param_dirs = []
    for base in sorted(output_bases):
        if not os.path.isdir(base):
            continue
        for entry in os.scandir(base):
            if entry.is_dir() and os.path.normpath(entry.path) not in expected_output_parent_dirs:
                stale_param_dirs.append(entry.path)
    if stale_param_dirs:
        raise RuntimeError(
            "Found parameter folder(s) under the output directory that this config did not produce - likely "
            "leftovers from a different config setup. Clear them (or rerun with the matching config) so "
            "results aren't mixed:\n" + "\n".join(sorted(stale_param_dirs)))


if __name__ == "__main__":
    main()
