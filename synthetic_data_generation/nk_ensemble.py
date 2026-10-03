"""The Kauffman NK ensemble and the measurements the make_random_nk_*_models scripts report on it.

Shared by the NK graphs_dirs rather than copied into each, so the ensembles differ only in the grid they
sweep and in the function family they draw, not in how a network is put together:
`make_random_nk_regime_sweep_models` holds p at 0.5 and varies K, which moves the regime, and
`make_random_nk_critical_models` moves the function parameter with K to hold the regime at criticality. A
model already written by either script is reproducible only as long as the generator draws in the order it
does here, so change the RNG usage only together with regenerating those directories.

The regime is set by the annealed approximation: a network is critical when 2*p*(1-p)*K = 1, ordered below
and chaotic above. That formula is about functions whose truth tables are i.i.d. coin flips; it does not
hold for the symmetric threshold family, which has its own closed form in `threshold_sensitivity`.
`sensitivity` measures the quantity itself, for either family, rather than assuming a formula.

Two function families are generated here:

* general Boolean, an i.i.d. bias-p truth table per node - `nk_network`, in one of two representations
  according to in-degree (see DENSE_MAX_ROWS);
* signed symmetric threshold, the class the symmetric / symmetric_topology inference methods search -
  `nk_threshold_network`.
"""
import math
import random

import numpy

# Truth-table rows above which a general Boolean function is stored as its minority rows rather than as a
# full table (SparseBooleanFunc rather than BooleanSymbolicFunc). 1024 rows is K = 10.
#
# The cap is a representation choice, but not only that: the two representations are also sampled
# differently - below it, one coin flip per row in table order; above it, a row count drawn from the
# binomial and then that many distinct row indices (see _sample_minority_rows). Both give exactly the same
# distribution over functions, but they consume the random stream differently, so a model at a given seed
# is NOT the same network on either side of this constant. Moving it regenerates every K that crosses it
# into different networks, which is why it sits above every K the existing directories were built with.
DENSE_MAX_ROWS = 1024


def critical_bias(k):
    """The larger of the two biases that put in-degree k exactly at criticality for the general Boolean
    family, i.e. the root of 2*p*(1-p)*k = 1 above 0.5: p = (1 + sqrt(1 - 2/k)) / 2. Undefined below k = 2,
    since p(1-p) <= 0.25 caps 2*p*(1-p)*k at k/2 - a k = 1 network cannot reach sensitivity 1 at any bias,
    which is why the ordered end of the regime sweep has no counterpart in the critical sweep.

    The symmetric root 1 - critical_bias(k) is equally critical; the upper one is used throughout so a
    single convention fixes which of the two a given K means."""
    if k < 2:
        raise ValueError("no bias makes in-degree {} critical: 2p(1-p)k <= k/2 < 1".format(k))
    return (1 + math.sqrt(1 - 2.0 / k)) / 2


def threshold_bias(k, threshold):
    """P(output = 1) over uniform random inputs, for a signed symmetric threshold function of in-degree k.

    Note what this does not depend on: the signs. SymmetricThresholdFunction counts the inputs that agree
    with their sign (a positive sign agreeing with 1, a negative one with 0) and fires at `threshold`, so
    over uniform random inputs that count is Binomial(k, 1/2) whatever the sign composition is. The bias is
    therefore a property of the threshold alone, and the achievable biases are the discrete upper tails of
    Binomial(k, 1/2) - a threshold function generally cannot hit critical_bias(k) exactly, only bracket
    it."""
    return sum(math.comb(k, i) for i in range(threshold, k + 1)) / 2.0 ** k


def threshold_sensitivity(k, threshold):
    """Exact average sensitivity (expected nodes flipped by a one-bit perturbation) for a signed symmetric
    threshold function of in-degree k - the threshold family's analogue of 2*p*(1-p)*k, which does NOT
    apply to it.

    Flipping one input moves the agreement count by one, so it flips the output exactly when the other
    k - 1 inputs contribute threshold - 1 agreements: every input has influence C(k-1, threshold-1) /
    2**(k-1), and the sensitivity is k times that. Like the bias, it is independent of the signs."""
    if not 1 <= threshold <= k:
        return 0.0        # constant function: no input has any influence
    return k * math.comb(k - 1, threshold - 1) / 2.0 ** (k - 1)


def critical_threshold(k):
    """The threshold in [1, k] whose average sensitivity comes closest to 1: the closest a signed symmetric
    threshold function of in-degree k can be brought to criticality. Ties (k = 2, where t = 1 and t = 2 are
    both exactly 1) go to the higher bias, matching critical_bias's upper-root convention.

    Criticality is a statement about sensitivity, and for this family that is the only thing it can be
    fitted to. The alternative of fitting the bias to critical_bias(k) - the bias a critical i.i.d. table
    would carry - does not deliver criticality, because the two families disagree at equal bias: at the
    thresholds in play here a threshold function has about half to two thirds of the sensitivity an i.i.d.
    table of the same bias would, and the gap widens with k, so fitting the bias undershoots. The two
    rules agree at k = 2, 5 and 10 and part at k = 15, where fitting the bias gives t = 4 (sensitivity
    0.33, firmly ordered) against t = 5 here (0.92).

    The family still cannot reach criticality exactly in general, since the achievable sensitivities are
    the discrete values k*C(k-1,t-1)/2**(k-1): 1.00 at k = 2, but 1.25, 0.70 and 0.92 at k = 5, 10 and 15.
    Each model's README reports the sensitivity it actually has."""
    return min(range(1, k + 1),
               key=lambda t: (abs(threshold_sensitivity(k, t) - 1.0), -threshold_bias(k, t)))


def mirrored_threshold(k, threshold):
    """The threshold giving the complementary bias: threshold_bias(k, t) + threshold_bias(k, k+1-t) = 1,
    and the two have exactly equal sensitivity, since C(k-1, t-1) = C(k-1, k-t). So a node can be put on
    either side of 0.5 without moving the network off the regime the threshold was chosen for."""
    return k + 1 - threshold


def _sample_minority_rows(n_rows, minority_prob, rng):
    """The set of truth-table rows disagreeing with the majority output, drawn without ever walking the
    table: the number of them is Binomial(n_rows, minority_prob), and given that number the rows are a
    uniform subset. That is exactly an i.i.d. coin flip per row, but it costs O(minority rows) instead of
    O(n_rows), which is what makes a large in-degree reachable at all.

    The rows are drawn by rejection into a set rather than by random.sample, which materializes the whole
    population when the requested count is a large fraction of it - the one thing this must not do."""
    count = numpy.random.RandomState(rng.randrange(2 ** 32)).binomial(n_rows, minority_prob)
    rows = set()
    while len(rows) < count:
        rows.add(rng.randrange(n_rows))
    return rows


def _nk_regulators(n, k, target, rng):
    """The k inputs of one node: uniform without replacement from the other nodes (no self-loops, matching
    the scale-free generator), sorted ascending so the first input is the truth table's most significant
    bit."""
    regulators = rng.sample([i for i in range(n) if i != target], k)
    regulators.sort()
    return regulators


def nk_network(n, k, bias, seed, dense_max_rows=DENSE_MAX_ROWS, mixed_polarity=False):
    """One NK network with general Boolean functions: exactly k inputs per node and an i.i.d. bias-`bias`
    truth table per node.

    With `mixed_polarity`, each node instead tosses a coin for `bias` or 1 - `bias`. Criticality is
    untouched by that choice - the annealed sensitivity 2*p*(1-p)*k is invariant under p -> 1-p, so every
    mixture of the two is exactly as critical as either alone - but the trajectories are not: a set drawn
    entirely at the upper root has nearly every truth table full of ones, and collapses onto the all-ones
    corner, while a mixture leaves roughly half the nodes biased each way and settles near density 0.5.
    The coin is tossed per node BEFORE that node's table is drawn, which is what fixes the random stream;
    off by default so the regime sweep keeps drawing exactly the sequence it always has (at p = 0.5 the
    mirror is a no-op in value, but not in RNG consumption).

    The table is held as a full BooleanSymbolicFunc while 2**k is at most `dense_max_rows`, and as a
    SparseBooleanFunc - only the rows disagreeing with the majority output - above that. The sparse form is
    what lifts the in-degree ceiling: the dense form costs 2**k in construction time (a sympy DNF clause
    per true row), in stored size and in load time, which is minutes and megabytes per model by k = 15,
    while the sparse form costs only the disagreeing rows, about 1100 of 32768 at k = 15. The ceiling it
    leaves is the minority count itself, so it does not reach an arbitrary k: at the critical bias the
    minority fraction falls only like 1/(2k), so k = 25 is about 680k rows per node and k = 50 about 1e13 -
    still out of reach."""
    from attractor_learning.graphs import Network
    from attractor_learning.logic import BooleanSymbolicFunc, SparseBooleanFunc

    rng = random.Random(seed)
    names = ["v{}".format(i) for i in range(n)]
    n_rows = 2 ** k
    edges, functions = [], []
    for target in range(n):
        regulators = _nk_regulators(n, k, target, rng)
        edges.extend((names[source], names[target]) for source in regulators)
        input_names = [names[source] for source in regulators]
        node_bias = 1.0 - bias if (mixed_polarity and rng.random() < 0.5) else bias
        if n_rows <= dense_max_rows:
            functions.append(BooleanSymbolicFunc(
                input_names=input_names,
                boolean_outputs=[rng.random() < node_bias for _ in range(n_rows)]))
        else:
            majority = node_bias >= 0.5
            minority_rows = _sample_minority_rows(n_rows, 1.0 - node_bias if majority else node_bias, rng)
            functions.append(SparseBooleanFunc(input_names=input_names, minority_rows=minority_rows,
                                               default_output=majority))
    return Network(vertex_names=names, edges=edges, vertex_functions=functions)


def nk_threshold_network(n, k, threshold, seed, mixed_polarity=False):
    """One NK network whose nodes are signed symmetric threshold functions: the same topology as
    `nk_network` - exactly k inputs per node, uniform without replacement from the other nodes - with a
    uniformly random sign per input (a fair coin: rng.choice over [1, -1]) and a threshold.

    `threshold` None draws it per node from Uniform[1, k] instead of sharing one, which is how
    make_random_scale_free_models draws its functions. That ensemble is critical at every k without
    anything being fitted, since summing k*C(k-1,t-1)/2**(k-1) over t = 1..k gives k*2**(k-1)/2**(k-1) and
    the mean over t is exactly 1. What a shared fitted threshold buys instead is zero variance: it is off
    criticality at most k (0.70 at k = 10) but identically so for every network, where uniform thresholds
    centre on 1 and scatter, with a per-node standard deviation of 0.6 to 1.1 and hence about 0.09 to 0.16
    for a 50-node network. Bias for variance, in other words. Uniform thresholds also spread the node
    biases symmetrically about 0.5 on their own, which is why `mixed_polarity` is meaningless with them -
    mirroring a Uniform[1, k] threshold gives back a Uniform[1, k] threshold.

    With a shared `threshold` and `mixed_polarity`, each node tosses a coin for `threshold` or its mirror
    k + 1 - threshold, the threshold family's version of the choice `nk_network` makes between p and 1 - p:
    the two have complementary biases and identical sensitivity (see `mirrored_threshold`), so the coin
    moves half the nodes below 0.5 without moving the network off its regime. Either draw happens per node
    before that node's signs are.

    Stored as signs plus a threshold, so a node costs a few hundred bytes whatever k is; there is no truth
    table to sample or to write, and so no in-degree ceiling at all. This is the class the symmetric and
    symmetric_topology inference methods search, so unlike the general Boolean sets these models sit inside
    the searched class - and the ILP method does not accept them at all (see SymmetricThresholdFunction's
    note about threshold not being Boolean).

    The signs are free: they decide which inputs activate and which inhibit, which is what topology
    inference has to recover, but the bias and the sensitivity depend only on the threshold (see
    `threshold_bias`)."""
    from attractor_learning.graphs import Network
    from attractor_learning.logic import SymmetricThresholdFunction

    rng = random.Random(seed)
    names = ["v{}".format(i) for i in range(n)]
    edges, functions = [], []
    for target in range(n):
        regulators = _nk_regulators(n, k, target, rng)
        edges.extend((names[source], names[target]) for source in regulators)
        if threshold is None:
            node_threshold = rng.randint(1, k)
        elif mixed_polarity and rng.random() < 0.5:
            node_threshold = mirrored_threshold(k, threshold)
        else:
            node_threshold = threshold
        functions.append(SymmetricThresholdFunction(
            signs=[rng.choice([1, -1]) for _ in range(k)], threshold=node_threshold))
    return Network(vertex_names=names, edges=edges, vertex_functions=functions)


def exact_threshold_sensitivity(network):
    """The network's Derrida slope, computed rather than measured: for threshold functions every input of
    a node has the same influence, so the expected number of nodes flipped by a one-bit perturbation is the
    mean of the nodes' own average sensitivities. The in-degree weighting cancels, so this holds whatever
    the degrees are.

    Constant across a set built on one shared threshold, and the whole point of reporting it when the
    thresholds are drawn per node, since then it is what separates one network from another."""
    from attractor_learning.logic import SymmetricThresholdFunction

    values = []
    for vertex in network.vertices:
        function = vertex.function
        if not isinstance(function, SymmetricThresholdFunction):
            raise ValueError("{} is not a threshold function".format(vertex.name))
        values.append(threshold_sensitivity(len(function.signs), function.threshold))
    return sum(values) / float(len(values))


def high_polarity_count(network):
    """How many of the network's nodes output 1 on more than half their truth-table rows - the realized
    result of the polarity coin, read back off the functions rather than trusted from the draw.

    Worth reporting per model because the split is Binomial(n, 1/2), so individual models run well away
    from even (about 18/32 either way at n = 50), and that is what explains a model's density sitting off
    0.5. At k = 2 the two critical biases coincide at 0.5, so there the count is just the tables' own
    fluctuation and carries no polarity meaning."""
    from attractor_learning.logic import SparseBooleanFunc, SymmetricThresholdFunction

    high = 0
    for vertex in network.vertices:
        function = vertex.function
        if isinstance(function, SparseBooleanFunc):
            high += 1 if function.default_output else 0
        elif isinstance(function, SymmetricThresholdFunction):
            high += 1 if threshold_bias(len(function.signs), function.threshold) > 0.5 else 0
        else:
            outputs = function.boolean_outputs
            high += 1 if 2 * sum(bool(out) for out in outputs) > len(outputs) else 0
    return high


def sensitivity(network, seed, samples=200):
    """Average number of nodes that differ one step after a single-bit perturbation - the Derrida slope at
    distance 1. Below 1 is ordered, 1 is critical, above 1 is chaotic; the annealed approximation puts it
    at 2*p*(1-p)*K for the general Boolean family, and threshold_sensitivity gives it exactly for the
    threshold family.

    A mean of `samples` draws, so it is noisy per model: measured across the n = 50 critical sets, the
    standard deviation from model to model runs from 0.06 (K = 5 upwards) to 0.11 (K = 2), which puts
    individual values as much as 0.2 from the ensemble mean. That is the estimator, not a network that
    missed its regime - a cell is only readable across its ten models."""
    rng = random.Random(seed)
    n = len(network)
    total = 0
    for _ in range(samples):
        state = [rng.randint(0, 1) for _ in range(n)]
        perturbed = list(state)
        flipped = rng.randrange(n)
        perturbed[flipped] = 1 - perturbed[flipped]
        after = [int(bool(v)) for v in network.next_state(state)]
        after_perturbed = [int(bool(v)) for v in network.next_state(perturbed)]
        total += sum(a != b for a, b in zip(after, after_perturbed))
    return total / float(samples)


def mean_density(network, seed, trajectories=10, burn_in=50, measured_steps=50):
    """Mean fraction of nodes holding the value 1 along a trajectory, measured after `burn_in` steps from a
    uniform random state have been run and discarded - so what is reported is where the network settles,
    not the arbitrary state it was started from.

    Reported because a biased ensemble reaches criticality by making nodes nearly constant: at p = 0.95 the
    truth tables are mostly ones, so trajectories collapse toward the all-ones corner and a large frozen
    core is what an inference method actually sees. A density near p is the expected behaviour, not a bug,
    but it is the number that says whether a graphs_dir exercises the methods or hands them a near-static
    signal that `all_constants` already fits."""
    rng = random.Random(seed)
    n = len(network)
    total = 0.0
    for _ in range(trajectories):
        state = [rng.randint(0, 1) for _ in range(n)]
        for _ in range(burn_in):
            state = [int(bool(v)) for v in network.next_state(state)]
        for _ in range(measured_steps):
            state = [int(bool(v)) for v in network.next_state(state)]
            total += sum(state) / float(n)
    return total / float(trajectories * measured_steps)


def attractor_summary(network, enumerate_up_to=10):
    """(number of attractors, longest period) by enumerating the whole state space, or None when too big."""
    n = len(network)
    if n > enumerate_up_to:
        return None
    successor = {}
    for bits in range(2 ** n):
        state = tuple((bits >> (n - 1 - i)) & 1 for i in range(n))
        successor[state] = tuple(int(bool(v)) for v in network.next_state(list(state)))
    attractors, seen = [], {}
    for state in successor:
        path, current = [], state
        while current not in seen and current not in path:
            path.append(current)
            current = successor[current]
        if current in path:                     # closed a new cycle
            cycle = path[path.index(current):]
            attractors.append(cycle)
            for member in cycle:
                seen[member] = len(attractors) - 1
        for member in path:
            seen.setdefault(member, seen.get(current, None))
    return len(attractors), max(len(cycle) for cycle in attractors)


def check_round_trip(network, model_dir, seed, samples=20):
    """Assert the written folder loads back with the same dynamics - a model directory is only useful if
    the pipeline's readers can read it, and a silently lossy write would surface much later as an inference
    result rather than as an error here."""
    from attractor_learning.graphs import Network

    n = len(network)
    assert Network.model_dir_size(model_dir) == n
    reloaded = Network.parse_model_dir(model_dir)
    rng = random.Random(seed)
    for _ in range(samples):
        state = [rng.randint(0, 1) for _ in range(n)]
        assert [int(bool(v)) for v in reloaded.next_state(state)] == \
               [int(bool(v)) for v in network.next_state(state)], \
            "{}: json round trip changed dynamics".format(model_dir)
