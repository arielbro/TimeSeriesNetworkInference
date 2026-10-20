"""Build the two folders of hand-made reference models under data/, both usable as a `graphs_dir`
(the role data/cellcollective_models plays for the real runs).

data/dynamic_toy_models - models whose state space keeps moving: several attractors, none of them (or
    almost none of them) fixed points. Three families:

      ring   repressilator rings of odd length, plus feed-forward readouts. The cycle IS the structure
             here; kept as the simple baseline, and the only family with long periods (6, 10, 30).
      motif  textbook transcription motifs: negative autoregulation (NAR), activator-inhibitor pairs
             (p53-Mdm2 / NF-kB-IkB shape), toggle switches, coherent and incoherent feed-forward loops,
             single-input modules.
      dense  the same requirement (multiple moving attractors) met without any part of the model being a
             ring: cliques, and feedback loops carrying chords. Every cycle here is a proper subgraph of
             a denser neighbourhood, so no node's regulation is just "the previous node in a cycle".

data/classic_toy_models - the classic feed-forward circuits, which have only fixed-point attractors and
    are meant for testing function and topology inference rather than dynamics:

      ffl    all eight non-equivalent three-node feed-forward loops (coherent C1-C4, incoherent I1-I4),
             each with an AND gate and with an OR gate at the target - sixteen models.
      gates  multi-level AND/OR combinations, majority and k-of-n thresholds, cascades, a single-input
             module and a dense overlapping regulon.

Every function in both folders is a symmetric threshold function (its inputs carry +/- signs and the node
fires when at least `threshold` of them agree), the model class the symmetric / symmetric_topology
inference methods search.

The whole state space is enumerated per model (n <= 10, so at most 1024 states), so attractors and basin
sizes are exact, and each family's requirements are asserted before anything is written.

Each model directory holds:
  SPECIES_KEY.csv + <index>.csv   cellcollective truth-table format, what generate.py's graphs_dir reads
  true_network.cnet               the same model as legible text
  network.png, network.svg        interaction graph + state graph (every state, or above
                                  MAX_FULL_STATE_GRAPH states the attractors with their basin sizes)
  README.md                       functions in words, and the attractor or truth-table analysis
"""
import itertools
import os
import shutil
import sys

import matplotlib
matplotlib.use("Agg")
# SVG text kept as text (editable, searchable) rather than outlined, and element ids made reproducible, so a
# regeneration that changes nothing leaves the .svg files unchanged too
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['svg.hashsalt'] = 'make_toy_models'
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from attractor_learning.graphs import Network
from attractor_learning.logic import SymmetricThresholdFunction

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
DYNAMIC_DIR = os.path.join(DATA_DIR, "dynamic_toy_models")
CLASSIC_DIR = os.path.join(DATA_DIR, "classic_toy_models")

MIN_ATTRACTORS = 2          # dynamic models: a single-attractor model is not what these are for
MIN_LONG_CYCLE = 3          # dynamic models: at least one attractor of at least this period
MIN_MEAN_HAMMING = 1.0      # dynamic models: mean nodes flipping per step along the attractors

MAX_LISTED_STATES = 10      # per attractor in the README, before eliding
MAX_TRUTH_ROWS = 16         # input patterns listed in a classic model's truth table


# ---------------------------------------------------------------------------------------------------
# vocabulary: every builder returns a list of (name, [(input name, sign)], threshold, formula) nodes
# ---------------------------------------------------------------------------------------------------
def node(name, inputs, threshold, text):
    return [(name, inputs, threshold, text)]


def gate_text(inputs, threshold):
    """Readable formula for a threshold function over signed inputs."""
    terms = ["{}{}".format("" if sign > 0 else "NOT ", name) for name, sign in inputs]
    if threshold == len(inputs):
        return " AND ".join(terms)
    if threshold == 1:
        return " OR ".join(terms)
    return "at least {} of ({})".format(threshold, ", ".join(terms))


def ring(prefix, length):
    names = ["{}{}".format(prefix, i) for i in range(length)]
    return [(name, [(names[i - 1], -1)], 1, "NOT {}".format(names[i - 1]))
            for i, name in enumerate(names)]


def nar(name):
    return node(name, [(name, -1)], 1, "NOT {}".format(name))


def activator_inhibitor(activator, inhibitor):
    return (node(activator, [(inhibitor, -1)], 1, "NOT {}".format(inhibitor)) +
            node(inhibitor, [(activator, 1)], 1, activator))


def toggle(first, second):
    return (node(first, [(second, -1)], 1, "NOT {}".format(second)) +
            node(second, [(first, -1)], 1, "NOT {}".format(first)))


def c1_ffl(source, intermediate, target):
    return (node(intermediate, [(source, 1)], 1, source) +
            node(target, [(source, 1), (intermediate, 1)], 2,
                 "{} AND {}".format(source, intermediate)))


def i1_ffl(source, intermediate, target):
    return (node(intermediate, [(source, 1)], 1, source) +
            node(target, [(source, 1), (intermediate, -1)], 2,
                 "{} AND NOT {}".format(source, intermediate)))


def sim(master, targets):
    return [(target, [(master, 1)], 1, master) for target in targets]


def clique(prefix, size, signs, threshold):
    """A fully connected group: every node is regulated by all the others (no self-loops). `signs` is the
    rule applied uniformly - signs[k] is the sign of the k-th neighbour clockwise from the node - so the
    whole group is described by one sentence rather than n unrelated functions."""
    names = ["{}{}".format(prefix, i) for i in range(size)]
    nodes = []
    for i, name in enumerate(names):
        neighbours = [names[(i + 1 + k) % size] for k in range(size - 1)]
        inputs = list(zip(neighbours, signs))
        nodes.append((name, inputs, threshold, gate_text(inputs, threshold)))
    return nodes


def chorded_loop(prefix, size, signs, threshold):
    """A feedback loop of `size` nodes carrying chords: each node is regulated by the three nodes before
    it, with `signs` giving their signs in that order. The long cycle is then a proper subgraph of a
    denser neighbourhood rather than the whole structure."""
    names = ["{}{}".format(prefix, i) for i in range(size)]
    nodes = []
    for i, name in enumerate(names):
        inputs = [(names[i - 1 - k], sign) for k, sign in enumerate(signs)]
        nodes.append((name, inputs, threshold, gate_text(inputs, threshold)))
    return nodes


def ffl(kind, gate, source="x", intermediate="y", target="z"):
    """One of the eight feed-forward loops, with an AND or OR gate at the target. `kind` names the type
    the way Alon's table does: coherent when sign(x->z) == sign(x->y) * sign(y->z), incoherent otherwise."""
    s_xy, s_yz, s_xz = FFL_SIGNS[kind]
    threshold = 2 if gate == "AND" else 1
    inputs = [(source, s_xz), (intermediate, s_yz)]
    return (node(source, [], 0, "input") +
            node(intermediate, [(source, s_xy)], 1,
                 source if s_xy > 0 else "NOT {}".format(source)) +
            node(target, inputs, threshold, gate_text(inputs, threshold)))


FFL_SIGNS = {          # (x->y, y->z, x->z), numbered as in Mangan & Alon, PNAS 2003
    "c1": (1, 1, 1), "c2": (-1, 1, -1), "c3": (1, -1, -1), "c4": (-1, -1, 1),
    "i1": (1, -1, 1), "i2": (-1, -1, -1), "i3": (1, 1, -1), "i4": (-1, 1, 1),
}


def build(parts):
    """parts: list of node lists. Returns the Network and {name: (formula text, threshold)}."""
    nodes = [n for part in parts for n in part]
    names = [name for name, _, _, _ in nodes]
    assert len(names) == len(set(names)), "duplicate node names: {}".format(names)
    order = {name: i for i, name in enumerate(names)}
    edges, functions, descriptions = [], [], {}
    for name, inputs, threshold, text in nodes:
        # predecessors() is ordered by vertex index, and SymmetricThresholdFunction reads its signs in
        # that same order, so sort each node's inputs by their position in `names` before taking signs
        ordered = sorted(inputs, key=lambda pair: order[pair[0]])
        edges.extend((source, name) for source, _ in ordered)
        functions.append(SymmetricThresholdFunction(signs=[sign > 0 for _, sign in ordered],
                                                    threshold=threshold))
        descriptions[name] = (text, threshold)
    network = Network(vertex_names=names, edges=edges, vertex_functions=functions)
    verify_functions(network, nodes)
    return network, descriptions


def verify_functions(network, nodes):
    """Each node's realized function must equal the formula it was declared with - the guard against a
    sign landing on the wrong input when a node mixes activation and repression."""
    specs = {name: (inputs, threshold) for name, inputs, threshold, _ in nodes}
    for v in network.vertices:
        inputs, threshold = specs[v.name]
        predecessors = [u.name for u in v.predecessors()]
        signs = dict(inputs)
        for values in itertools.product([False, True], repeat=len(predecessors)):
            expected = sum(1 for name, value in zip(predecessors, values)
                           if (signs[name] > 0) == value) >= threshold
            assert bool(v.function(*values)) == expected, \
                "{}: realized function disagrees with its formula at {}".format(v.name, values)


# ---------------------------------------------------------------------------------------------------
# state space
# ---------------------------------------------------------------------------------------------------
def state_space(network):
    n = len(network)
    transitions = {}
    for bits in range(2 ** n):
        state = tuple((bits >> (n - 1 - i)) & 1 for i in range(n))
        transitions[state] = tuple(int(bool(v)) for v in network.next_state(list(state)))
    return transitions


def attractors_and_basins(transitions):
    """Exact attractors and basin sizes, by walking every state to its cycle - deterministic dynamics on
    a finite state space, so this is complete, not an estimate."""
    cycle_id, cycles = {}, []
    for start in transitions:
        path, state = [], start
        while state not in cycle_id and state not in path:
            path.append(state)
            state = transitions[state]
        if state in path:
            cycle = path[path.index(state):]
            cycles.append(cycle)
            for member in cycle:
                cycle_id[member] = len(cycles) - 1
        target = cycle_id[state]
        for member in path:
            cycle_id.setdefault(member, target)
    basins = [0] * len(cycles)
    for state in transitions:
        basins[cycle_id[state]] += 1
    order = sorted(range(len(cycles)), key=lambda k: (-basins[k], cycles[k]))
    return [cycles[k] for k in order], [basins[k] for k in order]


def mean_hamming(cycles):
    distances = []
    for cycle in cycles:
        for i, state in enumerate(cycle):
            nxt = cycle[(i + 1) % len(cycle)]
            distances.append(sum(a != b for a, b in zip(state, nxt)))
    return float(np.mean(distances))


def dynamic_stats(network):
    cycles, basins = attractors_and_basins(state_space(network))
    return cycles, basins, mean_hamming(cycles)


def choose_uniform_rule(make_nodes, candidates):
    """Pick the (signs, threshold) candidate that gives the liveliest model, so a dense group stays
    describable by a single uniform rule instead of n hand-tuned functions. Among the candidates meeting
    the dynamic criteria, prefer the fewest fixed-point attractors (a clique that mostly settles is not
    what this folder is for), then the most movement per step, then the longest cycle. Deterministic."""
    scored = []
    for signs, threshold in candidates:
        network, _ = build([make_nodes(signs, threshold)])
        cycles, _, hamming = dynamic_stats(network)
        longest = max(len(c) for c in cycles)
        if len(cycles) < MIN_ATTRACTORS or longest < MIN_LONG_CYCLE or hamming < MIN_MEAN_HAMMING:
            continue
        fixed = sum(1 for c in cycles if len(c) == 1)
        scored.append(((fixed, -hamming, -longest), (signs, threshold)))
    if not scored:
        raise AssertionError("no uniform rule in the candidate set gives a moving, multi-attractor model")
    return min(scored, key=lambda pair: pair[0])[1]


def uniform_candidates(degree):
    """All (sign pattern, threshold) rules for a node of the given in-degree, ordered so that balanced
    thresholds and mixed signs come first (they are the ones that tend to oscillate)."""
    patterns = sorted(itertools.product([1, -1], repeat=degree),
                      key=lambda p: (abs(sum(p)), [-s for s in p]))
    thresholds = sorted(range(1, degree + 1), key=lambda t: abs(t - (degree + 1) / 2.0))
    return [(p, t) for t in thresholds for p in patterns]


# ---------------------------------------------------------------------------------------------------
# drawing
# ---------------------------------------------------------------------------------------------------
def layout(network):
    """Feedback groups each on their own circle, feed-forward nodes laid out in regulation order below."""
    graph = nx.DiGraph((u.name, v.name) for u, v in network.edges)
    graph.add_nodes_from(v.name for v in network.vertices)
    loops = [c for c in nx.strongly_connected_components(graph)
             if len(c) > 1 or any(u == v for u, v in graph.edges(c))]
    in_loop = {name for component in loops for name in component}
    positions = {}
    for k, component in enumerate(sorted(loops, key=lambda c: sorted(c))):
        members = sorted(component)
        radius = 0.55 if len(members) <= 2 else 1.0
        for i, name in enumerate(members):
            angle = 2 * np.pi * i / len(members) + np.pi / 2
            positions[name] = np.array([3.0 * k + radius * np.cos(angle), radius * np.sin(angle)])
    rest = [v.name for v in network.vertices if v.name not in in_loop]
    if not loops and rest:                       # purely feed-forward: layer by longest path from a source
        depth = {}
        for name in nx.topological_sort(graph):
            preds = list(graph.predecessors(name))
            depth[name] = 0 if not preds else 1 + max(depth[p] for p in preds)
        by_level = {}
        for name in rest:
            by_level.setdefault(depth[name], []).append(name)
        for level, members in by_level.items():
            # stagger odd levels sideways, so a skipping edge (x -> z over y, the FFL shape) doesn't run
            # straight through the node it skips
            offset = 0.75 if level % 2 else 0.0
            for i, name in enumerate(sorted(members)):
                positions[name] = np.array([offset + 1.5 * (i - (len(members) - 1) / 2.0), -1.2 * level])
        return positions
    columns = max(1, int(np.ceil(np.sqrt(len(rest)))))
    for i, name in enumerate(rest):
        positions[name] = np.array([1.6 * (i % columns) - 0.4, -2.2 - 1.0 * (i // columns)])
    return positions


NODE_RADIUS = 0.24               # interaction graph: node circle radius, in data units
PARALLEL_OFFSET = 0.10           # interaction graph: sideways shift of each edge of a reciprocal pair

MAX_FULL_STATE_GRAPH = 64        # draw every state up to this many; above it, attractors + basin sizes only
ATTRACTOR_COLOUR = '#58B7CA'
STATE_FILL = '#D9D9D9'
STATE_INCH_PER_UNIT = 0.6        # state graph scale; with STATE_FONT_SIZE it sizes the ellipses to the labels
STATE_FONT_SIZE = 8
STATE_CHAR_WIDTH = 0.105         # data units per label character at that scale and font size
STATE_HEIGHT = 0.34
STATE_ARROW_ROOM = 0.45          # space left between neighbouring ellipses for an arrow


def _clip_to_circle(start, end, radius):
    """The segment start -> end, shortened by `radius` at both ends so it runs between circle rims."""
    start, end = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
    unit = (end - start) / np.hypot(*(end - start))
    return start + unit * radius, end - unit * radius


def draw_graph(ax, network, descriptions):
    """Interaction graph with straight edges. A reciprocal pair (a -> b and b -> a) is drawn as two straight
    lines shifted to either side of the line between the nodes, so neither hides the other."""
    signs = {}
    for v in network.vertices:
        if v.function is None:
            continue
        for u, sign in zip(v.predecessors(), v.function.signs):
            signs[(u.name, v.name)] = 1 if sign else -1
    positions = layout(network)
    for name, (x, y) in positions.items():
        ax.add_patch(matplotlib.patches.Circle((x, y), NODE_RADIUS, facecolor='white', edgecolor='black',
                                               lw=1.0, zorder=2))
        ax.text(x, y, name, ha='center', va='center', fontsize=7 if len(name) <= 4 else 5.5, zorder=3)
    for (source, target), sign in signs.items():
        colour = 'tab:green' if sign > 0 else 'tab:red'
        style = '-|>' if sign > 0 else '-['
        if source == target:
            # a self-loop cannot be a straight line; a small loop over the node keeps autoregulation visible
            x, y = positions[source]
            ax.add_patch(matplotlib.patches.FancyArrowPatch(
                (x - 0.12, y + NODE_RADIUS * 0.9), (x + 0.12, y + NODE_RADIUS * 0.9),
                connectionstyle='arc3,rad=-2.2', arrowstyle=style, mutation_scale=9, lw=1.2, color=colour,
                zorder=1))
            continue
        start, end = np.asarray(positions[source]), np.asarray(positions[target])
        if (target, source) in signs:
            direction = (end - start) / np.hypot(*(end - start))
            normal = np.array([-direction[1], direction[0]]) * PARALLEL_OFFSET
            start, end = start + normal, end + normal
        start, end = _clip_to_circle(start, end, NODE_RADIUS + 0.02)
        ax.add_patch(matplotlib.patches.FancyArrowPatch(
            start, end, arrowstyle=style, mutation_scale=12 if sign > 0 else 6, lw=1.2, color=colour,
            shrinkA=0, shrinkB=0, zorder=1))
    xs = [x for x, _ in positions.values()]
    ys = [y for _, y in positions.values()]
    span = max(max(xs) - min(xs), max(ys) - min(ys)) / 2 + 0.7
    centre = ((max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2)
    ax.set_xlim(centre[0] - span, centre[0] + span)
    ax.set_ylim(centre[1] - span, centre[1] + span)
    ax.set_aspect('equal')
    ax.axis('off')


def _state_label(state):
    return "".join(str(bit) for bit in state)


def _cycle_layout(cycle, width):
    """An attractor on its own, laid out as a racetrack: the first half of the cycle left to right along a
    top row, the rest back right to left beneath it. A long cycle then takes two rows of labels rather than
    a circle whose middle is empty space."""
    step = width + STATE_ARROW_ROOM
    columns = (len(cycle) + 1) // 2
    rows_apart = STATE_HEIGHT + 2 * STATE_ARROW_ROOM
    positions = {}
    for i, state in enumerate(cycle):
        if i < columns:
            positions[state] = np.array([i * step, rows_apart / 2])
        else:
            positions[state] = np.array([(len(cycle) - 1 - i) * step, -rows_apart / 2])
    if len(cycle) == 2:
        # one above the other, a row's gap apart, which leaves the pair of bowed arrows between them room
        # to show (side by side, the gap would be one arrow's length)
        positions = {cycle[0]: np.array([0.0, rows_apart / 2]), cycle[1]: np.array([0.0, -rows_apart / 2])}
    return positions


def _basin_layout(cycle, predecessors, full, width):
    """Positions for one attractor and, when `full`, every transient state of its basin.

    An attractor with no transient states to draw is a racetrack (see _cycle_layout). Otherwise the basin -
    a cycle with in-trees hanging off it - goes through a Kamada-Kawai layout, which spreads a tree out
    from its root, and is then scaled up uniformly until no two labels' boxes come closer than an arrow's
    length."""
    members = set(cycle)
    if full:
        frontier = list(cycle)
        while frontier:
            for state in predecessors.get(frontier.pop(), ()):
                if state not in members:
                    members.add(state)
                    frontier.append(state)
    if len(members) == len(cycle):
        return _cycle_layout(cycle, width)
    graph = nx.Graph()
    graph.add_nodes_from(sorted(members))
    for state in sorted(members):
        for predecessor in sorted(predecessors.get(state, ())):
            if predecessor != state:
                graph.add_edge(predecessor, state)
    raw = nx.kamada_kawai_layout(graph)
    states = sorted(raw)
    points = np.array([raw[state] for state in states])
    # the smallest scale at which every pair is separated, horizontally or vertically, by a label plus an
    # arrow - the labels are wide and flat, so the two directions need different clearances
    dx = np.abs(points[:, None, 0] - points[None, :, 0]) / (width + STATE_ARROW_ROOM)
    dy = np.abs(points[:, None, 1] - points[None, :, 1]) / (STATE_HEIGHT + STATE_ARROW_ROOM)
    separation = np.maximum(dx, dy)
    np.fill_diagonal(separation, np.inf)
    scale = 1.0 / separation.min()
    return {state: point * scale for state, point in zip(states, points)}


def _clip_to_ellipse(start, end, half_width, half_height):
    """The segment start -> end, shortened at both ends to run between the rims of two equal ellipses."""
    start, end = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
    unit = (end - start) / np.hypot(*(end - start))
    reach = 1.0 / np.hypot(unit[0] / half_width, unit[1] / half_height)
    return start + unit * reach, end - unit * reach


def draw_state_graph(ax, network, transitions, cycles, basins):
    """The state transition graph: every state an ellipse labelled by its bits, an arrow to its successor,
    attractor states and their arrows in blue. Laid out one basin at a time.

    Above MAX_FULL_STATE_GRAPH states only the attractors are drawn, each captioned with its basin size in
    place of its transient states. Returns the extent (width, height) in data units, so the caller can size
    the figure to it."""
    n_states = len(transitions)
    full = n_states <= MAX_FULL_STATE_GRAPH
    predecessors = {}
    for state, successor in transitions.items():
        predecessors.setdefault(successor, []).append(state)
    width = STATE_CHAR_WIDTH * len(network) + 0.25
    half_width, half_height = width / 2, STATE_HEIGHT / 2

    basin_layouts = []
    for cycle in cycles:
        positions = _basin_layout(cycle, predecessors, full, width)
        xs = [x for x, _ in positions.values()]
        ys = [y for _, y in positions.values()]
        # room for the ellipses, a fixed point's self-loop above, and the caption below
        box = (min(xs) - half_width, max(xs) + half_width, min(ys) - half_height - 0.45,
               max(ys) + half_height + (0.35 if len(cycle) == 1 else 0.0))
        box = (min(box[0], -0.9), max(box[1], 0.9), box[2], box[3])   # wide enough for the caption
        basin_layouts.append((positions, box))

    # pack the basins into rows, in the order given (largest basin first), aiming for a roughly 3:2 panel
    area = sum((b[1] - b[0] + 0.4) * (b[3] - b[2] + 0.4) for _, b in basin_layouts)
    row_limit = max(max(b[1] - b[0] for _, b in basin_layouts) + 0.4, np.sqrt(area * 1.5))
    placed, x, y, row_height = [], 0.0, 0.0, 0.0
    for positions, box in basin_layouts:
        box_width, box_height = box[1] - box[0] + 0.4, box[3] - box[2] + 0.4
        if x > 0 and x + box_width > row_limit:
            x, y, row_height = 0.0, y - row_height, 0.0
        offset = np.array([x - box[0], y - box[3]])
        placed.append(({state: position + offset for state, position in positions.items()}, box, offset))
        x += box_width
        row_height = max(row_height, box_height)
    extent_width = max(o[0] + b[1] for _, b, o in placed)
    extent_height = row_height - y

    attractor_states = {state for cycle in cycles for state in cycle}
    for (positions, box, offset), cycle, basin in zip(placed, cycles, basins):
        for state, (sx, sy) in positions.items():
            on_attractor = state in attractor_states
            ax.add_patch(matplotlib.patches.Ellipse(
                (sx, sy), width, STATE_HEIGHT, facecolor=STATE_FILL,
                edgecolor=ATTRACTOR_COLOUR if on_attractor else 'black', lw=2.2 if on_attractor else 1.0,
                zorder=2))
            ax.text(sx, sy, _state_label(state), ha='center', va='center', fontsize=STATE_FONT_SIZE,
                    zorder=3)
        for state, (sx, sy) in positions.items():
            successor = transitions[state]
            on_attractor = state in attractor_states
            colour = ATTRACTOR_COLOUR if on_attractor else 'black'
            if successor == state:      # a fixed point: a loop over the ellipse
                ax.add_patch(matplotlib.patches.FancyArrowPatch(
                    (sx - 0.25 * width, sy + half_height), (sx + 0.25 * width, sy + half_height),
                    connectionstyle='arc3,rad=-1.6', arrowstyle='-|>', mutation_scale=10, lw=1.6,
                    color=colour, zorder=1))
                continue
            start, end = _clip_to_ellipse((sx, sy), positions[successor], half_width, half_height)
            # a period-2 attractor's two arrows would coincide as straight lines, so they bow apart
            # - by a fixed sideways distance rather than a fixed fraction of their length, so a pair laid far
            # apart does not sweep a wide arc across the rest of the basin
            rad = min(0.35, 0.4 / np.hypot(*(end - start))) if (on_attractor and len(cycle) == 2) else 0
            bend = 'arc3,rad={}'.format(rad)
            ax.add_patch(matplotlib.patches.FancyArrowPatch(
                start, end, connectionstyle=bend, arrowstyle='-|>', mutation_scale=10,
                lw=1.6 if on_attractor else 0.9, color=colour, shrinkA=0, shrinkB=0, zorder=1))
        caption = "period {}, basin {}/{}".format(len(cycle), basin, n_states)
        ax.text(offset[0] + (box[0] + box[1]) / 2, offset[1] + box[2] + 0.1, caption, ha='center',
                va='bottom', fontsize=7, color='dimgray')
    # a margin all round, so outlines and self-loops at the edge are not clipped by the axes
    ax.set_xlim(-0.2, extent_width + 0.2)
    ax.set_ylim(-extent_height - 0.2, 0.2)
    ax.set_aspect('equal')
    ax.axis('off')
    return extent_width + 0.4, extent_height + 0.4


def draw_model(network, descriptions, transitions, cycles, basins, model_dir, title):
    """Interaction graph beside the state graph, saved as network.png and network.svg."""
    n = len(network)
    names = [v.name for v in network.vertices]
    full = len(transitions) <= MAX_FULL_STATE_GRAPH
    # the state panel is drawn at STATE_INCH_PER_UNIT, which the ellipse widths assume, so its size comes
    # from its extent - learnt by drawing it once on a scratch figure
    probe = plt.figure()
    state_width, state_height = draw_state_graph(probe.add_subplot(), network, transitions, cycles, basins)
    plt.close(probe)
    graph_inches = 3.8
    state_inches = (state_width * STATE_INCH_PER_UNIT, state_height * STATE_INCH_PER_UNIT)
    header, margin = 1.0, 0.3
    fig = plt.figure(figsize=(graph_inches + state_inches[0] + 3 * margin,
                              max(graph_inches, state_inches[1]) + header + margin))
    total_width, total_height = fig.get_size_inches()
    top = total_height - header
    graph_ax = fig.add_axes((margin / total_width, (top - graph_inches) / total_height,
                             graph_inches / total_width, graph_inches / total_height))
    state_ax = fig.add_axes(((graph_inches + 2 * margin) / total_width, (top - state_inches[1]) / total_height,
                             state_inches[0] / total_width, state_inches[1] / total_height))
    draw_graph(graph_ax, network, descriptions)
    draw_state_graph(state_ax, network, transitions, cycles, basins)
    state_note = ("all {} states" if full else
                  "attractors only, of {} states - each captioned with its basin size").format(2 ** n)
    fig.text(0.5, 1 - 0.15 / total_height,
             "{}\n{} nodes, {} edges; green arrow = activates, red bar = represses\n"
             "state graph ({}); bits in the order {}; blue = attractor states".format(
                 title, n, len(network.edges), state_note, " ".join(names)),
             ha='center', va='top', fontsize=9)
    fig.savefig(os.path.join(model_dir, "network.png"), dpi=140)
    fig.savefig(os.path.join(model_dir, "network.svg"), metadata={'Date': None})
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------
# readmes
# ---------------------------------------------------------------------------------------------------
FUNCTION_NOTE = ("Every function is a symmetric threshold function: it fires when at least `threshold` of "
                 "its signed inputs agree (a `+` input agreeing means it is on, a `-` input agreeing means "
                 "it is off).")


def function_table(network, descriptions):
    lines = ["## Functions", "", "| node | function | threshold |", "|---|---|---|"]
    for v in network.vertices:
        text, threshold = descriptions[v.name]
        lines.append("| `{}` | `{}` | {} |".format(v.name, text, threshold))
    return lines


def write_dynamic_readme(path, name, title, network, descriptions, cycles, basins, hamming):
    n = len(network)
    names = [v.name for v in network.vertices]
    fixed = sum(1 for cycle in cycles if len(cycle) == 1)
    lines = ["# {}".format(name), "", title, "",
             "{} nodes, {} edges. {}".format(n, len(network.edges), FUNCTION_NOTE), ""]
    lines += function_table(network, descriptions)
    lines += ["", "## State space", "",
              "All {} states enumerated: {} attractors ({} of them fixed points), mean {:.2f} nodes "
              "change per step along them.".format(2 ** n, len(cycles), fixed, hamming), "",
              "| attractor | period | basin | states ({}) |".format(", ".join(names)), "|---|---|---|---|"]
    for k, (cycle, basin) in enumerate(zip(cycles, basins)):
        shown = ["".join(str(b) for b in state) for state in cycle[:MAX_LISTED_STATES]]
        states = " -> ".join(shown) + (" -> ..." if len(cycle) > MAX_LISTED_STATES else "")
        lines.append("| {} | {} | {}/{} | {} |".format(k + 1, len(cycle), basin, 2 ** n, states))
    lines += ["", "![interaction graph and state graph](network.png)", ""]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def write_classic_readme(path, name, title, network, descriptions, cycles, basins):
    n = len(network)
    names = [v.name for v in network.vertices]
    inputs = [v.name for v in network.vertices if len(v.predecessors()) == 0]
    lines = ["# {}".format(name), "", title, "",
             "{} nodes, {} edges. {}".format(n, len(network.edges), FUNCTION_NOTE), "",
             "Feed-forward, so the dynamics are not the point: held inputs settle the circuit into a fixed "
             "point within a few steps, and all {} attractors are fixed points - one per input pattern "
             "({} input node(s): {}).".format(len(cycles), len(inputs), ", ".join(inputs)), ""]
    lines += function_table(network, descriptions)
    lines += ["", "## Truth table", "",
              "| {} | {} |".format(" | ".join(inputs),
                                   " | ".join(v for v in names if v not in inputs)),
              "|" + "---|" * len(names)]
    combos = list(itertools.product([0, 1], repeat=len(inputs)))
    for combo in combos[:MAX_TRUTH_ROWS]:
        state = [0] * n
        for name_, value in zip(inputs, combo):
            state[names.index(name_)] = value
        for _ in range(n):                      # settle: depth is at most n for an acyclic circuit
            state = [int(bool(v)) for v in network.next_state(state)]
            for name_, value in zip(inputs, combo):
                state[names.index(name_)] = value
        lines.append("| {} | {} |".format(" | ".join(str(c) for c in combo),
                                          " | ".join(str(state[names.index(v)])
                                                     for v in names if v not in inputs)))
    if len(combos) > MAX_TRUTH_ROWS:
        lines.append("| ... | ... |")
        lines.append("")
        lines.append("(first {} of {} input patterns)".format(MAX_TRUTH_ROWS, len(combos)))
    lines += ["", "![interaction graph and state graph](network.png)", ""]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


# ---------------------------------------------------------------------------------------------------
# the models
# ---------------------------------------------------------------------------------------------------
CLIQUE4_SIGNS, CLIQUE4_THRESHOLD = None, None      # chosen by search at build time, see main()

RING_MODELS = [
    ("size03_ring3", "repressilator ring of 3", lambda: [ring("a", 3)]),
    ("size04_ring3_and", "ring of 3 + AND readout",
     lambda: [ring("a", 3), node("and_ab", [("a0", 1), ("a1", 1)], 2, "a0 AND a1")]),
    ("size05_ring5", "repressilator ring of 5", lambda: [ring("a", 5)]),
    ("size06_ring5_and", "ring of 5 + AND readout",
     lambda: [ring("a", 5), node("and_ab", [("a0", 1), ("a2", 1)], 2, "a0 AND a2")]),
    ("size07_two_ring3_and", "two rings of 3, coupled by an AND readout",
     lambda: [ring("a", 3), ring("b", 3), node("and_ab", [("a0", 1), ("b0", 1)], 2, "a0 AND b0")]),
    ("size08_two_ring3_and_or", "two rings of 3, coupled by AND and OR readouts",
     lambda: [ring("a", 3), ring("b", 3), node("and_ab", [("a0", 1), ("b0", 1)], 2, "a0 AND b0"),
              node("or_ab", [("a1", 1), ("b1", 1)], 1, "a1 OR b1")]),
    ("size09_ring5_ring3_and", "rings of 5 and 3, coupled by an AND readout",
     lambda: [ring("a", 5), ring("b", 3), node("and_ab", [("a0", 1), ("b0", 1)], 2, "a0 AND b0")]),
    ("size10_ring5_ring3_and_or", "rings of 5 and 3, coupled by AND and OR readouts",
     lambda: [ring("a", 5), ring("b", 3), node("and_ab", [("a0", 1), ("b0", 1)], 2, "a0 AND b0"),
              node("or_ab", [("a1", 1), ("b1", 1)], 1, "a1 OR b1")]),
]

MOTIF_MODELS = [
    ("size03_ai_nar", "activator-inhibitor pair + a self-repressing gene",
     lambda: [activator_inhibitor("p53", "mdm2"), nar("rep")]),
    ("size04_ai_toggle", "activator-inhibitor pair + a toggle switch",
     lambda: [activator_inhibitor("p53", "mdm2"), toggle("lacI", "tetR")]),
    ("size05_ai_i1ffl_nar", "activator-inhibitor pair driving an incoherent FFL, + a self-repressing gene",
     lambda: [activator_inhibitor("p53", "mdm2"), i1_ffl("p53", "ffl_y", "pulse"), nar("rep")]),
    ("size06_toggle_ai_i1ffl", "toggle switch + activator-inhibitor pair driving an incoherent FFL",
     lambda: [toggle("lacI", "tetR"), activator_inhibitor("p53", "mdm2"),
              i1_ffl("p53", "ffl_y", "pulse")]),
    ("size07_ai_toggle_c1ffl_nar",
     "activator-inhibitor pair + toggle switch + coherent FFL + a self-repressing gene",
     lambda: [activator_inhibitor("p53", "mdm2"), toggle("lacI", "tetR"),
              c1_ffl("p53", "ffl_y", "delayed"), nar("rep")]),
    ("size08_two_ai_toggle_i1ffl", "two activator-inhibitor pairs + toggle switch + incoherent FFL",
     lambda: [activator_inhibitor("p53", "mdm2"), activator_inhibitor("nfkb", "ikb"),
              toggle("lacI", "tetR"), i1_ffl("nfkb", "ffl_y", "pulse")]),
    ("size09_ai_toggle_sim_i1ffl",
     "activator-inhibitor pair + toggle switch + a single-input module + incoherent FFL",
     lambda: [activator_inhibitor("p53", "mdm2"), toggle("lacI", "tetR"),
              sim("p53", ["sim_t1", "sim_t2", "sim_t3"]), i1_ffl("mdm2", "ffl_y", "pulse")]),
    ("size10_two_ai_toggle_i1ffl_nar",
     "two activator-inhibitor pairs + toggle switch + incoherent FFL + a self-repressing gene, "
     "with a SIM target",
     lambda: [activator_inhibitor("p53", "mdm2"), activator_inhibitor("nfkb", "ikb"),
              toggle("lacI", "tetR"), i1_ffl("nfkb", "ffl_y", "pulse"), nar("rep"),
              sim("p53", ["sim_t1"])]),
]


def dense_models():
    """Built at call time, because the clique and chorded-loop rules are chosen by a small deterministic
    search (see choose_uniform_rule) rather than hand-picked."""
    c3_signs, c3_threshold = choose_uniform_rule(
        lambda signs, threshold: clique("c", 3, signs, threshold), uniform_candidates(2))
    c4_signs, c4_threshold = choose_uniform_rule(
        lambda signs, threshold: clique("c", 4, signs, threshold), uniform_candidates(3))
    k5_signs, k5_threshold = choose_uniform_rule(
        lambda signs, threshold: chorded_loop("k", 5, signs, threshold), uniform_candidates(3))
    k6_signs, k6_threshold = choose_uniform_rule(
        lambda signs, threshold: chorded_loop("k", 6, signs, threshold), uniform_candidates(3))

    def clique3(prefix="c"):
        return clique(prefix, 3, c3_signs, c3_threshold)

    def clique4(prefix="c"):
        return clique(prefix, 4, c4_signs, c4_threshold)

    def loop5():
        return chorded_loop("k", 5, k5_signs, k5_threshold)

    def loop6():
        return chorded_loop("k", 6, k6_signs, k6_threshold)

    return [
        ("size03_clique3", "3-clique: every gene regulated by both others",
         lambda: [clique3()]),
        ("size04_clique4", "4-clique: every gene regulated by all three others",
         lambda: [clique4()]),
        ("size05_chorded_loop5", "5-node feedback loop carrying every chord",
         lambda: [loop5()]),
        ("size06_chorded_loop6", "6-node feedback loop carrying every chord",
         lambda: [loop6()]),
        ("size07_clique4_clique3", "a 4-clique and a 3-clique, side by side",
         lambda: [clique4(), clique3("e")]),
        ("size08_clique4_clique3_readout", "a 4-clique and a 3-clique, coupled by an AND readout",
         lambda: [clique4(), clique3("e"), node("out", [("c0", 1), ("e0", 1)], 2, "c0 AND e0")]),
        ("size09_two_clique4_readout", "two 4-cliques + an AND readout coupling them",
         lambda: [clique4(), clique4("d"), node("out", [("c0", 1), ("d0", 1)], 2, "c0 AND d0")]),
        ("size10_clique4_chorded6", "4-clique + a 6-node chorded loop",
         lambda: [clique4(), loop6()]),
    ]


FFL_MODELS = [
    ("ffl_{}_{}".format(kind, gate.lower()),
     "{} feed-forward loop, {} gate at z".format(
         ("coherent type-" + kind[1]) if kind[0] == "c" else ("incoherent type-" + kind[1]), gate),
     (lambda kind=kind, gate=gate: [ffl(kind, gate)]))
    for kind in ["c1", "c2", "c3", "c4", "i1", "i2", "i3", "i4"] for gate in ["AND", "OR"]
]

def sources(*names):
    """Declare input nodes (no regulators, so they hold whatever they are driven with)."""
    return [n for name in names for n in node(name, [], 0, "input")]


GATE_MODELS = [
    ("gate_and2", "two-input AND",
     lambda: [sources("a", "b"), node("z", [("a", 1), ("b", 1)], 2, "a AND b")]),
    ("gate_or2", "two-input OR",
     lambda: [sources("a", "b"), node("z", [("a", 1), ("b", 1)], 1, "a OR b")]),
    ("gate_and3", "three-input AND",
     lambda: [sources("a", "b", "c"), node("z", [("a", 1), ("b", 1), ("c", 1)], 3, "a AND b AND c")]),
    ("gate_or3", "three-input OR",
     lambda: [sources("a", "b", "c"), node("z", [("a", 1), ("b", 1), ("c", 1)], 1, "a OR b OR c")]),
    ("gate_majority3", "majority of three (the middle threshold)",
     lambda: [sources("a", "b", "c"),
              node("z", [("a", 1), ("b", 1), ("c", 1)], 2, "at least 2 of (a, b, c)")]),
    ("gate_andnot", "AND NOT: one activator, one repressor",
     lambda: [sources("a", "b"), node("z", [("a", 1), ("b", -1)], 2, "a AND NOT b")]),
    ("gate_and_of_or", "two-level: AND of two ORs",
     lambda: [sources("a", "b", "c", "d"),
              node("u", [("a", 1), ("b", 1)], 1, "a OR b"),
              node("v", [("c", 1), ("d", 1)], 1, "c OR d"),
              node("z", [("u", 1), ("v", 1)], 2, "u AND v")]),
    ("gate_or_of_and", "two-level: OR of two ANDs",
     lambda: [sources("a", "b", "c", "d"),
              node("u", [("a", 1), ("b", 1)], 2, "a AND b"),
              node("v", [("c", 1), ("d", 1)], 2, "c AND d"),
              node("z", [("u", 1), ("v", 1)], 1, "u OR v")]),
    ("gate_k_of_5", "k-of-5 thresholds side by side (k = 1, 3, 5)",
     lambda: [sources(*"abcde"),
              node("any5", [(x, 1) for x in "abcde"], 1, "at least 1 of (a, b, c, d, e)"),
              node("mid5", [(x, 1) for x in "abcde"], 3, "at least 3 of (a, b, c, d, e)"),
              node("all5", [(x, 1) for x in "abcde"], 5, "at least 5 of (a, b, c, d, e)")]),
    ("gate_cascade3", "three-step cascade, alternating sign",
     lambda: [sources("a"), node("u", [("a", 1)], 1, "a"), node("v", [("u", -1)], 1, "NOT u"),
              node("z", [("v", 1)], 1, "v")]),
    ("gate_sim", "single-input module: one regulator, three targets",
     lambda: [sources("a"), sim("a", ["t1", "t2", "t3"])]),
    ("gate_dor", "dense overlapping regulon: two regulators, three targets with different gates",
     lambda: [sources("a", "b"),
              node("t_and", [("a", 1), ("b", 1)], 2, "a AND b"),
              node("t_or", [("a", 1), ("b", 1)], 1, "a OR b"),
              node("t_andnot", [("a", 1), ("b", -1)], 2, "a AND NOT b")]),
]


def main():
    dynamic_families = [
        ("ring", RING_MODELS,
         "One construction repeated: odd repressilator rings, which give the long cycles (period 6, 10, "
         "30). The AND/OR readouts are feed-forward, so they add churn and connect the rings without "
         "changing the attractor count. Here the cycle IS the structure - see the dense family for the "
         "same behaviour without that."),
        ("motif", MOTIF_MODELS,
         "Textbook transcription motifs. The movement comes from the activator-inhibitor pairs and the "
         "NAR nodes, which is why every attractor here has period 4; the toggles supply multiplicity "
         "rather than churn, sitting at one of their fixed points in most attractors while the AI pair "
         "oscillates around them. A NAR node removes fixed-point attractors from the whole model, since "
         "it flips every step. FFLs and SIMs are feed-forward."),
        ("dense", dense_models(),
         "Multiple moving attractors with no part of the model being a ring: in a clique every node is "
         "regulated by all the others, and in a chorded loop each node is regulated by the three before "
         "it, so every cycle is a proper subgraph of a denser neighbourhood. Each clique/loop uses one "
         "uniform rule for all its nodes (chosen by a small deterministic search over sign patterns and "
         "thresholds), so the whole group is still describable in a sentence."),
    ]
    classic_families = [
        ("ffl", FFL_MODELS,
         "All eight non-equivalent three-node feed-forward loops, each with an AND and an OR gate at the "
         "target. x is an input node; y = f(x); z = g(x, y). Coherent types (C1-C4) have "
         "sign(x->z) = sign(x->y) * sign(y->z), incoherent ones (I1-I4) do not - which is what makes I1 "
         "a pulse generator and C1 a sign-sensitive delay. In the state graph x holds its value, so each "
         "basin is one value of x, and its transients show z's path to the fixed point - through the "
         "pulse or the delay."),
        ("gates", GATE_MODELS,
         "Multi-level AND/OR combinations, the middle (majority) threshold, k-of-n thresholds, an "
         "AND NOT, a cascade, a single-input module and a dense overlapping regulon."),
    ]

    for directory, collections, dynamic in ((DYNAMIC_DIR, dynamic_families, True),
                                            (CLASSIC_DIR, classic_families, False)):
        if os.path.isdir(directory):
            shutil.rmtree(directory)
        os.makedirs(directory)
        header = ("# Dynamic toy models" if dynamic else "# Classic toy models")
        blurb = ("Small hand-built models whose state space keeps moving: several attractors, essentially "
                 "none of them fixed points, so time series generated from them are not near-static."
                 if dynamic else
                 "Small hand-built feed-forward circuits - the classic ones. Their attractors are all "
                 "fixed points (one per input pattern), so they are for testing function and topology "
                 "inference rather than dynamics. Generate data from them with `only_attractors = False`: "
                 "walking to an attractor first would start every trajectory at a fixed point and make "
                 "every row identical, whereas a random start captures the transient - which is where the "
                 "pulse of an incoherent FFL and the delay of a coherent one actually show up.")
        index = [header, "",
                 "A `graphs_dir` (the role `data/cellcollective_models` plays for the real runs). {} "
                 "Every function is a symmetric threshold function.".format(blurb), "",
                 "Regenerate with `python synthetic_data_generation/make_toy_models.py`.", ""]
        for family, models, note in collections:
            index += ["## {} family".format(family), "", note, ""]
            if dynamic:
                index += ["| model | nodes | edges | attractors | fixed points | periods | "
                          "mean nodes changing per step |", "|---|---|---|---|---|---|---|"]
            else:
                index += ["| model | nodes | edges | inputs | depth | fixed points |",
                          "|---|---|---|---|---|---|"]
            for name, title, spec in models:
                network, descriptions = build(spec())
                transitions = state_space(network)
                cycles, basins = attractors_and_basins(transitions)
                hamming = mean_hamming(cycles)
                model_dir = os.path.join(directory, name)
                graph = nx.DiGraph((u.name, v.name) for u, v in network.edges)
                graph.add_nodes_from(v.name for v in network.vertices)
                if dynamic:
                    assert len(cycles) >= MIN_ATTRACTORS, "{}: only {} attractor(s)".format(name, len(cycles))
                    assert max(len(c) for c in cycles) >= MIN_LONG_CYCLE, \
                        "{}: longest attractor has period {}".format(name, max(len(c) for c in cycles))
                    assert hamming >= MIN_MEAN_HAMMING, \
                        "{}: only {:.2f} nodes change per step".format(name, hamming)
                else:
                    assert nx.is_directed_acyclic_graph(graph), "{}: classic models must be feed-forward".format(name)
                    assert all(len(c) == 1 for c in cycles), \
                        "{}: a feed-forward circuit should only have fixed points".format(name)
                network.export_to_boolean_tables(directory, name)
                network.export_to_cnet(os.path.join(model_dir, "true_network.cnet"))
                if dynamic:
                    draw_model(network, descriptions, transitions, cycles, basins, model_dir, title)
                    write_dynamic_readme(os.path.join(model_dir, "README.md"), name, title, network,
                                         descriptions, cycles, basins, hamming)
                    index.append("| [{}]({}/README.md) | {} | {} | {} | {} | {} | {:.2f} |".format(
                        name, name, len(network), len(network.edges), len(cycles),
                        sum(1 for c in cycles if len(c) == 1),
                        ", ".join(str(p) for p in sorted(set(len(c) for c in cycles))), hamming))
                    print("{:<32} {:>2} nodes {:>2} edges  {:>2} attractors, periods {}, hamming {:.2f}"
                          .format(name, len(network), len(network.edges), len(cycles),
                                  sorted(set(len(c) for c in cycles)), hamming))
                else:
                    draw_model(network, descriptions, transitions, cycles, basins, model_dir, title)
                    write_classic_readme(os.path.join(model_dir, "README.md"), name, title, network,
                                         descriptions, cycles, basins)
                    inputs = [v.name for v in network.vertices if len(v.predecessors()) == 0]
                    depth = nx.dag_longest_path_length(graph)
                    index.append("| [{}]({}/README.md) | {} | {} | {} | {} | {} |".format(
                        name, name, len(network), len(network.edges), len(inputs), depth, len(cycles)))
                    print("{:<32} {:>2} nodes {:>2} edges  {} input(s), depth {}, {} fixed points".format(
                        name, len(network), len(network.edges), len(inputs), depth, len(cycles)))
                reloaded = Network.parse_boolean_tables(model_dir)
                assert state_space(reloaded) == state_space(network), \
                    "{}: boolean-tables round trip changed dynamics".format(name)
            index.append("")
        with open(os.path.join(directory, "README.md"), "w", encoding="utf-8") as f:
            f.write("\n".join(index) + "\n")
        print("-> {} models in {}\n".format(sum(len(m) for _, m, _ in collections), directory))


if __name__ == "__main__":
    main()
