import logic
import itertools
import random
import re
import time
import numpy as np
from enum import Enum
import sympy
import os
import csv

from logic import BooleanSymbolicFunc, SymmetricThresholdFunction
from utility import list_repr


class FunctionTypeRestriction(Enum):
    NONE = 1
    SYMMETRIC_THRESHOLD = 2
    SIMPLE_GATES = 3


class Network:
    def __init__(self, vertex_names=None, edges=None, vertex_functions=None):
        assert vertex_names
        assert len(set(vertex_names)) == len(vertex_names)
        if vertex_functions:
            assert len(vertex_names) == len(vertex_functions)
        for (x, y) in edges:
            assert (x in vertex_names and y in vertex_names)
        if not vertex_functions:
            vertex_functions = [None] * len(vertex_names)
        self.vertices = [Vertex(self, str(name), func, i) for
                         # TODO: switch to sets? Need to sort out function mutability
                         name, func, i in zip(vertex_names, vertex_functions, range(len(vertex_names)))] \
            if vertex_names else []  # order important!
        self.edges = [(self.get_vertex(str(a)), self.get_vertex(str(b))) for a, b in edges]
        for v in self.vertices:
            if len(v.predecessors()) == 0 and v.function not in [None, lambda _: False, lambda _:True, False, True]:
                print "warning, input node {} created with non-constant function {}. Removing function".\
                      format(v.name, v.function)
                v.function = None

    def get_vertex(self, name):
        matches = [vertex for vertex in self.vertices if vertex.name == name]

        assert len(matches) == 1
        return matches[0]

    def __len__(self):
        return len(self.vertices)

    def __str__(self):
        res = "Graph: \n\tV=" + list_repr(self.vertices) + \
              "\n\tE=" + list_repr([(a.name, b.name) for a, b in self.edges]) + \
              "\n\tfunctions:"
        for v in self.vertices:
            res += "\n\t\t f_{}: ".format(v.name)
            if len(v.predecessors()) > 0:
                res += str(v.function)
            else:
                res += "input node, value={}".format(v.function)
        return res

    def __repr__(self):
        return str(self)

    def __eq__(self, other):  # TODO: improve from name based matching?
        """
        Checks equality between two networks.
        Networks are equal if they have the same vertex names, the same edges between vertices,
        and the same functions on non-input nodes.
        :param other:
        :return:
        """
        if not isinstance(other, Network):
            return False
        if set(self.vertices) != set(other.vertices) or set(self.edges) != set(other.edges):
            return False
        return True

    def __ne__(self, other):
        return not self == other

    def randomize_edges(self, include_self_loops=False):
        """
        Replaces the ingoing edges of nodes by choosing uniformly (and without replacement) the same number
        of in-neighbors. Operates on graph inplace.
        :return:
        """
        # TODO: allow outgoing edges randomization, maybe a hybrid between the two.

        in_degrees = []
        for v in self.vertices:
            in_degrees.append(len(v.predecessors()))
            v.precomputed_predecessors = None
        new_edges = []
        for in_degree, v in zip(in_degrees, self.vertices):
            possible_neighbors = self.vertices if include_self_loops else \
                (self.vertices[:v.index] + self.vertices[v.index + 1:])
            in_neighbors = np.random.choice(possible_neighbors, in_degree, replace=False)
            new_edges.extend([(neighbor, v) for neighbor in in_neighbors])

            # BooleanSymbolicFunctions hold input names, so we need to recreate them
            if isinstance(v.function, BooleanSymbolicFunc):
                v.function = BooleanSymbolicFunc(input_names=[neighbor.name for neighbor in in_neighbors],
                                                 boolean_outputs=v.function.get_truth_table_outputs())

        self.edges = new_edges

    def randomize_functions(self, function_type_restriction=FunctionTypeRestriction.NONE, mutate_inputs=False):
        for v in self.vertices:
            n = len(v.predecessors())
            if n == 0 and mutate_inputs:
                v.function = random.choice([False, True])
            elif n == 0:
                v.function = None
            elif function_type_restriction == FunctionTypeRestriction.NONE:
                boolean_outputs = [random.choice([False, True]) for _ in range(2 ** n)]
                input_names = [u.name for u in v.predecessors()]
                v.function = BooleanSymbolicFunc(input_names=input_names, boolean_outputs=boolean_outputs)
            elif function_type_restriction == FunctionTypeRestriction.SYMMETRIC_THRESHOLD:
                signs = [random.choice([False, True]) for _ in range(n)]
                threshold = random.randint(1, n) if \
                    function_type_restriction != FunctionTypeRestriction.SIMPLE_GATES else \
                    random.choice([1, n])
                v.function = SymmetricThresholdFunction(signs, threshold)
            elif function_type_restriction == FunctionTypeRestriction.SIMPLE_GATES:
                raise NotImplementedError("randomizing simple gates funntions not yet implemented")
            else:
                raise ValueError("unknown function type restriction: {}".format(function_type_restriction))

    def __add__(self, other):
        return self.union(self, other)

    def next_state(self, vertex_states, return_as_string=False):
        v_states = []
        assert len(vertex_states) == len(self.vertices)
        if isinstance(vertex_states, str):
            vertex_states = [c for c in vertex_states]
        for v, v_state in zip(self.vertices, vertex_states):
            assert isinstance(v_state, bool) or v_state in [0, 1, "0", "1", sympy.false, sympy.true]
            v_states.append(True if v_state in [True, 1, "1", sympy.true] else False)
        v_next_states = []
        for v_index, v in enumerate(self.vertices):
            input_values = [v_states[u.index] for u in v.predecessors()]
            if len(input_values) == 0:
                v_next_states.append(v_states[v_index])
            else:
                v_next_states.append(v.function(*input_values))
        if return_as_string:
            res = ""
            for state in v_next_states:
                assert state in [False, True, sympy.false, sympy.true]
                res += '1' if state in [True, sympy.true] else '0'
            return res
        else:
            for state in v_next_states:
                assert state in [False, True, sympy.false, sympy.true, 0, 1]
            return tuple(1 if s else 0 for s in v_next_states)

    @staticmethod
    def union(a, b):
        assert isinstance(a, Network)
        assert isinstance(b, Network)

        sorted_a_vertices = sorted(a.vertices, key=lambda vertex: vertex.name)
        sorted_b_vertices = sorted(b.vertices, key=lambda vertex: vertex.name)
        union_vertex_names = ["a_{}".format(v.name) for v in sorted_a_vertices] + \
                             ["b_{}".format(v.name) for v in sorted_b_vertices]
        assert len(union_vertex_names) == len(set(union_vertex_names))
        union_edges = [("a_{}".format(u.name), ("a_{}".format(v.name))) for (u, v) in a.edges] + \
                      [("b_{}".format(u.name), "b_{}".format(v.name)) for (u, v) in b.edges]
        union_functions = [v.function for v in sorted_a_vertices] + [v.function for v in sorted_b_vertices]
        return Network(vertex_names=union_vertex_names, edges=union_edges, vertex_functions=union_functions)

    # TODO: generate scale-free graphs
    @staticmethod
    def generate_random(n_vertices, indegree_bounds=(1, 5), function_type_restriction=FunctionTypeRestriction.NONE,
                        indegree_geometric_p=None):
        vertices = list(range(n_vertices))
        edges = set()
        for v in vertices:
            if indegree_bounds is not None:
                indegree = random.randint(min(indegree_bounds[0], n_vertices), min(indegree_bounds[1], n_vertices))
            elif indegree_geometric_p is not None:
                indegree = min([np.random.geometric(indegree_geometric_p), len(vertices)])
            else:
                raise ValueError("either indegree_bounds or indegree_exponent need to be set")
            predecessors = random.sample(vertices, indegree)
            for u in predecessors:
                edges.add((u, v))

        G = Network(vertices, edges)
        G.randomize_functions(function_type_restriction=function_type_restriction)
        return G

    def export_to_cnet(self, path, verbose=False):
        """
        Given a network G, export G's structure (including the boolean model)
        to a cnet file, as described in BNS's user manual -
        https://people.kth.se/~dubrova/BNS/user_manual.html
        :param path: Path to the desired cnet file
        :return: G, a graph with a boolean model.
        """
        start = time.time()
        with open(path, 'w') as cnet_file:
            cnet_file.write("# Boolean network model exported from "
                            "https://github.com/arielbro/attractor_learning\n\n")
            cnet_file.write("# total number of nodes\n.v {}\n\n".format(len(self.vertices)))
            cnet_file.write("# labels of nodes and names of corresponding components\n")
            for v in sorted(self.vertices, key=lambda vertex: vertex.index):
                cnet_file.write("# {} = {}\n".format(v.index + 1, v.name))
            cnet_file.write("\n")
            for v in sorted(self.vertices, key=lambda vertex: vertex.index):
                cnet_file.write("# {} = {}\n".format(v.index + 1, v.name))
                cnet_file.write(".n {} {}".format(v.index + 1, len(v.predecessors())))
                for u in sorted(v.predecessors(), key=lambda vertex: vertex.index):
                    cnet_file.write(" {}".format(u.index + 1))
                cnet_file.write("\n")
                if (len(v.predecessors()) != 0) and v.function is None:
                    raise NotImplementedError("Can't export a graph with a non-fixed function.")
                if (len(v.predecessors()) == 0) and v.function is not None:
                    is_true = v.function == True or (isinstance(v.function, type(lambda _:_)) and v.function() == True)
                    out = 1 if is_true else 0
                    cnet_file.write("{}\n".format(out))
                elif len(v.predecessors()) != 0:
                    # Assumption - v.function inputs are ordered by index (not e.g. name).
                    rows = list(itertools.product([False, True], repeat=len(v.predecessors())))
                    if isinstance(v.function, BooleanSymbolicFunc):
                        row_outputs = v.function.get_truth_table_outputs()
                    else:
                        row_outputs = [v.function(*row) for row in rows]
                    for row, out in zip(rows, row_outputs):
                        comb_str = "".join(map(lambda t_val: "1" if t_val else "0", row))
                        assert (out == True) or (out == False)
                        out = 1 if out else 0
                        cnet_file.write("{} {}\n".format(comb_str, out))
                cnet_file.write("\n")
            cnet_file.write("\n\n")
        # print "time taken for graph export: {:.2f}".format(time.time() - start)

    @staticmethod
    def parse_cnet(path):
        """
        Given a cnet file, as described in BNS's user manual -
        https://people.kth.se/~dubrova/BNS/user_manual.html
        parse the file and produce a Network instance graph & boolean model.
        :param path: Path to the cnet file
        :return: G, a graph with a boolean model.
        """
        start = time.time()
        with open(path, 'r') as cnet_file:
            text = cnet_file.read()
        # TODO: if runtime is long, stop searching each section for each pattern and use your memory
        sections = re.split("\n\s*\n", text)
        bool_funcs = []
        edges = []
        last_v_index = -1
        for section in sections:
            if ".v" in section:
                # graph size
                n = int(re.search(r"[0-9]+", section).group())
                continue
            if "labels of nodes" in section:
                # vertex names
                indices = [int(re.search(r"([0-9]+)[ \t]*=", line).group(1)) for
                           line in section.split("\n")[1:]]
                # noinspection PyUnboundLocalVariable
                assert set(indices) == set(range(1, n + 1))
                names = [re.search(r"=[ \t]*([0-9a-zA-Z_\-\\ ]+)", line).group(1) for
                         line in section.split("\n")[1:]]
                # names with spaces cause troubles when making sympy vars
                names = [name.replace(" ", "_") for name in names]
                continue
            if ".n" in section:
                # boolean functions
                v_index_str, v_n_args_str, v_args_str = re.search(
                    r"\.n[ \t]*([0-9]+)[ \t]+([0-9]+)((?:[ \t][0-9]+)*)", section).groups()
                v_index = int(v_index_str) - 1  # cnet indices are 1 based, convert.
                assert v_index == last_v_index + 1
                last_v_index = v_index
                v_n_args = int(v_n_args_str)
                v_args = [int(arg_str) - 1 for arg_str in v_args_str.split()]
                assert len(v_args) == v_n_args
                # for i in range(len(v_args) - 1):  # the logic module assumes inputs are ordered by index
                #     assert v_args[i] < v_args[i + 1]
                edges.extend([(names[arg], names[v_index]) for arg in v_args])
                if v_n_args == 0:
                    if len(section.split("\n")[-1]) == 1:
                        val = section.split("\n")[-1]
                        assert val == '0' or val == '1'
                        bool_funcs.append(bool(int(val)))
                    else:
                        bool_funcs.append(None)
                    continue
                truth_table_dict = dict()
                # input is stated in bits, with - representing wildcards (/dontcares)
                for bool_rule_str in re.findall(r"[0-9\-]+[ \t]+[01]", section):
                    output = bool(int(bool_rule_str.split()[1]))
                    ordered_input_bits = [bit for arg_index, bit in zip(v_args, bool_rule_str.split()[0])]
                    input_value_lists = [[False, True] if bit == '-' else [bool(int(bit))]
                                         for bit in ordered_input_bits]
                    for input_combination in itertools.product(*input_value_lists):
                        truth_table_dict[tuple(input_combination)] = output
                ordered_bool_outputs = [truth_table_dict[tuple(input_value)] for input_value in
                                        itertools.product([False, True], repeat=v_n_args)]
                input_names = [names[arg] for arg in v_args]
                func = BooleanSymbolicFunc(input_names=input_names, boolean_outputs=ordered_bool_outputs)
                bool_funcs.append(func)
        assert len(bool_funcs) == n
        G = Network(vertex_names=names, edges=edges, vertex_functions=bool_funcs)
        print "time taken for graph import: {:.2f} seconds".format(time.time() - start)
        return G

    def export_to_boolean_tables(self, base_path, model_name):
        """
        Exports a network to (only possibly compatible) cellcollective's truth tables format.
        Given a base directory and a name for the model, exports it to a set of boolean table files in the following
        format - a directory is created with the model's name.
        Within, a csv (with tab separators) file is created for each vertex. The file's name is the vertex index
        (1 based). First row lists indices of predecessors and the node's index, others are truth table rows
        with predecessors' values and output value.
        Additionally, creates an empty external_components.ALL.txt file (we don't distinguish nodes) and a SPECIES_KEY
        csv file. This is a tab separated file, where each row has a node index and name.
        :param base_path:
        :param model_name:
        :return:
        """
        os.mkdir(os.path.join(base_path, model_name))
        open(os.path.join(base_path, model_name, "external_components.ALL.txt"), 'a').close()

        # create name mapping
        row_tuples = [(v.index + 1, v.name) for v in self.vertices]
        with open(os.path.join(base_path, model_name, "SPECIES_KEY.csv"), 'w') as mapping_file:
            writer = csv.writer(mapping_file, delimiter='\t')
            for row in row_tuples:
                writer.writerow(row)

        # write truth tables
        for v in self.vertices:
            row_tuples = [tup + (int(v.function(*tup)),) for tup in itertools.product(
                [0, 1], repeat=len(v.predecessors()))]
            with open(os.path.join(base_path, model_name, "{}.csv".format(v.index + 1)), 'w') as table_file:
                writer = csv.writer(table_file, delimiter='\t')
                writer.writerow([u.index + 1 for u in v.predecessors()] + [v.index + 1])
                for row in row_tuples:
                    writer.writerow(row)

    @staticmethod
    def parse_boolean_tables(path):
        """
        Parses cellcollective's truth tables format directory and returns a corresponding network.
        See export_toboolean_tables for format information. Note that this loses information about
        external components (we don't differentiate nodes).
        :param path:
        :return:
        """
        with(open(os.path.join(path, "SPECIES_KEY.csv"), 'r')) as mapping_file:
            reader = csv.reader(mapping_file, delimiter="\t")
            lines = list(reader)
            vertex_indices = [tup[0] for tup in lines]
            vertex_names = [tup[1] for tup in lines]
            indices_to_names = {tup[0]: tup[1] for tup in lines}

        edges = []
        functions = []
        for v_index, v_name in zip(vertex_indices, vertex_names):
            truth_table_path = os.path.join(path, "{}.csv".format(v_index))
            if not os.path.exists(truth_table_path):
                # an external species. For us this eans input node.
                functions.append(None)
            else:
                with(open(truth_table_path, 'r')) as truth_table_file:
                    reader = csv.reader(truth_table_file, delimiter="\t")
                    lines = list(reader)
                    # add edges
                    predecessor_indices = lines[0][:-1]
                    for predecessor_index in predecessor_indices:
                        edges.append((indices_to_names[predecessor_index], indices_to_names[v_index]))

                    # build function
                    # The order of predecessors might not be sorted same as in the SPECIES_KEY file, so for the
                    # function to be defined correctly (in ascending model index of variables)
                    # we need to permute the truth table rows.
                    n_vars = len(predecessor_indices)
                    sorted_indices = sorted(range(n_vars),
                                            key=lambda i: vertex_indices.index(predecessor_indices[i]))
                    outputs = [None] * (2 ** n_vars)
                    for truth_table_line in lines[1:]:
                        # MSB is first var
                        permuted_row_index = sum(2**(n_vars - i - 1) for i in range(n_vars)
                                                 if bool(int(truth_table_line[sorted_indices[i]])))
                        assert outputs[permuted_row_index] is None
                        outputs[permuted_row_index] = bool(int(truth_table_line[-1]))
                    for out in outputs:
                        assert out is not None
                    predecessor_names = [indices_to_names[predecessor_indices[sorted_indices[i]]] for
                                         i in range(n_vars)]
                    functions.append(BooleanSymbolicFunc(input_names=predecessor_names, boolean_outputs=outputs))
        return Network(vertex_names=vertex_names, edges=edges, vertex_functions=functions)

    def contract_vertex(self, vertex):
        """
        Assumes the vertex has no self loop, and currently assumes it has no successors.
        Removes all records for its existance
        :param vertex: A vertex (assumed to be contained in self.vertices) with no self loop or successors
        :return:
        """
        # Removes a vertex from the graph, and update its predecessors and successors.
        # Before that, if the vertex is an output node, the precomputed_successors field of its predecessors is reset.
        # All edges touching it are either removed (if an output vertex), or replaced in the following way -
        # if (u, vertex) and (vertex, v) are both in self.edges, an edge (u, v) is created.

        assert vertex in self.vertices
        # assert (vertex, vertex) not in self.edges
        assert vertex not in (u for (u, v) in self.edges)
        for u in vertex.predecessors():
            u.precomputed_successors = None
        self.edges = [(u, v) for (u, v) in self.edges if v is not vertex]
        self.vertices.remove(vertex)
        return

    def convert_inputs_to_loops(self):
        """
        For each input node, adds a self loop from it to itself, with the id function, so Dubrova will be able to
        run on it. Note that this will break any ILP where the node functions are allowed to be variables.
        :return:
        """
        # TODO: somehow enforce checking of this when variable functions are used in ILP
        for v in self.vertices:
            if len(v.predecessors()) == 0:
                v.precomputed_predecessors = None
                self.edges.append((v, v))
                v.function = sympy.And

    def copy(self):
        """
        returns a copy of self, assuming functions are immutable.
        :return:
        """
        return Network(vertex_names=[v.name for v in self.vertices], edges=[(u.name, v.name) for (u, v) in self.edges],
                       vertex_functions=[v.function for v in self.vertices])

    def __mul__(self, other):
        """
        Computes a composition of self's functions with other's functions.
        For each vertex, transforms its function to be taking values from its predecessors' predessesors.
        Only defined if the graphs share nodes.
        :param other:
        :return:
        """
        vertex_names = []
        edges = []
        functions = []
        source_functions = {v.name: v.function for v in self.vertices}
        any_converted = False
        for v in self.vertices:
            if not isinstance(v.function, logic.BooleanSymbolicFunc):
                if not isinstance(v.function, sympy.FunctionClass):
                    raise NotImplementedError("Multiplication of graphs with generic functions not yet implemented")
                predecessors_names = [u.name for u in v.predecessors()]
                source_functions[v.name] = logic.BooleanSymbolicFunc.from_sympy_func(v.function, predecessors_names)
                any_converted = True
        if any_converted:
            print "warning - graph with generic function types passed to __mul__, converting to BooleanSymbolicFunc"

        for v in self.vertices:
            vertex_names.append(v.name)
            predecessors_funcs = [source_functions[u.name] for u in v.predecessors()]
            functions.append(source_functions[v.name].compose(input_funcs=predecessors_funcs, simplify=True))
            new_predecessors = sorted([x.name for x in functions[-1].formula.free_symbols])
            edges.extend([(u_name, v.name) for u_name in new_predecessors])
        return Network(vertex_names, edges, functions)

    def __pow__(self, power, modulo=None):
        """
        Exponential by squaring, using the preivously defined __mul__ operation.
        :param power:
        :param modulo:
        :return:
        """
        if modulo is not None:
            raise NotImplementedError("can't modulo a graph")
        if power == 1:
            return self.copy()
        if power == 2:
            return self * self
        if power % 2 == 0:
            return (self * self) ** (power / 2)
        else:
            return self * ((self * self) ** ((power - 1) /2))


class Vertex:
    def __init__(self, graph, name, func, index=None):
        self.graph = graph
        self.name = name
        self.function = func
        self.index = index if index is not None else self.graph.vertices.index(self)
        self.precomputed_predecessors = None
        self.precomputed_successors = None

    def predecessors(self):
        if not self.precomputed_predecessors:
            # search using names (asserted to be unique during init) to avoid circular dependencies predecessors <> key
            name_based_edges = [(u.name, v.name) for (u, v) in self.graph.edges]
            predecessors = sorted([u for u in self.graph.vertices if (u.name, self.name) in name_based_edges],
                                  key=lambda vertex: vertex.index)
            self.precomputed_predecessors = predecessors
        return self.precomputed_predecessors

    def successors(self):
        if not self.precomputed_successors:
            # search using names (asserted to be unique during init) to avoid circular dependencies predecessors <> key
            name_based_edges = [(u.name, v.name) for (u, v) in self.graph.edges]
            successors = [v for v in self.graph.vertices if (self.name, v.name) in name_based_edges]
            self.precomputed_successors = successors
        return self.precomputed_successors

    def _key(self):
        if self.function is None:
            return self.name, None
        if len(self.predecessors()) == 0:
            if callable(self.function):
                return self.name, self.function()
            elif self.function in (False, True, sympy.false, sympy.true, 0 ,1):
                return bool(self.function)
            else:
                raise ValueError("An input node should not have a non function, "
                                 "non boolean and non None function field.")
        outputs = tuple(self.function(*row) for row in itertools.product([False, True],
                                                                         repeat=len(self.predecessors())))
        return self.name, outputs

    def __eq__(self, other):
        if not isinstance(other, Vertex):
            return False
        return self._key() == other._key()

    def __hash__(self):
        try:
            return hash(self._key())
        except Exception as e:
            raise e

    def __str__(self):
        # return "name:{}, function:{}".format(self.name, self.function) TODO: add dummy variables for printing
        return str(self.name)

    def __repr__(self):
        return str(self)
