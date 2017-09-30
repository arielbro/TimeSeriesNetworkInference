import itertools
import random
import re
import time

from logic import BooleanSymbolicFunc, SymmetricThresholdFunction
from utility import list_repr


class Network:

    def __init__(self, vertex_names=None, edges=None, vertex_functions=None):
        assert vertex_names
        assert len(set(vertex_names)) == len(vertex_names)
        for (x, y) in edges:
            assert (x in vertex_names and y in vertex_names)
        if not vertex_functions:
            vertex_functions = [None] * len(vertex_names)
        self.vertices = [Vertex(self, str(name), func, i) for # TODO: switch to sets? Need to sort out function mutability
                         name, func, i in zip(vertex_names, vertex_functions, range(len(vertex_names)))] \
            if vertex_names else []  # order important!
        self.edges = [(self.get_vertex(str(a)), self.get_vertex(str(b))) for a, b in edges]

    def get_vertex(self, name):
        matches = [vertex for vertex in self.vertices if vertex.name == name]
        assert len(matches) == 1
        return matches[0]

    def __str__(self):
        res = "Graph: \n\tV=" + list_repr(self.vertices) + \
               "\n\tE=" + list_repr([(a.name, b.name) for a, b in self.edges]) + \
               "\n\tfunctions:"
        for v in self.vertices:
            res += "\n\t\t f_{}: ".format(v.name)
            if len(v.predecessors()) > 0:
                res += str(v.function)
            else:
                res += "input node"
        return res

    def __repr__(self):
        return str(self)

    def __eq__(self, other): # TODO: improve from name based matching?
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

    def randomize_functions(self, restrict_signed_symmetric_threshold=False,
                            restrict_and_or_gates=False):
        for v in self.vertices:
            n = len(v.predecessors())
            if (not restrict_signed_symmetric_threshold) and (not restrict_and_or_gates):
                boolean_outputs = [random.choice([False, True]) for _ in range(2**n)]
                v.function = BooleanSymbolicFunc(boolean_outputs)
            else:
                if n == 0:
                    # either no function, or an AND and add node as its own predecessor, to assert stability.
                    v.function = None
                else:
                    signs = [random.choice([False, True]) for _ in range(n)]
                    threshold = random.randint(1, n) if not restrict_and_or_gates else random.choice([1, n])
                    v.function = SymmetricThresholdFunction(signs, threshold)

    def __add__(self, other):
        return self.union(self, other)

    @staticmethod
    def union(a, b):
        assert isinstance(a, Network)
        assert isinstance(b, Network)

        sorted_a_vertices = sorted(a.vertices, key=lambda v: v.name)
        sorted_b_vertices = sorted(b.vertices, key=lambda v: v.name)
        union_vertex_names = ["a_{}".format(v.name) for v in sorted_a_vertices] + \
                             ["b_{}".format(v.name) for v in sorted_b_vertices]
        assert len(union_vertex_names) == len(set(union_vertex_names))
        union_edges = [("a_{}".format(u.name), ("a_{}".format(v.name))) for (u, v) in a.edges] + \
                      [("b_{}".format(u.name), "b_{}".format(v.name)) for (u, v) in b.edges]
        union_functions = [v.function for v in sorted_a_vertices] + [v.function for v in sorted_b_vertices]
        return Network(vertex_names=union_vertex_names, edges=union_edges, vertex_functions=union_functions)

    @staticmethod
    def generate_random(n_vertices, indegree_bounds=[1, 5], restrict_signed_symmetric_threshold=False,
                        restrict_and_or_gates=False):
        vertices = list(range(n_vertices))
        edges = set()
        for v in vertices:
            indegree = random.randint(indegree_bounds[0], min(indegree_bounds[1], n_vertices))
            predecessors = random.sample(vertices, indegree)
            for u in predecessors:
                edges.add((u, v))

        G = Network(vertices, edges)
        G.randomize_functions(restrict_signed_symmetric_threshold, restrict_and_or_gates)
        return G

    def export_to_cnet(self, path):
        """
        Given a network G, export G's structure (including the boolean model)
        to a cnet file, as described in BNS's user manual -
        https://people.kth.se/~dubrova/BNS/user_manual.html
        :param G: A network model (assumes boolean valued functions defined for each non-input vertex
        :param path: Path to the desired cnet file
        :return: G, a graph with a boolean model.
        """
        start = time.time()
        with open(path, 'w') as cnet_file:
            cnet_file.write("# Boolean network model exported from "
                            "https://github.com/arielbro/attractor_learning\n\n")
            cnet_file.write("# total number of nodes\n.v {}\n\n".format(len(self.vertices)))
            cnet_file.write("# labels of nodes and names of corresponding components\n")
            for v in sorted(self.vertices, key=lambda v: v.index):
                cnet_file.write("# {} = {}\n".format(v.index + 1, v.name))
            cnet_file.write("\n")
            for v in sorted(self.vertices, key=lambda v: v.index):
                cnet_file.write("# {} = {}\n".format(v.index + 1, v.name))
                cnet_file.write(".n {} {}".format(v.index + 1, len(v.predecessors())))
                for u in v.predecessors():
                    cnet_file.write(" {}".format(u.index + 1))
                cnet_file.write("\n")
                for combination in itertools.product([False, True], repeat=len(v.predecessors())):
                    comb_str = "".join(map(lambda t_val: "1" if t_val else "0", combination))
                    out = v.function(*combination)
                    assert (out == True) or (out == False)
                    out = 1 if out else 0
                    cnet_file.write("{} {}\n".format(comb_str, out))
                cnet_file.write("\n")
            cnet_file.write("\n\n")
        print "time taken for graph export: {:.2f}".format(time.time() - start)

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
                assert set(indices) == set(range(1, n + 1))
                names = [re.search(r"=[ \t]*([0-9a-zA-Z_\-\\ ]+)", line).group(1) for
                         line in section.split("\n")[1:]]
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
                    bool_funcs.append(None)
                    continue
                truth_table_dict = dict()
                # input is stated in bits, with - representing wildcards (/dontcares)
                for bool_rule_str in re.findall(r"[0-9\-]+[ \t]+[01]", section):
                    output = bool(int(bool_rule_str.split()[1]))
                    ordered_input_bits = [bit for arg_index, bit in sorted(zip(v_args, bool_rule_str.split()[0]))]
                    input_value_lists = [[False, True] if bit == '-' else [bool(int(bit))]
                                         for bit in ordered_input_bits]
                    for input_combination in itertools.product(*input_value_lists):
                        truth_table_dict[tuple(input_combination)] = output
                ordered_bool_outputs = [truth_table_dict[tuple(input)] for input in
                                        itertools.product([False, True], repeat=v_n_args)]
                func = BooleanSymbolicFunc(ordered_bool_outputs)
                bool_funcs.append(func)
        assert len(bool_funcs) == n
        G = Network(vertex_names=names, edges=edges, vertex_functions=bool_funcs)
        print "time taken for graph import: {:.2f}".format(time.time() - start)
        return G


class Vertex:
    def __init__(self, graph, name, func, index=None):
        self.graph = graph
        self.name = name
        self.function = func
        self.index = index if index is not None else self.graph.vertices.index(self)
        self.precomputed_predecessors = None

    def predecessors(self):
        if not self.precomputed_predecessors:
            # search using names (asserted to be unique during init) to avoid cyrcular dependencies predecessors <> key
            name_based_edges  = [(u.name, v.name) for (u, v) in self.graph.edges]
            predecessors = [u for u in self.graph.vertices if (u.name, self.name) in name_based_edges]
            self.precomputed_predecessors = predecessors
        return self.precomputed_predecessors

    def __key(self):
        return self.name, (self.function if len(self.predecessors()) > 0 else None)

    def __eq__(self, other):
        if not isinstance(other, Vertex):
            return False
        return self.__key() == other.__key()

    def __hash__(self):
        try:
            return hash(self.__key())
        except Exception as e:
            raise e

    def __str__(self):
        # return "name:{}, function:{}".format(self.name, self.function) TODO: add dummy variables for printing
        return str(self.name)

    def __repr__(self):
        return str(self)
