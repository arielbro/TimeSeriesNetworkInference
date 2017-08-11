import itertools
import random
import functools
import sympy
from logic import pre_randomized_boolean_symbolic_func
from utility import list_repr



class Network:

    def __init__(self, vertex_names=None, edges=None, vertex_functions=None):
        if not vertex_names:
            assert not edges
            assert not vertex_functions
        for (x,y) in edges:
            assert x in vertex_names and y in vertex_names
        if not vertex_functions:
            vertex_functions = [None] * len(vertex_names)
        self.vertices = [Vertex(self, name, function, i) for
                         name, function, i in zip(vertex_names, vertex_functions, range(len(vertex_names)))] \
            if vertex_names else []  # order important!
        self.edges = [(self.get_vertex(a), self.get_vertex(b)) for a, b in edges]

    def get_vertex(self, name):
        matches = [vertex for vertex in self.vertices if vertex.name == name]
        assert len(matches) == 1
        return matches[0]

    def __str__(self):
        return "Graph: \n\tV=" + list_repr(self.vertices) + \
               "]\n\tE=" + list_repr(self.edges) + \
               "\n\tfunctions=" + list_repr([v.function.func.func_name for v in self.vertices])

    def __repr__(self):
        return str(self)

    @staticmethod
    def generate_random(n_vertices, edge_ratio=0.5):
        vertices = list(range(n_vertices))
        edges = []
        for couple in itertools.product(range(n_vertices), range(n_vertices)):
            if random.random() < edge_ratio:
                edges.append(couple)


        G = Network(vertices, edges)
        for i in range(len(vertices)):
            # randomize boolean using a truth table
            N = len(G.vertices[i].predecessors())
            boolean_outputs = [random.choice([False, True]) for k in range(2**N)]
            G.vertices[i].function = functools.partial(pre_randomized_boolean_symbolic_func, boolean_outputs)
        return G


class Vertex:
    def __init__(self, graph, name, function, index=None):
        self.graph = graph
        self.name = name
        self.function = function
        self.index = index if index is not None else self.graph.vertices.index(self)

    def predecessors(self):
        return [u for u in self.graph.vertices if (u, self) in self.graph.edges]

    def __eq__(self, other):
        if not isinstance(other, Vertex):
            return False
        return self.name == other.name

    def __str__(self):
        # return "name:{}, function:{}".format(self.name, self.function) TODO: add dummy variables for printing
        return str(self.name)

    def __repr__(self):
        return str(self)

