import itertools
import random
import functools
import sympy
from logic import PreRandomizedBooleanSymbolicFunc, PreRandomizedSymmetricThresholdFunction
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
        res = "Graph: \n\tV=" + list_repr(self.vertices) + \
               "]\n\tE=" + list_repr(self.edges) + \
               "\n\tfunctions:"
        for v in self.vertices:
            res += "\n\t\t"
            if v.function is not None:
                res += str(v.function)
            else:
                res += "input node"
        return res

    def __repr__(self):
        return str(self)

    def randomize_functions(self, restrict_signed_symmetric_threshold=False):
        for v in self.vertices:
            n = len(v.predecessors())
            if not restrict_signed_symmetric_threshold:
                boolean_outputs = [random.choice([False, True]) for _ in range(2**n)]
                v.function = PreRandomizedBooleanSymbolicFunc(boolean_outputs)
            else:
                if n == 0:
                    # either no function, or an AND and add node as its own predecessor, to assert stability.
                    v.function = None
                else:
                    signs = [random.choice([False, True]) for _ in range(2**n)]
                    threshold = random.randint(1, n)
                    v.function = PreRandomizedSymmetricThresholdFunction(signs, threshold)

    @staticmethod
    def generate_random(n_vertices, indegree_bounds=[1, 5], restrict_signed_symmetric_threshold=False):
        vertices = list(range(n_vertices))
        edges = []
        for v in vertices:
            indegree = random.randint(indegree_bounds[0], min(indegree_bounds[1], n_vertices))
            predecessors = random.sample(vertices, indegree)
            for u in predecessors:
                edges.append((u, v))

        G = Network(vertices, edges)
        G.randomize_functions(restrict_signed_symmetric_threshold)
        return G


class Vertex:
    def __init__(self, graph, name, function, index=None):
        self.graph = graph
        self.name = name
        self.function = function
        self.index = index if index is not None else self.graph.vertices.index(self)
        self.precomputed_predecessors = None

    def predecessors(self):
        if not self.precomputed_predecessors:
            predecessors = [u for u in self.graph.vertices if (u, self) in self.graph.edges]
            self.precomputed_predecessors = predecessors
        return self.precomputed_predecessors

    def __eq__(self, other):
        if not isinstance(other, Vertex):
            return False
        return self.name == other.name

    def __str__(self):
        # return "name:{}, function:{}".format(self.name, self.function) TODO: add dummy variables for printing
        return str(self.name)

    def __repr__(self):
        return str(self)

