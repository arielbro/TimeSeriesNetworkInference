import attractors, graphs, logic, stochastic, cnet_parser, ilp, sympy

# G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
#                    vertex_functions=[sympy.And])
#
G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                   vertex_functions=[sympy.Nand, sympy.And])
#
# G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
#                    vertex_functions=[sympy.Nand, sympy.Nand])
#
# G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
#                    vertex_functions=[lambda x: True, lambda x: False])
#
# G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "A"), ("C", "A")],
#                    vertex_functions=[sympy.Nand]*3)
#
# G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
#                    vertex_functions=[sympy.Nand]*3)
#
# G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("A", "C"), ("B", "C"), ("B", "A"),
#                                                         ("C", "A"), ("C", "B")],
#                    vertex_functions=[sympy.Nor]*3)
#
# G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("A", "C"), ("B", "C"), ("B", "A"),
#                                                         ("C", "A"), ("C", "B")],
#                    vertex_functions=[sympy.And] + [sympy.Nor] * 2)
# G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("A", "C"), ("B", "C"), ("B", "A"),
#                                                         ("C", "A"), ("C", "B")],
#                    vertex_functions=[lambda y, z: z, lambda x, z: sympy.logic.Not(x),
#                                      lambda x, y: sympy.logic.Or(sympy.logic.Not(x),  y)])
#
# # acyclic, should have 2**#input_nodes attractors of length 1
# G = graphs.Network(vertex_names=["v1", "v2", "v3", "v4", "v5", "v6"],
#                    edges=[("v1", "v4"), ("v2", "v4"), ("v1", "v5"), ("v4","v6")],
#                    vertex_functions=[lambda *args: sympy.Nand(*args)]*6)
#
# G.randomize_functions(restrict_signed_symmetric_threshold=True)
# for experiment in range(20):
#     G.randomize_functions()
#     stochastic.estimate_attractors(G, n_walks=100, max_walk_len=100)
#
# G = graphs.Network.generate_random(10, indegree_bounds=[1, 5], restrict_signed_symmetric_threshold=True)
# print G
# attractors.find_num_attractors_multistage(G, use_ilp=True)
# attractors.find_min_attractors_model(G)
# attractors.find_num_attractors_onestage(G, max_len=10, max_num=10, use_sat=False, verbose=False)
# attractors.stochastic_attractor_estimation(G, n_walks=100, max_walk_len=100)
# attractors.write_sat_sampling_analysis_table(10, 7, "C:/Users/Ariel/Downloads/graph_sampling.csv")
# attractors.write_random_graph_estimations_sampling(n_graphs=400, vertices_bounds=[3, 100],
#                                         indegree_bounds=[0, 20], restrict_symmetric_threshold=True,
#                                         n_walks=300, max_walk_len=300,
#                                         path="C:/Users/Ariel/Downloads/graph_sampling_symmetric_with_input_nodes.csv")
# attractors.write_random_graph_estimations_sampling(n_graphs=400, vertices_bounds=[3, 100],
#                                         indegree_bounds=[0, 5], restrict_symmetric_threshold=False,
#                                         n_walks=300, max_walk_len=300,
#                                         path="C:/Users/Ariel/Downloads/graph_sampling.csv")
# attractors.write_random_graph_estimations_sampling(n_graphs=400, vertices_bounds=[3, 100],
#                                         indegree_bounds=[0, 5], restrict_symmetric_threshold=True,
#                                         restrict_and_or_gates=True,
#                                         n_walks=300, max_walk_len=300,
#                                         path="C:/Users/Ariel/Downloads/graph_sampling_only_simple_gates.csv")

# G = cnet_parser.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
#                "\\Attractors - for Ariel\\BNS_Dubrova_2011\\arabidopsis.cnet.txt")
input_nodes = [u for u in G.vertices if len(u.predecessors()) == 0]
# stochastic_estimation = stochastic.estimate_attractors(G, n_walks=300, max_walk_len=30)
ilp_estimation = attractors.find_num_attractors_onestage(G, max_len=1, max_num=1, use_sat=False, verbose=True)
pass