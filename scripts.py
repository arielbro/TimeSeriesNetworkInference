import time
import attractors, graphs, logic, stochastic, utility, ilp, sympy


def estimate_size(n, m, T, P):
    # assumes find_symmetric_threshold model, and no input vertices
    vars = P*(T + 1)*(n + 2) + n + 2*m - 2*P
    # constr = 4*m + 2*(n+1)*(T+1)*P + 2*n*T*P + 2*(2*n+1)*(T-1)*P + 2*P + 1
    constr = 8*n*T*P - 2*n*P + 4*T*P + 4*m + 2*P + 1
    print "vars={}, constrs={}".format(vars, constr)


# G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
#                    vertex_functions=[sympy.Nand])
#
# G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
#                    vertex_functions=[sympy.Nand, sympy.And])

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
#
# start = time.time()
# conversions = 0
# for v in G.vertices:
#     try:
#         v.function = logic.SymmetricThresholdFunction.from_function(v.function, len(v.predecessors()))
#         print v.function
#         conversions += 1
#     except ValueError as e:
#         pass
# print "converted {} out of {} functions. Time taken: {:.2f}".format(
#     conversions, len(G.vertices), time.time() - start)

# for experiment in range(20):
#     G.randomize_functions()
#     stochastic.estimate_attractors(G, n_walks=100, max_walk_len=100)
#
G = graphs.Network.generate_random(20, indegree_bounds=[1, 10], restrict_signed_symmetric_threshold=True)
# print G
# attractors.find_num_attractors_multistage(G, use_ilp=False)
# attractors.find_num_attractors_onestage(G, use_sat=False, max_len=10, max_num=10)
# attractors.find_min_attractors_model(G)

# times = []
# for i in range(20):
#     start_time = time.time()
#     G = graphs.Network.generate_random(7, indegree_bounds=[1, 5], restrict_signed_symmetric_threshold=True)
#     attractors.find_num_attractors_onestage(G, max_len=5, max_num=5, use_sat=True, verbose=False)
#     times.append(time.time() - start_time)
# ave_sat = sum(times) / float(len(times))
# times = []
# for i in range(20):
#     start_time = time.time()
#     G = graphs.Network.generate_random(7, indegree_bounds=[1, 5], restrict_signed_symmetric_threshold=True)
#     attractors.find_num_attractors_onestage(G, max_len=5, max_num=5, use_sat=False, verbose=False)
#     times.append(time.time() - start_time)
# print "average run time for sat-based:{:.2f}".format(ave_sat)
# print "average run time for direct ILP:{:.2f}".format(sum(times) / float(len(times)))

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

# G = utility.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
#                "\\Attractors - for Ariel\\BNS_Dubrova_2011\\MAPK_large2.cnet")
# input_nodes = [u for u in G.vertices if len(u.predecessors()) == 0]
# stochastic_estimation = stochastic.estimate_attractors(G, n_walks=300, max_walk_len=30)
ilp_estimation = attractors.find_num_attractors_onestage(G, max_len=20, max_num=20, use_sat=False, verbose=True)
# attractors.write_random_fixed_graph_estimations_sampling(G=G, n_iter=400, restrict_symmetric_threshold=True,
#                                                          restrict_and_or_gates=True,
#                                                          n_walks=1500, max_walk_len=1000,
#                                                          path="C:/Users/Ariel/Downloads/MAPK_large2_sampling_gates.csv")

# print attractors.find_num_attractors_dubrova(G, "")
# print estimate_size(len(G.vertices), len(G.edges), 4, 200)
