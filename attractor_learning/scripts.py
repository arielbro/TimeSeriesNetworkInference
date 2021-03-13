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

# G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
#                    vertex_functions=[sympy.And])

# G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
#                    vertex_functions=[lambda _: True])

# G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
#                    vertex_functions=[sympy.Nand, sympy.And])

# G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
#                    vertex_functions=[sympy.Nand, sympy.Nand])
#
# G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
#                    vertex_functions=[lambda x: True, lambda x: False])

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
# acyclic, should have 2**#input_nodes attractors of length 1
# G = graphs.Network(vertex_names=["v1", "v2", "v3", "v4", "v5", "v6"],
#                    edges=[("v1", "v4"), ("v2", "v4"), ("v1", "v5"), ("v4","v6")],
#                    vertex_functions=[sympy.Nand]*6)
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
# G = graphs.Network.generate_random(5, indegree_bounds=[1, 4],
#                                    function_type_restriction=graphs.FunctionTypeRestriction.NONE)

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

# G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
#                "\\Attractors - for Ariel\\BNS_Dubrova_2011\\MAPK_large2.cnet")
# G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
#                "\\Attractors - for Ariel\\BNS_Dubrova_2011\\arabidopsis.cnet.txt")
# G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
#                "\\Attractors - for Ariel\\BNS_Dubrova_2011\\EGFR_man_with_inputs_all_zero.cnet")
# G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
#                "\\Attractors - for Ariel\\BNS_Dubrova_2011\\MAPK_large.cnet")
# G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
#                "\\Attractors - for Ariel\\BNS_Dubrova_2011\\thelper.cnet.txt")
# G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
#                "\\Attractors - for Ariel\\BNS_Dubrova_2011\\tcr.cnet")

# print G
# print len(G.vertices)
# input_nodes = [u for u in G.vertices if len(u.predecessors()) == 0]
# stochastic_estimation = stochastic.estimate_attractors(G, n_walks=300, max_walk_len=30)
# ilp_estimation = attractors.find_num_attractors_onestage(G, max_len=2, max_num=2, use_sat=False, verbose=True)
# attractors.write_random_fixed_graph_estimations_sampling(G=G, n_iter=400, restrict_symmetric_threshold=True,
#                                                          restrict_and_or_gates=True,
#                                                          n_walks=1500, max_walk_len=1000,
#                                                          path="C:/Users/Ariel/Downloads/MAPK_large2_sampling_gates.csv")

# print attractors.find_num_attractors_dubrova(G, "")
# print estimate_size(len(G.vertices), len(G.edges), 4, 200)

# attractors.find_max_attractor_model(G, verbose=False,
#                                     model_type_restriction=graphs.FunctionTypeRestriction.NONE,
#                                     attractor_length_threshold=100, attractor_num_threshold=100,
#                                     use_state_keys=True)
# print G
# attractors.find_num_attractors_multistage(G, use_ilp=False)
# attractors.find_num_attractors_onestage(G, use_sat=False, max_len=1, max_num=30,
#                                         verbose=True)
# attractors.find_min_attractors_model(G, max_len=12, min_attractors=2)
#
# G = graphs.Network(
#     vertex_names=['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16',
#                   '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31',
#                   '32', '33', '34'],
#     edges=[('1', '2'), ('2', '16'), ('3', '17'), ('5', '15'), ('6', '29'), ('7', '28'), ('8', '22'),
#            ('9', '28'), ('10', '18'), ('11', '15'), ('12', '24'), ('13', '14'), ('15', '18'), ('16', '26'),
#            ('17', '27'), ('18', '20'), ('19', '23'), ('20', '27'), ('23', '26'), ('24', '29'), ('25', '33'),
#            ('26', '30'), ('27', '32'), ('28', '32'), ('30', '32'), ('31', '34'), ('32', '33'), ('33', '34')],
#     vertex_functions=[None, None, sympy.Nand, None, None, None, None, None, None, None, None, None, None, None,
#                       sympy.Or, sympy.Nand,
#                       sympy.Nand, sympy.Nand, sympy.Nand, None, sympy.Xor, None, sympy.And, sympy.Nand,
#                       sympy.Xor, None, sympy.And, sympy.Nand, sympy.And, sympy.Xor, sympy.Or, None, sympy.Or,
#                       sympy.And, sympy.And])

# G = graphs.Network.generate_random(15, indegree_bounds=[1, 5])


# n_attractors = attractors.find_num_attractors_onestage(G=G, max_len=13, max_num=11, verbose=True, sample_mip_start=True,
#                                                        simplify_general_boolean=True)

# G = graphs.Network.generate_random(50, indegree_bounds=[1, 5])
# # attractors.find_num_steady_states(G, verbose=True, simplify_general_boolean=True)

# attractors.find_attractors_onestage_enumeration(G, max_len=10, verbose=True, simplify_general_boolean=True,
#                                                 key_slice_size=10)

# print attractors.find_bitchange_for_new_attractor(G=G, max_len=10, verbose=True, key_slice_size=5, use_dubrova=False)

# print attractors.find_model_bitchange_probability_for_different_attractors(G, n_iter=20)

#
# G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
#                    vertex_functions=[logic.SymmetricThresholdFunction(signs=[-1], threshold=1)])
# print attractors.find_num_attractors_onestage(G, max_len=2, max_num=3, verbose=True, sampling_bounds=(10, 10),
#                                               use_sampling_for_mip_start=False, simplify_general_boolean=False,
#                                               key_slice_size=5)
#

# import random
# from graphs import FunctionTypeRestriction, Network
# for _ in range(20):
#     n = random.randint(1, 20)
#     for restriction in [FunctionTypeRestriction.NONE, FunctionTypeRestriction.SYMMETRIC_THRESHOLD]:
#         # TODO: implement and test for FunctionTypeRestriction.SIMPLE_GATES
#         if random.choice([False, True]):
#             indegree_bounds = [5, 9]
#             indegree_geometric_p = None
#         else:
#             indegree_bounds = None
#             indegree_geometric_p = 0.3
#         G = Network.generate_random(n_vertices=n, indegree_bounds=indegree_bounds,
#                                     function_type_restriction=restriction,
#                                     indegree_geometric_p=indegree_geometric_p)
#         G.export_to_cnet("./temp_test_network.cnet")
#         G_tag = Network.parse_cnet("./temp_test_network.cnet")
#         assert G == G_tag