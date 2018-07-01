from collections import namedtuple
import attractors, graphs, time, random, itertools, csv

TimingResult = namedtuple("TimingResult", "graph_name algorithm n P T simplify sample slice_size time")
timings = []

graph_pool = []
# 'EGFR_man_with_inputs_all_zero.cnet'
# 'MAPK_large.cnet' 'arabidopsis.cnet.txt'
# for model_name in ['thelper_zeroed_inputs_ariel.cnet', 'tcr.cnet']:
#     G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
#                     "\\Attractors - for Ariel\\BNS_Dubrova_2011\\" + model_name)
#     G.convert_inputs_to_loops()
#     graph_pool.append((G, model_name))

for i in range(5000):
    n = random.randint(1, 10)
    G = graphs.Network.generate_random(n_vertices=n, indegree_bounds=(0, 4))
    G.convert_inputs_to_loops()
    graph_pool.append((G, "random"))

for i, (G, name) in enumerate(graph_pool):
    print "Iteration #{}".format(i)
    print "working on {}".format(name)
    print "running Dubrova"
    start = time.time()
    num_attractors, max_length = attractors.find_num_attractors_dubrova(G, "C:/Users/ariel/Downloads/Attractors - "
                                                                           "for Ariel/Attractors - for Ariel/BNS_Dubrova_2011",
                                                                        return_max_length=True)
    timings.append(TimingResult(algorithm="Dubrova", graph_name=name, n=len(G.vertices),
                                P=num_attractors, T=max_length, simplify=None, sample=None,
                                slice_size=None, time=time.time() - start))
    print "n={}, P={}, T={}".format(len(G.vertices), num_attractors, max_length)
    print "running ILPs"
    # take slice size 1 with 50% chance, otherwise uniform over 2-25
    slice_size = random.randint(0, 1)
    if not slice_size:
        slice_size = 1
    else:
        slice_size = random.randint(2, 25)
    for simplify, sample in itertools.product([False, True], repeat=2):
        try:
            start = time.time()
            sample_params = (num_attractors, max_length + 2) if sample else None
            found_attractors = attractors.find_num_attractors_onestage(G, max_len=max_length,
                                                                       max_num=num_attractors + 1, verbose=True,
                                                                       simplify_general_boolean=simplify,
                                                                       sample_mip_start_bounds=sample_params,
                                                                       key_slice_size=slice_size)
            if found_attractors != num_attractors:
                raise ValueError("Did not find the right amount of attractors in ILP")
            timings.append(TimingResult(algorithm="ILP", graph_name=name, n=len(G.vertices), P=num_attractors,
                                        T=max_length, simplify=simplify, sample=sample, slice_size=slice_size,
                                        time=time.time() - start))
        except attractors.TimeoutError:
            print "warning, timeout encountered in ILP"
            timings.append(TimingResult(algorithm="ILP", graph_name=name, n=len(G.vertices), P=num_attractors,
                                        T=max_length, simplify=simplify, sample=sample, slice_size=slice_size,
                                        time="inf"))
        except ValueError:
            print "warning, value error encountered in ILP"
            timings.append(TimingResult(algorithm="ILP", graph_name=name, n=len(G.vertices), P=num_attractors,
                                        T=max_length, simplify=simplify, sample=sample, slice_size=slice_size,
                                        time="invalid"))

    # save on each iteration, why not
    with open('temp_runtime_dict.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["graph_name", "algorithm", "n", "P", "T", "simplify", "sample", "slice_size", "time"])
        for timing_result in timings:
            writer.writerow(timing_result)


