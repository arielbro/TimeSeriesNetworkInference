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

for i in range(10000):
    n = random.randint(3, 15)
    G = graphs.Network.generate_random(n_vertices=n, indegree_bounds=(0, 4))
    G.convert_inputs_to_loops()
    name = "random"
    print "Iteration #{}".format(i)
    print "working on {}".format(name)
    print "running Dubrova"
    start = time.time()
    attractors = attractors.find_attractors_dubrova(G, "C:/Users/ariel/Downloads/Attractors - "
                                                                           "for Ariel/Attractors - for Ariel/BNS_Dubrova_2011",
                                                                        return_max_length=True)
    num_attractors, max_length = len(attractors), max(len(att) for att in attractors)
    timings.append(TimingResult(algorithm="Dubrova", graph_name=name, n=len(G.vertices),
                                P=num_attractors, T=max_length, simplify=None, sample=None,
                                slice_size=None, time=time.time() - start))
    print "n={}, P={}, T={}".format(len(G.vertices), num_attractors, max_length)
    print "running ILPs"
    slice_size = random.randint(1, 15)
    for simplify, sample in itertools.product([False, True], [False, "mip", "unique"]):
        try:
            start = time.time()
            sample_params = (num_attractors, max_length + 2) if sample else None
            mip_start = True if sample=="mip" else False
            found_attractors = attractors.find_num_attractors_onestage(G, max_len=max_length,
                                                                       max_num=num_attractors + 1, verbose=True,
                                                                       simplify_general_boolean=simplify,
                                                                       sampling_bounds=sample_params,
                                                                       use_sampling_for_mip_start=mip_start,
                                                                       key_slice_size=slice_size)
            if found_attractors != num_attractors:
                print "Warning - did not find the right amount of attractors in ILP"
                timings.append(TimingResult(algorithm="ILP", graph_name=name, n=len(G.vertices), P=num_attractors,
                                            T=max_length, simplify=simplify, sample=sample, slice_size=slice_size,
                                            time='invalid'))
            else:
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
    # enumeration
    for simplify in [False, True]:
        try:
            start = time.time()
            found_attractors = attractors.find_attractors_onestage_enumeration(G, max_len=max_length, verbose=True,
                                                                               simplify_general_boolean=simplify,
                                                                               key_slice_size=slice_size)
            if len(found_attractors) != num_attractors:
                print "Warning - did not find the right amount of attractors in ILP"
                timings.append(TimingResult(algorithm="ILP_enum", graph_name=name, n=len(G.vertices), P=num_attractors,
                                            T=max_length, simplify=simplify, sample=None, slice_size=slice_size,
                                            time='invalid'))
            else:
                timings.append(TimingResult(algorithm="ILP_enum", graph_name=name, n=len(G.vertices), P=num_attractors,
                                            T=max_length, simplify=simplify, sample=None, slice_size=slice_size,
                                            time=time.time() - start))

        except attractors.TimeoutError:
            print "warning, timeout encountered in ILP"
            timings.append(TimingResult(algorithm="ILP_enum", graph_name=name, n=len(G.vertices), P=num_attractors,
                                        T=max_length, simplify=simplify, sample=None, slice_size=slice_size,
                                        time="inf"))
        except ValueError:
            print "warning, value error encountered in ILP"
            timings.append(TimingResult(algorithm="ILP_enum", graph_name=name, n=len(G.vertices), P=num_attractors,
                                        T=max_length, simplify=simplify, sample=None, slice_size=slice_size,
                                        time="invalid"))

    # save on each iteration, why not
    with open('temp_runtime_dict.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["graph_name", "algorithm", "n", "P", "T", "simplify", "sample", "slice_size", "time"])
        for timing_result in timings:
            writer.writerow(timing_result)


