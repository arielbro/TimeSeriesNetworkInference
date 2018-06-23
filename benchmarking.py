from collections import namedtuple
import attractors, graphs, time, random, csv

TimingResult = namedtuple("TimingResult", "algorithm n P T simplify sample time")
timings = []

graph_pool = []
# 'EGFR_man_with_inputs_all_zero.cnet'
for model_name in ['thelper.cnet.txt', 'arabidopsis.cnet.txt', 'MAPK_large.cnet', 'tcr.cnet']:
    G = graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
                    "\\Attractors - for Ariel\\BNS_Dubrova_2011\\" + model_name)
    G.convert_inputs_to_loops()
    graph_pool.append(G)

for i in range(30):
    n = random.randint(1, 30)
    G = graphs.Network.generate_random(n_vertices=n)
    G.convert_inputs_to_loops()
    graph_pool.append(G)

for graph in graph_pool:
    start = time.time()
    num_attractors, max_length = attractors.find_num_attractors_dubrova(G, "C:/Users/ariel/Downloads/Attractors - "
                                                                        "for Ariel/Attractors - for Ariel/BNS_Dubrova_2011",
                                                                        return_max_length=True)
    timings.append(TimingResult(algorithm="Dubrova", n=len(G.vertices),
                                P=num_attractors, T=max_length, simplify=None, sample=None, time=time.time() - start))
    start = time.time()
    attractors.find_num_attractors_onestage(G, max_len=max_length, max_num=num_attractors + 1,
                                            verbose=True, simplify_general_boolean=False)
    timings.append(TimingResult(algorithm="ILP", n=len(G.vertices), P=num_attractors,
                                T=max_length, simplify=False, sample=False, time=time.time() - start))
    start = time.time()
    attractors.find_num_attractors_onestage(G, max_len=max_length, max_num=num_attractors + 1,
                                            verbose=True, simplify_general_boolean=True)
    timings.append(TimingResult(algorithm="ILP", n=len(G.vertices), P=num_attractors,
                                T=max_length, simplify=True, sample=False, time=time.time() - start))
    start = time.time()
    attractors.find_num_attractors_onestage(G, max_len=max_length, max_num=num_attractors + 1,
                                            verbose=True, simplify_general_boolean=False,
                                            sample_mip_start_bounds=(num_attractors, max_length + 2))
    timings.append(TimingResult(algorithm="ILP", n=len(G.vertices), P=num_attractors,
                                T=max_length, simplify=False, sample=True, time=time.time() - start))
    start = time.time()
    attractors.find_num_attractors_onestage(G, max_len=max_length, max_num=num_attractors + 1,
                                            verbose=True, simplify_general_boolean=True,
                                            sample_mip_start_bounds=(num_attractors, max_length + 2))
    timings.append(TimingResult(algorithm="ILP", n=len(G.vertices), P=num_attractors,
                                T=max_length, simplify=True, sample=True, time=time.time() - start))

    # save on each iteration, in case run gets stuck.
    with open('temp_runtime_dict.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for timing_result in timings:
            writer.writerow(timing_result)


