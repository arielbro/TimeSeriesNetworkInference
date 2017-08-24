import time
import re
import itertools
from logic import PreRandomizedBooleanSymbolicFunc
from graphs import Network


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
            assert indices == list(range(1, n + 1))
            names = [re.search(r"=[ \t]*([0-9a-zA-Z_]+)", line).group(1) for
                     line in section.split("\n")[1:]]
            continue
        if ".n" in section:
            # boolean functions
            v_index_str, v_n_args_str, v_args_str = re.search(
                r"\.n\s+([0-9]+)[ \t]+([0-9]+)((?:[ \t][0-9]+)*)\n", section).groups()
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
            func = PreRandomizedBooleanSymbolicFunc(ordered_bool_outputs)
            bool_funcs.append(func)
    assert len(bool_funcs) == n
    G = Network(vertex_names=names, edges=edges, vertex_functions=bool_funcs)
    return G


if __name__ == "__main__":
    # G = parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
    #            "\\Attractors - for Ariel\\BNS_Dubrova_2011\\MAPK_large2.cnet")
    G = parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel"
               "\\Attractors - for Ariel\\BNS_Dubrova_2011\\ariel_test_network.txt")
    import stochastic
    initial_state = (0, 0, 0)
    visited_states = set()
    cur = initial_state
    while True:
        visited_states.add(cur)
        next = stochastic.next_network_state(G, cur)
        print next
        if next in visited_states:
            break
        cur = next
    pass
