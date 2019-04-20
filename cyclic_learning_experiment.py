"""
A test for the accuracy of model learning under steady-state and attractor assumptions.
Done on specific models, IL-6 and MAPK, for which there are specific changes to the model listed in the
corresponding papers for removing cycles.
"""

import graphs
import logic
import sympy
import itertools

IL6_cyclic = graphs.Network.parse_boolean_tables("cellcollective_models/IL-6 Signalling")
# I found a discrepancy between the description in the paper and the cellcollective version in one vertex.
# This is the correction.
vertex_names = [v.name for v in IL6_cyclic.vertices]
edges = [(u.name, v.name) if v.name != "gp130s" else (u.name, "gp130m") for (u, v) in IL6_cyclic.edges]
gp130_func = IL6_cyclic.get_vertex("gp130s")
# Intentionally leave the new input node with its function, it will be deleted with a warning in model construction.
functions = [v.function if v.name != "gp130m" else gp130_func for v in IL6_cyclic.vertices]
IL6_cyclic = graphs.Network(vertex_names, edges, functions)

# Create the acyclic version - some nodes need to be made inputs, other need to have some inputs removed.
# Justification for the changes in plan_for_making_IL6_acyclic.png
vertex_names = [v.name for v in IL6_cyclic.vertices]
new_input_vertex_names = ["gp130m", "ras_gap", "jak1"]
removed_edges = [("gab1_mem_p", "pi3k"), ("shp2_a", "pi3k"),
                 ("shp2", "il6rc_p"), ("shp2_a", "il6rc_p"),
                 ("socs3", "shp2"),
                 ("socs1", "vav")]
ras, jak1, il6rc, il6rc_p, ros, sirp1a, dum = \
    sympy.symbols(["ras", "jak1", "il6rc", "il6rc_p", "ros", "sirp1a", "dum_il6rc_p_or_grb2_vav"])
new_functions = {"pi3k": logic.BooleanSymbolicFunc(input_names=["ras"], formula=ras),
                 "il6rc_p": logic.BooleanSymbolicFunc(input_names=["jak1", "il6rc"], formula=jak1 & il6rc),
                 "shp2": logic.BooleanSymbolicFunc(input_names=["il6rc_p", "jak1", "ros", "sirp1a"],
                                                   formula=il6rc_p & jak1 & ~ros & ~sirp1a),
                 "vav": logic.BooleanSymbolicFunc(input_names=["dum_il6rc_p_or_grb2_vav"], formula=dum)}
edges = [(u.name, v.name) for (u, v) in IL6_cyclic.edges if (v.name not in new_input_vertex_names) and
                                                            ((u.name, v.name) not in removed_edges)]
functions = []
for v in IL6_cyclic.vertices:
    if v.name in new_input_vertex_names:
        functions.append(None)
    elif v.name in new_functions:
        functions.append(new_functions[v.name])
    else:
        functions.append(v.function)
IL6_acyclic = graphs.Network(vertex_names, edges, functions)

# Make sure no cycles remain?
print IL6_cyclic
print "\n\n"
print IL6_acyclic

# Now we take gene expression data for IL-6
