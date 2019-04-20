"""
A script for the specific conversion (and a small model bug fix) for models IL-1, IL-6, EFGR and MAPK,
creating a cyclic and acyclic version that I could later test model building on.
"""

import graphs
import logic
import sympy
import networkx


def to_acyclic(G, cycle_edges):
    """
    Creates and returns an acyclic version of G, given a set of edges that, when removed, remove all cycles.
    The changes to Boolean functions on edge removals are as defined in graphs.Network.remove_edge_dependency.
    :param G:
    :param cycle_edges:
    :return:
    """
    G_copy = G.copy()
    for edge in cycle_edges:
        G_copy.remove_edge_dependency(edge)
    return G_copy


networks = dict()

"""
IL-1
"""
print "IL-1"
IL1_cyclic = graphs.Network.parse_boolean_tables("cellcollective_models/IL-1 Signalling")

# Create the acyclic version - most nodes here should be labeled "late" and have all their outgoing edges removed,
# but there are dummies depending on them.
# See justification for the modifications in figure 1 of
# "Large-scale network models of IL-1 and IL-6 signalling and their hepatocellular specification"
late_nodes_list = [IL1_cyclic.get_vertex("inos")]
removed_edges = [(u, v) for (u, v) in IL1_cyclic.edges if u in late_nodes_list]
removed_edges.append((IL1_cyclic.get_vertex("il1ra"), IL1_cyclic.get_vertex("il1r1")))
removed_edges.append((IL1_cyclic.get_vertex("il1ra"), IL1_cyclic.get_vertex("il1r2")))
removed_edges.append((IL1_cyclic.get_vertex("hsp27_ps"), IL1_cyclic.get_vertex("traf6_ub")))
removed_edges.append((IL1_cyclic.get_vertex("cyt_p38"), IL1_cyclic.get_vertex("tak1_tab")))
removed_edges.append((IL1_cyclic.get_vertex("mkp1"), IL1_cyclic.get_vertex("erk12")))
removed_edges.append((IL1_cyclic.get_vertex("mkp1"), IL1_cyclic.get_vertex("jnk")))
removed_edges.append((IL1_cyclic.get_vertex("mkp1"), IL1_cyclic.get_vertex("nuc_p38")))
removed_edges.append((IL1_cyclic.get_vertex("tpl2_degr"), IL1_cyclic.get_vertex("tpl2")))
removed_edges.append((IL1_cyclic.get_vertex("nfkb"), IL1_cyclic.get_vertex("ikba")))
removed_edges.append((IL1_cyclic.get_vertex("nfkb"), IL1_cyclic.get_vertex("ros")))

IL1_acyclic = to_acyclic(IL1_cyclic, removed_edges)
networks["IL1_acyclic"] = IL1_acyclic

"""
IL-6
"""
print "IL-6"
IL6_cyclic = graphs.Network.parse_boolean_tables("cellcollective_models/IL-6 Signalling")
# I found a discrepancy between the description in the paper and the cellcollective version in one vertex.
# This is the correction.
vertex_names = [v.name for v in IL6_cyclic.vertices]
edges = [(u.name, v.name) if v.name != "gp130s" else (u.name, "gp130m") for (u, v) in IL6_cyclic.edges]
correction_dict = {"gp130s": None, "gp130m": IL6_cyclic.get_vertex("gp130s").function}
functions = [v.function if v.name not in correction_dict else correction_dict[v.name] for v in IL6_cyclic.vertices]
IL6_cyclic = graphs.Network(vertex_names, edges, functions)
networks["IL6_cyclic"] = IL6_cyclic

# Create the acyclic version - some nodes need to be made inputs, other need to have some inputs removed.
# See justification for the modifications in figure 2 of
# "Large-scale network models of IL-1 and IL-6 signalling and their hepatocellular specification"
new_input_vertices = [IL6_cyclic.get_vertex("gp130m"), IL6_cyclic.get_vertex("ras_gap"), IL6_cyclic.get_vertex("jak1")]
removed_edges = [(u, v) for (u, v) in IL6_cyclic.edges if v in new_input_vertices]
removed_edges.append((IL6_cyclic.get_vertex("gab1_mem_p"), IL6_cyclic.get_vertex("pi3k")))
removed_edges.append((IL6_cyclic.get_vertex("shp2_a"), IL6_cyclic.get_vertex("pi3k")))
removed_edges.append((IL6_cyclic.get_vertex("shp2"), IL6_cyclic.get_vertex("il6rc_p")))
removed_edges.append((IL6_cyclic.get_vertex("shp2_a"), IL6_cyclic.get_vertex("il6rc_p")))
removed_edges.append((IL6_cyclic.get_vertex("socs3"), IL6_cyclic.get_vertex("shp2")))
removed_edges.append((IL6_cyclic.get_vertex("socs1"), IL6_cyclic.get_vertex("vav")))

IL6_acyclic = to_acyclic(IL6_cyclic, removed_edges)
networks["IL6_acyclic"] = IL6_acyclic

"""
EGFR
"""
print "EGFR"
egfr_cyclic = graphs.Network.parse_boolean_tables("cellcollective_models/EGFR & ErbB Signaling")
networks["EGFR_cyclic"] = egfr_cyclic

# In this model most of the late reactions (edges we need to remove) originate
# in nodes for which all outgoing reactions are late. Where they aren't (sometimes when they are too),
# there's usually a dummy node whose only outgoing reactions are late. We can then mark nodes as "late" and remove
# all of their successors for most of the changes.
# See justification in figure 1 of
# "The Logic of EGFR/ErbB Signaling: Theoretical Properties and Analysis of High-Throughput Data"
late_nodes_list = [egfr_cyclic.get_vertex("aktd"), egfr_cyclic.get_vertex("endocyt_degrad"),
                   egfr_cyclic.get_vertex("shp1d"), egfr_cyclic.get_vertex("p90rskerk12d")]
print "late nodes done"
removed_edges = []
for u, v in egfr_cyclic.edges:
    print u, v
    print len(u.predecessors())
    if u in late_nodes_list:
        print "u in there"
        removed_edges.append((u, v))
    print "check done"
removed_edges = [(u, v) for (u, v) in egfr_cyclic.edges if u in late_nodes_list]
print "removed edge"
removed_edges.append((egfr_cyclic.get_vertex("ptend"), egfr_cyclic.get_vertex("pip3")))
print "removed edge"
removed_edges.append((egfr_cyclic.get_vertex("ship2d"), egfr_cyclic.get_vertex("pip3")))
print "removed edge"
removed_edges.append((egfr_cyclic.get_vertex("ptend"), egfr_cyclic.get_vertex("pi34p2")))
print "removed edge"
removed_edges.append((egfr_cyclic.get_vertex("pip3"), egfr_cyclic.get_vertex("gab1")))
print "removed edge"

egfr_acyclic = to_acyclic(egfr_cyclic, removed_edges)
print "acyclic now"
networks["EGFR_acyclic"] = egfr_acyclic

print "running tests"

# Check if transformation worked
for name, G in networks.items():
    Gx = networkx.DiGraph()
    for v in G.vertices:
        Gx.add_node(v.name)
    for u, v in G.edges:
        Gx.add_edge(u.name, v.name)
    if "_cyclic" in name:
        try:
            networkx.algorithms.cycles.find_cycle(Gx)
        except networkx.NetworkXNoCycle:
            raise ValueError("graph {} has no cycles".format(name))
    elif "_acyclic" in name:
        try:
            cycles = networkx.algorithms.cycles.find_cycle(Gx)
            raise ValueError("graph {} has cycles: {}".format(name, cycles))
        except networkx.NetworkXNoCycle:
            pass
    else:
        raise ValueError("unfamiliar graph name pattern")

print "exporting"

for name, G in networks.items():
    print "exporting {}".format(name)
    graphs.Network.export_to_boolean_tables(G, "cyclic_acyclic_models", name)
