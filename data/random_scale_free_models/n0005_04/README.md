# n0005_04

A random scale-free model: topology from networkx's directed scale-free generator (seed 5004), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 5 |
| edges | 6 |
| input nodes (in-degree 0) | 2 |
| mean in-degree | 1.20 |
| max in-degree | 3 |
| max out-degree | 2 |
| attractors (full enumeration) | 7 |
| attractor periods | 1, 3 |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**3 rows for this model's largest node.

![topology and degree distribution](network.png)
