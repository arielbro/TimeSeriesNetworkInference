# n0005_05

A random scale-free model: topology from networkx's directed scale-free generator (seed 5005), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 5 |
| edges | 5 |
| input nodes (in-degree 0) | 2 |
| mean in-degree | 1.00 |
| max in-degree | 3 |
| max out-degree | 2 |
| attractors (full enumeration) | 4 |
| attractor periods | 5 |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**3 rows for this model's largest node.

![topology and degree distribution](network.png)
