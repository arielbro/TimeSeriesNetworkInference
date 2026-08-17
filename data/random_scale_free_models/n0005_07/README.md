# n0005_07

A random scale-free model: topology from networkx's directed scale-free generator (seed 5007), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 5 |
| edges | 7 |
| input nodes (in-degree 0) | 2 |
| mean in-degree | 1.40 |
| max in-degree | 4 |
| max out-degree | 2 |
| attractors (full enumeration) | 6 |
| attractor periods | 1 |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**4 rows for this model's largest node.

![topology and degree distribution](network.png)
