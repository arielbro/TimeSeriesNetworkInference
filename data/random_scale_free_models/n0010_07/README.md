# n0010_07

A random scale-free model: topology from networkx's directed scale-free generator (seed 10007), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 10 |
| edges | 14 |
| input nodes (in-degree 0) | 4 |
| mean in-degree | 1.40 |
| max in-degree | 7 |
| max out-degree | 4 |
| attractors (full enumeration) | 19 |
| attractor periods | 1, 2 |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**7 rows for this model's largest node.

![topology and degree distribution](network.png)
