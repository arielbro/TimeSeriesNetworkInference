# n0010_02

A random scale-free model: topology from networkx's directed scale-free generator (seed 10002), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 10 |
| edges | 9 |
| input nodes (in-degree 0) | 6 |
| mean in-degree | 0.90 |
| max in-degree | 4 |
| max out-degree | 1 |
| attractors (full enumeration) | 70 |
| attractor periods | 1, 2, 6 |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**4 rows for this model's largest node.

![topology and degree distribution](network.png)
