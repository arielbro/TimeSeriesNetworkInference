# n0010_06

A random scale-free model: topology from networkx's directed scale-free generator (seed 10006), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 10 |
| edges | 14 |
| input nodes (in-degree 0) | 4 |
| mean in-degree | 1.40 |
| max in-degree | 4 |
| max out-degree | 4 |
| attractors (full enumeration) | 17 |
| attractor periods | 1, 4 |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**4 rows for this model's largest node.

![topology and degree distribution](network.png)
