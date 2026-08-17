# n0010_09

A random scale-free model: topology from networkx's directed scale-free generator (seed 10009), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 10 |
| edges | 15 |
| input nodes (in-degree 0) | 6 |
| mean in-degree | 1.50 |
| max in-degree | 6 |
| max out-degree | 3 |
| attractors (full enumeration) | 64 |
| attractor periods | 1 |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**6 rows for this model's largest node.

![topology and degree distribution](network.png)
