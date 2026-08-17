# n0005_09

A random scale-free model: topology from networkx's directed scale-free generator (seed 5009), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 5 |
| edges | 5 |
| input nodes (in-degree 0) | 1 |
| mean in-degree | 1.00 |
| max in-degree | 2 |
| max out-degree | 2 |
| attractors (full enumeration) | 5 |
| attractor periods | 1, 3 |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**2 rows for this model's largest node.

![topology and degree distribution](network.png)
