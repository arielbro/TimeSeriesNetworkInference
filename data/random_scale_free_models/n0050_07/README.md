# n0050_07

A random scale-free model: topology from networkx's directed scale-free generator (seed 50007), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 50 |
| edges | 68 |
| input nodes (in-degree 0) | 36 |
| mean in-degree | 1.36 |
| max in-degree | 32 |
| max out-degree | 5 |
| attractors | not enumerated (2**50 states) |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**32 rows for this model's largest node.

![topology and degree distribution](network.png)
