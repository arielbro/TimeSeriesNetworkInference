# n0050_09

A random scale-free model: topology from networkx's directed scale-free generator (seed 50009), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 50 |
| edges | 82 |
| input nodes (in-degree 0) | 35 |
| mean in-degree | 1.64 |
| max in-degree | 26 |
| max out-degree | 7 |
| attractors | not enumerated (2**50 states) |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**26 rows for this model's largest node.

![topology and degree distribution](network.png)
