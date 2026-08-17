# n0050_02

A random scale-free model: topology from networkx's directed scale-free generator (seed 50002), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 50 |
| edges | 79 |
| input nodes (in-degree 0) | 32 |
| mean in-degree | 1.58 |
| max in-degree | 16 |
| max out-degree | 8 |
| attractors | not enumerated (2**50 states) |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**16 rows for this model's largest node.

![topology and degree distribution](network.png)
