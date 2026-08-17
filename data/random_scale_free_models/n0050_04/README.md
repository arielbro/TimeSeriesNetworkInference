# n0050_04

A random scale-free model: topology from networkx's directed scale-free generator (seed 50004), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 50 |
| edges | 71 |
| input nodes (in-degree 0) | 41 |
| mean in-degree | 1.42 |
| max in-degree | 34 |
| max out-degree | 5 |
| attractors | not enumerated (2**50 states) |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**34 rows for this model's largest node.

![topology and degree distribution](network.png)
