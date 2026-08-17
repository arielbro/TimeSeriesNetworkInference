# n1000_02

A random scale-free model: topology from networkx's directed scale-free generator (seed 1000002), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 1000 |
| edges | 1680 |
| input nodes (in-degree 0) | 747 |
| mean in-degree | 1.68 |
| max in-degree | 239 |
| max out-degree | 38 |
| attractors | not enumerated (2**1000 states) |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**239 rows for this model's largest node.

![topology and degree distribution](network.png)
