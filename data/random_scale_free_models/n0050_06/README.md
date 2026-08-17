# n0050_06

A random scale-free model: topology from networkx's directed scale-free generator (seed 50006), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 50 |
| edges | 60 |
| input nodes (in-degree 0) | 37 |
| mean in-degree | 1.20 |
| max in-degree | 36 |
| max out-degree | 5 |
| attractors | not enumerated (2**50 states) |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**36 rows for this model's largest node.

![topology and degree distribution](network.png)
