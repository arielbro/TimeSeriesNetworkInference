# n0250_05

A random scale-free model: topology from networkx's directed scale-free generator (seed 250005), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 250 |
| edges | 385 |
| input nodes (in-degree 0) | 187 |
| mean in-degree | 1.54 |
| max in-degree | 91 |
| max out-degree | 14 |
| attractors | not enumerated (2**250 states) |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**91 rows for this model's largest node.

![topology and degree distribution](network.png)
