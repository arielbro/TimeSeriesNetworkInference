# n0250_08

A random scale-free model: topology from networkx's directed scale-free generator (seed 250008), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 250 |
| edges | 395 |
| input nodes (in-degree 0) | 188 |
| mean in-degree | 1.58 |
| max in-degree | 130 |
| max out-degree | 12 |
| attractors | not enumerated (2**250 states) |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**130 rows for this model's largest node.

![topology and degree distribution](network.png)
