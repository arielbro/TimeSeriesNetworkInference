# n0250_09

A random scale-free model: topology from networkx's directed scale-free generator (seed 250009), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 250 |
| edges | 414 |
| input nodes (in-degree 0) | 201 |
| mean in-degree | 1.66 |
| max in-degree | 95 |
| max out-degree | 17 |
| attractors | not enumerated (2**250 states) |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**95 rows for this model's largest node.

![topology and degree distribution](network.png)
