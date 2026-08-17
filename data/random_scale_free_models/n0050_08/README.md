# n0050_08

A random scale-free model: topology from networkx's directed scale-free generator (seed 50008), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 50 |
| edges | 75 |
| input nodes (in-degree 0) | 34 |
| mean in-degree | 1.50 |
| max in-degree | 16 |
| max out-degree | 9 |
| attractors | not enumerated (2**50 states) |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**16 rows for this model's largest node.

![topology and degree distribution](network.png)
