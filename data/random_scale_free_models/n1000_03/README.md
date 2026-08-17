# n1000_03

A random scale-free model: topology from networkx's directed scale-free generator (seed 1000003), with a uniformly random symmetric threshold function per non-input node - random signs, and a threshold drawn uniformly from [1, in-degree].

| property | value |
|---|---|
| nodes | 1000 |
| edges | 1803 |
| input nodes (in-degree 0) | 733 |
| mean in-degree | 1.80 |
| max in-degree | 291 |
| max out-degree | 54 |
| attractors | not enumerated (2**1000 states) |

Stored as `network.json` (signs + threshold per node) rather than truth tables, which would need 2**291 rows for this model's largest node.

![topology and degree distribution](network.png)
