# nk_k2_n0050_00

Kauffman NK network: every node has exactly 2 inputs, chosen uniformly without replacement from the other nodes, and a truth table whose 4 rows are independent coin flips with P(1) = 0.5.

| property | value |
|---|---|
| nodes | 50 |
| edges | 100 |
| K | 2 |
| bias | 0.5 |
| regime | critical |
| sensitivity (Derrida slope at distance 1) | 1.00 |
| expected sensitivity, 2p(1-p)K | 1.00 |
| seed | 2050000 |
| attractors | not enumerated (2**50 states) |

Stored as `network.json`. `Network.parse_model_dir` reads it, so this folder works anywhere a graphs_dir is expected.
