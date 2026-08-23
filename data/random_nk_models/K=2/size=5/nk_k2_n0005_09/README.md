# nk_k2_n0005_09

Kauffman NK network: every node has exactly 2 inputs, chosen uniformly without replacement from the other nodes, and a truth table whose 4 rows are independent coin flips with P(1) = 0.5.

| property | value |
|---|---|
| nodes | 5 |
| edges | 10 |
| K | 2 |
| bias | 0.5 |
| regime | critical |
| sensitivity (Derrida slope at distance 1) | 1.02 |
| expected sensitivity, 2p(1-p)K | 1.00 |
| seed | 2005009 |
| attractors (exact) | 4 |
| longest period | 3 |

Stored as `network.json`. `Network.parse_model_dir` reads it, so this folder works anywhere a graphs_dir is expected.
