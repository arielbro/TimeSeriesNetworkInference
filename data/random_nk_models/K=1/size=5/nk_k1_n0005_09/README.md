# nk_k1_n0005_09

Kauffman NK network: every node has exactly 1 input, chosen uniformly without replacement from the other nodes, and a truth table whose 2 rows are independent coin flips with P(1) = 0.5.

| property | value |
|---|---|
| nodes | 5 |
| edges | 5 |
| K | 1 |
| bias | 0.5 |
| regime | ordered |
| sensitivity (Derrida slope at distance 1) | 0.56 |
| expected sensitivity, 2p(1-p)K | 0.50 |
| seed | 1005009 |
| attractors (exact) | 2 |
| longest period | 6 |

Stored as `network.json`. `Network.parse_model_dir` reads it, so this folder works anywhere a graphs_dir is expected.
