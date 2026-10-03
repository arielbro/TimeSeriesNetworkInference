# nk_k1_n0010_01

Kauffman NK network: every node has exactly 1 input, chosen uniformly without replacement from the other nodes, and a truth table whose 2 rows are independent coin flips with P(1) = 0.5.

| property | value |
|---|---|
| nodes | 10 |
| edges | 10 |
| K | 1 |
| bias | 0.5 |
| regime | ordered |
| sensitivity (Derrida slope at distance 1) | 0.75 |
| expected sensitivity, 2p(1-p)K | 0.50 |
| seed | 1010001 |
| attractors (exact) | 1 |
| longest period | 1 |

Stored as `network.json`. `Network.parse_model_dir` reads it, so this folder works anywhere a graphs_dir is expected.
