# nk_k3_n0010_00

Kauffman NK network: every node has exactly 3 inputs, chosen uniformly without replacement from the other nodes, and a truth table whose 8 rows are independent coin flips with P(1) = 0.5.

| property | value |
|---|---|
| nodes | 10 |
| edges | 30 |
| K | 3 |
| bias | 0.5 |
| regime | chaotic |
| sensitivity (Derrida slope at distance 1) | 1.41 |
| expected sensitivity, 2p(1-p)K | 1.50 |
| seed | 3010000 |
| attractors (exact) | 1 |
| longest period | 6 |

Stored as `network.json`. `Network.parse_model_dir` reads it, so this folder works anywhere a graphs_dir is expected.
