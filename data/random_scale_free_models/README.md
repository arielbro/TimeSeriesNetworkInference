# Random scale-free models

A `graphs_dir` of 50 completely random models: 10 at each of the sizes 5, 10, 50, 250, 1000. Topology from networkx's directed scale-free generator; every non-input node gets a uniformly random symmetric threshold function (random signs, threshold uniform in [1, in-degree]).

Stored as `network.json` (each function as signs + threshold), not as cellcollective truth tables: hub nodes here have in-degrees in the dozens, and a truth table needs 2**in-degree rows. `Network.parse_model_dir` / `Network.model_dir_size` read either format, so this folder works anywhere a graphs_dir is expected.

Regenerate with `python synthetic_data_generation/make_random_scale_free_models.py`.

| model | nodes | edges | input nodes | mean in-degree | max in-degree | attractors |
|---|---|---|---|---|---|---|
| [n0005_00](n0005_00/README.md) | 5 | 5 | 2 | 1.00 | 2 | 5 |
| [n0005_01](n0005_01/README.md) | 5 | 7 | 2 | 1.40 | 4 | 5 |
| [n0005_02](n0005_02/README.md) | 5 | 7 | 2 | 1.40 | 3 | 4 |
| [n0005_03](n0005_03/README.md) | 5 | 5 | 2 | 1.00 | 2 | 7 |
| [n0005_04](n0005_04/README.md) | 5 | 6 | 2 | 1.20 | 3 | 7 |
| [n0005_05](n0005_05/README.md) | 5 | 5 | 2 | 1.00 | 3 | 4 |
| [n0005_06](n0005_06/README.md) | 5 | 5 | 2 | 1.00 | 2 | 5 |
| [n0005_07](n0005_07/README.md) | 5 | 7 | 2 | 1.40 | 4 | 6 |
| [n0005_08](n0005_08/README.md) | 5 | 5 | 2 | 1.00 | 2 | 7 |
| [n0005_09](n0005_09/README.md) | 5 | 5 | 1 | 1.00 | 2 | 5 |
| [n0010_00](n0010_00/README.md) | 10 | 13 | 5 | 1.30 | 5 | 32 |
| [n0010_01](n0010_01/README.md) | 10 | 12 | 7 | 1.20 | 6 | 134 |
| [n0010_02](n0010_02/README.md) | 10 | 9 | 6 | 0.90 | 4 | 70 |
| [n0010_03](n0010_03/README.md) | 10 | 12 | 6 | 1.20 | 5 | 66 |
| [n0010_04](n0010_04/README.md) | 10 | 11 | 6 | 1.10 | 6 | 66 |
| [n0010_05](n0010_05/README.md) | 10 | 11 | 6 | 1.10 | 8 | 85 |
| [n0010_06](n0010_06/README.md) | 10 | 14 | 4 | 1.40 | 4 | 17 |
| [n0010_07](n0010_07/README.md) | 10 | 14 | 4 | 1.40 | 7 | 19 |
| [n0010_08](n0010_08/README.md) | 10 | 12 | 7 | 1.20 | 6 | 128 |
| [n0010_09](n0010_09/README.md) | 10 | 15 | 6 | 1.50 | 6 | 64 |
| [n0050_00](n0050_00/README.md) | 50 | 71 | 39 | 1.42 | 24 | not enumerated |
| [n0050_01](n0050_01/README.md) | 50 | 73 | 34 | 1.46 | 27 | not enumerated |
| [n0050_02](n0050_02/README.md) | 50 | 79 | 32 | 1.58 | 16 | not enumerated |
| [n0050_03](n0050_03/README.md) | 50 | 70 | 37 | 1.40 | 19 | not enumerated |
| [n0050_04](n0050_04/README.md) | 50 | 71 | 41 | 1.42 | 34 | not enumerated |
| [n0050_05](n0050_05/README.md) | 50 | 71 | 32 | 1.42 | 23 | not enumerated |
| [n0050_06](n0050_06/README.md) | 50 | 60 | 37 | 1.20 | 36 | not enumerated |
| [n0050_07](n0050_07/README.md) | 50 | 68 | 36 | 1.36 | 32 | not enumerated |
| [n0050_08](n0050_08/README.md) | 50 | 75 | 34 | 1.50 | 16 | not enumerated |
| [n0050_09](n0050_09/README.md) | 50 | 82 | 35 | 1.64 | 26 | not enumerated |
| [n0250_00](n0250_00/README.md) | 250 | 420 | 188 | 1.68 | 103 | not enumerated |
| [n0250_01](n0250_01/README.md) | 250 | 376 | 189 | 1.50 | 155 | not enumerated |
| [n0250_02](n0250_02/README.md) | 250 | 426 | 192 | 1.70 | 73 | not enumerated |
| [n0250_03](n0250_03/README.md) | 250 | 379 | 202 | 1.52 | 118 | not enumerated |
| [n0250_04](n0250_04/README.md) | 250 | 399 | 185 | 1.60 | 108 | not enumerated |
| [n0250_05](n0250_05/README.md) | 250 | 385 | 187 | 1.54 | 91 | not enumerated |
| [n0250_06](n0250_06/README.md) | 250 | 424 | 197 | 1.70 | 88 | not enumerated |
| [n0250_07](n0250_07/README.md) | 250 | 403 | 188 | 1.61 | 113 | not enumerated |
| [n0250_08](n0250_08/README.md) | 250 | 395 | 188 | 1.58 | 130 | not enumerated |
| [n0250_09](n0250_09/README.md) | 250 | 414 | 201 | 1.66 | 95 | not enumerated |
| [n1000_00](n1000_00/README.md) | 1000 | 1712 | 760 | 1.71 | 308 | not enumerated |
| [n1000_01](n1000_01/README.md) | 1000 | 1778 | 758 | 1.78 | 336 | not enumerated |
| [n1000_02](n1000_02/README.md) | 1000 | 1680 | 747 | 1.68 | 239 | not enumerated |
| [n1000_03](n1000_03/README.md) | 1000 | 1803 | 733 | 1.80 | 291 | not enumerated |
| [n1000_04](n1000_04/README.md) | 1000 | 1736 | 752 | 1.74 | 380 | not enumerated |
| [n1000_05](n1000_05/README.md) | 1000 | 1724 | 772 | 1.72 | 299 | not enumerated |
| [n1000_06](n1000_06/README.md) | 1000 | 1700 | 742 | 1.70 | 311 | not enumerated |
| [n1000_07](n1000_07/README.md) | 1000 | 1782 | 763 | 1.78 | 230 | not enumerated |
| [n1000_08](n1000_08/README.md) | 1000 | 1641 | 756 | 1.64 | 487 | not enumerated |
| [n1000_09](n1000_09/README.md) | 1000 | 1647 | 765 | 1.65 | 506 | not enumerated |
