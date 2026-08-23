# Random NK models

A `graphs_dir` of 120 Kauffman NK networks: 10 at each combination of K in 1, 2, 3 and size in 5, 10, 50, 250. Every node has exactly K inputs drawn uniformly without replacement from the other nodes, and a truth table whose 2**K rows are independent coin flips with P(1) = 0.5.

The regimes follow the annealed approximation, critical at 2*p*(1-p)*K = 1: at p = 0.5 that is K/2, so K = 1, 2, 3 are ordered, critical and chaotic. The sensitivity column is the measured Derrida slope at distance 1 (mean nodes differing one step after a one-bit perturbation), which should sit near K/2.

Unlike `data/random_scale_free_models`, the functions here are general Boolean rather than symmetric threshold, so these models are outside the class the symmetric / symmetric_topology methods search.

Grouped `K=<k>/size=<n>`, one directory level per parameter, so both reach the analysis notebook as swept numeric parameters.

Regenerate with `python synthetic_data_generation/make_random_nk_models.py`.

## K = 1 (ordered)

| model | nodes | edges | sensitivity | attractors | longest period |
|---|---|---|---|---|---|
| [nk_k1_n0005_00](K=1/size=5/nk_k1_n0005_00/README.md) | 5 | 5 | 0.56 | 3 | 2 |
| [nk_k1_n0005_01](K=1/size=5/nk_k1_n0005_01/README.md) | 5 | 5 | 0.60 | 1 | 1 |
| [nk_k1_n0005_02](K=1/size=5/nk_k1_n0005_02/README.md) | 5 | 5 | 0.00 | 1 | 1 |
| [nk_k1_n0005_03](K=1/size=5/nk_k1_n0005_03/README.md) | 5 | 5 | 0.69 | 1 | 1 |
| [nk_k1_n0005_04](K=1/size=5/nk_k1_n0005_04/README.md) | 5 | 5 | 0.41 | 1 | 1 |
| [nk_k1_n0005_05](K=1/size=5/nk_k1_n0005_05/README.md) | 5 | 5 | 0.41 | 1 | 1 |
| [nk_k1_n0005_06](K=1/size=5/nk_k1_n0005_06/README.md) | 5 | 5 | 0.84 | 1 | 1 |
| [nk_k1_n0005_07](K=1/size=5/nk_k1_n0005_07/README.md) | 5 | 5 | 0.68 | 1 | 4 |
| [nk_k1_n0005_08](K=1/size=5/nk_k1_n0005_08/README.md) | 5 | 5 | 0.55 | 1 | 1 |
| [nk_k1_n0005_09](K=1/size=5/nk_k1_n0005_09/README.md) | 5 | 5 | 0.56 | 2 | 6 |
| [nk_k1_n0010_00](K=1/size=10/nk_k1_n0010_00/README.md) | 10 | 10 | 0.72 | 1 | 1 |
| [nk_k1_n0010_01](K=1/size=10/nk_k1_n0010_01/README.md) | 10 | 10 | 0.75 | 1 | 1 |
| [nk_k1_n0010_02](K=1/size=10/nk_k1_n0010_02/README.md) | 10 | 10 | 0.39 | 1 | 1 |
| [nk_k1_n0010_03](K=1/size=10/nk_k1_n0010_03/README.md) | 10 | 10 | 0.41 | 1 | 1 |
| [nk_k1_n0010_04](K=1/size=10/nk_k1_n0010_04/README.md) | 10 | 10 | 0.64 | 1 | 1 |
| [nk_k1_n0010_05](K=1/size=10/nk_k1_n0010_05/README.md) | 10 | 10 | 0.28 | 1 | 1 |
| [nk_k1_n0010_06](K=1/size=10/nk_k1_n0010_06/README.md) | 10 | 10 | 0.64 | 1 | 1 |
| [nk_k1_n0010_07](K=1/size=10/nk_k1_n0010_07/README.md) | 10 | 10 | 0.44 | 1 | 1 |
| [nk_k1_n0010_08](K=1/size=10/nk_k1_n0010_08/README.md) | 10 | 10 | 0.56 | 1 | 1 |
| [nk_k1_n0010_09](K=1/size=10/nk_k1_n0010_09/README.md) | 10 | 10 | 0.31 | 1 | 1 |
| [nk_k1_n0050_00](K=1/size=50/nk_k1_n0050_00/README.md) | 50 | 50 | 0.40 | not enumerated | - |
| [nk_k1_n0050_01](K=1/size=50/nk_k1_n0050_01/README.md) | 50 | 50 | 0.54 | not enumerated | - |
| [nk_k1_n0050_02](K=1/size=50/nk_k1_n0050_02/README.md) | 50 | 50 | 0.47 | not enumerated | - |
| [nk_k1_n0050_03](K=1/size=50/nk_k1_n0050_03/README.md) | 50 | 50 | 0.43 | not enumerated | - |
| [nk_k1_n0050_04](K=1/size=50/nk_k1_n0050_04/README.md) | 50 | 50 | 0.55 | not enumerated | - |
| [nk_k1_n0050_05](K=1/size=50/nk_k1_n0050_05/README.md) | 50 | 50 | 0.57 | not enumerated | - |
| [nk_k1_n0050_06](K=1/size=50/nk_k1_n0050_06/README.md) | 50 | 50 | 0.60 | not enumerated | - |
| [nk_k1_n0050_07](K=1/size=50/nk_k1_n0050_07/README.md) | 50 | 50 | 0.55 | not enumerated | - |
| [nk_k1_n0050_08](K=1/size=50/nk_k1_n0050_08/README.md) | 50 | 50 | 0.50 | not enumerated | - |
| [nk_k1_n0050_09](K=1/size=50/nk_k1_n0050_09/README.md) | 50 | 50 | 0.69 | not enumerated | - |
| [nk_k1_n0250_00](K=1/size=250/nk_k1_n0250_00/README.md) | 250 | 250 | 0.51 | not enumerated | - |
| [nk_k1_n0250_01](K=1/size=250/nk_k1_n0250_01/README.md) | 250 | 250 | 0.46 | not enumerated | - |
| [nk_k1_n0250_02](K=1/size=250/nk_k1_n0250_02/README.md) | 250 | 250 | 0.50 | not enumerated | - |
| [nk_k1_n0250_03](K=1/size=250/nk_k1_n0250_03/README.md) | 250 | 250 | 0.49 | not enumerated | - |
| [nk_k1_n0250_04](K=1/size=250/nk_k1_n0250_04/README.md) | 250 | 250 | 0.57 | not enumerated | - |
| [nk_k1_n0250_05](K=1/size=250/nk_k1_n0250_05/README.md) | 250 | 250 | 0.59 | not enumerated | - |
| [nk_k1_n0250_06](K=1/size=250/nk_k1_n0250_06/README.md) | 250 | 250 | 0.47 | not enumerated | - |
| [nk_k1_n0250_07](K=1/size=250/nk_k1_n0250_07/README.md) | 250 | 250 | 0.62 | not enumerated | - |
| [nk_k1_n0250_08](K=1/size=250/nk_k1_n0250_08/README.md) | 250 | 250 | 0.51 | not enumerated | - |
| [nk_k1_n0250_09](K=1/size=250/nk_k1_n0250_09/README.md) | 250 | 250 | 0.54 | not enumerated | - |

## K = 2 (critical)

| model | nodes | edges | sensitivity | attractors | longest period |
|---|---|---|---|---|---|
| [nk_k2_n0005_00](K=2/size=5/nk_k2_n0005_00/README.md) | 5 | 10 | 1.24 | 1 | 4 |
| [nk_k2_n0005_01](K=2/size=5/nk_k2_n0005_01/README.md) | 5 | 10 | 0.94 | 2 | 2 |
| [nk_k2_n0005_02](K=2/size=5/nk_k2_n0005_02/README.md) | 5 | 10 | 1.28 | 2 | 2 |
| [nk_k2_n0005_03](K=2/size=5/nk_k2_n0005_03/README.md) | 5 | 10 | 0.63 | 1 | 1 |
| [nk_k2_n0005_04](K=2/size=5/nk_k2_n0005_04/README.md) | 5 | 10 | 0.99 | 1 | 1 |
| [nk_k2_n0005_05](K=2/size=5/nk_k2_n0005_05/README.md) | 5 | 10 | 0.99 | 1 | 1 |
| [nk_k2_n0005_06](K=2/size=5/nk_k2_n0005_06/README.md) | 5 | 10 | 0.83 | 1 | 1 |
| [nk_k2_n0005_07](K=2/size=5/nk_k2_n0005_07/README.md) | 5 | 10 | 0.77 | 3 | 2 |
| [nk_k2_n0005_08](K=2/size=5/nk_k2_n0005_08/README.md) | 5 | 10 | 1.00 | 1 | 4 |
| [nk_k2_n0005_09](K=2/size=5/nk_k2_n0005_09/README.md) | 5 | 10 | 1.02 | 4 | 3 |
| [nk_k2_n0010_00](K=2/size=10/nk_k2_n0010_00/README.md) | 10 | 20 | 1.22 | 4 | 3 |
| [nk_k2_n0010_01](K=2/size=10/nk_k2_n0010_01/README.md) | 10 | 20 | 0.94 | 3 | 2 |
| [nk_k2_n0010_02](K=2/size=10/nk_k2_n0010_02/README.md) | 10 | 20 | 0.86 | 1 | 1 |
| [nk_k2_n0010_03](K=2/size=10/nk_k2_n0010_03/README.md) | 10 | 20 | 0.84 | 1 | 4 |
| [nk_k2_n0010_04](K=2/size=10/nk_k2_n0010_04/README.md) | 10 | 20 | 1.18 | 3 | 6 |
| [nk_k2_n0010_05](K=2/size=10/nk_k2_n0010_05/README.md) | 10 | 20 | 1.11 | 1 | 1 |
| [nk_k2_n0010_06](K=2/size=10/nk_k2_n0010_06/README.md) | 10 | 20 | 0.85 | 1 | 4 |
| [nk_k2_n0010_07](K=2/size=10/nk_k2_n0010_07/README.md) | 10 | 20 | 1.10 | 3 | 2 |
| [nk_k2_n0010_08](K=2/size=10/nk_k2_n0010_08/README.md) | 10 | 20 | 0.94 | 2 | 2 |
| [nk_k2_n0010_09](K=2/size=10/nk_k2_n0010_09/README.md) | 10 | 20 | 0.61 | 1 | 1 |
| [nk_k2_n0050_00](K=2/size=50/nk_k2_n0050_00/README.md) | 50 | 100 | 1.00 | not enumerated | - |
| [nk_k2_n0050_01](K=2/size=50/nk_k2_n0050_01/README.md) | 50 | 100 | 1.21 | not enumerated | - |
| [nk_k2_n0050_02](K=2/size=50/nk_k2_n0050_02/README.md) | 50 | 100 | 1.07 | not enumerated | - |
| [nk_k2_n0050_03](K=2/size=50/nk_k2_n0050_03/README.md) | 50 | 100 | 1.10 | not enumerated | - |
| [nk_k2_n0050_04](K=2/size=50/nk_k2_n0050_04/README.md) | 50 | 100 | 0.96 | not enumerated | - |
| [nk_k2_n0050_05](K=2/size=50/nk_k2_n0050_05/README.md) | 50 | 100 | 1.06 | not enumerated | - |
| [nk_k2_n0050_06](K=2/size=50/nk_k2_n0050_06/README.md) | 50 | 100 | 0.86 | not enumerated | - |
| [nk_k2_n0050_07](K=2/size=50/nk_k2_n0050_07/README.md) | 50 | 100 | 0.91 | not enumerated | - |
| [nk_k2_n0050_08](K=2/size=50/nk_k2_n0050_08/README.md) | 50 | 100 | 1.07 | not enumerated | - |
| [nk_k2_n0050_09](K=2/size=50/nk_k2_n0050_09/README.md) | 50 | 100 | 0.90 | not enumerated | - |
| [nk_k2_n0250_00](K=2/size=250/nk_k2_n0250_00/README.md) | 250 | 500 | 0.99 | not enumerated | - |
| [nk_k2_n0250_01](K=2/size=250/nk_k2_n0250_01/README.md) | 250 | 500 | 1.09 | not enumerated | - |
| [nk_k2_n0250_02](K=2/size=250/nk_k2_n0250_02/README.md) | 250 | 500 | 1.02 | not enumerated | - |
| [nk_k2_n0250_03](K=2/size=250/nk_k2_n0250_03/README.md) | 250 | 500 | 0.95 | not enumerated | - |
| [nk_k2_n0250_04](K=2/size=250/nk_k2_n0250_04/README.md) | 250 | 500 | 0.94 | not enumerated | - |
| [nk_k2_n0250_05](K=2/size=250/nk_k2_n0250_05/README.md) | 250 | 500 | 0.99 | not enumerated | - |
| [nk_k2_n0250_06](K=2/size=250/nk_k2_n0250_06/README.md) | 250 | 500 | 0.97 | not enumerated | - |
| [nk_k2_n0250_07](K=2/size=250/nk_k2_n0250_07/README.md) | 250 | 500 | 0.95 | not enumerated | - |
| [nk_k2_n0250_08](K=2/size=250/nk_k2_n0250_08/README.md) | 250 | 500 | 1.01 | not enumerated | - |
| [nk_k2_n0250_09](K=2/size=250/nk_k2_n0250_09/README.md) | 250 | 500 | 1.08 | not enumerated | - |

## K = 3 (chaotic)

| model | nodes | edges | sensitivity | attractors | longest period |
|---|---|---|---|---|---|
| [nk_k3_n0005_00](K=3/size=5/nk_k3_n0005_00/README.md) | 5 | 15 | 1.39 | 2 | 3 |
| [nk_k3_n0005_01](K=3/size=5/nk_k3_n0005_01/README.md) | 5 | 15 | 1.26 | 1 | 1 |
| [nk_k3_n0005_02](K=3/size=5/nk_k3_n0005_02/README.md) | 5 | 15 | 1.34 | 1 | 4 |
| [nk_k3_n0005_03](K=3/size=5/nk_k3_n0005_03/README.md) | 5 | 15 | 1.35 | 2 | 1 |
| [nk_k3_n0005_04](K=3/size=5/nk_k3_n0005_04/README.md) | 5 | 15 | 1.22 | 1 | 1 |
| [nk_k3_n0005_05](K=3/size=5/nk_k3_n0005_05/README.md) | 5 | 15 | 1.37 | 1 | 7 |
| [nk_k3_n0005_06](K=3/size=5/nk_k3_n0005_06/README.md) | 5 | 15 | 1.41 | 2 | 4 |
| [nk_k3_n0005_07](K=3/size=5/nk_k3_n0005_07/README.md) | 5 | 15 | 1.57 | 2 | 3 |
| [nk_k3_n0005_08](K=3/size=5/nk_k3_n0005_08/README.md) | 5 | 15 | 1.45 | 3 | 5 |
| [nk_k3_n0005_09](K=3/size=5/nk_k3_n0005_09/README.md) | 5 | 15 | 1.41 | 1 | 5 |
| [nk_k3_n0010_00](K=3/size=10/nk_k3_n0010_00/README.md) | 10 | 30 | 1.41 | 1 | 6 |
| [nk_k3_n0010_01](K=3/size=10/nk_k3_n0010_01/README.md) | 10 | 30 | 1.56 | 3 | 8 |
| [nk_k3_n0010_02](K=3/size=10/nk_k3_n0010_02/README.md) | 10 | 30 | 1.47 | 2 | 20 |
| [nk_k3_n0010_03](K=3/size=10/nk_k3_n0010_03/README.md) | 10 | 30 | 1.49 | 2 | 4 |
| [nk_k3_n0010_04](K=3/size=10/nk_k3_n0010_04/README.md) | 10 | 30 | 1.50 | 2 | 9 |
| [nk_k3_n0010_05](K=3/size=10/nk_k3_n0010_05/README.md) | 10 | 30 | 1.41 | 5 | 8 |
| [nk_k3_n0010_06](K=3/size=10/nk_k3_n0010_06/README.md) | 10 | 30 | 1.60 | 1 | 5 |
| [nk_k3_n0010_07](K=3/size=10/nk_k3_n0010_07/README.md) | 10 | 30 | 1.53 | 3 | 6 |
| [nk_k3_n0010_08](K=3/size=10/nk_k3_n0010_08/README.md) | 10 | 30 | 1.23 | 1 | 5 |
| [nk_k3_n0010_09](K=3/size=10/nk_k3_n0010_09/README.md) | 10 | 30 | 1.49 | 2 | 13 |
| [nk_k3_n0050_00](K=3/size=50/nk_k3_n0050_00/README.md) | 50 | 150 | 1.58 | not enumerated | - |
| [nk_k3_n0050_01](K=3/size=50/nk_k3_n0050_01/README.md) | 50 | 150 | 1.43 | not enumerated | - |
| [nk_k3_n0050_02](K=3/size=50/nk_k3_n0050_02/README.md) | 50 | 150 | 1.38 | not enumerated | - |
| [nk_k3_n0050_03](K=3/size=50/nk_k3_n0050_03/README.md) | 50 | 150 | 1.36 | not enumerated | - |
| [nk_k3_n0050_04](K=3/size=50/nk_k3_n0050_04/README.md) | 50 | 150 | 1.54 | not enumerated | - |
| [nk_k3_n0050_05](K=3/size=50/nk_k3_n0050_05/README.md) | 50 | 150 | 1.61 | not enumerated | - |
| [nk_k3_n0050_06](K=3/size=50/nk_k3_n0050_06/README.md) | 50 | 150 | 1.52 | not enumerated | - |
| [nk_k3_n0050_07](K=3/size=50/nk_k3_n0050_07/README.md) | 50 | 150 | 1.58 | not enumerated | - |
| [nk_k3_n0050_08](K=3/size=50/nk_k3_n0050_08/README.md) | 50 | 150 | 1.68 | not enumerated | - |
| [nk_k3_n0050_09](K=3/size=50/nk_k3_n0050_09/README.md) | 50 | 150 | 1.47 | not enumerated | - |
| [nk_k3_n0250_00](K=3/size=250/nk_k3_n0250_00/README.md) | 250 | 750 | 1.42 | not enumerated | - |
| [nk_k3_n0250_01](K=3/size=250/nk_k3_n0250_01/README.md) | 250 | 750 | 1.50 | not enumerated | - |
| [nk_k3_n0250_02](K=3/size=250/nk_k3_n0250_02/README.md) | 250 | 750 | 1.54 | not enumerated | - |
| [nk_k3_n0250_03](K=3/size=250/nk_k3_n0250_03/README.md) | 250 | 750 | 1.40 | not enumerated | - |
| [nk_k3_n0250_04](K=3/size=250/nk_k3_n0250_04/README.md) | 250 | 750 | 1.46 | not enumerated | - |
| [nk_k3_n0250_05](K=3/size=250/nk_k3_n0250_05/README.md) | 250 | 750 | 1.47 | not enumerated | - |
| [nk_k3_n0250_06](K=3/size=250/nk_k3_n0250_06/README.md) | 250 | 750 | 1.50 | not enumerated | - |
| [nk_k3_n0250_07](K=3/size=250/nk_k3_n0250_07/README.md) | 250 | 750 | 1.34 | not enumerated | - |
| [nk_k3_n0250_08](K=3/size=250/nk_k3_n0250_08/README.md) | 250 | 750 | 1.44 | not enumerated | - |
| [nk_k3_n0250_09](K=3/size=250/nk_k3_n0250_09/README.md) | 250 | 750 | 1.47 | not enumerated | - |

