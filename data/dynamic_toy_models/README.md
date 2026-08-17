# Dynamic toy models

A `graphs_dir` (the role `data/cellcollective_models` plays for the real runs). Small hand-built models whose state space keeps moving: several attractors, essentially none of them fixed points, so time series generated from them are not near-static. Every function is a symmetric threshold function.

Regenerate with `python synthetic_data_generation/make_toy_models.py`.

## ring family

One construction repeated: odd repressilator rings, which give the long cycles (period 6, 10, 30). The AND/OR readouts are feed-forward, so they add churn and connect the rings without changing the attractor count. Here the cycle IS the structure - see the dense family for the same behaviour without that.

| model | nodes | edges | attractors | fixed points | periods | mean nodes changing per step |
|---|---|---|---|---|---|---|
| [size03_ring3](size03_ring3/README.md) | 3 | 3 | 2 | 0 | 2, 6 | 1.50 |
| [size04_ring3_and](size04_ring3_and/README.md) | 4 | 5 | 2 | 0 | 2, 6 | 2.00 |
| [size05_ring5](size05_ring5/README.md) | 5 | 5 | 4 | 0 | 2, 10 | 2.50 |
| [size06_ring5_and](size06_ring5_and/README.md) | 6 | 7 | 4 | 0 | 2, 10 | 2.88 |
| [size07_two_ring3_and](size07_two_ring3_and/README.md) | 7 | 8 | 12 | 0 | 2, 6 | 3.38 |
| [size08_two_ring3_and_or](size08_two_ring3_and_or/README.md) | 8 | 10 | 12 | 0 | 2, 6 | 3.75 |
| [size09_ring5_ring3_and](size09_ring5_ring3_and/README.md) | 9 | 10 | 16 | 0 | 2, 6, 10, 30 | 4.38 |
| [size10_ring5_ring3_and_or](size10_ring5_ring3_and_or/README.md) | 10 | 12 | 16 | 0 | 2, 6, 10, 30 | 4.75 |

## motif family

Textbook transcription motifs. The movement comes from the activator-inhibitor pairs and the NAR nodes, which is why every attractor here has period 4; the toggles supply multiplicity rather than churn, sitting at one of their fixed points in most attractors while the AI pair oscillates around them. A NAR node removes fixed-point attractors from the whole model, since it flips every step. FFLs and SIMs are feed-forward.

| model | nodes | edges | attractors | fixed points | periods | mean nodes changing per step |
|---|---|---|---|---|---|---|
| [size03_ai_nar](size03_ai_nar/README.md) | 3 | 3 | 2 | 0 | 4 | 2.00 |
| [size04_ai_toggle](size04_ai_toggle/README.md) | 4 | 4 | 4 | 0 | 4 | 2.00 |
| [size05_ai_i1ffl_nar](size05_ai_i1ffl_nar/README.md) | 5 | 6 | 2 | 0 | 4 | 3.00 |
| [size06_toggle_ai_i1ffl](size06_toggle_ai_i1ffl/README.md) | 6 | 7 | 4 | 0 | 4 | 3.00 |
| [size07_ai_toggle_c1ffl_nar](size07_ai_toggle_c1ffl_nar/README.md) | 7 | 8 | 8 | 0 | 4 | 4.00 |
| [size08_two_ai_toggle_i1ffl](size08_two_ai_toggle_i1ffl/README.md) | 8 | 9 | 16 | 0 | 4 | 4.00 |
| [size09_ai_toggle_sim_i1ffl](size09_ai_toggle_sim_i1ffl/README.md) | 9 | 10 | 4 | 0 | 4 | 4.50 |
| [size10_two_ai_toggle_i1ffl_nar](size10_two_ai_toggle_i1ffl_nar/README.md) | 10 | 11 | 32 | 0 | 4 | 5.50 |

## dense family

Multiple moving attractors with no part of the model being a ring: in a clique every node is regulated by all the others, and in a chorded loop each node is regulated by the three before it, so every cycle is a proper subgraph of a denser neighbourhood. Each clique/loop uses one uniform rule for all its nodes (chosen by a small deterministic search over sign patterns and thresholds), so the whole group is still describable in a sentence.

| model | nodes | edges | attractors | fixed points | periods | mean nodes changing per step |
|---|---|---|---|---|---|---|
| [size03_clique3](size03_clique3/README.md) | 3 | 6 | 2 | 1 | 1, 3 | 1.50 |
| [size04_clique4](size04_clique4/README.md) | 4 | 12 | 3 | 0 | 2, 4 | 3.00 |
| [size05_chorded_loop5](size05_chorded_loop5/README.md) | 5 | 15 | 4 | 0 | 2, 5, 10 | 2.73 |
| [size06_chorded_loop6](size06_chorded_loop6/README.md) | 6 | 18 | 4 | 0 | 2, 6 | 2.40 |
| [size07_clique4_clique3](size07_clique4_clique3/README.md) | 7 | 18 | 6 | 0 | 2, 4, 6, 12 | 4.50 |
| [size08_clique4_clique3_readout](size08_clique4_clique3_readout/README.md) | 8 | 20 | 6 | 0 | 2, 4, 6, 12 | 5.12 |
| [size09_two_clique4_readout](size09_two_clique4_readout/README.md) | 9 | 26 | 20 | 0 | 2, 4 | 6.47 |
| [size10_clique4_chorded6](size10_clique4_chorded6/README.md) | 10 | 30 | 24 | 0 | 2, 4, 6, 12 | 5.40 |

