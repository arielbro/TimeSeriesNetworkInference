# Classic toy models

A `graphs_dir` (the role `data/cellcollective_models` plays for the real runs). Small hand-built feed-forward circuits - the classic ones. Their attractors are all fixed points (one per input pattern), so they are for testing function and topology inference rather than dynamics. Generate data from them with `only_attractors = False`: walking to an attractor first would start every trajectory at a fixed point and make every row identical, whereas a random start captures the transient - which is where the pulse of an incoherent FFL and the delay of a coherent one actually show up. Every function is a symmetric threshold function.

Regenerate with `python synthetic_data_generation/make_toy_models.py`.

## ffl family

All eight non-equivalent three-node feed-forward loops, each with an AND and an OR gate at the target. x is an input node; y = f(x); z = g(x, y). Coherent types (C1-C4) have sign(x->z) = sign(x->y) * sign(y->z), incoherent ones (I1-I4) do not - which is what makes I1 a pulse generator and C1 a sign-sensitive delay. The step-response panel drives x up and back down, so the delay or pulse is visible.

| model | nodes | edges | inputs | depth | fixed points |
|---|---|---|---|---|---|
| [ffl_c1_and](ffl_c1_and/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_c1_or](ffl_c1_or/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_c2_and](ffl_c2_and/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_c2_or](ffl_c2_or/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_c3_and](ffl_c3_and/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_c3_or](ffl_c3_or/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_c4_and](ffl_c4_and/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_c4_or](ffl_c4_or/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_i1_and](ffl_i1_and/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_i1_or](ffl_i1_or/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_i2_and](ffl_i2_and/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_i2_or](ffl_i2_or/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_i3_and](ffl_i3_and/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_i3_or](ffl_i3_or/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_i4_and](ffl_i4_and/README.md) | 3 | 3 | 1 | 2 | 2 |
| [ffl_i4_or](ffl_i4_or/README.md) | 3 | 3 | 1 | 2 | 2 |

## gates family

Multi-level AND/OR combinations, the middle (majority) threshold, k-of-n thresholds, an AND NOT, a cascade, a single-input module and a dense overlapping regulon.

| model | nodes | edges | inputs | depth | fixed points |
|---|---|---|---|---|---|
| [gate_and2](gate_and2/README.md) | 3 | 2 | 2 | 1 | 4 |
| [gate_or2](gate_or2/README.md) | 3 | 2 | 2 | 1 | 4 |
| [gate_and3](gate_and3/README.md) | 4 | 3 | 3 | 1 | 8 |
| [gate_or3](gate_or3/README.md) | 4 | 3 | 3 | 1 | 8 |
| [gate_majority3](gate_majority3/README.md) | 4 | 3 | 3 | 1 | 8 |
| [gate_andnot](gate_andnot/README.md) | 3 | 2 | 2 | 1 | 4 |
| [gate_and_of_or](gate_and_of_or/README.md) | 7 | 6 | 4 | 2 | 16 |
| [gate_or_of_and](gate_or_of_and/README.md) | 7 | 6 | 4 | 2 | 16 |
| [gate_k_of_5](gate_k_of_5/README.md) | 8 | 15 | 5 | 1 | 32 |
| [gate_cascade3](gate_cascade3/README.md) | 4 | 3 | 1 | 3 | 2 |
| [gate_sim](gate_sim/README.md) | 4 | 3 | 1 | 1 | 2 |
| [gate_dor](gate_dor/README.md) | 5 | 6 | 2 | 1 | 4 |

