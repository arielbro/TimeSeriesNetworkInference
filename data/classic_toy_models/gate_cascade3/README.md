# gate_cascade3

three-step cascade, alternating sign

4 nodes, 3 edges. Every function is a symmetric threshold function: it fires when at least `threshold` of its signed inputs agree (a `+` input agreeing means it is on, a `-` input agreeing means it is off).

Feed-forward, so the dynamics are not the point: held inputs settle the circuit into a fixed point within a few steps, and all 2 attractors are fixed points - one per input pattern (1 input node(s): a).

## Functions

| node | function | threshold |
|---|---|---|
| `a` | `input` | 0 |
| `u` | `a` | 1 |
| `v` | `NOT u` | 1 |
| `z` | `v` | 1 |

## Truth table

| a | u | v | z |
|---|---|---|---|
| 0 | 0 | 1 | 1 |
| 1 | 1 | 0 | 0 |

![interaction graph and step response](network.png)
