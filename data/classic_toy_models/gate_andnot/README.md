# gate_andnot

AND NOT: one activator, one repressor

3 nodes, 2 edges. Every function is a symmetric threshold function: it fires when at least `threshold` of its signed inputs agree (a `+` input agreeing means it is on, a `-` input agreeing means it is off).

Feed-forward, so the dynamics are not the point: held inputs settle the circuit into a fixed point within a few steps, and all 4 attractors are fixed points - one per input pattern (2 input node(s): a, b).

## Functions

| node | function | threshold |
|---|---|---|
| `a` | `input` | 0 |
| `b` | `input` | 0 |
| `z` | `a AND NOT b` | 2 |

## Truth table

| a | b | z |
|---|---|---|
| 0 | 0 | 0 |
| 0 | 1 | 0 |
| 1 | 0 | 1 |
| 1 | 1 | 0 |

![interaction graph and step response](network.png)
