# Thomas Algorithm (TDMA) for Tridiagonal Systems

Teaching material for a CFD course: the Thomas algorithm for tridiagonal linear systems, derived
step by step and then run on example systems from CFD.

## Interactive page

**Open it here: https://computationalgasdynamicslab.github.io/ME464/thomas-algorithm/**

It runs in any web browser, with nothing to install. It has two parts.

**1 · Derive the algorithm.** A playable symbolic derivation on a general 5×5 system. Each step
shows one row operation on the augmented matrix $`[A \mid d]`$, with the changed entries highlighted,
and the derived formulas build up into the complete algorithm:

- the tridiagonal system and the plan (reduce to a unit upper-bidiagonal matrix)
- row 1: divide by $`b_1`$
- row 2: eliminate $`x_1`$ with $`R_2 \leftarrow R_2 - a_2 R_1`$, then divide by the pivot $`m_2`$
- any row $`i`$: the forward-sweep recurrence
- back substitution from $`x_n`$ up to $`x_1`$
- operation count, and why diagonal dominance keeps every pivot nonzero

**2 · Run an example.** The same algorithm with real numbers. You can:

- step forward and back, or play, one row operation at a time (keys: ← → and Space)
- watch the full matrix, the working vectors $`a, b, c, d, m, c', d', x`$ and the pseudocode update
  together; each step marks the values it reads and the values it writes
- pick an example system: a warm-up with exact answer $`x = 1`$, 1D steady conduction,
  one implicit (BTCS) step of the heat equation, convection–diffusion with cell Péclet number 1,
  or a random diagonally dominant system
- choose the number of unknowns $`n`$ (3 to 10) and edit any coefficient
- turn on **Pause to predict** to hide each new value until you reveal it
- see the zero-pivot failure when the system is not diagonally dominant

## The algorithm

For $`a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i`$, $`i = 1, \dots, n`$, with $`a_1 = c_n = 0`$:

```math
c'_1 = \frac{c_1}{b_1}, \qquad d'_1 = \frac{d_1}{b_1}
```

```math
m_i = b_i - a_i c'_{i-1}, \qquad c'_i = \frac{c_i}{m_i}, \qquad d'_i = \frac{d_i - a_i d'_{i-1}}{m_i}, \qquad i = 2, \dots, n
```

```math
x_n = d'_n, \qquad x_i = d'_i - c'_i x_{i+1}, \qquad i = n-1, \dots, 1
```

The cost is about $`8n`$ operations, compared with about $`\tfrac{2}{3}n^3`$ for dense Gaussian elimination.
