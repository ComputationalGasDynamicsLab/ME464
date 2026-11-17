# Thomas Algorithm (TDMA) for Tridiagonal Systems

Teaching material for a CFD course: the Thomas algorithm for tridiagonal linear systems, derived
step by step and then run on example systems from CFD.

## Interactive page

**Open it here: https://computationalgasdynamicslab.github.io/ME464/thomas-algorithm/**

It runs in any web browser, with nothing to install. A new illustrated introduction connects a rod to the linear system before the two algorithm views.

**Start with a physical problem.** Select one of three interior nodes, adjust r and the wall/old temperatures, and watch its BTCS equation and matrix row highlight together. There are three unknown temperatures and four spatial intervals; the known wall values move to the right-hand side. **Solve these three equations below** loads that exact assembled system into the numerical stepper.

The two-pass overview explains the direction and purpose of each sweep. In the normalized row xᵢ + c′ᵢxᵢ₊₁ = d′ᵢ, the primes mean *modified coefficients*, not derivatives. A repeated-BTCS example factors the rod matrix once and solves three successive right-hand sides, using each new temperature profile as the next old profile.

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
  a random diagonally dominant system, or a diagnostic failure/scaling example
- choose the number of unknowns $`n`$ (3 to 10) and edit any coefficient
- turn on **Pause to predict**, enter the current d′ or x, and get sign/index/pivot feedback before revealing the step
- check an individual original equation after solving, and read the normalized infinity-norm residual
- distinguish a nonsingular system needing pivoting, a singular disconnected block, a tiny relative pivot, and a uniformly scaled warm-up

The default warm-up has three equations and the hand-checkable solution x = [1,1,1]. Numerical examples still support 3–10 unknowns.

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

## Failure cases and numerical checks

The solver rejects an exactly zero pivot or nonfinite arithmetic. It does **not** reject a uniformly small but well-scaled matrix using an absolute 10⁻¹² cutoff. It flags pivots smaller than 10⁻¹⁰ times their local coefficient/elimination scale, while allowing the demonstration to continue so students can inspect cancellation and the resulting residual. Small numbers use scientific notation rather than being rounded to zero on screen.

The input badge distinguishes strict diagonal dominance, irreducible weak diagonal dominance with at least one strict row, and symmetric positive definiteness. These are sufficient non-breakdown conditions in exact arithmetic, not requirements for every solvable tridiagonal system. A missing certificate does not mean singularity. In particular, one strict row does not protect a disconnected singular block. The derivation and error explanation now state this explicitly.

The final verification uses the original coefficients, not the transformed identity matrix. Select a row to inspect its three contributions, right-hand side and residual. The scale-aware check is

```math
\eta = \frac{\|A x-d\|_\infty}{\|A\|_\infty\|x\|_\infty+\|d\|_\infty}.
```

A small scaled residual is evidence of a small backward error; forward solution error also depends on conditioning. The tiny-pivot example illustrates that elimination without pivoting can lose accuracy even when the matrix is nonsingular. General-purpose tridiagonal solvers can use [partial pivoting, as in LAPACK DGTSV](https://netlib.org/lapack/explore-html/d1/db3/dgtsv_8f_source.html).

## Verification

```bash
node verify.cjs
```

The dependency-free checks compare the actual solver with a pivoted dense reference, verify the original-equation residuals, preserve the uniformly scaled solution, distinguish nonsingular zero-pivot and singular examples, exercise the positive-definite certificate, and compare factor reuse with fresh solves. A minimal DOM smoke test executes the real page initialization, rod loading, prediction feedback, row checks and failure-case rendering. MathJax rendering and visual browser layout are outside the DOM adapter's coverage.
