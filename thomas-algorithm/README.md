# Thomas Algorithm (TDMA) for Tridiagonal Systems

Teaching material for a CFD course: the Thomas algorithm for tridiagonal linear systems, derived
step by step and then run on example systems from CFD.

## Interactive page

**Open it here: https://computationalgasdynamicslab.github.io/ME464/thomas-algorithm/**

It runs in any web browser, with nothing to install. An illustrated introduction connects a rod to the linear system, and three tabs then move from general Gaussian elimination to the Thomas algorithm and a worked example.

**Notation.** The right-hand side of $`A\mathbf{x} = \mathbf{f}`$ is written $`f`$ throughout the page and this README, matching the course lecture notes; the modified right-hand side produced by the forward sweep is $`f'_i`$. Many textbooks call the same vector $`d`$ or $`r`$.

**0 · Gauss elimination first.** The new opening tab runs ordinary, structure-blind Gaussian elimination one row operation at a time, so the Thomas algorithm arrives as a specialization rather than a recipe. You can:

- step through a **dense 4×4** system (exact answer $`x = (1, 2, -1, 1)`$) and watch every elimination touch every remaining column
- run the identical procedure on **tridiagonal** systems (4×4, the 5×5 of part 1, and the 3-unknown BTCS rod) and see that most multipliers are zero and that no elimination ever creates a new nonzero — tridiagonal elimination has no fill-in
- read the multipliers as they accumulate into the $`L`$ of $`A = LU`$, with $`U`$ the matrix on screen
- watch an honest operation counter: arithmetic actually performed, what the textbook triple loop would execute, and how much of it was spent on entries already known to be zero
- hide the zero-multiplier rows once the point is made, and follow the **From Gauss to Thomas** panel, which shows that the Gauss pivot $`u_{ii} = b_i - \ell_i c_{i-1}`$ and the Thomas pivot $`m_i = b_i - a_i c'_{i-1}`$ are the same number, since $`c'_{i-1} = c_{i-1}/m_{i-1}`$

The summary step states the three differences that turn $`\tfrac{2}{3}n^3`$ into $`8n`$ — store three diagonals, do one elimination per column, update one coefficient and one right-hand side per elimination — and the one thing Thomas gives up, the row exchange that general Gaussian elimination uses for stability.

**Start with a physical problem.** Select one of three interior nodes, adjust r and the wall/old temperatures, and watch its BTCS equation and matrix row highlight together. There are three unknown temperatures and four spatial intervals; the known wall values move to the right-hand side. **Solve these three equations below** loads that exact assembled system into the numerical stepper.

The two-pass overview explains the direction and purpose of each sweep. In the normalized row xᵢ + c′ᵢxᵢ₊₁ = f′ᵢ, the primes mean *modified coefficients*, not derivatives. A repeated-BTCS example factors the rod matrix once and solves three successive right-hand sides, using each new temperature profile as the next old profile.

**1 · Derive the algorithm.** A playable symbolic derivation on a general 5×5 system. Each step
shows one row operation on the augmented matrix $`[A \mid f]`$, with the changed entries highlighted,
and the derived formulas build up into the complete algorithm:

- the tridiagonal system and the plan (reduce to a unit upper-bidiagonal matrix)
- row 1: divide by $`b_1`$
- row 2: eliminate $`x_1`$ with $`R_2 \leftarrow R_2 - a_2 R_1`$, then divide by the pivot $`m_2`$
- any row $`i`$: the forward-sweep recurrence
- back substitution from $`x_n`$ up to $`x_1`$
- operation count, and why diagonal dominance keeps every pivot nonzero

**2 · Run an example.** The same algorithm with real numbers. You can:

- step forward and back, or play, one row operation at a time (keys: ← → and Space)
- watch the full matrix, the working vectors $`a, b, c, f, m, c', f', x`$ and the pseudocode update
  together; each step marks the values it reads and the values it writes
- pick an example system: a warm-up with exact answer $`x = 1`$, 1D steady conduction,
  one implicit (BTCS) step of the heat equation, convection–diffusion with cell Péclet number 1,
  a random diagonally dominant system, or a diagnostic failure/scaling example
- choose the number of unknowns $`n`$ (3 to 10) and edit any coefficient
- turn on **Pause to predict**, enter the current f′ or x, and get sign/index/pivot feedback before revealing the step
- check an individual original equation after solving, and read the normalized infinity-norm residual
- distinguish a nonsingular system needing pivoting, a singular disconnected block, a tiny relative pivot, and a uniformly scaled warm-up

The default warm-up has three equations and the hand-checkable solution x = [1,1,1]. Numerical examples still support 3–10 unknowns.

## Old temperature, boundary contribution, and time-step figures

The node equation now explicitly shows its old-temperature RHS before substituting a known wall. With the default values, node 1 is

```math
-0.5T_0^{k+1}+2T_1^{k+1}-0.5T_2^{k+1}=T_1^k=20.
```

Here 20 is the **old interior temperature**, independently specified from the walls. Substituting the known left wall T₀ = 100 and moving its term to the right produces

```math
2T_1^{k+1}-0.5T_2^{k+1}=f_1=20+0.5(100)=70.
```

Both stages are shown before the matrix. Selecting a different node or editing the old/wall temperatures updates the explanation.

**Watch the temperature after each BTCS time step** shows a separate figure for the initial profile and each completed solve. **Compute next BTCS time step** advances the latest solution using the previous step's temperatures; **Take 3 BTCS steps** starts a three-step history, and **Restart** returns to the initial profile. The matrix factors are reused throughout a history. Up to 20 steps can be displayed.

Select a saved figure to inspect its actual old temperatures, boundary contributions, assembled f, new interior temperatures and scaled residual. Selecting an older figure does not change which state the next step advances. Selecting another rod node preserves the history; changing rod parameters resets it.

All figures use shared temperature axes and include both fixed wall nodes. Blue identifies the current time level, dashed gray the previous one, and orange squares the prescribed walls. Lines join the selected grid nodes; they are not exact continuum curves. The uniform old interior and imposed walls can disagree initially, which is stated beside the initial figure. For the illustrative L = 1 and α = 0.01, Δx = 1/N and r = αΔt/Δx². The default N = 4 and Δt = 3.125 match the three-node example at r = 0.5.

**Refine the temperature grid.** The history has independent controls for N (4–400 spatial intervals) and Δt (default 3.125). N intervals give N + 1 nodes, including the two walls, and N − 1 interior unknowns. Try N = 40, 100 or 400 for more spatial detail. The Thomas matrix is assembled and factored for the actual selected grid; the three-unknown equation demonstration and row-operation stepper stay small enough to follow by hand. Wall and initial temperatures are shared. The top r control belongs to the three-unknown algebra example; the history computes its own r from N and Δt.

Changing N holds Δt fixed, so the same number of steps reaches the same physical time across grids. Changing N or Δt restarts the history and refactors the matrix. Larger N reduces spatial discretization error, while temporal error still depends on Δt; unconditional BTCS stability does not imply accuracy. Dense-grid figures use continuous node-joining lines and omit interior markers above N = 80. Long vectors show their first/last three entries and the full entry count.

A Thomas elimination/back-substitution operation is a stage of solving **one** linear system. It is not another physical time step. The original partial-x plot and new BTCS history explicitly distinguish these two meanings of “step.”

## The algorithm

For $`a_i x_{i-1} + b_i x_i + c_i x_{i+1} = f_i`$, $`i = 1, \dots, n`$, with $`a_1 = c_n = 0`$:

```math
c'_1 = \frac{c_1}{b_1}, \qquad f'_1 = \frac{f_1}{b_1}
```

```math
m_i = b_i - a_i c'_{i-1}, \qquad c'_i = \frac{c_i}{m_i}, \qquad f'_i = \frac{f_i - a_i f'_{i-1}}{m_i}, \qquad i = 2, \dots, n
```

```math
x_n = f'_n, \qquad x_i = f'_i - c'_i x_{i+1}, \qquad i = n-1, \dots, 1
```

The cost is about $`8n`$ operations, compared with about $`\tfrac{2}{3}n^3`$ for dense Gaussian elimination.

## Failure cases and numerical checks

The solver rejects an exactly zero pivot or nonfinite arithmetic. It does **not** reject a uniformly small but well-scaled matrix using an absolute 10⁻¹² cutoff. It flags pivots smaller than 10⁻¹⁰ times their local coefficient/elimination scale, while allowing the demonstration to continue so students can inspect cancellation and the resulting residual. Small numbers use scientific notation rather than being rounded to zero on screen.

The input badge distinguishes strict diagonal dominance, irreducible weak diagonal dominance with at least one strict row, and symmetric positive definiteness. These are sufficient non-breakdown conditions in exact arithmetic, not requirements for every solvable tridiagonal system. A missing certificate does not mean singularity. In particular, one strict row does not protect a disconnected singular block. The derivation and error explanation now state this explicitly.

The final verification uses the original coefficients, not the transformed identity matrix. Select a row to inspect its three contributions, right-hand side and residual. The scale-aware check is

```math
\eta = \frac{\|A x-f\|_\infty}{\|A\|_\infty\|x\|_\infty+\|f\|_\infty}.
```

A small scaled residual is evidence of a small backward error; forward solution error also depends on conditioning. The tiny-pivot example illustrates that elimination without pivoting can lose accuracy even when the matrix is nonsingular. General-purpose tridiagonal solvers can use [partial pivoting, as in LAPACK DGTSV](https://netlib.org/lapack/explore-html/d1/db3/dgtsv_8f_source.html).

## Verification

```bash
node verify.cjs
```

The dependency-free checks verify the Gaussian-elimination trace against a pivoted dense reference, confirm $`A = LU`$ from the stored multipliers, assert that tridiagonal elimination produces no fill-in and no spurious multipliers, check that each Gauss pivot equals the corresponding Thomas pivot, and confirm that the two operation counters differ by exactly the work spent on known zeros. A further check renders every step of every Gauss system and fails on any unresolved value, and one more asserts that the right-hand side appears as $`f`$, never $`d`$, in the rendered panels. The remaining checks compare the actual solver with a pivoted dense reference, verify the original-equation residuals, preserve the uniformly scaled solution, distinguish nonsingular zero-pivot and singular examples, exercise the positive-definite certificate, and compare factor reuse with fresh solves. A minimal DOM smoke test executes the real page initialization, rod loading, prediction feedback, row checks and failure-case rendering. Variable-grid checks through N = 400 verify a discrete sine eigenmode, spatial convergence for a smooth solution, maximum-principle bounds at large diffusion numbers, and fixed-time comparisons when N changes. MathJax rendering and visual browser layout are outside the DOM adapter's coverage.
