# ME464 · 1D Burgers Equation Lab

An interactive finite-volume laboratory for Burgers’ equation. The page introduces the **non-conservative form** first:

```math
u_t + u u_x = \nu u_{xx},\qquad 0<x<1,
```

then the **conservative form**:

```math
u_t + [f(u)]_x = \nu u_{xx},\qquad f(u)=u^2/2.
```

The forms are equivalent for differentiable solutions. Inviscid shocks require the conservative weak formulation and entropy selection; set ν = 0 for the inviscid model.

Open `dist/index.html` directly; all calculations run in the browser with no installation or build step. To serve locally:

```bash
python3 -m http.server 8766 --directory dist
```

The Sites identity and static directory are recorded in `.openai/hosting.json`. The hosted version is initially owner-private.

Live owner-private Site: https://me464-burgers-lab.zhangchonglin.chatgpt.site

## Illustrated introduction for students

The page now begins with **Before you press Play**, before the live experiment. It explains the velocity-like meaning of u, local wave speed, nonlinear transport, viscosity, and the distinction between this scalar model and a full gas-dynamics system. Three diagrams compare compression, rarefaction, and a viscous layer.

A selected-problem initial-condition plot and formula follow the main controls, along with boundary illustrations that switch between periodic ghost-cell copies and exact exterior data. The notes distinguish cell averages from pointwise curves and explain why boundary fluxes matter to conservation.

The finite-volume illustration shows neighboring cell averages, face fluxes, and the incoming-minus-outgoing update. Four algorithm summaries explain Godunov, Rusanov, MUSCL/SSP-RK2, and MacCormack. An independent five-cell worked example runs the actual solver kernels for one selected method; it reveals the new averages, center-cell update and conserved total without advancing the main experiment. A Godunov hand calculation yields center value 0.92.

The final introductory reading guide explains profiles, signed errors, norms, conservation residuals, time-step limits and the post-breaking reference cutoff. The introductory diagrams remain visible on small screens and include text descriptions for accessibility.

## Learning objectives

1. Predict nonlinear characteristic speeds, shock speed, rarefaction structure and shock-formation time.
2. Compare first-order Godunov and Rusanov with MUSCL/SSP-RK2 and MacCormack.
3. Distinguish physical viscosity, numerical diffusion, bounded oscillations and growing instability.
4. Verify a conservative update using periodic mass conservation or an integrated boundary-flux budget.
5. Measure errors against exact **cell averages** and observe refinement at a common physical time.

## Experiments and references

| Case | Initial data / exact reference | Boundaries |
|---|---|---|
| Moving shock | 1 → 0 at x₀ = 0.3; shock speed s = 0.5 | Exact whole-line ghost-cell averages |
| Rarefaction | 0 → 1 at x₀ = 0.3; u = (x−x₀)/t inside the fan | Exact whole-line ghost-cell averages |
| Transonic rarefaction | −1 → 1 at x₀ = 0.5; entropy fan through u = 0 | Exact whole-line ghost-cell averages |
| Stationary shock | 1 → −1 at x₀ = 0.5; s = 0 | Exact whole-line ghost-cell averages |
| Smooth sine (default) | u₀ = ½ + A sin(2πx); tₛ = 1/(2πA) | Periodic |
| Two smooth sine waves | u₀ = ½ + A sin(4πx); tₛ = 1/(4πA) | Periodic |
| Smooth periodic pulse | u₀ = 0.2 + A exp{4[cos(2π(x−0.3))−1]}; analytic gradient sets tₛ | Periodic |
| Viscous traveling shock | s − (Δu/2)tanh[Δu(x−x₀−st)/(4ν)], with 1 → 0 | Exact traveling-wave ghost-cell averages |

The inviscid Riemann references are entropy solutions. The smooth-profile references solve x = ξ + u₀(ξ)t by bisection, uniquely before the classical shock time. After tₛ, the numerical solutions continue, but the exact line and error norms are explicitly unavailable. This is not a post-shock entropy reference.

Exact averages for shocks, rarefactions and viscous shocks are evaluated analytically. The pre-shock smooth profiles use eight-point Gauss–Legendre quadrature of the characteristic solution. The plotted exact line is pointwise; the error comparison uses cell averages. The traveling shock uses a stable log-cosh antiderivative, including at small viscosity. Initial sine averages are analytic; the smooth pulse initial averages use the same eight-point quadrature.

For nonperiodic cases, the ghost values are supplied by the exact solution at each numerical stage. They are not fixed wall or generic outflow boundaries. This keeps a controlled reference experiment even if a wave reaches an endpoint. Two ghost layers support MUSCL reconstruction.

## Numerical methods

All methods use the conservative update

```math
\bar u_i^{n+1}=\bar u_i^n-\frac{\Delta t}{\Delta x}(F_{i+1/2}-F_{i-1/2}),\qquad
F=F_{adv}-\nu\frac{\bar u_{i+1}-\bar u_i}{\Delta x}.
```

- **Godunov:** exact scalar Riemann flux, piecewise constant states, forward Euler.
- **Rusanov:** local Lax–Friedrichs flux with interface speed max(|uL|,|uR|), forward Euler.
- **MUSCL / SSP-RK2:** MC-limited linear reconstruction with Godunov flux; two-stage SSP-RK2. Centered diffusion uses the cell means, independently of advective reconstruction.
- **MacCormack:** conservative forward-flux predictor and backward-flux corrector, with centered explicit diffusion in both stages. It is second order for smooth resolved solutions but has no nonlinear monotonicity guarantee.

The same Δt is used for every active run. The default adaptive policy takes

```math
\Delta t=\sigma\frac{\Delta x}{\max|u|+2\nu/\Delta x},\qquad \sigma=0.45.
```

The maximum includes all active numerical methods and boundary ghost values. Fixed mode uses the initial maximum; nonlinear evolution may subsequently violate the intended CFL. The last step is shortened to reach the requested final time exactly.

The cards report actual C + 2D using the largest speed over both stages, where C = max|u|Δt/Δx and D = νΔt/Δx². C + 2D ≤ 1 is the standard first-order working bound; ≤ 0.5 is used as a conservative guideline for the implemented MUSCL update. The MacCormack C ≤ 1 criterion is a linearized inviscid constant-state result, not a guarantee of entropy selection or nonlinear stability. The default 0.45 also allows for explicit diffusion.

## Read the figures

The lab starts from a continuous sine profile. The **Smooth initial data** group contains a single sine (A = 0.5), a two-wave sine (A = 0.5), and a Gaussian-like periodic pulse (A = 0.8). All three are smooth at the periodic seam. The amplitude slider spans 0.1–1.0; changing it resets all runs and sets the default final time to 1.6 times the predicted breaking time (final times up to 5 are supported).

The initial-profile panel shows the formula, initial curve, predicted tₛ and current stage (initial, steepening, or after breaking). The **Start smooth → form a shock** preset restores the single-sine experiment. The characteristic diagram uses the chosen profile and marks tₛ. The numerical curves continue after breaking, while the classical reference and corresponding error norms stop.

For inviscid Burgers, tₛ = −1/min u₀′. For the pulse with κ = 4, write c = (sqrt(1+4κ²)−1)/(2κ). Then tₛ = 1/[2πκA sqrt(1−c²) exp(κ(c−1))]. Doubling A halves the predicted breaking time. At default amplitudes, the sine breaks at approximately 0.3183, the double sine at 0.1592, and the pulse at approximately 0.1692 (the page computes the exact formula).

Every method has its own solution and signed-error figure, L₁/RMS errors, min/max, CFL load, total variation and conservation residual. Solution axes can be shared or fitted independently; error axes always fit each method and are labeled accordingly. Characteristic diagrams are guides; straight continuations beyond a crossing or shock are not physical multivalued solutions.

The conservation residual is

```math
M^n-M^0-\sum_{steps}\Delta t\,(F_{left}-F_{right}),\qquad M=\sum_i\bar u_i\Delta x.
```

For two-stage methods, boundary fluxes receive the same stage weights as the interior update. Thus nonperiodic mass change is not incorrectly labeled a conservation error. Total variation includes the seam for the periodic case, and interior cell differences for other cases.

A run freezes at its first nonfinite value or max |u| > 10³; its recorded time, error and reference remain synchronized. Other methods continue. A step limit of 20,000 and a time-step underflow guard keep the demonstration responsive. The error-history display is decimated after 1,000 samples; the numerical steps are unaffected.

CSV export includes one row per method/cell, its actual time and stopping step, cell bounds, numerical mean, exact mean and signed error. Reference fields are empty where the exact reference is unavailable.

The webpage contains six worksheet questions covering shock speed, entropy, grid convergence, viscosity, CFL and smooth-wave breaking.

## Verification

```bash
node verify.cjs
node verify-ui.cjs
node --check dist/app.js
```

The dependency-free checks exercise Riemann fluxes, constant preservation, periodic conservation (inviscid and viscous), time-dependent boundary-flux budgets, transonic entropy fans, shock propagation, exact viscous averages and reference validity. Additional smooth-profile checks verify periodicity, analytic gradients, shock times, characteristic solutions, inverse-amplitude scaling, post-breaking reference cutoff, and actual numerical steepening. A resolved traveling shock with ν = 0.04 yields L₁ orders near 1 for Godunov/Rusanov and 2 for MUSCL/MacCormack on N = 80,160,320.

## Sources

The [Clawpack Riemann Book: Burgers](https://www.clawpack.org/riemann_book/html/Burgers.html) explains the exact shock and rarefaction solutions and entropy selection. Its [approximate Burgers solvers chapter](https://www.clawpack.org/riemann_book/html/Burgers_approximate.html) discusses numerical fluxes and entropy issues. The worksheet and implementation are designed for this ME464 lab.
