# Linear convection: six schemes and their stability

An interactive companion to `ftbs-linear-convection` and `ftbs-stability`.
Open `dist/index.html` directly in a browser; no dependencies or build step are needed.
For a local HTTP preview, run `python3 -m http.server 8765 --directory dist` in this folder.

The page compares FTBS, FTCS, leapfrog, implicit BTCS (backward Euler / centered space), Lax–Wendroff, and MacCormack for `u_t + a u_x = 0` on a periodic domain `[0,2)`.

- All six numerical solutions advance together against the exact translated profile.
- A dedicated initial-condition graph and formulas explain each profile, including the extra periodic seam jump for the step. An optional initial-profile overlay remains on every solution plot.
- A connected-domain diagram tracks a characteristic through wraparound. Live ghost-cell values show how the selected scheme copies periodic neighbors at its current time level.
- Signed Courant number `C = a dt/dx` ranges from −2 to 2; grid and profile are adjustable.
- At nonzero C, `a = sign(C)` and `dt = |C| dx`. At C = 0, `a = 0` and `dt = dx`.
- Inspect amplification factors across wavelengths and evolve a selected small Fourier perturbation directly, without subtraction cancellation.
- See update formulas, stability derivations for the demonstrated methods, and the leapfrog repeated-root endpoint caveat.
- The implicit method remains available for experimentation, but its analytical amplification factor, stability verdict, prediction curve, and stated stability condition are withheld for homework. The forward-time/forward-space method is omitted from this student lab.
- Download the current solution arrays as CSV. A frozen run includes its actual stopping step.

## Stability summary

| Method | Requirement (constant coefficients, periodic grid) |
|---|---|
| FTBS | 0 ≤ C ≤ 1 |
| FTCS | Unstable for nonzero C |
| Leapfrog | \|C\| < 1 robust root condition; usual CFL bound \|C\| ≤ 1 has repeated-root endpoints |
| BTCS | Reserved for student derivation |
| Lax–Wendroff | \|C\| ≤ 1 |
| MacCormack | \|C\| ≤ 1; identical to LW for this linear problem |

Leapfrog solution runs start with one Lax–Wendroff step. The perturbation experiment also offers an independent zero second time level to expose repeated-root growth at the CFL endpoint. The solution runs use no time filter. BTCS solves its periodic circulant system by direct discrete Fourier transform. Solution runs stop individually at max |u| > 10^6, perturbations at amplitude > 10^12, and the global run stops at 600 steps. Per-plot vertical scaling is explicit. The perturbation experiment models propagation of an initial error, not continuous injection of machine round-off. Its amplitude is measured by projecting onto its initial Fourier mode; for leapfrog the spectral-radius curve is a guide rather than an exact two-root transient solution.

## Verification

Run `node verify.cjs`. It checks Fourier-mode response, periodic conservation, the BTCS linear-system residual, LW/MacCormack equivalence, stability bounds, and leapfrog endpoint growth for incompatible two-level data. Run `node --check dist/app.js` for UI-script syntax.

The deployable static source lives in `dist/`; `.openai/hosting.json` records the private hosted Site.
