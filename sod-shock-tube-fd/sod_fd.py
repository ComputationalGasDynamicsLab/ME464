"""Sod shock tube with classic finite-difference schemes for the 1D Euler equations.

    U_t + F(U)_x = 0,  U = (rho, rho u, E),  F = (rho u, rho u^2 + p, u (E + p)),
    p = (gamma - 1)(E - rho u^2 / 2),  gamma = 1.4

Values live at grid points x_i = i dx (i = 0..nx-1) on 0 <= x <= 1, with lam = dt / dx:

    lf   Lax-Friedrichs (1st order)
         U_i^{n+1} = (U_{i+1}^n + U_{i-1}^n)/2 - (lam/2)(F_{i+1}^n - F_{i-1}^n)
    lw   Richtmyer two-step Lax-Wendroff (2nd order)
         U_{i+1/2}^{n+1/2} = (U_i^n + U_{i+1}^n)/2 - (lam/2)(F_{i+1}^n - F_i^n)
         U_i^{n+1} = U_i^n - lam (F_{i+1/2}^{n+1/2} - F_{i-1/2}^{n+1/2})
    mc   MacCormack (2nd order), forward predictor and backward corrector
         U*_i = U_i^n - lam (F_{i+1}^n - F_i^n)
         U_i^{n+1} = (U_i^n + U*_i - lam (F*_i - F*_{i-1}))/2

No artificial viscosity is added, so the 2nd-order schemes oscillate near the shock and contact.
The numerical solutions are compared with the exact Riemann solution at t = 0.2.

Examples:
    python sod_fd.py                      # all three schemes, 101 points, CFL 0.8
    python sod_fd.py --cfl 0.3            # smaller time step: more smearing and more wiggles
    python sod_fd.py --scheme mc --nx 401 # one scheme, finer grid
    python sod_fd.py --animate            # also save an animated GIF
"""
import argparse

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

G = 1.4  # ratio of specific heats
SOD = dict(rL=1.0, uL=0.0, pL=1.0, rR=0.125, uR=0.0, pR=0.1, x0=0.5)


# --------------------------------------------------------------------------------------
# Exact Riemann solver (Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, ch. 4)
# --------------------------------------------------------------------------------------
def _f(p, r, pk, a):
    """Pressure function f_K(p) and its derivative for the left or right state."""
    if p > pk:  # shock
        A, B = 2 / ((G + 1) * r), (G - 1) / (G + 1) * pk
        q = np.sqrt(A / (p + B))
        return (p - pk) * q, q * (1 - (p - pk) / (2 * (B + p)))
    # rarefaction
    return (2 * a / (G - 1) * ((p / pk) ** ((G - 1) / (2 * G)) - 1),
            1 / (r * a) * (p / pk) ** (-(G + 1) / (2 * G)))


def star_state(rL, uL, pL, rR, uR, pR):
    """Pressure and velocity in the star region between the left and right waves."""
    aL, aR = np.sqrt(G * pL / rL), np.sqrt(G * pR / rR)
    p = max(1e-8, 0.5 * (pL + pR) - 0.125 * (uR - uL) * (rL + rR) * (aL + aR))  # PVRS guess
    for _ in range(100):
        fL, dL = _f(p, rL, pL, aL)
        fR, dR = _f(p, rR, pR, aR)
        p_new = max(1e-10, p - (fL + fR + uR - uL) / (dL + dR))
        if abs(p_new - p) < 1e-14 * p:
            p = p_new
            break
        p = p_new
    fL, _ = _f(p, rL, pL, aL)
    fR, _ = _f(p, rR, pR, aR)
    return p, 0.5 * (uL + uR) + 0.5 * (fR - fL)


def exact_solution(x, t, rL, uL, pL, rR, uR, pR, x0):
    """Exact density, velocity and pressure at positions x and time t > 0."""
    ps, us = star_state(rL, uL, pL, rR, uR, pR)
    aL, aR = np.sqrt(G * pL / rL), np.sqrt(G * pR / rR)
    gm = (G - 1) / (G + 1)
    out = np.zeros((3, len(x)))
    for k, xi in enumerate((np.asarray(x) - x0) / t):
        if xi <= us:  # left of the contact
            if ps > pL:  # left shock
                S = uL - aL * np.sqrt((G + 1) / (2 * G) * ps / pL + (G - 1) / (2 * G))
                w = (rL, uL, pL) if xi <= S else (rL * (ps / pL + gm) / (gm * ps / pL + 1), us, ps)
            else:  # left rarefaction
                a_s = aL * (ps / pL) ** ((G - 1) / (2 * G))
                if xi <= uL - aL:
                    w = (rL, uL, pL)
                elif xi > us - a_s:
                    w = (rL * (ps / pL) ** (1 / G), us, ps)
                else:  # inside the fan
                    a = 2 / (G + 1) * (aL + (G - 1) / 2 * (uL - xi))
                    w = (rL * (a / aL) ** (2 / (G - 1)), 2 / (G + 1) * (aL + (G - 1) / 2 * uL + xi),
                         pL * (a / aL) ** (2 * G / (G - 1)))
        else:  # right of the contact
            if ps > pR:  # right shock
                S = uR + aR * np.sqrt((G + 1) / (2 * G) * ps / pR + (G - 1) / (2 * G))
                w = (rR, uR, pR) if xi >= S else (rR * (ps / pR + gm) / (gm * ps / pR + 1), us, ps)
            else:  # right rarefaction
                a_s = aR * (ps / pR) ** ((G - 1) / (2 * G))
                if xi >= uR + aR:
                    w = (rR, uR, pR)
                elif xi <= us + a_s:
                    w = (rR * (ps / pR) ** (1 / G), us, ps)
                else:
                    a = 2 / (G + 1) * (aR - (G - 1) / 2 * (uR - xi))
                    w = (rR * (a / aR) ** (2 / (G - 1)), 2 / (G + 1) * (-aR + (G - 1) / 2 * uR + xi),
                         pR * (a / aR) ** (2 * G / (G - 1)))
        out[:, k] = w
    return out


# --------------------------------------------------------------------------------------
# Finite-difference schemes
# --------------------------------------------------------------------------------------
def conserved(r, u, p):
    return np.array([r, r * u, p / (G - 1) + 0.5 * r * u * u])


def primitives(U):
    r = U[0]
    u = U[1] / r
    return r, u, (G - 1) * (U[2] - 0.5 * r * u * u)


def flux(U):
    r, u, p = primitives(U)
    return np.array([r * u, r * u * u + p, u * (U[2] + p)])


def step_lf(U, lam):
    F = flux(U)
    Un = U.copy()
    Un[:, 1:-1] = 0.5 * (U[:, 2:] + U[:, :-2]) - 0.5 * lam * (F[:, 2:] - F[:, :-2])
    return Un


def step_lw(U, lam):
    F = flux(U)
    Uh = 0.5 * (U[:, 1:] + U[:, :-1]) - 0.5 * lam * (F[:, 1:] - F[:, :-1])   # at i+1/2, t+dt/2
    Fh = flux(Uh)
    Un = U.copy()
    Un[:, 1:-1] = U[:, 1:-1] - lam * (Fh[:, 1:] - Fh[:, :-1])
    return Un


def step_mc(U, lam):
    F = flux(U)
    Us = U.copy()
    Us[:, :-1] = U[:, :-1] - lam * (F[:, 1:] - F[:, :-1])                   # predictor: forward
    Fs = flux(Us)
    Un = U.copy()
    Un[:, 1:-1] = 0.5 * (U[:, 1:-1] + Us[:, 1:-1] - lam * (Fs[:, 1:-1] - Fs[:, :-2]))  # corrector: backward
    return Un


SCHEMES = {"lf": step_lf, "lw": step_lw, "mc": step_mc}
NAMES = {"lf": "Lax-Friedrichs", "lw": "Lax-Wendroff (Richtmyer)", "mc": "MacCormack"}


def solve(name, nx, cfl, t_end, n_frames=0):
    x = np.linspace(0.0, 1.0, nx)
    dx = x[1] - x[0]
    s = SOD
    UL, UR = conserved(s["rL"], s["uL"], s["pL"]), conserved(s["rR"], s["uR"], s["pR"])
    U = np.where(x < s["x0"], UL[:, None], UR[:, None])
    U[:, np.isclose(x, s["x0"])] = (0.5 * (UL + UR))[:, None]   # a point on the diaphragm gets the average
    t, steps, frames, failed = 0.0, 0, [], None
    next_frame = 0.0
    while t < t_end - 1e-12:
        r, u, p = primitives(U)
        if not (np.all(np.isfinite(U)) and r.min() > 0 and p.min() > 0):
            failed = t
            break
        dt = min(cfl * dx / np.max(np.abs(u) + np.sqrt(G * p / r)), t_end - t)
        U = SCHEMES[name](U, dt / dx)
        U[:, 0], U[:, -1] = U[:, 1], U[:, -2]                        # transmissive boundaries
        t += dt
        steps += 1
        if n_frames and t >= next_frame - 1e-12:
            frames.append((t, np.array(primitives(U))))
            next_frame += t_end / n_frames
    return x, np.array(primitives(U)), steps, frames, failed


# --------------------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scheme", nargs="+", choices=list(SCHEMES), default=list(SCHEMES),
                   help="one or more of lf, lw, mc (default: all three)")
    p.add_argument("--nx", type=int, default=101, help="number of grid points (default 101)")
    p.add_argument("--cfl", type=float, default=0.8, help="CFL number (default 0.8)")
    p.add_argument("--t-end", type=float, default=0.2, help="final time (default 0.2)")
    p.add_argument("--animate", action="store_true", help="also save sod_fd.gif")
    args = p.parse_args()

    s = SOD
    exact = lambda x, t: exact_solution(x, t, s["rL"], s["uL"], s["pL"], s["rR"], s["uR"], s["pR"], s["x0"])
    xf = np.linspace(0, 1, 2001)
    results = {}
    for name in args.scheme:
        x, W, steps, frames, failed = solve(name, args.nx, args.cfl, args.t_end, n_frames=60 if args.animate else 0)
        if failed is not None:
            print(f"{NAMES[name]:25s} FAILED: negative pressure or density at t = {failed:.4f} (step {steps}). "
                  f"Try a smaller --cfl.")
            continue
        ex = exact(x, args.t_end)
        l1 = np.sum(np.abs(W[0] - ex[0])) * (x[1] - x[0])
        print(f"{NAMES[name]:25s} steps = {steps:4d}   L1 error in density = {l1:.5f}   max u = {W[1].max():.3f}")
        results[name] = (x, W, frames)
    if not results:
        return

    titles = ["Density ρ", "Velocity u", "Pressure p", "Internal energy e = p / ((γ−1)ρ)"]
    colors = {"lf": "tab:blue", "lw": "tab:orange", "mc": "tab:green"}
    to_plot = lambda W: [W[0], W[1], W[2], W[2] / ((G - 1) * W[0])]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    ex_f = to_plot(exact(xf, args.t_end))
    for k, ax in enumerate(axes.flat):
        ax.plot(xf, ex_f[k], "k-", lw=1.5, label="Exact")
        for name, (x, W, _) in results.items():
            ax.plot(x, to_plot(W)[k], "o-", ms=2.5, lw=1, color=colors[name], label=NAMES[name])
        ax.set(title=titles[k], xlabel="x")
        ax.grid(alpha=0.3)
    axes[0, 0].legend()
    fig.suptitle(f"Sod shock tube, finite differences, t = {args.t_end}, {args.nx} points, CFL = {args.cfl}")
    fig.tight_layout()
    fig.savefig("sod_fd.png", dpi=120)
    print("saved sod_fd.png")

    if args.animate:
        fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
        n = min(len(f) for (_, _, f) in results.values())
        lines = []
        for k, ax in enumerate(axes.flat):
            ex_line, = ax.plot([], [], "k-", lw=1.5, label="Exact")
            num = {name: ax.plot([], [], "o-", ms=2.5, lw=1, color=colors[name], label=NAMES[name])[0] for name in results}
            lines.append((ex_line, num))
            ax.set(title=titles[k], xlabel="x", xlim=(0, 1))
            ax.grid(alpha=0.3)
        for ax, lo, hi in zip(axes.flat, [0, -0.1, 0, 1.6], [1.1, 1.3, 1.1, 3.2]):
            ax.set_ylim(lo, hi)
        axes[0, 0].legend(loc="lower left")

        def update(i):
            t = next(iter(results.values()))[2][i][0]
            ex = to_plot(exact(xf, t))
            for k, (ex_line, num) in enumerate(lines):
                ex_line.set_data(xf, ex[k])
                for name, (x, _, frames) in results.items():
                    num[name].set_data(x, to_plot(frames[i][1])[k])
            fig.suptitle(f"Sod shock tube, finite differences, {args.nx} points, CFL = {args.cfl}, t = {t:.3f}")

        FuncAnimation(fig, update, frames=n).save("sod_fd.gif", writer=PillowWriter(fps=12))
        print("saved sod_fd.gif")


if __name__ == "__main__":
    main()
