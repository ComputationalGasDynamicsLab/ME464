"""Sod shock tube: the 1D Euler equations solved with Roe, HLLC and AUSM+ fluxes.

    U_t + F(U)_x = 0,  U = (rho, rho u, E),  F = (rho u, rho u^2 + p, u (E + p)),
    p = (gamma - 1)(E - rho u^2 / 2),  gamma = 1.4

Initial condition (Sod 1978) on 0 <= x <= 1 with the diaphragm at x = 0.5:
    left:  rho = 1,     u = 0, p = 1
    right: rho = 0.125, u = 0, p = 0.1
The numerical solutions are compared with the exact Riemann solution at t = 0.2.

Finite-volume method: first order, or second order with MUSCL reconstruction of the
primitive variables (with a slope limiter) and two-stage SSP Runge-Kutta time stepping.

Examples:
    python sod_shock_tube.py                         # all three fluxes, 1st order, 100 cells, CFL 0.5
    python sod_shock_tube.py --order 2               # MUSCL + van Leer limiter
    python sod_shock_tube.py --flux ausm+ --nx 400   # one flux, finer grid
    python sod_shock_tube.py --order 2 --animate     # also save an animated GIF
    python sod_shock_tube.py --cfl 0.8               # AUSM+ fails in the first steps (see README)
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
# Numerical fluxes. Inputs are left/right primitive states (arrays over all interfaces).
# --------------------------------------------------------------------------------------
def energy(r, u, p):
    return p / (G - 1) + 0.5 * r * u * u


def physical_flux(r, u, p):
    E = energy(r, u, p)
    return np.array([r * u, r * u * u + p, u * (E + p)])


def flux_roe(rL, uL, pL, rR, uR, pR):
    """Roe's approximate Riemann solver with Harten's entropy fix."""
    HL, HR = (energy(rL, uL, pL) + pL) / rL, (energy(rR, uR, pR) + pR) / rR
    sL, sR = np.sqrt(rL), np.sqrt(rR)
    u = (sL * uL + sR * uR) / (sL + sR)          # Roe averages
    H = (sL * HL + sR * HR) / (sL + sR)
    a = np.sqrt((G - 1) * (H - 0.5 * u * u))
    r = sL * sR
    dr, du, dp = rR - rL, uR - uL, pR - pL
    alpha = [(dp - r * a * du) / (2 * a * a), dr - dp / (a * a), (dp + r * a * du) / (2 * a * a)]
    lam = [np.abs(u - a), np.abs(u), np.abs(u + a)]
    delta = 0.1 * (np.abs(u) + a)                # Harten's entropy fix on the acoustic waves
    for k in (0, 2):
        lam[k] = np.where(lam[k] < delta, (lam[k] ** 2 + delta ** 2) / (2 * delta), lam[k])
    one = np.ones_like(u)
    K = [np.array([one, u - a, H - u * a]), np.array([one, u, 0.5 * u * u]), np.array([one, u + a, H + u * a])]
    F = 0.5 * (physical_flux(rL, uL, pL) + physical_flux(rR, uR, pR))
    for k in range(3):
        F -= 0.5 * lam[k] * alpha[k] * K[k]
    return F


def flux_hllc(rL, uL, pL, rR, uR, pR):
    """HLLC (Toro, Spruce and Speares 1994): HLL with the contact wave restored."""
    aL, aR = np.sqrt(G * pL / rL), np.sqrt(G * pR / rR)
    SL = np.minimum(uL - aL, uR - aR)            # Davis wave-speed estimates
    SR = np.maximum(uL + aL, uR + aR)
    Ss = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / (rL * (SL - uL) - rR * (SR - uR))
    EL, ER = energy(rL, uL, pL), energy(rR, uR, pR)
    FL, FR = physical_flux(rL, uL, pL), physical_flux(rR, uR, pR)
    UL, UR = np.array([rL, rL * uL, EL]), np.array([rR, rR * uR, ER])

    def u_star(r, u, p, E, S):
        c = r * (S - u) / (S - Ss)
        return np.array([c, c * Ss, c * (E / r + (Ss - u) * (Ss + p / (r * (S - u))))])

    FsL = FL + SL * (u_star(rL, uL, pL, EL, SL) - UL)
    FsR = FR + SR * (u_star(rR, uR, pR, ER, SR) - UR)
    return np.where(SL >= 0, FL, np.where(Ss >= 0, FsL, np.where(SR > 0, FsR, FR)))


def flux_ausm_plus(rL, uL, pL, rR, uR, pR, alpha=3 / 16, beta=1 / 8):
    """AUSM+ (Liou 1996): separate convective (mass flux) and pressure parts."""
    HL, HR = (energy(rL, uL, pL) + pL) / rL, (energy(rR, uR, pR) + pR) / rR
    acL = np.sqrt(2 * (G - 1) / (G + 1) * HL)    # critical speeds of sound
    acR = np.sqrt(2 * (G - 1) / (G + 1) * HR)
    a12 = np.minimum(acL ** 2 / np.maximum(acL, uL), acR ** 2 / np.maximum(acR, -uR))
    ML, MR = uL / a12, uR / a12

    def M4(M, s):   # split Mach number, s = +1 or -1
        return np.where(np.abs(M) >= 1, 0.5 * (M + s * np.abs(M)),
                        s * 0.25 * (M + s) ** 2 + s * beta * (M * M - 1) ** 2)

    def P5(M, s):   # split pressure
        return np.where(np.abs(M) >= 1, 0.5 * (1 + s * np.sign(M)),
                        0.25 * (M + s) ** 2 * (2 - s * M) + s * alpha * M * (M * M - 1) ** 2)

    m12 = M4(ML, 1) + M4(MR, -1)
    p12 = P5(ML, 1) * pL + P5(MR, -1) * pR
    mdot = a12 * np.where(m12 > 0, m12 * rL, m12 * rR)
    up = mdot >= 0
    return np.array([mdot,
                     mdot * np.where(up, uL, uR) + p12,
                     mdot * np.where(up, HL, HR)])


FLUXES = {"roe": flux_roe, "hllc": flux_hllc, "ausm+": flux_ausm_plus}
NAMES = {"roe": "Roe", "hllc": "HLLC", "ausm+": "AUSM+"}


# --------------------------------------------------------------------------------------
# Finite-volume solver
# --------------------------------------------------------------------------------------
def slope(dm, dp, limiter):
    if limiter == "minmod":
        return np.where(dm * dp > 0, np.sign(dm) * np.minimum(np.abs(dm), np.abs(dp)), 0.0)
    if limiter == "vanleer":
        return np.where(dm * dp > 0, 2 * dm * dp / (dm + dp + 1e-300), 0.0)
    if limiter == "superbee":
        mm = lambda a, b: np.where(a * b > 0, np.sign(a) * np.minimum(np.abs(a), np.abs(b)), 0.0)
        s1, s2 = mm(2 * dm, dp), mm(dm, 2 * dp)
        return np.where(np.abs(s1) > np.abs(s2), s1, s2)
    raise ValueError(limiter)


def primitives(U):
    r = U[0]
    u = U[1] / r
    return r, u, (G - 1) * (U[2] - 0.5 * r * u * u)


def residual(U, dx, flux, order, limiter):
    """-(F_{i+1/2} - F_{i-1/2}) / dx for the interior cells. U has 2 ghost cells per side."""
    U = U.copy()
    U[:, :2] = U[:, 2:3]                         # transmissive boundaries
    U[:, -2:] = U[:, -3:-2]
    W = np.array(primitives(U))
    S = np.zeros_like(W)
    if order == 2:
        S[:, 1:-1] = slope(W[:, 1:-1] - W[:, :-2], W[:, 2:] - W[:, 1:-1], limiter)
    WL = W[:, 1:-2] + 0.5 * S[:, 1:-2]           # left state at interfaces 1+1/2 .. N+1+1/2
    WR = W[:, 2:-1] - 0.5 * S[:, 2:-1]
    F = flux(*WL, *WR)
    R = np.zeros_like(U)
    R[:, 2:-2] = -(F[:, 1:] - F[:, :-1]) / dx
    return R


def solve(flux_name, nx, cfl, t_end, order, limiter, n_frames=0):
    dx = 1.0 / nx
    xc = (np.arange(nx) + 0.5) * dx
    s = SOD
    r = np.where(xc < s["x0"], s["rL"], s["rR"])
    u = np.where(xc < s["x0"], s["uL"], s["uR"])
    p = np.where(xc < s["x0"], s["pL"], s["pR"])
    U = np.zeros((3, nx + 4))
    U[:, 2:-2] = [r, r * u, energy(r, u, p)]
    flux = FLUXES[flux_name]
    t, steps, frames = 0.0, 0, []
    next_frame = 0.0
    failed = None
    while t < t_end - 1e-12:
        rr, uu, pp = primitives(U[:, 2:-2])
        if not (np.all(np.isfinite(U)) and rr.min() > 0 and pp.min() > 0):
            failed = t                           # negative density or pressure: the run cannot continue
            break
        dt = min(cfl * dx / np.max(np.abs(uu) + np.sqrt(G * pp / rr)), t_end - t)
        if order == 1:
            U = U + dt * residual(U, dx, flux, 1, limiter)
        else:  # SSP-RK2 (Heun)
            U1 = U + dt * residual(U, dx, flux, 2, limiter)
            U = 0.5 * (U + U1 + dt * residual(U1, dx, flux, 2, limiter))
        t += dt
        steps += 1
        if n_frames and t >= next_frame - 1e-12:
            frames.append((t, np.array(primitives(U[:, 2:-2]))))
            next_frame += t_end / n_frames
    return xc, np.array(primitives(U[:, 2:-2])), steps, frames, failed


# --------------------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--flux", nargs="+", choices=list(FLUXES), default=list(FLUXES),
                   help="one or more of roe, hllc, ausm+ (default: all three)")
    p.add_argument("--nx", type=int, default=100, help="number of cells (default 100)")
    p.add_argument("--cfl", type=float, default=0.5, help="CFL number (default 0.5)")
    p.add_argument("--order", type=int, choices=[1, 2], default=1, help="1 (default) or 2 (MUSCL + RK2)")
    p.add_argument("--limiter", choices=["minmod", "vanleer", "superbee"], default="vanleer",
                   help="slope limiter for --order 2 (default vanleer)")
    p.add_argument("--t-end", type=float, default=0.2, help="final time (default 0.2)")
    p.add_argument("--animate", action="store_true", help="also save sod_shock_tube.gif")
    args = p.parse_args()

    s = SOD
    ps, us = star_state(s["rL"], s["uL"], s["pL"], s["rR"], s["uR"], s["pR"])
    print(f"exact star state: p* = {ps:.5f}, u* = {us:.5f}")

    xf = np.linspace(0, 1, 2001)
    exact = lambda x, t: exact_solution(x, t, s["rL"], s["uL"], s["pL"], s["rR"], s["uR"], s["pR"], s["x0"])
    ex_f = exact(xf, args.t_end)
    results = {}
    for name in args.flux:
        xc, W, steps, frames, failed = solve(name, args.nx, args.cfl, args.t_end, args.order, args.limiter,
                                             n_frames=60 if args.animate else 0)
        if failed is not None:
            print(f"{NAMES[name]:6s} FAILED: negative pressure or density at t = {failed:.4f} "
                  f"(step {steps}). Try a smaller --cfl.")
            continue
        ex_c = exact(xc, args.t_end)
        l1 = np.sum(np.abs(W[0] - ex_c[0])) / args.nx
        print(f"{NAMES[name]:6s} steps = {steps:4d}   L1 error in density = {l1:.5f}")
        results[name] = (xc, W, frames)

    if not results:
        return
    titles = ["Density ρ", "Velocity u", "Pressure p", "Internal energy e = p / ((γ−1)ρ)"]
    colors = {"roe": "tab:blue", "hllc": "tab:orange", "ausm+": "tab:green"}
    to_plot = lambda W: [W[0], W[1], W[2], W[2] / ((G - 1) * W[0])]
    order_txt = "1st order" if args.order == 1 else f"2nd order MUSCL ({args.limiter})"

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    for k, ax in enumerate(axes.flat):
        ax.plot(xf, to_plot(ex_f)[k], "k-", lw=1.5, label="Exact")
        for name, (xc, W, _) in results.items():
            ax.plot(xc, to_plot(W)[k], "o-", ms=2.5, lw=1, color=colors[name], label=NAMES[name])
        ax.set(title=titles[k], xlabel="x")
        ax.grid(alpha=0.3)
    axes[0, 0].legend()
    fig.suptitle(f"Sod shock tube, t = {args.t_end}, {args.nx} cells, CFL = {args.cfl}, {order_txt}")
    fig.tight_layout()
    fig.savefig("sod_shock_tube.png", dpi=120)
    print("saved sod_shock_tube.png")

    if args.animate:
        fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
        n = min(len(f) for (_, _, f) in results.values())
        lines = []
        for k, ax in enumerate(axes.flat):
            ex_line, = ax.plot([], [], "k-", lw=1.5, label="Exact")
            num = {name: ax.plot([], [], "o-", ms=2.5, lw=1, color=colors[name], label=NAMES[name])[0]
                   for name in results}
            lines.append((ex_line, num))
            ax.set(title=titles[k], xlabel="x", xlim=(0, 1))
            ax.grid(alpha=0.3)
        for ax, lo, hi in zip(axes.flat, [0, -0.1, 0, 1.6], [1.1, 1.1, 1.1, 3.0]):
            ax.set_ylim(lo, hi)
        axes[0, 0].legend(loc="lower left")

        def update(i):
            t = next(iter(results.values()))[2][i][0]
            ex = to_plot(exact(xf, t))
            for k, (ex_line, num) in enumerate(lines):
                ex_line.set_data(xf, ex[k])
                for name, (xc, _, frames) in results.items():
                    num[name].set_data(xc, to_plot(frames[i][1])[k])
            fig.suptitle(f"Sod shock tube, {order_txt}, {args.nx} cells, t = {t:.3f}")

        FuncAnimation(fig, update, frames=n).save("sod_shock_tube.gif", writer=PillowWriter(fps=12))
        print("saved sod_shock_tube.gif")


if __name__ == "__main__":
    main()
