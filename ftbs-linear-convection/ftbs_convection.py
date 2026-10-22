"""Solve the 1D linear convection equation u_t + a u_x = 0.

Schemes (c = a dt / dx is the CFL number):
    ftbs  FTBS / first-order upwind (numerical diffusion)
          u_i^{n+1} = u_i^n - c (u_i^n - u_{i-1}^n)
    lw    Lax-Wendroff, second order (numerical dispersion)
          u_i^{n+1} = u_i^n - (c/2)(u_{i+1}^n - u_{i-1}^n) + (c^2/2)(u_{i+1}^n - 2 u_i^n + u_{i-1}^n)
    tvd   flux-limited (TVD) scheme: sharp front, no wiggles
          u_i^{n+1} = u_i^n - c (u_i^n - u_{i-1}^n)
                      - (c(1-c)/2) [phi(r_i)(u_{i+1}^n - u_i^n) - phi(r_{i-1})(u_i^n - u_{i-1}^n)]
          r_i = (u_i^n - u_{i-1}^n) / (u_{i+1}^n - u_i^n),  phi = minmod, vanleer or superbee

Examples:
    python ftbs_convection.py                 # c = 0.5, step initial condition
    python ftbs_convection.py --cfl 1.0       # exact shift
    python ftbs_convection.py --cfl 1.1       # unstable
    python ftbs_convection.py --scheme lw     # Lax-Wendroff: wiggles behind the jump
    python ftbs_convection.py --scheme tvd --limiter superbee
    python ftbs_convection.py --animate       # also save an animated GIF
"""
import argparse

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


def initial_condition(x, shape):
    if shape == "step":      # shock-like jump: u = 1 left of x = 0.5, u = 2 right of it
        return np.where(x < 0.5, 1.0, 2.0)
    if shape == "square":    # u = 2 on [0.5, 1], u = 1 elsewhere
        return np.where((x >= 0.5) & (x <= 1.0), 2.0, 1.0)
    if shape == "gauss":
        return 1.0 + np.exp(-(((x - 0.75) / 0.12) ** 2))
    raise ValueError(f"unknown shape: {shape}")


LIMITERS = {
    "minmod": lambda r: np.maximum(0.0, np.minimum(1.0, r)),
    "vanleer": lambda r: (r + np.abs(r)) / (1.0 + np.abs(r)),
    "superbee": lambda r: np.maximum.reduce([np.zeros_like(r), np.minimum(2 * r, 1.0), np.minimum(r, 2.0)]),
}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--a", type=float, default=1.0, help="wave speed (default 1)")
    p.add_argument("--cfl", type=float, default=0.5, help="CFL number c = a dt/dx (default 0.5)")
    p.add_argument("--nx", type=int, default=81, help="number of grid points (default 81)")
    p.add_argument("--t-end", type=float, default=1.0, help="final time (default 1)")
    p.add_argument("--shape", choices=["step", "square", "gauss"], default="step")
    p.add_argument("--scheme", choices=["ftbs", "lw", "tvd"], default="ftbs",
                   help="ftbs (default), lw (Lax-Wendroff) or tvd (flux limiter)")
    p.add_argument("--limiter", choices=list(LIMITERS), default="vanleer", help="limiter for --scheme tvd (default vanleer)")
    p.add_argument("--animate", action="store_true", help="save ftbs_convection.gif")
    args = p.parse_args()

    a, c, L = args.a, args.cfl, 2.0
    x = np.linspace(0, L, args.nx)
    dx = x[1] - x[0]
    dt = c * dx / a
    nt = int(round(args.t_end / dt))

    u0 = initial_condition(x, args.shape)
    u = u0.copy()
    history = [u.copy()]
    for n in range(1, nt + 1):
        un = u.copy()
        if args.scheme == "ftbs":
            u[1:] = un[1:] - c * (un[1:] - un[:-1])
        elif args.scheme == "lw":
            u[1:-1] = (un[1:-1] - 0.5 * c * (un[2:] - un[:-2])
                       + 0.5 * c**2 * (un[2:] - 2 * un[1:-1] + un[:-2]))
            u[-1] = un[-1] - c * (un[-1] - un[-2])  # outflow boundary: upwind
        else:  # TVD: Lax-Wendroff correction scaled by the limiter phi(r)
            dR = un[2:] - un[1:-1]                   # u_{j+1} - u_j for j = 1..nx-2
            dL = un[1:-1] - un[:-2]                  # u_j - u_{j-1}
            safe = np.abs(dR) > 1e-14
            r = np.where(safe, dL / np.where(safe, dR, 1.0), 0.0)  # r does not matter where dR = 0
            D = np.zeros_like(un)                    # limited correction at face j+1/2, zero at both ends
            D[1:-1] = 0.5 * (1 - c) * LIMITERS[args.limiter](r) * dR
            u[1:] = un[1:] - c * ((un[1:] - un[:-1]) + D[1:] - D[:-1])
        u[0] = initial_condition(np.array([-a * n * dt]), args.shape)[0]  # inflow boundary
        history.append(u.copy())

    t = nt * dt
    u_exact = initial_condition(x - a * t, args.shape)
    print(f"dx = {dx:.4f}, dt = {dt:.4f}, c = {c}, steps = {nt}, t = {t:.3f}")
    print(f"max |u - u_exact| = {np.max(np.abs(u - u_exact)):.4g}")
    print(f"L1 error = dx * sum|u - u_exact| = {dx * np.sum(np.abs(u - u_exact)):.4g}")

    xf = np.linspace(0, L, 1000)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(xf, initial_condition(xf, args.shape), "k--", lw=1.2, label="Initial, t = 0")
    ax.plot(xf, initial_condition(xf - a * t, args.shape), "g-", lw=2, label=f"Exact, t = {t:.2f}")
    name = {"ftbs": "FTBS", "lw": "Lax-Wendroff", "tvd": f"TVD ({args.limiter})"}[args.scheme]
    ax.plot(x, u, "bo-", ms=3, label=f"{name}, t = {t:.2f}")
    ax.set(xlabel="x", ylabel="u", ylim=(0.5, 2.5), title=f"1D linear convection, {name} (a = {a}, c = {c})")
    ax.legend(); ax.grid(alpha=0.3)
    fig.savefig("ftbs_convection.png", dpi=120, bbox_inches="tight")
    print("saved ftbs_convection.png")

    if args.animate:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(xf, initial_condition(xf, args.shape), "k--", lw=1.2, label="Initial")
        exact_line, = ax.plot([], [], "g-", lw=2, label="Exact")
        num_line, = ax.plot([], [], "bo-", ms=3, label=name)
        ax.set(xlabel="x", ylabel="u", xlim=(0, L), ylim=(0.5, 2.5))
        ax.legend(loc="upper left"); ax.grid(alpha=0.3)
        frames = range(0, nt + 1, max(1, nt // 80))

        def update(k):
            tk = k * dt
            exact_line.set_data(xf, initial_condition(xf - a * tk, args.shape))
            num_line.set_data(x, history[k])
            ax.set_title(f"{name}, a = {a}, c = {c}, t = {tk:.3f}")
            return exact_line, num_line

        FuncAnimation(fig, update, frames=frames).save("ftbs_convection.gif", writer=PillowWriter(fps=15))
        print("saved ftbs_convection.gif")


if __name__ == "__main__":
    main()
