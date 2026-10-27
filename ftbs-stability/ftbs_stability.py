"""Von Neumann stability of the FTBS scheme for u_t + a u_x = 0, shown with round-off growth.

    u_i^{n+1} = u_i^n - c (u_i^n - u_{i-1}^n),   c = a dt / dx
    amplification factor  G(theta) = 1 - c + c exp(-i theta),
    |G|^2 = 1 - 2 c (1 - c)(1 - cos theta)  ->  stable for 0 <= c <= 1

Two runs on a periodic grid start from the same smooth profile; one gets a tiny random
perturbation (default 1e-16, the size of double-precision round-off). Their difference delta
obeys the same linear scheme, so max|delta| grows or decays like |G|^n.

Examples:
    python ftbs_stability.py                      # c = 0.8, 1.0, 1.05, 1.1, 1.3
    python ftbs_stability.py --cfl 0.9 1.2        # your own CFL numbers
    python ftbs_stability.py --eps 1e-8 --nx 100  # larger perturbation, finer grid
"""
import argparse

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def run(c, nx, eps, n_max, seed=1):
    x = np.arange(nx) / nx
    u_clean = np.exp(-((x - 0.5) / 0.1) ** 2)
    rng = np.random.default_rng(seed)
    u_pert = u_clean + eps * (2 * rng.random(nx) - 1)
    hist = [np.max(np.abs(u_pert - u_clean))]
    for _ in range(n_max):
        u_clean = u_clean - c * (u_clean - np.roll(u_clean, 1))   # periodic FTBS
        u_pert = u_pert - c * (u_pert - np.roll(u_pert, 1))
        d = np.max(np.abs(u_pert - u_clean))
        hist.append(d)
        if not np.isfinite(d) or d > 1e6:
            break
    return np.array(hist)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cfl", type=float, nargs="+", default=[0.8, 1.0, 1.05, 1.1, 1.3], help="CFL numbers to compare")
    p.add_argument("--nx", type=int, default=50, help="number of grid points, periodic (default 50)")
    p.add_argument("--eps", type=float, default=1e-16, help="size of the initial perturbation (default 1e-16)")
    p.add_argument("--steps", type=int, default=500, help="maximum number of time steps (default 500)")
    args = p.parse_args()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1.6, 1]})
    theta = np.linspace(-np.pi, np.pi, 400)
    ax2.fill(np.cos(theta), np.sin(theta), color="tab:green", alpha=0.12)
    ax2.plot(np.cos(theta), np.sin(theta), "--", color="tab:green", label="unit circle |G| = 1")

    print(f"{'c':>6} {'|G(pi)|':>8} {'measured growth':>16} {'steps to |delta| = 1':>22}")
    for c in args.cfl:
        h = run(c, args.nx, args.eps, args.steps)
        n = np.arange(len(h))
        G = max(1.0, abs(1 - 2 * c))
        line, = ax1.semilogy(n, np.maximum(h, 1e-20), lw=2, label=f"c = {c}")
        ax1.semilogy(n, h[0] * G ** n, "--", lw=1, color=line.get_color())
        reached = np.argmax(h >= 1) if np.any(h >= 1) else None
        rate = (h[-1] / h[-11]) ** 0.1 if len(h) > 11 and h[-11] > 0 else float("nan")
        print(f"{c:6.2f} {abs(1 - 2 * c):8.3f} {rate:16.4f} {str(reached) if reached else 'never':>22}")
        Gz = 1 - c + c * np.exp(-1j * theta)
        ax2.plot(Gz.real, Gz.imag, lw=2, color=line.get_color(), label=f"c = {c}")
        ax2.plot(1 - 2 * c, 0, "o", color=line.get_color())

    ax1.axhline(1, color="gray", ls=":", lw=1)
    ax1.text(2, 2, "|δ| = 1: visible in the solution", color="gray", fontsize=9)
    ax1.set(xlabel="time step n", ylabel="max |δ|", ylim=(1e-20, 1e6),
            title=f"Growth of a {args.eps:g} perturbation (solid) and ε|G(π)|ⁿ (dashed), N = {args.nx}")
    ax1.grid(alpha=0.3, which="both")
    ax1.legend(loc="lower right")
    ax2.axhline(0, color="gray", lw=0.8)
    ax2.axvline(0, color="gray", lw=0.8)
    ax2.set(xlabel="Re G", ylabel="Im G", title="Amplification factor G(θ) = 1 − c + c e^(−iθ)", aspect="equal")
    ax2.legend(loc="upper left", fontsize=8)
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("ftbs_stability.png", dpi=120)
    print("saved ftbs_stability.png")


if __name__ == "__main__":
    main()
