"""
1D unsteady heat equation: explicit FTCS vs implicit BTCS and Crank-Nicolson.

    dT/dt = alpha * d2T/dx2,   0 < x < L,  t > 0,   T(0,t) = T(L,t) = 0

    FTCS:  T_i^{n+1} = T_i^n + r (T_{i+1}^n - 2 T_i^n + T_{i-1}^n)               (stable for r <= 1/2)
    BTCS:  -r T_{i-1}^{n+1} + (1 + 2r) T_i^{n+1} - r T_{i+1}^{n+1} = T_i^n          (always stable)
    CN:    average of the two                                                     (always stable)

    r = alpha dt / dx^2 is the diffusion number.

Examples:
    python heat1d.py                          # FTCS, BTCS and CN at r = 0.4, N = 50, t = 10
    python heat1d.py --r 0.6 --scheme ftcs    # FTCS unstable: round-off grows into a sawtooth
    python heat1d.py --r 5 --ic pulse         # implicit only; Crank-Nicolson oscillates
    python heat1d.py --animate                # also save an animated GIF
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L, ALPHA = 1.0, 0.01
NAMES = {"ftcs": "FTCS (explicit)", "btcs": "BTCS (implicit)", "cn": "Crank–Nicolson"}
COLORS = {"ftcs": "#2552c4", "btcs": "#c2410c", "cn": "#7c3aed"}
PULSE = (0.4, 0.6)        # square pulse T = 1 on 0.4 <= x <= 0.6


def initial(x, ic):
    if ic == "sine":
        return np.sin(np.pi * x / L)
    if ic == "modes":
        return np.sin(np.pi * x / L) + 0.4 * np.sin(6 * np.pi * x / L)
    a, b = PULSE                                     # pulse: half value exactly on a jump
    T = np.where((x > a) & (x < b), 1.0, 0.0)
    T[np.isclose(x, a) | np.isclose(x, b)] = 0.5
    return T


def exact(x, t, ic):
    decay = lambda k: np.exp(-ALPHA * (k * np.pi / L) ** 2 * t)
    if ic == "sine":
        return np.sin(np.pi * x / L) * decay(1)
    if ic == "modes":
        return np.sin(np.pi * x / L) * decay(1) + 0.4 * np.sin(6 * np.pi * x / L) * decay(6)
    if t == 0:
        return initial(x, ic)
    a, b = PULSE                                     # Fourier sine series of the pulse
    k = np.arange(1, 8001)
    bk = 2 / (k * np.pi) * (np.cos(k * np.pi * a / L) - np.cos(k * np.pi * b / L)) * decay(k)
    return np.sin(np.outer(x, k) * np.pi / L) @ bk


def thomas(sub, diag, rhs):
    """Solve the constant-coefficient tridiagonal system sub*x[i-1] + diag*x[i] + sub*x[i+1] = rhs[i]."""
    m = len(rhs)
    cp, dp = np.empty(m), np.empty(m)
    cp[0], dp[0] = sub / diag, rhs[0] / diag
    for i in range(1, m):
        den = diag - sub * cp[i - 1]
        cp[i], dp[i] = sub / den, (rhs[i] - sub * dp[i - 1]) / den
    out = np.empty(m)
    out[-1] = dp[-1]
    for i in range(m - 2, -1, -1):
        out[i] = dp[i] - cp[i] * out[i + 1]
    return out


def step(T, scheme, r):
    Tn = T.copy()
    lap = Tn[2:] - 2 * Tn[1:-1] + Tn[:-2]
    if scheme == "ftcs":
        T[1:-1] = Tn[1:-1] + r * lap
    elif scheme == "btcs":
        T[1:-1] = thomas(-r, 1 + 2 * r, Tn[1:-1])
    else:
        T[1:-1] = thomas(-r / 2, 1 + r, Tn[1:-1] + 0.5 * r * lap)
    T[0] = T[-1] = 0.0                               # Dirichlet boundaries
    return T


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scheme", nargs="+", choices=list(NAMES), default=list(NAMES), help="schemes to run")
    p.add_argument("--r", type=float, default=0.4, help="diffusion number r = alpha dt / dx^2 (default 0.4)")
    p.add_argument("--N", type=int, default=50, help="number of intervals, dx = L/N (default 50)")
    p.add_argument("--t-end", type=float, default=10.0, help="final time, at most 10 (default 10)")
    p.add_argument("--ic", choices=["sine", "modes", "pulse"], default="sine",
                   help="sine: sin(pi x); modes: sin(pi x) + 0.4 sin(6 pi x); pulse: T = 1 on 0.4..0.6")
    p.add_argument("--animate", action="store_true", help="also save heat1d.gif")
    a = p.parse_args()
    if not 0 < a.t_end <= 10:
        p.error("--t-end must be between 0 and 10")

    dx = L / a.N
    dt = a.r * dx**2 / ALPHA
    nt = int(round(a.t_end / dt)) if abs(a.t_end / dt - round(a.t_end / dt)) < 1e-9 else int(a.t_end / dt)
    x = np.linspace(0, L, a.N + 1)
    T0 = initial(x, a.ic)
    sol = {s: T0.copy() for s in a.scheme}
    blown = {s: None for s in a.scheme}
    frames, every = [], max(1, nt // 120)

    for n in range(1, nt + 1):
        for s in a.scheme:
            if blown[s] is None:
                step(sol[s], s, a.r)
                if not np.all(np.isfinite(sol[s])) or np.abs(sol[s]).max() > 1e6:
                    blown[s] = n
        if a.animate and (n % every == 0 or n == nt):
            frames.append((n * dt, {s: v.copy() for s, v in sol.items()}))
    t = nt * dt

    print(f"r = {a.r:g}, N = {a.N}, dx = {dx:g}, dt = {dt:g}, steps = {nt}, t = {t:g}")
    Te = exact(x, t, a.ic)
    for s in a.scheme:
        if blown[s]:
            print(f"  {NAMES[s]:18s} blew up at step {blown[s]} (t = {blown[s] * dt:g})")
        else:
            print(f"  {NAMES[s]:18s} max |T - T_exact| = {np.abs(sol[s] - Te).max():.3e}")

    def draw(ax, t, solns):
        xf = np.linspace(0, L, 801)
        ax.plot(xf, initial(xf, a.ic) if a.ic != "pulse" else exact(xf, 0, a.ic), "--", color="#8a97a4", lw=1.2, label="Initial, t = 0")
        ax.plot(xf, exact(xf, t, a.ic), color="#17795a", lw=2.5, label=f"Exact, t = {t:.3g}")
        for s, v in solns.items():
            ax.plot(x, np.clip(v, -5, 5), "o-", ms=3, lw=1.2, color=COLORS[s], label=NAMES[s])
        ax.set_xlim(0, L)
        ax.set_ylim(-0.5 if a.ic == "modes" else -0.15, 1.55 if a.ic == "modes" else 1.15)
        ax.set_xlabel("x"); ax.set_ylabel("T"); ax.grid(alpha=0.3); ax.legend(loc="upper right", fontsize=9)
        ax.set_title(f"1D heat equation, r = {a.r:g}, N = {a.N}, t = {t:.3g}")

    fig, ax = plt.subplots(figsize=(8, 4.8))
    draw(ax, t, sol)
    fig.tight_layout(); fig.savefig("heat1d.png", dpi=130)
    print("Saved heat1d.png")

    if a.animate:
        from matplotlib.animation import FuncAnimation, PillowWriter
        fig, ax = plt.subplots(figsize=(8, 4.8))
        def upd(k):
            ax.clear(); draw(ax, *frames[k])
        anim = FuncAnimation(fig, upd, frames=len(frames))
        anim.save("heat1d.gif", writer=PillowWriter(fps=20))
        print("Saved heat1d.gif")


if __name__ == "__main__":
    main()
