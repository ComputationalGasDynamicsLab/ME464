"""
1D unsteady heat equation: explicit FTCS vs implicit BTCS and Crank-Nicolson.

    dT/dt = alpha * d2T/dx2,   0 < x < L,  t > 0

    Boundary conditions at each end (--bc, values --v0 and --v1):
        Dirichlet  T = value            (fixed temperature)
        Neumann    dT/dx = value        (fixed gradient, i.e. heat flux; 0 = insulated)
    A Neumann end uses a ghost node from the central difference, e.g. T_{N+1} = T_{N-1} + 2 dx g.

    FTCS:  T_i^{n+1} = T_i^n + r (T_{i+1}^n - 2 T_i^n + T_{i-1}^n)               (stable for r <= 1/2)
    BTCS:  -r T_{i-1}^{n+1} + (1 + 2r) T_i^{n+1} - r T_{i+1}^{n+1} = T_i^n          (always stable)
    CN:    average of the two                                                     (always stable)

    r = alpha dt / dx^2 is the diffusion number.

Examples:
    python heat1d.py                              # FTCS, BTCS and CN at r = 0.4, N = 50, t = 10
    python heat1d.py --r 0.6 --scheme ftcs        # FTCS unstable: round-off grows into a sawtooth
    python heat1d.py --r 5 --ic pulse --t-end 1   # implicit only; Crank-Nicolson oscillates
    python heat1d.py --bc nn --ic hat             # insulated rod: relaxes to its mean
    python heat1d.py --bc dn --v1 1 --ic zero     # heat flux in at x = 1
    python heat1d.py --bc nn --v0 0 --v1 1 --ic zero   # net heating: no steady state
    python heat1d.py --animate                    # also save an animated GIF
    python heat1d.py --study 10 20 40 80 160 --r 0.25 --t-end 1     # grid convergence, fixed r (dt ~ dx^2)
    python heat1d.py --study 10 20 40 80 160 --lam 1 --t-end 1      # grid convergence, fixed dt/dx (dt ~ dx)
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L, ALPHA, NSERIES, PI = 1.0, 0.01, 8000, np.pi
NAMES = {"ftcs": "FTCS (explicit)", "btcs": "BTCS (implicit)", "cn": "Crank–Nicolson"}
COLORS = {"ftcs": "#2552c4", "btcs": "#c2410c", "cn": "#7c3aed"}
BCS = {"dd": "Dirichlet at both ends", "dn": "Dirichlet at x = 0, Neumann at x = 1",
       "nd": "Neumann at x = 0, Dirichlet at x = 1", "nn": "Neumann at both ends"}

# Initial conditions as pieces on [a, b]: ("poly", a, b, [c0, c1, c2]) = c0 + c1 x + c2 x^2, or ("sin"/"cos", A, nu) = A sin(nu x)
ICS = {
    "sine":     [("sin", 1, PI)],
    "zero":     [],
    "modes":    [("sin", 1, PI), ("sin", 0.4, 6 * PI)],
    "cosine":   [("cos", 1, PI)],
    "uniform":  [("poly", 0, 1, [1, 0, 0])],
    "ramp":     [("poly", 0, 1, [0, 1, 0])],
    "hat":      [("poly", 0, 0.5, [0, 2, 0]), ("poly", 0.5, 1, [2, -2, 0])],
    "parabola": [("poly", 0, 1, [0, 4, -4])],
    "pulse":    [("poly", 0.4, 0.6, [1, 0, 0])],
    "step":     [("poly", 0, 0.5, [1, 0, 0])],
}


def initial(x, ic):
    """Initial condition; a piece counts half on an interior jump, so a node on the jump gets the average."""
    T = np.zeros_like(x)
    for p in ICS[ic]:
        if p[0] != "poly":
            T += p[1] * (np.sin(p[2] * x) if p[0] == "sin" else np.cos(p[2] * x))
            continue
        _, a, b, c = p
        v = c[0] + c[1] * x + c[2] * x * x
        w = np.where((x > a + 1e-9) & (x < b - 1e-9), 1.0, 0.0)
        w[np.isclose(x, a) & (a > 1e-9)] = 0.5
        w[np.isclose(x, b) & (b < L - 1e-9)] = 0.5
        w[np.isclose(x, a) & (a <= 1e-9)] = 1.0
        w[np.isclose(x, b) & (b >= L - 1e-9)] = 1.0
        T += w * v
    return T


def _Ic(lam, a, b):   # integral of cos(lam x) over [a, b]
    lam = np.asarray(lam, float)
    safe = np.where(np.abs(lam) < 1e-12, 1.0, lam)
    return np.where(np.abs(lam) < 1e-12, b - a, (np.sin(lam * b) - np.sin(lam * a)) / safe)


def _Is(lam, a, b):   # integral of sin(lam x) over [a, b]
    lam = np.asarray(lam, float)
    safe = np.where(np.abs(lam) < 1e-12, 1.0, lam)
    return np.where(np.abs(lam) < 1e-12, 0.0, (np.cos(lam * a) - np.cos(lam * b)) / safe)


def _integ(p, fn, m):
    """Closed-form integral of an initial-condition piece times fn(m x) over the piece, for an array of m."""
    if p[0] != "poly":
        kind, A, nu = p
        a, b, A = 0.0, L, A / 2
        if kind == "sin" and fn == "sin":
            return A * (_Ic(nu - m, a, b) - _Ic(nu + m, a, b))
        if kind == "sin" and fn == "cos":
            return A * (_Is(nu + m, a, b) + _Is(nu - m, a, b))
        if kind == "cos" and fn == "sin":
            return A * (_Is(m + nu, a, b) + _Is(m - nu, a, b))
        return A * (_Ic(nu - m, a, b) + _Ic(nu + m, a, b))
    _, a, b, (c0, c1, c2) = p
    out = np.zeros_like(m)
    z = m == 0
    if fn == "cos":
        F0 = lambda X: c0 * X + c1 * X**2 / 2 + c2 * X**3 / 3
        out[z] = F0(b) - F0(a)
    mm = m[~z]
    if fn == "sin":
        F = lambda X: (c0 * (-np.cos(mm * X) / mm) + c1 * (-X * np.cos(mm * X) / mm + np.sin(mm * X) / mm**2)
                       + c2 * (-X**2 * np.cos(mm * X) / mm + 2 * X * np.sin(mm * X) / mm**2 + 2 * np.cos(mm * X) / mm**3))
    else:
        F = lambda X: (c0 * (np.sin(mm * X) / mm) + c1 * (X * np.sin(mm * X) / mm + np.cos(mm * X) / mm**2)
                       + c2 * (X**2 * np.sin(mm * X) / mm + 2 * X * np.cos(mm * X) / mm**2 - 2 * np.sin(mm * X) / mm**3))
    out[~z] = F(b) - F(a)
    return out


class Exact:
    """T = p(x, t) + sum_k c_k exp(-alpha mu_k^2 t) phi_k(x).

    p = q0 + q1 x + q2 x^2 + qt t satisfies the boundary conditions; phi_k = sin or cos(mu_k x) with
    dd: sin, k pi;  dn: sin, (k - 1/2) pi;  nd: cos, (k - 1/2) pi;  nn: cos, k pi (k >= 0)."""

    def __init__(self, ic, bc, v0, v1):
        self.ic, self.bc, self.v0, self.v1 = ic, bc, v0, v1
        if bc == "dd":
            self.q = [v0, v1 - v0, 0, 0]
        elif bc == "dn":
            self.q = [v0, v1, 0, 0]
        elif bc == "nd":
            self.q = [v1 - v0, v0, 0, 0]
        else:
            B = (v1 - v0) / 2
            self.q = [0, v0, B, 2 * ALPHA * B]
        self.fn = "sin" if bc in ("dd", "dn") else "cos"
        mu0 = PI if bc == "dd" else 0 if bc == "nn" else PI / 2
        self.mu = mu0 + PI * np.arange(NSERIES)
        norm = np.where(self.mu == 0, 1.0, 0.5)
        c = -_integ(("poly", 0, 1, self.q[:3]), self.fn, self.mu)
        for p in ICS[ic]:
            c = c + _integ(p, self.fn, self.mu)
        self.c = c / norm

    def __call__(self, x, t):
        if t == 0:
            T = initial(x, self.ic)
            if self.bc[0] == "d":
                T[np.isclose(x, 0)] = self.v0
            if self.bc[1] == "d":
                T[np.isclose(x, L)] = self.v1
            return T
        q = self.q
        ck = self.c * np.exp(-ALPHA * self.mu**2 * t)
        keep = np.abs(ck) > 1e-17
        f = np.sin if self.fn == "sin" else np.cos
        return q[0] + q[1] * x + q[2] * x**2 + q[3] * t + f(np.outer(x, self.mu[keep])) @ ck[keep]

    def steady(self, x):
        if self.q[3] != 0:
            return None                              # Neumann-Neumann with net heat input: no steady state
        return self.q[0] + self.q[1] * x + (self.c[0] if self.bc == "nn" else 0)


def thomas(a, b, c, d):
    """General tridiagonal solve: a[j] x[j-1] + b[j] x[j] + c[j] x[j+1] = d[j]."""
    m = len(d)
    cp, dp = np.empty(m), np.empty(m)
    cp[0], dp[0] = c[0] / b[0], d[0] / b[0]
    for j in range(1, m):
        den = b[j] - a[j] * cp[j - 1]
        cp[j], dp[j] = c[j] / den, (d[j] - a[j] * dp[j - 1]) / den
    out = np.empty(m)
    out[-1] = dp[-1]
    for j in range(m - 2, -1, -1):
        out[j] = dp[j] - cp[j] * out[j + 1]
    return out


def make_stepper(N, dx, r, bc, v0, v1):
    i0 = 1 if bc[0] == "d" else 0                      # unknown nodes i0..i1
    i1 = N - 1 if bc[1] == "d" else N

    def lap(T):                                      # T_{i-1} - 2 T_i + T_{i+1} at i0..i1, with ghost nodes at Neumann ends
        ext = np.concatenate(([T[1] - 2 * dx * v0], T, [T[N - 1] + 2 * dx * v1]))
        i = np.arange(i0, i1 + 1) + 1
        return ext[i - 1] - 2 * ext[i] + ext[i + 1]

    def step(T, scheme):
        Tn = T.copy()
        if scheme == "ftcs":
            T[i0:i1 + 1] = Tn[i0:i1 + 1] + r * lap(Tn)
            return T
        th = 1.0 if scheme == "btcs" else 0.5          # theta = 1 BTCS, 1/2 Crank-Nicolson
        idx = np.arange(i0, i1 + 1)
        a = np.where(idx == N, -2.0, -1.0) * th * r    # the ghost node doubles the inner neighbor at a Neumann end
        c = np.where(idx == 0, -2.0, -1.0) * th * r
        b = np.full(len(idx), 1 + 2 * th * r)
        d = Tn[i0:i1 + 1] + (1 - th) * r * lap(Tn)
        d[0] += th * r * Tn[0] if i0 == 1 else -th * r * 2 * dx * v0
        d[-1] += th * r * Tn[N] if i1 == N - 1 else th * r * 2 * dx * v1
        a[0] = c[-1] = 0.0
        T[i0:i1 + 1] = thomas(a, b, c, d)
        return T

    return step


def study(Ns, a):
    """Grid convergence: run every N to the same final time and print the error and observed order of each scheme."""
    exact = Exact(a.ic, a.bc, a.v0, a.v1)
    mode = f"fixed dt/dx = {a.lam:g}" if a.lam else f"fixed r = {a.r:g}"
    print(f"Grid convergence, {mode}, t = {a.t_end:g}, {BCS[a.bc]}, initial condition {a.ic}")
    print("dt is rounded down slightly so every run lands exactly on t; r is the value actually used.\n")
    hdr = f"{'N':>5} {'dx':>9} {'dt':>10} {'r':>7} {'steps':>7}" + "".join(f" {s.upper() + ' max err':>13} {'order':>6}" for s in a.scheme)
    print(hdr)
    prev = None
    for N in sorted(set(Ns)):
        dx = L / N
        dt0 = a.lam * dx if a.lam else a.r * dx**2 / ALPHA
        nt = max(1, int(np.ceil(a.t_end / dt0 - 1e-9)))
        dt = a.t_end / nt
        r = ALPHA * dt / dx**2
        x = np.linspace(0, L, N + 1)
        step = make_stepper(N, dx, r, a.bc, a.v0, a.v1)
        T0, Te = exact(x, 0), exact(x, a.t_end)
        errs = {}
        for s in a.scheme:
            T = T0.copy()
            with np.errstate(all="ignore"):
                for _ in range(nt):
                    step(T, s)
                    if not np.all(np.isfinite(T)) or np.abs(T).max() > 1e6:
                        break
            ok = np.all(np.isfinite(T)) and np.abs(T).max() <= 1e6
            errs[s] = np.abs(T - Te).max() if ok else None
        line = f"{N:5d} {dx:9.5f} {dt:10.3e} {r:7.3g} {nt:7d}"
        for s in a.scheme:
            e = errs[s]
            o = (np.log(prev[1][s] / e) / np.log(prev[0] / dx)
                 if prev and e and prev[1][s] else None)
            line += f" {('unstable' if e is None else f'{e:.4e}'):>13} {('' if o is None else f'{o:.3f}'):>6}"
        print(line)
        prev = (dx, errs)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scheme", nargs="+", choices=list(NAMES), default=list(NAMES), help="schemes to run")
    p.add_argument("--r", type=float, default=0.4, help="diffusion number r = alpha dt / dx^2 (default 0.4)")
    p.add_argument("--N", type=int, default=50, help="number of intervals, dx = L/N (default 50)")
    p.add_argument("--t-end", type=float, default=10.0, help="final time, at most 10 (default 10)")
    p.add_argument("--ic", choices=list(ICS), default="sine", help="initial condition (default sine)")
    p.add_argument("--bc", choices=list(BCS), default="dd",
                   help="d = Dirichlet (fixed T), n = Neumann (fixed dT/dx), at x = 0 then x = 1 (default dd)")
    p.add_argument("--v0", type=float, default=0.0, help="boundary value at x = 0: T or dT/dx (default 0)")
    p.add_argument("--v1", type=float, default=0.0, help="boundary value at x = 1: T or dT/dx (default 0)")
    p.add_argument("--animate", action="store_true", help="also save heat1d.gif")
    p.add_argument("--study", type=int, nargs="+", metavar="N", help="grid convergence study over these N values")
    p.add_argument("--lam", type=float, default=None, help="with --study: fix dt/dx = lam instead of r (dt ~ dx)")
    a = p.parse_args()
    if not 0 < a.t_end <= 10:
        p.error("--t-end must be between 0 and 10")
    if a.study:
        if len(set(a.study)) < 2 or min(a.study) < 4:
            p.error("--study needs at least two different N values, each 4 or more")
        study(a.study, a)
        return

    dx = L / a.N
    dt = a.r * dx**2 / ALPHA
    nt = max(1, int(np.floor(a.t_end / dt + 1e-9)))
    x = np.linspace(0, L, a.N + 1)
    exact = Exact(a.ic, a.bc, a.v0, a.v1)
    step = make_stepper(a.N, dx, a.r, a.bc, a.v0, a.v1)
    T0 = exact(x, 0)
    sol = {s: T0.copy() for s in a.scheme}
    blown = {s: None for s in a.scheme}
    frames, every = [], max(1, nt // 120)

    for n in range(1, nt + 1):
        for s in a.scheme:
            if blown[s] is None:
                step(sol[s], s)
                if not np.all(np.isfinite(sol[s])) or np.abs(sol[s]).max() > 1e6:
                    blown[s] = n
        if a.animate and (n % every == 0 or n == nt):
            frames.append((n * dt, {s: v.copy() for s, v in sol.items()}))
    t = nt * dt

    print(f"{BCS[a.bc]}, v0 = {a.v0:g}, v1 = {a.v1:g}, initial condition {a.ic}")
    print(f"r = {a.r:g}, N = {a.N}, dx = {dx:g}, dt = {dt:g}, steps = {nt}, t = {t:g}")
    Te = exact(x, t)
    for s in a.scheme:
        if blown[s]:
            print(f"  {NAMES[s]:18s} blew up at step {blown[s]} (t = {blown[s] * dt:g})")
        else:
            print(f"  {NAMES[s]:18s} max |T - T_exact| = {np.abs(sol[s] - Te).max():.3e}")

    xf = np.linspace(0, L, 801)
    ys = np.concatenate([initial(xf, a.ic), exact(xf, t), [0.0]])
    lo, hi = ys.min(), max(ys.max(), ys.min() + 0.5)
    pad = 0.12 * (hi - lo)

    def draw(ax, t, solns):
        ax.plot(xf, initial(xf, a.ic), "--", color="#8a97a4", lw=1.2, label="Initial, t = 0")
        ss = exact.steady(xf)
        if ss is not None:
            ax.plot(xf, ss, ":", color="#5a6876", lw=1.5, label="Steady state")
        ax.plot(xf, exact(xf, t), color="#17795a", lw=2.5, label=f"Exact, t = {t:.3g}")
        for s, v in solns.items():
            ax.plot(x, np.clip(v, lo - 5, hi + 5), "o-", ms=3, lw=1.2, color=COLORS[s], label=NAMES[s])
        ax.set_xlim(0, L)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlabel("x"); ax.set_ylabel("T"); ax.grid(alpha=0.3); ax.legend(loc="best", fontsize=9)
        ax.set_title(f"1D heat equation, {a.bc.upper()} BC, r = {a.r:g}, N = {a.N}, t = {t:.3g}")

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
