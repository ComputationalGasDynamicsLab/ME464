# FTBS Linear Convection

Teaching material for a CFD course: the 1D linear convection equation solved with the
forward-in-time, backward-in-space (FTBS, first-order upwind) scheme.

![FTBS animation](ftbs_convection.gif)

## The equation

$$
\frac{\partial u}{\partial t} + a\,\frac{\partial u}{\partial x} = 0, \qquad a = 1
$$

The exact solution is the initial profile shifted to the right: $u(x,t) = u_0(x - at)$.

## The scheme

With the CFL (Courant) number

$$
c = \frac{a\,\Delta t}{\Delta x}
$$

the FTBS update is

$$
u_i^{n+1} = u_i^n - c\left(u_i^n - u_{i-1}^n\right) = (1-c)\,u_i^n + c\,u_{i-1}^n
$$

It is stable for $a > 0$ when $0 \le c \le 1$.

| CFL number | What you see |
|---|---|
| $c < 1$ | Stable, but the sharp front smears into a ramp (numerical diffusion). Smaller $c$ smears more. |
| $c = 1$ | $u_i^{n+1} = u_{i-1}^n$: an exact one-cell shift, zero error. |
| $c > 1$ | Unstable. The weight $1 - c$ is negative and oscillations grow until the solution blows up. |

The default initial condition is a shock-like step: $u = 1$ for $x < 0.5$ and $u = 2$ for $x \ge 0.5$,
on the domain $x \in [0, 2]$.

## Interactive animation

Open [`index.html`](index.html) in a browser (or enable GitHub Pages for this repository).
You can play the solution step by step, change the CFL number and grid size, switch the initial
condition, and see on a stencil diagram why the scheme needs $c \le 1$.

## Python script

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

python ftbs_convection.py                # c = 0.5, saves ftbs_convection.png
python ftbs_convection.py --cfl 1.0      # exact shift
python ftbs_convection.py --cfl 1.1      # unstable
python ftbs_convection.py --animate      # also saves ftbs_convection.gif
python ftbs_convection.py --help         # all options (--nx, --t-end, --shape, --a)
```
