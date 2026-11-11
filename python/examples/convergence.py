"""Convergence of the Möbius-transformed trapezoidal rule.

This reproduces the example of the Matlab code cot_trap.m, the integral of
|x|^p with respect to the standard Gaussian weight. A convergence plot is shown
if matplotlib is installed.

Run from the python directory with

    python examples/convergence.py
"""

import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import moebius_trapezoidal as mt  # noqa: E402

p = 1  # parameter of the integrand |x|^p
exact = math.sqrt(2**p / math.pi) * math.gamma((p + 1) / 2)

ns = np.array([2**m + 1 for m in range(2, 17)])
errors = []
for n in ns:
    value = mt.integrate(lambda x: np.abs(x) ** p, int(n), "gaussian")
    errors.append(abs(value - exact))
    print(f"n = {n:6d}   estimate = {value:.10f}   error = {errors[-1]:.3e}")
errors = np.array(errors)

rate = np.polyfit(np.log(ns[2:]), np.log(errors[2:]), 1)[0]
print(f"estimated rate of convergence: n^{rate:.2f}")

# automatic refinement until a relative tolerance is met
result = mt.quad(lambda x: np.abs(x) ** p, rtol=1e-8)
print(f"quad: value = {result.value:.10f}, error estimate = {result.error:.1e}, n = {result.n}")

try:
    import matplotlib.pyplot as plt
except ImportError:
    pass
else:
    plt.loglog(ns, errors, "o-")
    plt.title("Convergence graph")
    plt.xlabel("Number of points")
    plt.ylabel("Absolute error")
    plt.show()
