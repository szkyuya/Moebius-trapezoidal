"""Root-mean-squared error of the two randomized rules.

integrate_randomized uses n points with a random shift. integrate_rmt also
draws the number of points at random from {⌊n/2⌋, ..., n} (Definition 4.1 of
arXiv:2407.13650), which improves the worst-case RMSE by half an order.

The comparison uses two integrands with the Gaussian weight:

* |x|, whose transformed Fourier coefficients decay regularly. Both rules
  converge at the same rate n^-2, and integrate_rmt has a larger constant
  since it uses fewer points on average.
* cos(n θ(x)), where θ(x) is the inverse of x = -cot(θ/2), i.e. a frequency
  equal to n after the transformation. This is the kind of integrand that
  determines the worst case: the n point rule aliases it completely, so its
  error does not decrease, while integrate_rmt still converges.

Run from the python directory with

    python examples/rmse_comparison.py
"""

import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import moebius_trapezoidal as mt  # noqa: E402
import moebius_trapezoidal_randomized as mtr  # noqa: E402


def theta(x):
    """Inverse of the transformation x = -cot(θ/2)."""
    return 2.0 * np.arctan2(1.0, -x)


def rmse(estimates, exact):
    return math.sqrt(np.mean((np.asarray(estimates) - exact) ** 2))


repetitions = 200
rng = np.random.default_rng(2025)
ns = [32, 64, 128, 256, 512]

print(f"RMSE over {repetitions} repetitions, Gaussian weight\n")
print(f"{'n':>5} | {'|x|, shift':>11} {'|x|, RMT':>11} | {'cos(nθ), shift':>15} {'cos(nθ), RMT':>13}")
rows = []
for n in ns:
    exact_abs = math.sqrt(2 / math.pi)
    f_osc = lambda x, n=n: np.cos(n * theta(x))
    exact_osc = mt.integrate(f_osc, 3**10)
    row = [
        rmse(mtr.integrate_randomized(np.abs, n, n_shifts=repetitions, rng=rng).estimates, exact_abs),
        rmse(mtr.integrate_rmt(np.abs, n, n_samples=repetitions, rng=rng).estimates, exact_abs),
        rmse(mtr.integrate_randomized(f_osc, n, n_shifts=repetitions, rng=rng).estimates, exact_osc),
        rmse(mtr.integrate_rmt(f_osc, n, n_samples=repetitions, rng=rng).estimates, exact_osc),
    ]
    rows.append(row)
    print(f"{n:5d} | {row[0]:11.2e} {row[1]:11.2e} | {row[2]:15.2e} {row[3]:13.2e}")

slopes = np.polyfit(np.log(ns), np.log(np.array(rows)), 1)[0]
print(f"{'rate':>5} | {slopes[0]:11.2f} {slopes[1]:11.2f} | {slopes[2]:15.2f} {slopes[3]:13.2f}")
