"""Randomly shifted Möbius-transformed trapezoidal rule with an error estimate.

The integral of |x| with respect to the standard Gaussian weight, as in
cot_trap.m, now with random shifts. The standard error is compared with the
actual error. Run from the python directory with

    python examples/randomized.py
"""

import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import moebius_trapezoidal_randomized as mtr  # noqa: E402

exact = math.sqrt(2 / math.pi)
rng = np.random.default_rng(2025)
n_shifts = 16

print(f"{'n':>6} {'estimate':>14} {'actual error':>13} {'std error':>11}")
for m in range(2, 13):
    n = 2**m
    r = mtr.integrate_randomized(np.abs, n, n_shifts=n_shifts, rng=rng)
    print(f"{n:6d} {r.value:14.10f} {abs(r.value - exact):13.3e} {r.std_error:11.3e}")

print(f"\nEach row uses {n_shifts} random shifts, i.e. n * {n_shifts} function values.")
print("value ± 2 * std error is an approximate 95% confidence interval,")
print("so the actual error is mostly, but not always, below two standard errors.")
