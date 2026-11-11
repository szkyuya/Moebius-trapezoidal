import math

import numpy as np
import pytest

import moebius_trapezoidal as mt
import moebius_trapezoidal_randomized as mtr


def test_half_shift_is_deterministic_rule():
    x, w = mtr.shifted_nodes_weights(81, 0.5)
    x0, w0 = mt.nodes_weights(81)
    np.testing.assert_allclose(x, x0, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(w, w0, rtol=1e-12, atol=1e-300)


def test_zero_shift_skips_infinity():
    x, w = mtr.shifted_nodes_weights(81, 0.0)
    assert np.all(np.isfinite(x))
    assert w.sum() == pytest.approx(1.0, rel=1e-6)


def test_every_shift_integrates_smooth_functions():
    for shift in [0.0, 0.1, 0.37, 0.9]:
        x, w = mtr.shifted_nodes_weights(243, shift)
        assert w @ x**2 == pytest.approx(1.0, rel=1e-12)


def test_unbiased_with_standard_error():
    exact = math.sqrt(2 / math.pi)
    r = mtr.integrate_randomized(np.abs, 33, n_shifts=64, rng=1)
    assert r.std_error > 0
    assert abs(r.value - exact) < 5 * r.std_error
    assert r.estimates.shape == (64,)


def test_standard_error_decreases_with_n():
    errs = [mtr.integrate_randomized(np.abs, n, rng=2).std_error for n in (16, 64, 256)]
    assert errs[0] > errs[1] > errs[2]


def test_reproducible_with_seed():
    a = mtr.integrate_randomized(np.abs, 17, rng=3)
    b = mtr.integrate_randomized(np.abs, 17, rng=np.random.default_rng(3))
    assert a.value == b.value and a.std_error == b.std_error


def test_vector_valued_integrand():
    r = mtr.integrate_randomized(lambda x: np.stack([x**2, np.abs(x)], axis=1), 65, rng=4)
    assert r.value.shape == (2,) and r.std_error.shape == (2,)
    assert r.value[0] == pytest.approx(1.0, rel=1e-6)
    assert abs(r.value[1] - math.sqrt(2 / math.pi)) < 5 * r.std_error[1]


def test_other_weights_and_options():
    r = mtr.integrate_randomized(lambda x: x**2, 243, "logistic", rng=5)
    assert r.value == pytest.approx(math.pi**2 / 3, rel=1e-6)
    r = mtr.integrate_randomized(lambda x: x, 243, scale=2.0, loc=1.0, rng=6)
    assert abs(r.value) < 1e-8


def test_bad_arguments():
    with pytest.raises(ValueError):
        mtr.integrate_randomized(np.abs, 9, n_shifts=1)
    with pytest.raises(ValueError):
        mtr.shifted_nodes_weights(9, 1.0)
