import math

import numpy as np
import pytest

import moebius_trapezoidal as mt


def test_gaussian_moments():
    for p, exact in [(0, 1.0), (2, 1.0), (4, 3.0), (6, 15.0)]:
        assert mt.integrate(lambda x: x**p, 243) == pytest.approx(exact, rel=1e-12)


def test_odd_moments_vanish():
    for p in [1, 3, 5]:
        assert abs(mt.integrate(lambda x: x**p, 81)) < 1e-12


def test_logistic_weight():
    assert mt.integrate(lambda x: np.ones_like(x), 243, "logistic") == pytest.approx(1.0, rel=1e-9)
    assert mt.integrate(lambda x: x**2, 729, "logistic") == pytest.approx(math.pi**2 / 3, rel=1e-9)


def test_lebesgue_weight():
    value = mt.integrate(lambda x: np.exp(-x * x), 81, None)
    assert value == pytest.approx(math.sqrt(math.pi), rel=1e-12)


def test_callable_weight():
    def normal_sd2(x):
        return np.exp(-x**2 / 8) / math.sqrt(8 * math.pi)  # N(0, 4) density
    assert mt.integrate(lambda x: x**2, 243, normal_sd2, scale=2.0) == pytest.approx(4.0, rel=1e-10)


def test_finite_smoothness_rate():
    """|x| with the Gaussian weight converges at the rate n^-2."""
    exact = math.sqrt(2 / math.pi)
    ns = np.array([2**m + 1 for m in range(4, 12)])
    errs = np.array([abs(mt.integrate(np.abs, int(n)) - exact) for n in ns])
    rate = np.polyfit(np.log(ns), np.log(errs), 1)[0]
    assert -2.2 < rate < -1.8


def test_scale_and_loc():
    # the integral does not depend on the change of variables
    for scale, loc in [(0.5, 0.0), (2.0, 1.0), (3.0, -1.0)]:
        assert mt.integrate(lambda x: x**2, 729, scale=scale, loc=loc) == pytest.approx(1.0, rel=1e-10)


def test_nodes_weights():
    x, w = mt.nodes_weights(81)
    assert np.all(np.diff(x) > 0)
    assert np.all(w > 0)
    assert w.sum() == pytest.approx(1.0, rel=1e-8)
    x_all, w_all = mt.nodes_weights(81, prune=False)
    assert x_all.size == 81 and w_all.size == 81


def test_vector_valued_integrand():
    value = mt.integrate(lambda x: np.stack([x**2, x**4], axis=1), 243)
    assert value.shape == (2,)
    np.testing.assert_allclose(value, [1.0, 3.0], rtol=1e-12)


def test_scalar_integrand():
    value = mt.integrate(lambda x: math.cos(x), 243, vectorized=False)
    assert value == pytest.approx(math.exp(-0.5), rel=1e-12)


def test_wrong_shape_raises():
    with pytest.raises(ValueError):
        mt.integrate(lambda x: 1.0, 9)


def test_quad_reuses_points():
    result = mt.quad(lambda x: x**2, rtol=1e-13)
    assert result.converged
    assert result.value == pytest.approx(1.0, rel=1e-13)
    assert result.value == pytest.approx(mt.integrate(lambda x: x**2, result.n), rel=1e-14)


def test_quad_finite_smoothness():
    result = mt.quad(np.abs, rtol=1e-8)
    assert result.converged
    assert abs(result.value - math.sqrt(2 / math.pi)) < 1e-8


def test_bad_arguments():
    with pytest.raises(ValueError):
        mt.integrate(np.abs, 0)
    with pytest.raises(ValueError):
        mt.integrate(np.abs, 9, scale=0.0)
    with pytest.raises(ValueError):
        mt.integrate(np.abs, 9, "unknown")
