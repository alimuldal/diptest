import pytest
import numpy as np

import diptest as dt
import diptest.lib._diptest as _dt


def _generator(N, distance=0.5, sigma=1.0):
    """Draw a sample from a mixture of Gaussians.

    Parameters
    ----------
    N : int
        number of samples to draw
    distance : float, default=0.5
        distance between the means, the means will be zero +/- (distance/2)
    sigma : float, tuple, default=1.0
        standard deviation(s) of the Gaussians

    Returns
    -------
    sample : np.ndarray
        sample drawn from the mixture distribution

    """
    hd = distance / 2
    if N % 2 == 0:
        lN = uN =N // 2
    else:
        lN = int(np.ceil(N / 2))
        uN = N - lN
    mu_a = 0.0 - hd
    mu_b = 0.0 + hd

    if isinstance(sigma, float):
        sigma_a = sigma_b = sigma
    elif isinstance(sigma, (tuple, list)):
        sigma_a, sigma_b = sigma
    else:
        raise TypeError("`sigma` must be a float or tuple of floats")

    sample = np.empty(N, dtype=np.float64)
    sample[:lN] = np.random.normal(loc=mu_a, scale=sigma_a, size=lN)
    sample[lN:] = np.random.normal(loc=mu_b, scale=sigma_b, size=uN)
    return sample


def test_dipstat_default():
    """Test diptest.dipstat with default settings."""
    sample = _generator(100)
    dip = dt.dipstat(sample)
    assert not np.isnan(dip)
    assert not np.isinf(dip)


def test_dipstat_small():
    """Test diptest.dipstat with small array."""
    sample = _generator(4)
    dip = dt.dipstat(sample)
    assert not np.isnan(dip)
    assert not np.isinf(dip)


def test_dipstat_large():
    """Test diptest.dipstat with large array."""
    sample = _generator(10000)
    dip = dt.dipstat(sample)
    assert not np.isnan(dip)
    assert not np.isinf(dip)


def test_dipstat_sort_x():
    """Test diptest.dipstat behaviour for sort_x parameter."""
    sample = _generator(100)
    with pytest.raises(RuntimeError):
        dip = dt.dipstat(sample, sort_x=False)
    sample.sort()
    dip = dt.dipstat(sample, sort_x=True)
    assert not np.isnan(dip)
    assert not np.isinf(dip)


def test_dipstat_allow_zero():
    """Test diptest.dipstat behaviour for allow_zero parameter."""
    sample = np.ones(100)
    assert not np.isclose(
        dt.dipstat(sample, allow_zero=False),
        0.0
    )
    assert np.isclose(
        dt.dipstat(sample, allow_zero=True),
        0.0
    )
    sample = _generator(2)
    assert not np.isclose(
        dt.dipstat(sample, allow_zero=False),
        0.0
    )
    assert np.isclose(
        dt.dipstat(sample, allow_zero=True),
        0.0
    )


def test_dipstat_full_output():
    """Test diptest.dipstat with full output."""
    sample = _generator(100)
    dip, res = dt.dipstat(sample, full_output=True)
    exp_keys = {'lo', 'hi', 'xl', 'xu', 'gcm', 'lcm'}
    obs_keys = set(res.keys())
    assert len(obs_keys - exp_keys) == 0
    assert len(exp_keys - obs_keys) == 0
    assert not np.isnan(dip)
    assert not np.isinf(dip)


def test_dipstat_non_contiguous():
    """Test diptest.dipstat with non contiguous arrays."""
    sample = _generator(20)
    sample.sort()
    farr = sample.copy(order='F')
    carr = sample.copy(order='C')

    dip = dt.dipstat(farr, sort_x=False)
    assert not np.isnan(dip)
    assert not np.isinf(dip)
    np.arange(0, 20, step=2)

    dip = dt.dipstat(farr[::2], sort_x=False)
    assert not np.isnan(dip)
    assert not np.isinf(dip)


    dip = dt.dipstat(carr, sort_x=False)
    assert not np.isnan(dip)
    assert not np.isinf(dip)

    dip = dt.dipstat(carr[::2], sort_x=False)
    assert not np.isnan(dip)
    assert not np.isinf(dip)


def test_dipstat_2d():
    """Test diptest.dipstat with 2d arrays."""
    sample = _generator(20)
    sample.sort()

    dip = dt.dipstat(sample[:, None], sort_x=False)
    assert not np.isnan(dip)
    assert not np.isinf(dip)
    dip = dt.dipstat(sample[None, :], sort_x=False)
    assert not np.isnan(dip)
    assert not np.isinf(dip)

    sample_2d = sample.reshape(5, 4).copy()
    with pytest.raises(TypeError):
        dip = dt.dipstat(sample_2d, sort_x=False)
