from os import getpid
from time import sleep

import numpy as np
import pytest

from ..progress import progress_map, progress_map_vectorized


def func(x):
    sleep(np.random.uniform(0, 0.1))
    return np.square(x)


def func2(x):
    sleep(np.random.uniform(0, 0.1))
    return np.cos(x), np.sin(x)


@pytest.mark.parametrize("jobs", [1, 8, None])
def test_map(jobs):
    x = np.arange(20)
    result = list(progress_map(func, x, jobs=jobs))
    np.testing.assert_array_equal(result, np.square(x))


def map0(_):
    return getpid()


def map1(_):
    return list(progress_map(map0, range(8), jobs=8))


def map2():
    return sum(progress_map(map1, range(8), jobs=8), [])


def test_no_nested_pools():
    """Test that parallelism is disabled in nested calls to progress_map."""
    assert len(set(map2())) <= 8


@pytest.mark.parametrize("jobs", [1, 8, None])
def test_indefinite(jobs):
    """Test iteration over a collection of indefinite length."""
    assert list(progress_map(np.square, (i for i in range(3)), jobs=jobs)) == [0, 1, 4]


@pytest.mark.parametrize("jobs", [1, 8, None])
def test_map_vectorized(jobs):
    x = np.arange(20)
    result = progress_map_vectorized(func, x, jobs=jobs)
    np.testing.assert_array_equal(result, np.square(x))


@pytest.mark.parametrize("jobs", [1, 8, None])
def test_map_vectorized_nout(jobs):
    x = np.arange(20)
    result1, result2 = progress_map_vectorized(func2, x, jobs=jobs, nout=2)
    np.testing.assert_array_equal(result1, np.cos(x))
    np.testing.assert_array_equal(result2, np.sin(x))
