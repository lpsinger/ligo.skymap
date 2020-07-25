from time import sleep

import numpy as np
import pytest

from ..progress import progress_map


def func(x):
    sleep(np.random.uniform(0, 0.1))
    return np.square(x)


@pytest.mark.parametrize('jobs', [1, 8])
def test_map(jobs):
    x = np.arange(20)
    result = list(progress_map(func, x, jobs=jobs))
    np.testing.assert_array_equal(result, np.square(x))
