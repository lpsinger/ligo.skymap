import numpy as np
import pytest

from ..interpolation import interpolate_max


@pytest.mark.parametrize('method,expected_imax,expected_ymax', [
    ['lanczos', 4.183314387190362, 9.415514264876181+6.056606011770283j],
    ['catmull-rom', 4.606828337716582, 8.69890509129855+7.057010747638154j],
    ['quadratic-fit', 4.047784375601079, 9.635253882349538+5.677856407731738j],
    ['nearest-neighbor', 4, 9.713133+5.5589147j]])
def test_interpolate_max(method, expected_imax, expected_ymax):
    y = np.asarray([ 9.135017 -2.8185585j,  9.995214 -1.1222992j,  # noqa
                    10.682851 +0.8188147j, 10.645139 +3.0268786j,  # noqa
                     9.713133 +5.5589147j,  7.9043484+7.9039335j,  # noqa
                     5.511646 +9.333084j ,  2.905198 +9.715742j ,  # noqa
                     0.5302934+9.544538j ])  # noqa
    i = np.argmax(np.abs(y))
    window = 4
    imax, ymax = interpolate_max(i, y, window, method)
    assert imax == pytest.approx(expected_imax)
    assert ymax == pytest.approx(expected_ymax)
