from ..bayestar_inject import get_decisive_snr


def test_get_decisive_snr():
    assert get_decisive_snr([]) == 0.0
    assert get_decisive_snr([1.0, 3.0, 2.0, 4.0]) == 3.0
    assert get_decisive_snr([4.0]) == 4.0
