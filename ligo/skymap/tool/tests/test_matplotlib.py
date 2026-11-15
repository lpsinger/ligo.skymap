import argparse
from unittest import mock

import pytest
from matplotlib import pyplot as plt

from ..matplotlib import MatplotlibFigureType


def test_output_file(tmpdir):
    """Specifying an output file should never require a terminal."""
    filename = str(tmpdir / "plot.png")
    figure_type = MatplotlibFigureType()
    save = figure_type(filename)
    assert callable(save)


def test_show_without_terminal(monkeypatch):
    """Asking to show the plot interactively should fail when standard
    input is not a terminal, rather than silently opening a GUI window.
    """
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    figure_type = MatplotlibFigureType()
    with pytest.raises(argparse.ArgumentTypeError):
        figure_type("-")


def test_show_with_terminal(monkeypatch):
    """Asking to show the plot interactively should succeed when standard
    input is a terminal.

    Patch out ``switch_backend`` so that this test does not actually change
    the active Matplotlib backend as a side effect.
    """
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    with mock.patch.object(plt, "switch_backend") as switch_backend:
        figure_type = MatplotlibFigureType()
        show = figure_type("-")
    assert callable(show)
    switch_backend.assert_called_once_with(plt.rcParamsOrig["backend"])
