#
# Copyright (C) 2014-2024  Leo Singer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
"""Construct a LIGO-LW XML power spectral density file for a network of
detectors by evaluating a model power noise sensitivity curve.
"""

from argparse import SUPPRESS
import inspect

from . import ArgumentParser, FileType, register_to_xmldoc, write_fileobj

psd_name_prefix = 'SimNoisePSD'


def parser():
    import lal
    import lalsimulation

    # Get names of PSD functions.
    psd_names = sorted(
        name[len(psd_name_prefix):]
        for name, func in inspect.getmembers(lalsimulation)
        if name.startswith(psd_name_prefix) and callable(func) and (
            '(double f) -> double' in func.__doc__ or
            '(REAL8FrequencySeries psd, double flow) -> int' in func.__doc__))

    parser = ArgumentParser()
    parser.add_argument(
        '-o', '--output', metavar='OUT.xml[.gz]', type=FileType('wb'),
        default='-', help='Name of output file [default: stdout]')
    parser.add_argument(
        '--df', metavar='Hz', type=float, default=1.0,
        help='Frequency step size [default: %(default)s]')
    parser.add_argument(
        '--f-max', metavar='Hz', type=float, default=2048.0,
        help='Maximum frequency [default: %(default)s]')

    detector_group = parser.add_argument_group(
        'detector noise curves',
        'Options to select noise curves for detectors.\n\n'
        'All detectors support the following options:\n\n' +
        '\n'.join(psd_names))

    scale_group = parser.add_argument_group(
        'detector scaling',
        'Options to apply scale factors to noise curves for detectors.\n'
        'For example, a scale factor of 2 means that the amplitude spectral\n'
        'density is multiplied by 1/2 so that the range is multiplied by a 2.')

    # Add options for individual detectors
    for detector in lal.CachedDetectors:
        name = detector.frDetector.name
        prefix = detector.frDetector.prefix
        detector_group.add_argument(
            '--' + prefix, choices=psd_names,
            metavar='func', default=SUPPRESS,
            help='PSD function for {0} detector'.format(name))
        scale_group.add_argument(
            '--' + prefix + '-scale', type=float, default=SUPPRESS,
            help='Scale range for {0} detector'.format(name))

    return parser


def main(args=None):
    p = parser()
    with p.parse_args(args) as opts:
        import lal.series
        import lalsimulation
        import numpy as np
        from ..bayestar.filter import vectorize_swig_psd_func

        # Add basic options.

        psds = {}

        n = int(opts.f_max // opts.df)
        f = np.arange(n) * opts.df

        detectors = [d.frDetector.prefix for d in lal.CachedDetectors]

        for detector in detectors:
            psd_name = getattr(opts, detector, None)
            if psd_name is None:
                continue
            scale = 1 / np.square(getattr(opts, detector + '_scale', 1.0))
            func = getattr(lalsimulation, psd_name_prefix + psd_name)
            series = lal.CreateREAL8FrequencySeries(
                psd_name, 0, 0, opts.df, lal.SecondUnit, n)
            if '(double f) -> double' in func.__doc__:
                series.data.data = vectorize_swig_psd_func(
                    psd_name_prefix + psd_name)(f)
            else:
                func(series, 0.0)

                # Find indices of first and last nonzero samples.
                nonzero = np.flatnonzero(series.data.data)
                # FIXME: int cast seems to be needed on old versions of Numpy
                first_nonzero = int(nonzero[0])
                last_nonzero = int(nonzero[-1])

                # Truncate
                series = lal.CutREAL8FrequencySeries(
                    series, first_nonzero, last_nonzero - first_nonzero + 1)
                series.f0 = first_nonzero * series.deltaF

                series.name = psd_name
            series.data.data *= scale
            psds[detector] = series

        xmldoc = lal.series.make_psd_xmldoc(psds)
        register_to_xmldoc(xmldoc, p, opts)
        write_fileobj(xmldoc, opts.output)
