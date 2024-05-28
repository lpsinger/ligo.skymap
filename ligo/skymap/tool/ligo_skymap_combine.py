#
# Copyright (C) 2018-2024  Tito Dal Canton, Eric Burns, Leo Singer
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""Combine different sky localizations of a common event observed by different
instruments in order to form a more constrained localization.

If one of the input maps contains distance information (for instance from
BAYESTAR or LALInference) then the marginal distance posterior in the output
map is updated according to the restriction in sky location imposed by the
other input map(s). Only one input map can currently have distance
information.
"""

from . import ArgumentParser, FileType


def parser():
    parser = ArgumentParser()
    parser.add_argument('input', metavar='INPUT.fits[.gz]',
                        type=FileType('rb'), nargs='+',
                        help='Input sky localizations')
    # FIXME the output option has type str because astropy.io.fits.writeto()
    # only honors the .gz extension when given a file name string (as of 3.0.1)
    parser.add_argument('output', metavar='OUTPUT.fits[.gz]', type=str,
                        help='Output combined sky localization')
    parser.add_argument('--origin', type=str,
                        help='Optional tag describing the organization'
                             ' responsible for the combined output')
    return parser


def main(args=None):
    with parser().parse_args(args) as args:
        from textwrap import wrap
        import numpy as np
        import astropy_healpix as ah
        from astropy.io import fits
        from astropy.time import Time
        import healpy as hp

        from ..distance import parameters_to_marginal_moments
        from ..io import read_sky_map, write_sky_map

        input_skymaps = []
        dist_mu = dist_sigma = dist_norm = None
        for input_file in args.input:
            with fits.open(input_file) as hdus:
                header = hdus[0].header.copy()
                header.extend(hdus[1].header)
                has_distance = 'DISTMU' in hdus[1].columns.names
                data, meta = read_sky_map(hdus, nest=True,
                                          distances=has_distance)

            if has_distance:
                if dist_mu is not None:
                    raise RuntimeError('only one input localization can have '
                                       'distance information')
                dist_mu = data[1]
                dist_sigma = data[2]
                dist_norm = data[3]
            else:
                data = (data,)

            nside = ah.npix_to_nside(len(data[0]))
            input_skymaps.append((nside, data[0], meta, header))

        max_nside = max(x[0] for x in input_skymaps)

        # upsample sky posteriors to maximum resolution and combine them
        combined_prob = None
        for nside, prob, _, _ in input_skymaps:
            if nside < max_nside:
                prob = hp.ud_grade(prob, max_nside, order_in='NESTED',
                                   order_out='NESTED')
            if combined_prob is None:
                combined_prob = np.ones_like(prob)
            combined_prob *= prob

        # normalize joint posterior
        norm = combined_prob.sum()
        if norm == 0:
            raise RuntimeError('input sky localizations are disjoint')
        combined_prob /= norm

        out_kwargs = {'gps_creation_time': Time.now().gps,
                      'nest': True}
        if args.origin is not None:
            out_kwargs['origin'] = args.origin

        # average the various input event times
        input_gps = [
            x[2]['gps_time'] for x in input_skymaps if 'gps_time' in x[2]]
        if input_gps:
            out_kwargs['gps_time'] = np.mean(input_gps)

        # combine instrument tags
        out_instruments = set()
        for x in input_skymaps:
            if 'instruments' in x[2]:
                out_instruments.update(x[2]['instruments'])
        out_kwargs['instruments'] = ','.join(out_instruments)

        # update marginal distance posterior, if available
        if dist_mu is not None:
            if ah.npix_to_nside(len(dist_mu)) < max_nside:
                dist_mu = hp.ud_grade(dist_mu, max_nside, order_in='NESTED',
                                      order_out='NESTED')
                dist_sigma = hp.ud_grade(dist_sigma, max_nside,
                                         order_in='NESTED', order_out='NESTED')
                dist_norm = hp.ud_grade(dist_norm, max_nside,
                                        order_in='NESTED', order_out='NESTED')
            distmean, diststd = parameters_to_marginal_moments(combined_prob,
                                                               dist_mu,
                                                               dist_sigma)
            out_data = (combined_prob, dist_mu, dist_sigma, dist_norm)
            out_kwargs['distmean'] = distmean
            out_kwargs['diststd'] = diststd
        else:
            out_data = combined_prob

        # save input headers in output history
        out_kwargs['HISTORY'] = []
        for i, x in enumerate(input_skymaps):
            out_kwargs['HISTORY'].append('')
            out_kwargs['HISTORY'].append(
                'Headers of HDUs 0 and 1 of input file {:d}:'.format(i))
            out_kwargs['HISTORY'].append('')
            for line in x[3].tostring(sep='\n',
                                      endcard=False,
                                      padding=False).split('\n'):
                out_kwargs['HISTORY'].extend(wrap(line, 72))

        write_sky_map(args.output, out_data, **out_kwargs)
