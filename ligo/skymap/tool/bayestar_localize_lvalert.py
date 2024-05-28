#
# Copyright (C) 2013-2024  Leo Singer
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
r"""
Listen for new events from IGWN Alert and perform sky localization.

`bayestar-localize-lvalert` supports two modes of operation. You can
explicitly specify the GraceDb ID on the command line, as in::

    $ bayestar-localize-lvalert T90713

Or, `bayetar-localize-lvalert` can read GraceDB IDs from stdin (e.g., from the
terminal, or redirected from a fifo)::

    $ mkfifo /var/run/bayestar
    $ tail -F /var/run/bayestar | bayestar_localize_lvalert &
    $ echo T90713 > /var/run/bayestar
"""

from . import (
    ArgumentParser, EnableAction, get_waveform_parser, get_posterior_parser,
    get_mcmc_parser, get_random_parser, iterlines)


def parser():
    parser = ArgumentParser(
        parents=[get_waveform_parser(), get_posterior_parser(),
                 get_mcmc_parser(), get_random_parser()])
    parser.add_argument(
        '-d', '--disable-detector', metavar='X1', type=str, nargs='+',
        help='disable certain detectors')
    parser.add_argument(
        '-N', '--dry-run', action='store_true',
        help='Dry run; do not update GraceDB entry')
    parser.add_argument(
        '--no-tag', action='store_true',
        help='Do not set lvem tag for GraceDB entry')
    parser.add_argument(
        '-o', '--output', metavar='FILE.fits[.gz]', default='bayestar.fits',
        help='Name for uploaded file')
    parser.add_argument(
        '--enable-multiresolution', action=EnableAction, default=True,
        help='generate a multiresolution HEALPix map')
    parser.add_argument(
        'graceid', metavar='G123456', nargs='*',
        help='Run on these GraceDB IDs. If no GraceDB IDs are listed on the '
        'command line, then read newline-separated GraceDB IDs from stdin.')
    return parser


def main(args=None):
    with parser().parse_args(args) as opts:
        import logging
        import os
        import re
        import sys
        import tempfile
        import urllib.parse
        from ..bayestar import localize, rasterize
        from ..io import fits
        from ..io import events
        from .. import omp
        from ..util.file import rename
        import ligo.gracedb.logging
        import ligo.gracedb.rest

        log = logging.getLogger('BAYESTAR')

        log.info('Using %d OpenMP thread(s)', omp.num_threads)

        # If no GraceDB IDs were specified on the command line, then read them
        # from stdin line-by-line.
        graceids = opts.graceid if opts.graceid else iterlines(sys.stdin)

        # Fire up a GraceDb client
        # FIXME: Mimic the behavior of the GraceDb command line client, where
        # the environment variable GRACEDB_SERVICE_URL overrides the default
        # service URL. It would be nice to get this behavior into the gracedb
        # package itself.
        gracedb = ligo.gracedb.rest.GraceDb(
            os.environ.get(
                'GRACEDB_SERVICE_URL', ligo.gracedb.rest.DEFAULT_SERVICE_URL))

        # Determine the base URL for event pages.
        scheme, netloc, *_ = urllib.parse.urlparse(gracedb._service_url)
        base_url = urllib.parse.urlunparse(
            (scheme, netloc, 'events', '', '', ''))

        if opts.chain_dump:
            chain_dump = re.sub(r'.fits(.gz)?$', r'.hdf5', opts.output)
        else:
            chain_dump = None

        tags = ("sky_loc",)
        if not opts.no_tag:
            tags += ("lvem",)

        event_source = events.gracedb.open(graceids, gracedb)

        if opts.disable_detector:
            event_source = events.detector_disabled.open(
                event_source, opts.disable_detector)

        for graceid in event_source.keys():

            try:
                event = event_source[graceid]
            except:  # noqa: E722
                log.exception('failed to read event %s from GraceDB', graceid)
                continue

            # Send log messages to GraceDb too
            if not opts.dry_run:
                handler = ligo.gracedb.logging.GraceDbLogHandler(
                    gracedb, graceid)
                handler.setLevel(logging.INFO)
                logging.root.addHandler(handler)

            # A little bit of Cylon humor
            log.info('by your command...')

            try:
                # perform sky localization
                log.info("starting sky localization")
                sky_map = localize(
                    event, opts.waveform, opts.f_low, opts.min_distance,
                    opts.max_distance, opts.prior_distance_power,
                    opts.cosmology, mcmc=opts.mcmc, chain_dump=chain_dump,
                    enable_snr_series=opts.enable_snr_series,
                    f_high_truncate=opts.f_high_truncate,
                    rescale_loglikelihood=opts.rescale_loglikelihood)
                if not opts.enable_multiresolution:
                    sky_map = rasterize(sky_map)
                sky_map.meta['objid'] = str(graceid)
                sky_map.meta['url'] = '{}/{}'.format(base_url, graceid)
                log.info("sky localization complete")

                # upload FITS file
                with tempfile.TemporaryDirectory() as fitsdir:
                    fitspath = os.path.join(fitsdir, opts.output)
                    fits.write_sky_map(fitspath, sky_map, nest=True)
                    log.debug('wrote FITS file: %s', opts.output)
                    if opts.dry_run:
                        rename(fitspath, os.path.join('.', opts.output))
                    else:
                        gracedb.writeLog(
                            graceid, "BAYESTAR rapid sky localization ready",
                            filename=fitspath, tagname=tags)
                    log.debug('uploaded FITS file')
            except KeyboardInterrupt:
                # Produce log message and then exit if we receive SIGINT
                # (ctrl-C).
                log.exception("sky localization failed")
                raise
            except:  # noqa: E722
                # Produce log message for any otherwise uncaught exception.
                # Unless we are in dry-run mode, keep going.
                log.exception("sky localization failed")
                if opts.dry_run:
                    # Then re-raise the exception if we are in dry-run mode
                    raise

            if not opts.dry_run:
                # Remove old log handler
                logging.root.removeHandler(handler)
                del handler
