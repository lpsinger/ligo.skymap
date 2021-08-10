.. highlight:: sh

Localizing Simulated Events with BAYESTAR
=========================================

This tutorial will show you how to create a collection of sky maps for
simulated binary neutron star mergers using BAYESTAR.

Some of the command-line tools used in this tutorial come from LALSuite or
python-ligo-lw. If you followed the :doc:`Quick Start installation instructions
<install>`, then these tools should have been installed automatically with the
other dependencies of ligo.skymap.

For each command below, we show some of the most useful options for customizing
their behavior. Inline comments look ```# like this` \`` and should be
ignored by your shell.

.. contents::
    :local:

1. Draw Injections
------------------

Draw a sample of injections, or simulated events, from an astrophysical
distribution using the ``lalapps_inspinj`` tool from LALSuite. This tool
has many options to control the distribution of each parameter. To see all
options, run ``lalapps_inspinj --help``.

Run the following command to generate a handful of BNS events::

    lalapps_inspinj \
    `# Write output to inj.xml.` \
    -o inj.xml \
    `# Mass distribution.` \
    `# In this example, the masses are pinned to 1.4 and 1.4 Msun.` \
    --m-distr fixMasses --fixed-mass1 1.4 --fixed-mass2 1.4 \
    `# Coalescence time distribution: adjust time step, start, and stop` \
    `# time to control the number of injections.` \
    --t-distr uniform --time-step 7200 \
    --gps-start-time 1000000000 \
    --gps-end-time 1000086400 \
    `# Distance distribution: uniform in Euclidean volume.` \
    `# WARNING: distances are in kpc.` \
    --d-distr volume \
    --min-distance 1 --max-distance 600e3 \
    `# Sky position and inclination distribution: isotropic.` \
    --l-distr random --i-distr uniform \
    `# Write a table of CBC injections to inj.xml.` \
    --f-lower 30 --disable-spin \
    --waveform TaylorF2threePointFivePN

The output is saved in LIGO-LW format to ``inj.xml``.

2. Get Noise Curves
-------------------

Create discretely sampled noise PSDs for all of your detectors using
:doc:`bayestar-sample-model-psd <../tool/bayestar_sample_model_psd>`::

    bayestar-sample-model-psd \
    `# Write output to psd.xml.` \
    -o psd.xml \
    `# Specify noise models for desired detectors.` \
    --H1=aLIGOZeroDetHighPower \
    --L1=aLIGOZeroDetHighPower \
    --I1=aLIGOZeroDetHighPower \
    --V1=AdvVirgo \
    --K1=KAGRA \
    `# Optional: apply scale factor to selected PSDs to increase or` \
    `# decrease their sensitivity. The PSD is multiplied by a factor of one ` \
    `# over scale squared; the horizon distance is multiplied by the scale.` \
    --I1-scale=0.75

The output is saved in LIGO-LW format to ``psd.xml``.

3. Run Matched Filter Pipeline
------------------------------

Run a simulated matched-filter search pipeline to add measurement error and
determine which injections are "detected" in coincidence by two or more of your
simulated GW observatories. Run the following command::

    bayestar-realize-coincs \
    `# Write output to coinc.xml.` \
    -o coinc.xml \
    `# Use the injections and noise PSDs that we generated.` \
    inj.xml --reference-psd psd.xml \
    `# Specify which detectors are in science mode.` \
    --detector H1 L1 V1 I1 K1 \
    `# Optionally, add Gaussian noise (rather than zero noise).` \
    --measurement-error gaussian-noise \
    `# Optionally, adjust the detection threshold: single-detector` \
    `# SNR, network SNR, and minimum number of detectors above` \
    `# threshold to form a coincidence.` \
    --snr-threshold 4.0 \
    --net-snr-threshold 12.0 \
    --min-triggers 2 \
    `# Optionally, save triggers that were below the single-detector` \
    `# threshold.` \
    --keep-subthreshold

The output is saved in LIGO-LW format to ``coinc.xml``.

4. Run BAYESTAR
---------------

Finally, make sky maps for your simulated events using
:doc:`bayestar-localize-coincs <../tool/bayestar_localize_coincs>`.
If you are working on a computing cluster that uses the
`HTCondor <https://research.cs.wisc.edu/htcondor/>`_ job scheduler, then you
can add the ``--condor-submit`` to automatically submit the BAYESTAR jobs to
the cluster. Run this command::

    # IMPORTANT: HIGHLY RECOMMENDED IF USING A SHARED WORKSTATION.
    # Explicitly set the number of OpenMP threads
    # instead of using all available cores.
    export OMP_NUM_THREADS=4

    # Run BAYESTAR on all coincident events in coinc.xml.
    bayestar-localize-coincs coinc.xml \
    `# Optional: submit jobs to Condor` \
    `# instead of running BAYESTAR locally.` \
    --condor-submit

The output is saved in the current working directory to FITS files named
``0.fits``, ``1.fits``, etc.

5. Find Injections (Optional)
-----------------------------

Optionally, if you want to generate P-P plots, you need to convert the
``coinc.xml`` file to SQLite using the ``ligolw_sqlite`` tool from
python-ligo-lw. Run the following command::

    ligolw_sqlite --preserve-ids --replace --database coinc.sqlite coinc.xml

The output is saved in SQLite format as ``coinc.sqlite``.

6. Analyze Sky Maps
-------------------

Use the :doc:`ligo-skymap-stats <../tool/ligo_skymap_stats>` tool
to gather summary statistics including credible areas for each sky map::

    ligo-skymap-stats \
    `# Write output to bayestar.tsv.` \
    -o bayestar.tsv \
    `# Include this option to enable P-P plots.` \
    --database coinc.sqlite \
    `# Read all sky maps in this directory.` \
    *.fits \
    `# Optional: calculate the 50% and 90% credible areas.` \
    --contour 50 90 \
    `# Optional: calculate the probability contained within the smallest` \
    `# credible regions of 10 and 100 deg2.` \
    --area 10 100 \
    `# Optional: count the number of disjoint patches on the sky.` \
    `# WARNING: this option makes the script very slow!` \
    --modes \
    `# Optional, but highly recommended: analyze sky maps using multiple` \
    `# threads. In this example, we use 8 worker processes.` \
    -j 8

The output is saved in tab-separated value format as ``bayestar.tsv``.

7. Make Summary Plots
---------------------

Lastly, make summary graphs including histograms and P-P plots using
:doc:`ligo-skymap-plot-stats <../tool/ligo_skymap_plot_stats>`::

    ligo-skymap-plot-stats bayestar.tsv
