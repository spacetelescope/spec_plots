"""
.. module:: make_fits
   :synopsis: Makes a binary FITS table out of the provided COS spectrum.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from __future__ import division
import os
import sys
#--------------------
# External Imports
#--------------------
from astropy.io import fits
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots.utils.specutils_cos.get_segment_names import get_segment_names
from spec_plots import __version__
#--------------------

#--------------------
def make_fits(cos_spectrum, output_file, plot_metrics, input_file):
    """
    Accepts a COSSpectrum object from the READSPEC function and produces a
    binary FITS table.

    :param cos_spectrum: COS spectrum as returned by READSPEC.

    :type cos_spectrum: COSSpectrum

    :param output_file: Name of output file (including full path).

    :type output_file: str

    :param plot_metrics: Collection of plot metrics (flux statistics, axis
        ranges, etc.) to use when making the plots.  These are computed using
        `utils.specutils.calc_plot_metrics()`.

    :type plot_metrics: list

    :param input_file: The full path and name of the input file we are making
        a preview spectrum for.

    :type input_file: str
    """

    # Make sure the output path exists, if not, create it.
    if (os.path.dirname(output_file) != "" and
            not os.path.isdir(os.path.dirname(output_file))):
        try:
            os.mkdir(os.path.dirname(output_file))
        except OSError as this_error:
            if this_error.errno == 13:
                sys.stderr.write("*** MAKE_HST_SPEC_PREVIEWS ERROR: Output"
                                 " directory could not be created, " +
                                 repr(this_error.strerror)+"\n")
                exit(1)
            else:
                raise

    # Get list of segment names.
    segment_names = get_segment_names(cos_spectrum)

    # Loop over each segment.
    segment_ylims = []
    for i, segment in enumerate(segment_names):
        if i == 0:
            all_wls = cos_spectrum.segments[segment].wavelengths
            n_wls = len(cos_spectrum.segments[segment].wavelengths)
            all_fls = cos_spectrum.segments[segment].fluxes
            all_flerrs = cos_spectrum.segments[segment].fluxerrs
            all_dqs = cos_spectrum.segments[segment].dqs
            all_segments = numpy.repeat(segment, n_wls)
        else:
            all_wls = numpy.append(
                all_wls, cos_spectrum.segments[segment].wavelengths)
            n_wls = len(cos_spectrum.segments[segment].wavelengths)
            all_fls = numpy.append(
                all_fls, cos_spectrum.segments[segment].fluxes)
            all_flerrs = numpy.append(
                all_flerrs, cos_spectrum.segments[segment].fluxerrs)
            all_dqs = numpy.append(
                all_dqs, cos_spectrum.segments[segment].dqs)
            all_segments = numpy.append(
                all_segments, numpy.repeat(segment, n_wls))
        # Add this segments suggested y-range to the list so we can add it
        # to the header as keywords later.
        segment_ylims.append({'segment' : segment,
                              'ymin' : plot_metrics[i]['y_axis_range'][0],
                              'ymax' : plot_metrics[i]['y_axis_range'][1]})

    # Mask any wavelengths that are within an Avoid Region.
    all_masks = numpy.repeat(0, len(all_wls))
    for pmetric in plot_metrics:
        for avregion in pmetric['avoid_regions']:
            mask_indices = numpy.where((all_wls >= avregion.minwl) &
                                       (all_wls <= avregion.maxwl))[0]
            all_masks[mask_indices] = 1

    # Output the FITS file.
    wl_col = fits.Column(name='wls', format='E', array=all_wls, unit="angstrom")
    fl_col = fits.Column(name='fls', format='D', array=all_fls,
                         unit="erg/s/cm**2/angstrom")
    fle_col = fits.Column(name='flerrs', format='D', array=all_flerrs,
                          unit="erg/s/cm**2/angstrom")
    dq_col = fits.Column(name='dq_wgts', format='I', array=all_dqs)
    seg_col = fits.Column(name='segment', format='4A', array=all_segments)
    mask_col = fits.Column(name='masked', format='I', array=all_masks)

    # Create the extension HDU.
    cols = fits.ColDefs([wl_col, fl_col, fle_col, dq_col, seg_col, mask_col])
    hdu1 = fits.BinTableHDU.from_columns(cols)

    # Add the suggested ymin and ymax values for each segment.
    for ylim in segment_ylims:
        hdu1.header.set(ylim['segment'].strip().upper()+'YMIN', ylim['ymin'],
                        "Suggested min. for " +
                        ylim['segment'].strip().upper() + " y-axis plot range.")
        hdu1.header.set(ylim['segment'].strip().upper()+'YMAX', ylim['ymax'],
                        "Suggested max. for " +
                        ylim['segment'].strip().upper() + " y-axis plot range.")

    # Create a barebones primary HDU.
    hdr = fits.Header()
    hdr.set('ORIGFILE', os.path.basename(input_file), "Input file name.")
    hdu0 = fits.PrimaryHDU(header=hdr)

    # Output the FITS file.
    hdulist = fits.HDUList([hdu0, hdu1])
    hdulist.writeto(output_file, overwrite=True)
#--------------------
