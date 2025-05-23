"""
.. module:: plotspec
   :synopsis: Creates preview plots for the provided JWST spectrum.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from __future__ import division
import copy
import os
import sys
from builtins import str
#--------------------
# External Imports
#--------------------
import matplotlib
from matplotlib.ticker import FormatStrFormatter
from matplotlib import pyplot
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots.utils.specutils.specutilserror import SpecUtilsError
from spec_plots.utils.specutils.debug_oplot import debug_oplot
from spec_plots.utils.specutils.calc_covering_fraction import (
    calc_covering_fraction)
from spec_plots import __version__

if matplotlib.get_backend().lower() != 'agg':
    pyplot.switch_backend('Agg')
#--------------------

#--------------------
def plotspec(jwst_spectrum, output_type, output_file, flux_scale_factor,
             fluxerr_scale_factor, plot_metrics, dpi_val=96., output_size=1024,
             debug=False, full_ylabels=False, optimize=True):
    """
    Accepts a JWSTSpectrum object from the READSPEC function and produces
    preview plots.

    :param jwst_spectrum: JWST spectrum as returned by READSPEC.

    :type jwst_spectrum: JWSTSpectrum

    :param output_type: What kind of output to make?

    :type output_type: str

    :param output_file: Name of output file (including full path).

    :type output_file: str

    :param flux_scale_factor: Max. allowed ratio between the flux and a median
        flux value, used in edge trimming.  Default = 10.

    :type flux_scale_factor: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty
        and a median flux uncertainty value, used in edge trimming.
        Default = 5.

    :type fluxerr_scale_factor: float

    :param plot_metrics: Collection of plot metrics (flux statistics, axis
        ranges, etc.) to use when making the plots.  These are computed using
        `utils.specutils.calc_plot_metrics()`.

    :type plot_metrics: dict

    :param dpi_val: The DPI value of your device's monitor.  Affects the size of
        the output plots.  Default = 96. (applicable to most modern monitors).

    :type dpi_val: float

    :param output_size: Size of plot in pixels (plots are square in dimensions).
        Defaults to 1024.

    :param output_size: int

    :param debug: Should the output plots include debugging information
        (color-coded data points based on rejection criteria, shaded exclude
        regions)?  Default = False.

    :type debug: bool

    :param full_ylabels: Should the y-labels contain the full values (including
        the power of 10 in scientific notation)?  Default = False.

    :type full_ylabels: bool

    :param optimize: If set to True, will use a slightly optimized version of
        determining the plot covering fraction.

    :type optimize: bool

    :raises: OSError, utils.specutils.SpecUtilsError

    .. note::

         This function assumes a screen resolution of 96 DPI in order to
         generate plots of the desired sizes.  This is because matplotlib works
         in units of inches and DPI rather than pixels.
    """

    # Make sure the plot size is set to an integer value.
    if output_size is not None:
        if not isinstance(output_size, int):
            output_size = int(round(output_size))

    # Make sure the output path exists, if not, create it.
    if output_type != 'screen':
        if (os.path.dirname(output_file) != "" and
                not os.path.isdir(os.path.dirname(output_file))):
            try:
                os.mkdir(os.path.dirname(output_file))
            except OSError as this_error:
                if this_error.errno == 13:
                    sys.exit("*** MAKE_JWST_SPEC_PREVIEWS ERROR:"
                                 " Output directory could not be created,"
                                 " " + repr(this_error.strerror) + "\n")
                else:
                    raise

    is_bigplot = output_size > 128

    # Start plot figure.
    this_figure, this_plotarea = pyplot.subplots(nrows=1, ncols=1,
                                                 figsize=(
                                                     output_size/dpi_val,
                                                     output_size/dpi_val),
                                                 dpi=dpi_val)

    # Adjust the plot geometry (margins, etc.) based on plot size.
    if is_bigplot:
        this_figure.subplots_adjust(hspace=0.3, top=0.915)
    else:
        this_figure.subplots_adjust(top=0.85, bottom=0.3, left=0.25, right=0.8)

    # Get the wavelengths, fluxes, flux uncertainties, and data quality
    # flags out of the spectrum.
    try:
        all_wls = jwst_spectrum.wavelengths
        all_fls = jwst_spectrum.fluxes
        all_flerrs = jwst_spectrum.fluxerrs
        all_dqs = jwst_spectrum.dqs
    except KeyError as the_error:
        raise SpecUtilsError("The provided stitched spectrum does not have"
                             " the expected format, missing key " +
                             str(the_error)+".")

    # Extract the optimal x-axis plot range from the plot_metrics dict,
    # since we use it a lot.
    optimal_xaxis_range = plot_metrics["optimal_xaxis_range"]

    # Plot the spectrum, but only if valid wavelength ranges for x-axis
    # are returned, otherwise plot a special "Fluxes Are All Zero or NaN" plot.
    if (
            all(numpy.isfinite(optimal_xaxis_range)) and
            not all(numpy.isnan(plot_metrics["y_axis_range"]))
        ):
        # We plot the spectrum as a regular line for use
        # in calc_covering_fraction, it will be removed later.
        this_line = this_plotarea.plot(all_wls, all_fls, 'b')

        # Update y-axis range, but only adjust the ranges if this isn't
        # an all-zero flux case (and not in debug mode, in which case I want
        # to see the entire y-axis range).
        if not debug:
            this_plotarea.set_ylim(plot_metrics["y_axis_range"])

        covering_fractions = calc_covering_fraction(
            this_figure, numpy.asarray([this_plotarea]), 0, optimize=optimize)
        # Note: here we remove the line we plotted before, it was only
        # so that calc_covering_fraction would have someting to draw on the
        # canvas and thereby determine which pixels were "blue" (i.e., part
        # of the plotted spectrum vs. background).
        this_line.pop(0).remove()
        # Now we plot the spectrum as a LineCollection so that the
        # transparency will have the desired effect, but, this is not
        # rendered on the canvas inside calc_covering_fraction, hence why we
        # need to plot it both as a regular line first.
        # Note: because we are re-using a LineCollection object in the
        # array of plot_metrics (specifically, when creating the
        # thumbnail-sized plot), we have to use a copy of the LineCollection
        # object, otherwise it will have Axes, Figure, etc. all defined and
        # resetting them to None does not work.  Since this is only an issue
        # with thumbnail-sizes, this is only relevant for the first
        # LineCollection.
        this_collection = this_plotarea.add_collection(copy.copy(
            plot_metrics["line_collection"]))

        if covering_fractions > 30.:
            this_collection.set_alpha(0.1)

        # Turn on plot grid lines.
        this_plotarea.grid(True, linestyle="dashed")

        if is_bigplot:
            this_figure.suptitle(os.path.basename(
                jwst_spectrum.orig_file), fontsize=18, color='r')

        if debug:
            # Overplot points color-coded based on rejection criteria.
            debug_oplot(this_plotarea, "jwst", all_wls, all_fls,
                        all_flerrs, all_dqs,
                        plot_metrics["median_flux"],
                        plot_metrics["median_fluxerr"],
                        flux_scale_factor, fluxerr_scale_factor,
                        plot_metrics["fluxerr_95th"])

            # Overplot the x-axis edges that are trimmed to define the
            # y-axis plot range as a shaded area.
            this_plotarea.axvspan(numpy.nanmin(all_wls),
                                  optimal_xaxis_range[0],
                                  facecolor="lightgrey", alpha=0.5)
            this_plotarea.axvspan(optimal_xaxis_range[1],
                                  numpy.nanmax(all_wls),
                                  facecolor="lightgrey", alpha=0.5)

            # Overplot the avoid regions in a light color as a shaded area.
            for region in plot_metrics["avoid_regions"]:
                this_plotarea.axvspan(region.minwl, region.maxwl,
                                      facecolor="lightgrey", alpha=0.5)

        # This is where we ensure the x-axis range is based on the full
        # x-axis range, rather than using the optimum x-axis range.  This is
        # done so that all the plots for a similar instrument setting will
        # have the same starting and ending plot values.
        min_wl = numpy.nanmin(all_wls)
        max_wl = numpy.nanmax(all_wls)
        xplot_buffer = (max_wl - min_wl) * 0.05
        this_plotarea.set_xlim([min_wl-xplot_buffer, max_wl+xplot_buffer])

        # Only use two tick labels (min and max wavelengths) for
        # thumbnails, because there isn't enough space otherwise.
        if not is_bigplot:
            this_plotarea.tick_params(labelsize=6., labelrotation=45.,
                                          axis='y', pad=0.)
            minwl = numpy.nanmin(all_wls)
            maxwl = numpy.nanmax(all_wls)
            this_plotarea.set_xticks([minwl, maxwl])
            this_plotarea.set_xticklabels(this_plotarea.get_xticks(),
                                          rotation=25.)
            this_plotarea.xaxis.set_major_formatter(FormatStrFormatter(
                "%6.1f"))
        else:
            # Make sure the font properties go back to normal.
            pyplot.rcdefaults()
            this_plotarea.tick_params(axis='x', labelsize=14)
            if jwst_spectrum.wavelengths_unit:
                this_plotarea.set_xlabel(
                    f"Wavelength ({jwst_spectrum.wavelengths_unit:latex})",
                    fontsize=16, color='k')
            else:
                this_plotarea.set_xlabel(
                    f"Wavelength (unknown unit)", fontsize=16, color='k')
            if jwst_spectrum.fluxes_unit:
                this_plotarea.set_ylabel(
                    f"Flux ({jwst_spectrum.fluxes_unit:latex})",
                    fontsize=16, color='k')
            else:
                this_plotarea.set_xlabel(
                    f"Wavelength (unknown unit)", fontsize=16, color='k')

            # If requested, include the powers of 10 part of the y-axis
            # tickmarks.
            if full_ylabels:
                this_plotarea.yaxis.set_major_formatter(FormatStrFormatter(
                    '%3.2E'))

    else:
        # Otherwise this is a spectrum that has all zero fluxes, all NaN fluxes
        # or some other problem, and we make a default plot. Define the optimal
        # x-axis range to span the original spectrum.
        optimal_xaxis_range = [numpy.nanmin(all_wls), numpy.nanmax(all_wls)]
        this_plotarea.set_xlim(optimal_xaxis_range)

        # Make the plot background grey to distinguish that this is a
        # `special` plot.  Turn off y-tick labels.
        this_plotarea.set_facecolor("lightgrey")
        this_plotarea.set_yticklabels([])

        # Configure the plot units, text size, and other markings based
        # on whether this is a large or thumbnail-sized plot.
        if not is_bigplot:
            this_plotarea.tick_params(labelsize=6., labelrotation=45.,
                                          axis='y', pad=0.)
            minwl = numpy.nanmin(all_wls)
            maxwl = numpy.nanmax(all_wls)
            this_plotarea.set_xticks([minwl, maxwl])
            this_plotarea.set_xticklabels(this_plotarea.get_xticks(),
                                          rotation=25.)
            this_plotarea.xaxis.set_major_formatter(FormatStrFormatter(
                "%6.1f"))
            textsize = "small"
            plottext = "Fluxes are \n all 0 or NaN."
        else:
            # Make sure the font properties go back to normal.
            pyplot.rcdefaults()
            this_plotarea.tick_params(axis='x', labelsize=14)
            if jwst_spectrum.wavelengths_unit:
                this_plotarea.set_xlabel(
                    f"Wavelength ({jwst_spectrum.wavelengths_unit:latex})",
                    fontsize=16, color='k')
            else:
                this_plotarea.set_xlabel(
                    f"Wavelength (unknown unit)", fontsize=16, color='k')
            if jwst_spectrum.fluxes_unit:
                this_plotarea.set_ylabel(
                    f"Flux ({jwst_spectrum.fluxes_unit:latex})",
                    fontsize=16, color='k')
            else:
                this_plotarea.set_xlabel(
                    f"Wavelength (unknown unit)", fontsize=16, color='k')

            # If requested, include the powers of 10 part of the y-axis
            # tickmarks.
            if full_ylabels:
                this_plotarea.yaxis.set_major_formatter(FormatStrFormatter(
                    '%3.2E'))

            textsize = "x-large"
            plottext = "Fluxes are all 0 or NaN."

        # Place the text with the informational message in the center of
        # the plot.
        this_plotarea.text(0.5, 0.5, plottext, horizontalalignment="center",
                           verticalalignment="center",
                           transform=this_plotarea.transAxes, size=textsize)

    # Display or plot to the desired format.
    if output_type != "screen":
        if output_size == 128:
            revised_output_file = output_file.split('.png')[0] + '_thumb.png'
        else:
            revised_output_file = output_file

        # Save figure.
        this_figure.savefig(revised_output_file, format=output_type,
                            dpi=dpi_val)

    elif output_type == "screen":
        pyplot.show()
#--------------------
