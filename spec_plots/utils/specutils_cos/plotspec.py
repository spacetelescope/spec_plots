"""
.. module:: plotspec
   :synopsis: Makes preview plots of the provided COS spectrum.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from __future__ import division
import os
import sys
from builtins import str
#--------------------
# External Imports
#--------------------
import matplotlib
from matplotlib.ticker import FormatStrFormatter
from matplotlib import pyplot
from matplotlib import rc
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots.utils.specutils_cos.get_segment_names import get_segment_names
from spec_plots.utils.specutils.specutilserror import SpecUtilsError
from spec_plots.utils.specutils.debug_oplot import debug_oplot
from spec_plots.utils.specutils.calc_covering_fraction import (
    calc_covering_fraction)
from spec_plots import __version__

if matplotlib.get_backend().lower() != 'agg':
    pyplot.switch_backend('Agg')
#--------------------

#--------------------
def plotspec(cos_spectrum, output_type, output_file,
             flux_scale_factor, fluxerr_scale_factor, plot_metrics, dpi_val=96.,
             output_size=1024, debug=False, full_ylabels=False,
             stitched_spectrum=None, title_addendum="", optimize=True):
    """
    Accepts a COSSpectrum object from the READSPEC function and produces preview
    plots.

    :param cos_spectrum: COS spectrum as returned by READSPEC.

    :type cos_spectrum: COSSpectrum

    :param output_type: What kind of output to make?

    :type output_type: str

    :param output_file: Name of output file (including full path).

    :type output_file: str

    :param flux_scale_factor: Max. allowed ratio between the flux and a median
        flux value, used in edge trimming.  Default = 10.

    :type flux_scale_factor: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty
        and a median flux uncertainty value, used in edge trimming.  Default =
        5.

    :type fluxerr_scale_factor: float

    :param plot_metrics: Collection of plot metrics (flux statistics, axis
        ranges, etc.) to use when making the plots.  These are computed using
        `utils.specutils.calc_plot_metrics()`.

    :type plot_metrics: list

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

    :param stitched_spectrum: The stitched version of the COS spectrum, where
        each segment has been stitched using
        `utils.specutils.stitch_components()`.
        This is required is making a small (thumb-sized) plot, but not used if
        making a large-sized plot.

    :type stitched_spectrum: dict

    :param title_addendum: A plot title addendum that contains possible warnings
        (e.g., if all DQ or DQ_WGT flags are defined as bad).  This is only used
        in the case of large-sized COS plots, otherwise this information is
        included in the stitched spectrum (COS large-sized plots do not use a
        stitched version of the spectrum).

    :type title_addendum: str

    :param optimize: If set to True, will use a slightly optimized version of
        determining the plot covering fraction.

    :type optimize: bool

    :raises: OSError
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
                    sys.stderr.write("*** MAKE_HST_SPEC_PREVIEWS ERROR: Output"
                                     " directory could not be created, " +
                                     repr(this_error.strerror)+"\n")
                    sys.exit()
                else:
                    raise

    # Get list of segment names.
    segment_names = get_segment_names(cos_spectrum)

    # Determine the number of subplots to make depending on output size.
    # Get a list of subplot segment names to iterate over when plotting the
    # figure.
    if output_size > 128:
        is_bigplot = True
        n_subplots = len(segment_names)
        subplot_segment_names = segment_names
    else:
        is_bigplot = False
        # Change segment_names to be a one-element array, because we need to
        # join all the segments into a single array if making a thumbnail plot.
        # Make sure it's still a list...
        n_subplots = 1
        subplot_segment_names = ["-".join(segment_names)]

    # Start the plot figure.
    this_figure, these_plotareas = pyplot.subplots(
        nrows=n_subplots, ncols=1, figsize=(output_size/dpi_val,
                                            output_size/dpi_val), dpi=dpi_val)

    # Make sure the subplots are in a numpy array (I think it's not if there
    # is only one).
    if not isinstance(these_plotareas, numpy.ndarray):
        these_plotareas = numpy.asarray([these_plotareas])

    # Adjust the plot geometry (margins, etc.) based on plot size.
    if is_bigplot:
        this_figure.subplots_adjust(hspace=0.3, top=0.915)
    else:
        this_figure.subplots_adjust(top=0.85, bottom=0.3, left=0.25, right=0.8)

    # Loop over each segment.
    covering_fractions = [0.] * len(subplot_segment_names)
    for i, segment in enumerate(subplot_segment_names):
        this_plotarea = these_plotareas[i]

        # If this is a large plot, then get the spectrum for this segment.
        # Otherwise, stitch all the segments together (there is only one subplot
        # in this case).  Note that the DQ values are the bits from DQ_WGT in
        # the header.
        if is_bigplot:
            all_wls = cos_spectrum.segments[segment].wavelengths
            all_fls = cos_spectrum.segments[segment].fluxes
            all_flerrs = cos_spectrum.segments[segment].fluxerrs
            all_dqs = cos_spectrum.segments[segment].dqs
            title_addendum = title_addendum
        else:
            if stitched_spectrum is not None:
                try:
                    all_wls = stitched_spectrum["wls"]
                    all_fls = stitched_spectrum["fls"]
                    all_flerrs = stitched_spectrum["flerrs"]
                    all_dqs = stitched_spectrum["dqs"]
                    title_addendum = stitched_spectrum["title"]
                except KeyError as the_error:
                    raise SpecUtilsError("The provided stitched spectrum does"
                                         " not have the expected format,"
                                         " missing key "+str(the_error)+".")

            else:
                raise SpecUtilsError("You must provide a stitched spectrum"
                                     " through the `stitched_spectrum`"
                                     " input argument if creating a small"
                                     " (thumbnail-sized) plot.")

        # Only plot information in the plot title if the plot is large (and
        # therefore sufficient space exists on the plot).
        if is_bigplot:
            this_plotarea.set_title(title_addendum, loc="right", size="small",
                                    color="red")

        # Extract the optimal x-axis plot range from the plot_metrics dict,
        # since we use it a lot.
        optimal_xaxis_range = plot_metrics[i]["optimal_xaxis_range"]

        # Plot the spectrum, but only if valid wavelength ranges for x-axis
        # are returned, otherwise plot a special "Fluxes Are All Zero" plot.
        if all(numpy.isfinite(optimal_xaxis_range)):
            # We plot the spectrum as a regular line for use in
            # calc_covering_fraction, it will be removed later.
            this_line = this_plotarea.plot(all_wls, all_fls, 'b')

            # Update y-axis range, but only adjust the ranges if this isn't
            # an all-zero flux case (and not in debug mode, in which case I want
            # to see the entire y-axis range).
            if not debug:
                this_plotarea.set_ylim(plot_metrics[i]["y_axis_range"])

            covering_fractions[i] = calc_covering_fraction(this_figure,
                                                           these_plotareas, i,
                                                           optimize=optimize)
            # Note: here we remove the line we plotted before, it was only so
            # that calc_covering_fraction would have someting to draw on the
            # canvas and thereby determine which pixels were "blue" (i.e., part
            # of the plotted spectrum vs. background).
            this_line.pop(0).remove()
            # Now we plot the spectrum as a LineCollection so that the
            # transparency will have the desired effect, but, this is not
            # rendered on the canvas inside calc_covering_fraction, hence why we
            # need to plot it both as a regular line first.
            this_collection = this_plotarea.add_collection(
                plot_metrics[i]["line_collection"])

            if covering_fractions[i] > 30.:
                this_collection.set_alpha(0.1)

            # Turn on plot grid lines.
            this_plotarea.grid(True, linestyle="dashed")

            # Add the super title AFTER determining plot transparency (to
            # minimize number of colored pixels).
            if is_bigplot:
                if i == len(subplot_segment_names)-1:
                    this_figure.suptitle(os.path.basename(
                        cos_spectrum.orig_file), fontsize=18, color='r')
                else:
                    this_figure.suptitle(os.path.basename(
                        cos_spectrum.orig_file), fontsize=18, color='white')
                # Uncomment the lines below to include the covering fraction
                # as part of the suptitle.
#                if i == len(subplot_segment_names)-1:
#                    this_figure.suptitle(os.path.basename(
#                            cos_spectrum.orig_file) + ": " +
#                                         ','.join(['{0:6.2f}'.format(y)
#                                                   for y in covering_fractions]
#                                                  ), fontsize=18, color='r')
#                else:
#                    this_figure.suptitle(os.path.basename(
#                            cos_spectrum.orig_file) + ": " +
#                                         ','.join(['{0:6.2f}'.format(y)
#                                                   for y in covering_fractions]
#                                                  ), fontsize=18,
#                                         color='white')
            if debug:
                # Overplot points color-coded based on rejection criteria.
                debug_oplot(this_plotarea, "cos", all_wls, all_fls, all_flerrs,
                            all_dqs, plot_metrics[i]["median_flux"],
                            plot_metrics[i]["median_fluxerr"],
                            flux_scale_factor, fluxerr_scale_factor,
                            plot_metrics[i]["fluxerr_95th"])

                # Overplot regions excluded by the optimal x-axis range as a
                # shaded area.
                this_plotarea.axvspan(numpy.nanmin(all_wls),
                                      optimal_xaxis_range[0],
                                      facecolor="lightgrey", alpha=0.5)
                this_plotarea.axvspan(optimal_xaxis_range[1],
                                      numpy.nanmax(all_wls),
                                      facecolor="lightgrey", alpha=0.5)

                # Overplot the Avoid Regions as a shaded area.
                for region in plot_metrics[i]["avoid_regions"]:
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
                rc('font', size=10)
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
                this_plotarea.set_xlabel(r"Wavelength $(\AA)$", fontsize=16,
                                         color='k')
                this_plotarea.set_ylabel(r"Flux $\mathrm{(erg/s/cm^2\!/\AA)}$",
                                         fontsize=16, color='k')

                # If requested, include the powers of 10 part of the y-axis
                # tickmarks.
                if full_ylabels:
                    this_plotarea.yaxis.set_major_formatter(FormatStrFormatter(
                        '%3.2E'))

        else:
            # Otherwise this is a spectrum that has all zero fluxes, or some
            # other problem, and we make a default plot.  Define the optimal
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
                rc('font', size=10)
                minwl = numpy.nanmin(all_wls)
                maxwl = numpy.nanmax(all_wls)
                this_plotarea.set_xticks([minwl, maxwl])
                this_plotarea.set_xticklabels(this_plotarea.get_xticks(),
                                              rotation=25.)
                this_plotarea.xaxis.set_major_formatter(FormatStrFormatter(
                    "%6.1f"))
                textsize = "small"
                plottext = "Fluxes are \n all 0."
            else:
                # Make sure the font properties go back to normal.
                pyplot.rcdefaults()
                this_plotarea.tick_params(axis='x', labelsize=14)
                this_plotarea.set_xlabel(r"Wavelength $(\AA)$", fontsize=16,
                                         color='k')
                this_plotarea.set_ylabel(r"Flux $\mathrm{(erg/s/cm^2\!/\AA)}$",
                                         fontsize=16, color='k')

                # If requested, include the powers of 10 part of the y-axis
                # tickmarks.
                if full_ylabels:
                    this_plotarea.yaxis.set_major_formatter(FormatStrFormatter(
                        '%3.2E'))

                textsize = "x-large"
                plottext = "Fluxes are all 0."

            # Place the text with the informational message in the center of
            # the plot.
            this_plotarea.text(0.5, 0.5, plottext, horizontalalignment="center",
                               verticalalignment="center",
                               transform=this_plotarea.transAxes, size=textsize)

    # Display or plot to the desired device.
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
