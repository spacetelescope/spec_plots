__version__ = '1.31'

"""
.. module:: plotspec

   :synopsis: Creates preview plots for the provided STIS spectrum.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as pyplot
from matplotlib import rc
import numpy
import os
import sys

""" <DEVEL> Note that this hack to make it so that the user can import `plotspec` directly as a module or run it from the command line as __main__ has the side effect of importing this module twice, despite my best efforts to work around it.  I don't think it will be a major issue, but worth thinking about in the future. </DEVEL> """
if __package__ is None:
    specutils_cos_dir = os.path.dirname(os.path.abspath(__file__))
    utils_dir = os.path.dirname(specutils_cos_dir)
    parent_dir = os.path.dirname(utils_dir)
    sys.path.insert(1, parent_dir)
    import utils
    __package__ = str("utils.specutils")
    __name__ = str(__package__+"."+__name__)

from ..specutils.specutilserror import SpecUtilsError
from ..specutils.debug_oplot import debug_oplot

#--------------------

def plotspec(stis_spectrum, association_indices, stitched_spectra, output_type, output_file, n_consecutive, flux_scale_factor, fluxerr_scale_factor, plot_metrics, dpi_val=96., output_size=1024, debug=False, full_ylabels=False):
    """
    Accepts a STIS1DSpectrum object from the READSPEC function and produces preview plots.

    :param stis_spectrum: STIS spectrum as returned by READSPEC.

    :type stis_spectrum: STIS1DSpectrum

    :param association_indices: The indices of the original list of associations for this spectrum that will be plotted.  Used in the informational plot title, e.g., `Association <x> of <y>.`, where <x> is the association index and <y> is the total number of associations in this spectrum.

    :type association_indices: list

    :param stitched_spectra: List of order-stitched spectra for each association, as created in `utils.specutils.stitch_components()`.

    :type stitched_spectra: list

    :param output_type: What kind of output to make?

    :type output_type: str

    :param output_file: Name of output file (including full path).

    :type output_file: str

    :param n_consecutive: How many consecutive points must pass the test for the index to count as the valid start/end of the spectrum?  Default = 20.

    :type n_consecutive: int

    :param flux_scale_factor: Max. allowed ratio between the flux and a median flux value, used in edge trimming.  Default = 10.

    :type flux_scale_factor: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty and a median flux uncertainty value, used in edge trimming.  Default = 5.

    :type fluxerr_scale_factor: float

    :param plot_metrics: Collection of plot metrics (flux statistics, axis ranges, etc.) to use when making the plots.  These are computed using `utils.specutils.calc_plot_metrics()`.
    
    :type plot_metrics: list

    :param dpi_val: The DPI value of your device's monitor.  Affects the size of the output plots.  Default = 96. (applicable to most modern monitors).
    
    :type dpi_val: float

    :param output_size: Size of plot in pixels (plots are square in dimensions).  Defaults to 1024.

    :param output_size: int

    :param debug: Should the output plots include debugging information (color-coded data points based on rejection criteria, shaded exclude regions)?  Default = False.

    :type debug: bool

    :param full_ylabels: Should the y-labels contain the full values (including the power of 10 in scientific notation)?  Default = False.

    :type full_ylabels: bool

    :raises: OSError, utils.specutils.SpecUtilsError

    .. note::

         This function assumes a screen resolution of 96 DPI in order to generate plots of the desired sizes.  This is because matplotlib works in units of inches and DPI rather than pixels.
    """

    """ Make sure the plot size is set to an integer value. """
    if output_size is not None:
        if not isinstance(output_size, int):
            output_size = int(round(output_size))

    """ Make sure the output path exists, if not, create it. """
    if output_type != 'screen':
        if not os.path.isdir(os.path.dirname(output_file)):
            try:
                os.mkdir(os.path.dirname(output_file))
            except OSError as this_error:
                if this_error.errno == 13: 
                    sys.stderr.write("*** MAKE_HST_SPEC_PREVIEWS ERROR: Output directory could not be created, "+repr(this_error.strerror)+"\n")
                    exit(1)
                else:
                    raise

    """ If the figure size is large, then plot up to three associations, otherwise force only one association on the plot. """
    n_associations = len(stis_spectrum.associations)
    n_subplots = len(stitched_spectra)
    if output_size > 128:
        is_bigplot = True
    else:
        is_bigplot = False

    """ Start plot figure. """
    this_figure, these_plotareas = pyplot.subplots(nrows=n_subplots, ncols=1, figsize=(output_size/dpi_val, output_size/dpi_val), dpi=dpi_val)

    """ Make sure the subplots are in a numpy array (I think by default it is not if there is only one). """
    if not isinstance(these_plotareas, numpy.ndarray):
        these_plotareas = numpy.asarray([these_plotareas])

    """ Adjust the plot geometry (margins, etc.) based on plot size. """
    if is_bigplot:
        this_figure.subplots_adjust(hspace=0.3,top=0.915)
        this_figure.suptitle(os.path.basename(stis_spectrum.orig_file), fontsize=18, color=r'r')
    else:
        this_figure.subplots_adjust(top=0.85,bottom=0.3,left=0.25,right=0.8)

    """ Loop over each association. """
    for i in xrange(len(stitched_spectra)):
        this_plotarea = these_plotareas[i]

        """ Get the wavelengths, fluxes, flux uncertainties, and data quality flags out of the stitched spectrum for this association. """
        try:
            all_wls = stitched_spectra[i]["wls"]
            all_fls = stitched_spectra[i]["fls"]
            all_flerrs = stitched_spectra[i]["flerrs"]
            all_dqs = stitched_spectra[i]["dqs"]
            title_addendum = stitched_spectra[i]["title"]
        except KeyError as the_error:
            raise SpecUtilsError("The provided stitched spectrum does not have the expected format, missing key "+str(the_error)+".")

        if is_bigplot:
            this_plotarea.set_title(title_addendum, loc="right", size="small", color="red")
        if n_associations > 1 and is_bigplot:
            this_plotarea.set_title("Association "+str(association_indices[i]+1)+"/"+str(n_associations), loc="center", size="small", color="black")

        """ Extract the optimal x-axis plot range from the plot_metrics dict, since we use it a lot. """
        optimal_xaxis_range = plot_metrics[i]["optimal_xaxis_range"]

        """ Plot the spectrum, but only if valid wavelength ranges for x-axis are returned, otherwise plot a special "Fluxes Are All Zero" plot. """
        if all(numpy.isfinite(optimal_xaxis_range)):
            """ Plot the spectrum, turn on plot grid lines. """
            this_plotarea.plot(all_wls, all_fls, 'b')
            this_plotarea.grid(True)

            """ Get the values of the x- and y-axes plot limits that *would* have been used by pyplot automatically.  We still use these x-axis plot ranges so that plots from the same instrument configuration have the same x-axis range. """
            pyplot_xrange = this_plotarea.get_xlim()
            pyplot_yrange = this_plotarea.get_ylim()

            if debug:
                """ Overplot points color-coded based on rejection criteria. """
                debug_oplot(this_plotarea, "stis", all_wls, all_fls, all_flerrs, all_dqs, plot_metrics[i]["median_flux"], plot_metrics[i]["median_fluxerr"], flux_scale_factor, fluxerr_scale_factor, plot_metrics[i]["fluxerr_95th"])

                """ Overplot the x-axis edges that are trimmed to define the y-axis plot range as a shaded area. """
                this_plotarea.axvspan(numpy.nanmin(all_wls), optimal_xaxis_range[0],facecolor="lightgrey",alpha=0.5)
                this_plotarea.axvspan(optimal_xaxis_range[1], numpy.nanmax(all_wls),facecolor="lightgrey",alpha=0.5)

                """ Overplot the avoid regions in a light color as a shaded area. """
                for ar in plot_metrics[i]["avoid_regions"]:
                    this_plotarea.axvspan(ar.minwl,ar.maxwl,facecolor="lightgrey",alpha=0.5)

            """ This is where we ensure the x-axis range is set to the pyplot-determined x-axis range, rather than using the optimum x-axis range.  This is done so that all the plots for a similar instrument setting will have the same starting and ending plot values. """
            this_plotarea.set_xlim(pyplot_xrange)

            """ Only use two tick labels (min and max wavelengths) for thumbnails, because there isn't enough space otherwise. """
            if not is_bigplot:
                rc('font', size=10)
                minwl = numpy.nanmin(all_wls) ; maxwl = numpy.nanmax(all_wls)
                this_plotarea.set_xticks([minwl, maxwl])
                this_plotarea.set_xticklabels(this_plotarea.get_xticks(),rotation=25.)
                this_plotarea.xaxis.set_major_formatter(FormatStrFormatter("%6.1f"))
            else:
                """ Make sure the font properties go back to normal. """
                pyplot.rcdefaults()
                this_plotarea.tick_params(axis='x', labelsize=14)
                this_plotarea.set_xlabel(r"Wavelength $(\AA)$", fontsize=16, color='k')
                this_plotarea.set_ylabel(r"Flux $\mathrm{(erg/s/cm^2\!/\AA)}$", fontsize=16, color='k')

                """ If requested, include the powers of 10 part of the y-axis tickmarks. """
                if full_ylabels:
                    this_plotarea.yaxis.set_major_formatter(FormatStrFormatter('%3.2E'))

            """ Update y-axis range, but only adjust the ranges if this isn't an all-zero flux case (and not in debug mode, in which case I want to see the entire y-axis range). """
            if not debug:
                this_plotarea.set_ylim(plot_metrics[i]["y_axis_range"])

        else:
            """ Otherwise this is a spectrum that has all zero fluxes, or some other problem, and we make a default plot.  Define the optimal x-axis range to span the original spectrum. """
            optimal_xaxis_range = [numpy.nanmin(all_wls),numpy.nanmax(all_wls)]
            this_plotarea.set_xlim(optimal_xaxis_range)

            """ Make the plot background grey to distinguish that this is a `special` plot.  Turn off y-tick labels. """
            this_plotarea.set_axis_bgcolor("lightgrey")
            this_plotarea.set_yticklabels([])

            """ Configure the plot units, text size, and other markings based on whether this is a large or thumbnail-sized plot. """
            if not is_bigplot:
                rc('font', size=10)
                minwl = numpy.nanmin(all_wls) ; maxwl = numpy.nanmax(all_wls)
                this_plotarea.set_xticks([minwl, maxwl])
                this_plotarea.set_xticklabels(this_plotarea.get_xticks(),rotation=25.)
                this_plotarea.xaxis.set_major_formatter(FormatStrFormatter("%6.1f"))
                textsize = "small"
                plottext = "Fluxes are \n all 0."
            else:
                """ Make sure the font properties go back to normal. """
                pyplot.rcdefaults()
                this_plotarea.tick_params(axis='x', labelsize=14)
                this_plotarea.set_xlabel(r"Wavelength $(\AA)$", fontsize=16, color='k')
                this_plotarea.set_ylabel(r"Flux $\mathrm{(erg/s/cm^2\!/\AA)}$", fontsize=16, color='k')

                """ If requested, include the powers of 10 part of the y-axis tickmarks. """
                if full_ylabels:
                    this_plotarea.yaxis.set_major_formatter(FormatStrFormatter('%3.2E'))

                textsize = "x-large"
                plottext = "Fluxes are all 0."

            """ Place the text with the informational message in the center of the plot. """
            this_plotarea.text(0.5,0.5,plottext,horizontalalignment="center",verticalalignment="center",transform=this_plotarea.transAxes,size=textsize)

    """ Display or plot to the desired format. """
    if output_type != "screen":
        """ Deconstruct output file to include plot size information. """
        output_splits = os.path.split(output_file)
        file_splits = os.path.splitext(output_splits[1])
        revised_output_file = output_splits[0]+os.path.sep+file_splits[0]+'_{0:04d}'.format(output_size)+file_splits[1]

        """ Save figure. """
        this_figure.savefig(revised_output_file, format=output_type, dpi=dpi_val)

    elif output_type == "screen":
        pyplot.show()

#--------------------
