"""
.. module:: stitch_components
   :synopsis: Given a COS or STIS spectrum, will stitch each segment/order,
       respectively, into a contiguous array.  Does not handle any overlap in
       wavelength between sections.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from __future__ import division
import os
import sys
from builtins import range
from builtins import str
#--------------------
# External Imports
#--------------------
import numpy
#--------------------
# Package Imports
#--------------------
from spec_plots.utils.specutils.is_bad_dq import is_bad_dq
from spec_plots.utils.specutils.get_flux_stats import get_flux_stats
from spec_plots.utils.specutils.edge_trim import edge_trim
from spec_plots.utils.specutils_cos.cosspectrum import COSSpectrum
from spec_plots.utils.specutils_stis.stis1dspectrum import STISExposureSpectrum
from spec_plots import __version__
#--------------------

#--------------------
def stitch_components(input_exposure, n_consecutive, flux_scale_factor,
                      fluxerr_scale_factor, segment_names=None):
    """
    Given a COSSpectrum or STISExposureSpectrum object, stitches each
    segment/order, respectively, into a contiguous array.

    :param input_exposure: The COS segment or STIS exposure spectrum to stitch.

    :type input_exposure: COSSpectrum or STISExposureSpectrum

    :param n_consecutive: How many consecutive points must pass the test for the
        index to count as the valid start/end of the spectrum?

    :type n_consecutive: int

    :param flux_scale_factor: Max. allowed ratio between the flux and a median
        flux value, used in edge trimming.

    :type flux_scale_factor: float

    :param fluxerr_scale_factor: Max. allowed ratio between the flux uncertainty
        and a median flux uncertainty value, used in edge trimming.

    :type fluxerr_scale_factor: float

    :param segment_names: List of segment names if input_exposure is a
        COSSpectrum object.

    :type segment_names: list

    :returns: dict{numpy array, numpy array, numpy array, str} -- The stitched
        wavelengths, fluxes, flux errors, and an informational plot title in the
        event that all the fluxes had the DQ flag set.  These are packaged in a
        dict.

    :raises: ValueError
    """

    # These lits will contain the stitched spectrum.
    all_wls = []
    all_fls = []
    all_flerrs = []
    all_dqs = []

    # Determine how many pieces there are to stitch, and how to loop through
    # them (different depending on instrument type).
    if isinstance(input_exposure, COSSpectrum):
        if segment_names is not None:
            n_components = len(segment_names)
            loop_iterable = segment_names
            inst_type = "cos"
        else:
            raise ValueError("Must provide a list of segment names for COS"
                             " spectra.")

    elif isinstance(input_exposure, STISExposureSpectrum):
        n_components = len(input_exposure.orders)
        loop_iterable = list(range(n_components))
        inst_type = "stis"

    else:
        raise ValueError("Input must be either a COSSpectrum or"
                         " STISExposureSpectrum object.")

    # Create an array that will keep track whether each component to stitch
    # is filled with bad DQ flags.  The return title will be set to a non-empty
    # string if every component to stitch is filled with bad DQ flags.  Note
    # that for COS the DQ flags are the DQ_WGT bits from the header.
    all_dq_flags = numpy.zeros(n_components)
    return_title = ""

    # Get the wavelengths, fluxes, flux uncertainties, and DQ flags for this
    # component.
    for j_j, j in enumerate(loop_iterable):
        if inst_type == "stis":
            these_wls = input_exposure.orders[j].wavelengths
            these_fls = input_exposure.orders[j].fluxes
            these_flerrs = input_exposure.orders[j].fluxerrs
            these_dqs = input_exposure.orders[j].dqs
            where_bad_dq = numpy.where(is_bad_dq(inst_type, these_dqs))[0]
        elif inst_type == "cos":
            these_wls = input_exposure.segments[j].wavelengths
            these_fls = input_exposure.segments[j].fluxes
            these_flerrs = input_exposure.segments[j].fluxerrs
            these_dqs = input_exposure.segments[j].dqs
            where_bad_dq = numpy.where(is_bad_dq(inst_type, these_dqs))[0]
        n_elems = len(these_wls)

        # Calculate some statistics for this component.
        median_flux, median_fluxerr, fluxerr_95th = get_flux_stats(these_fls,
                                                                   these_flerrs)

        # Find the start and end indices using the edge trimming function.
        (start_index_nodq, end_index_nodq, start_index_withdq,
         end_index_withdq) = edge_trim(inst_type, these_fls, these_flerrs,
                                       these_dqs, n_consecutive, median_flux,
                                       flux_scale_factor, median_fluxerr,
                                       fluxerr_scale_factor, fluxerr_95th)

        # Check if this component has all bad DQ flags.
        if len(where_bad_dq) == n_elems:
            all_dq_flags[j_j] = 1

        # Only append the parts of this order's spectrum that pass the edge
        # trimming, but do NOT trim too much.
        if (start_index_withdq+1) / float(n_elems) <= 0.1:
            start_index_to_use = start_index_withdq
        elif (start_index_nodq+1) / float(n_elems) <= 0.1:
            start_index_to_use = start_index_nodq
        else:
            start_index_to_use = 0

        if abs(end_index_withdq) / float(n_elems) <= 0.1:
            end_index_to_use = end_index_withdq
        elif abs(end_index_nodq) / float(n_elems) <= 0.1:
            end_index_to_use = end_index_nodq
        else:
            end_index_to_use = -1

        # Append the portion of this component to the spectrum we are
        # building up.
        if end_index_to_use == -1:
            all_wls += list(these_wls[start_index_to_use:])
            all_fls += list(these_fls[start_index_to_use:])
            all_flerrs += list(these_flerrs[start_index_to_use:])
            all_dqs += list(these_dqs[start_index_to_use:])
        else:
            all_wls += list(these_wls[start_index_to_use:end_index_to_use+1])
            all_fls += list(these_fls[start_index_to_use:end_index_to_use+1])
            all_flerrs += list(these_flerrs[start_index_to_use:end_index_to_use+
                                            1])
            all_dqs += list(these_dqs[start_index_to_use:end_index_to_use+1])

    # If every single order had all DQ flags, then we want to print a warning
    # on the plot. """
    if sum(all_dq_flags) == n_components:
        if inst_type == 'stis':
            return_title = "Warning: All fluxes have bad DQ."
        elif inst_type == 'cos':
            return_title = "Warning: All fluxes have bad DQ_WGT."

    # Convert these to numpy arrays.
    # <DEVEL> Should these just be numpy arrays to begin with? </DEVEL>
    all_wls = numpy.asarray(all_wls)
    all_fls = numpy.asarray(all_fls)
    all_flerrs = numpy.asarray(all_flerrs)
    all_dqs = numpy.asarray(all_dqs)

    # Make sure the spectrum is sorted in wavelength.
    # <DEVEL> Not sure how to handle components that might overlap in
    # wavelength space, but are not monotonically increasing/decreasing as each
    # component is stitched. </DEVEL>
    sorted_indexes = numpy.argsort(all_wls)
    all_wls = all_wls[sorted_indexes]
    all_fls = all_fls[sorted_indexes]
    all_flerrs = all_flerrs[sorted_indexes]
    all_dqs = all_dqs[sorted_indexes]

    return {"wls":all_wls, "fls":all_fls, "flerrs":all_flerrs, "dqs":all_dqs,
            "title":return_title}
#--------------------
