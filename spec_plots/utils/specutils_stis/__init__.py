"""
.. module:: __init__

   :synopsis: Used to treat "specutils_stis" as a directory and package.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
#--------------------
# Package Imports
#--------------------
from spec_plots import __version__
from spec_plots.utils.specutils_stis.stis1dspectrum import (
    STIS1DSpectrum, STISExposureSpectrum)
from spec_plots.utils.specutils_stis.get_association_indices import (
    get_association_indices)
from spec_plots.utils.specutils_stis.plotspec import plotspec
from spec_plots.utils.specutils_stis.readspec import readspec
