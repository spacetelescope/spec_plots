"""
.. module:: __init__

   :synopsis: Used to treat "specutils_stis" as a directory and package.

.. moduleauthor:: Rob Swaters <rswaters@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
#--------------------
# Package Imports
#--------------------
from spec_plots import __version__
from spec_plots.utils.specutils_hasp.haspspectrum import HASPSpectrum
from spec_plots.utils.specutils_hasp.plotspec import plotspec
from spec_plots.utils.specutils_hasp.readspec import readspec
