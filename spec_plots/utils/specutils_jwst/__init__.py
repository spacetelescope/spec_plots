"""
.. module:: __init__

   :synopsis: Used to treat "specutils_jwst" as a directory and package.

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
from spec_plots.utils.specutils_jwst.plotspec import plotspec
from spec_plots.utils.specutils_jwst.readspec import readspec
