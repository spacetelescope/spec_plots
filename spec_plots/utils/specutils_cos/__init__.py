"""
.. module:: __init__

   :synopsis: Used to treat "specutils_cos" as a directory and package.

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
from spec_plots.utils.specutils_cos.cosspectrum import COSSpectrum
from spec_plots.utils.specutils_cos.check_segments import check_segments
from spec_plots.utils.specutils_cos.extract_subspec import extract_subspec
from spec_plots.utils.specutils_cos.plotspec import plotspec
from spec_plots.utils.specutils_cos.readspec import readspec
from spec_plots.utils.specutils_cos.get_segment_names import get_segment_names
from spec_plots.utils.specutils_cos.make_fits import make_fits
