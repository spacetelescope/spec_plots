"""
.. module:: __init__

   :synopsis: Used to treat "specutils_stis" as a directory and package.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

from __future__ import unicode_literals
from .stis1dspectrum import STIS1DSpectrum, STISExposureSpectrum
from .get_association_indices import get_association_indices
from .plotspec import plotspec
from .readspec import readspec
