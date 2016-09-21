"""
.. module:: __init__

   :synopsis: Used to treat "specutils_cos" as a directory and package.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""
from __future__ import unicode_literals
from .cosspectrum import COSSpectrum
from .check_segments import check_segments
from .extract_subspec import extract_subspec
from .plotspec import plotspec
from .readspec import readspec
from .get_segment_names import get_segment_names
