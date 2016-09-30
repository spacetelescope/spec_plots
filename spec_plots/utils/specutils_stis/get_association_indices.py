"""
.. module:: get_association_indices
   :synopsis: Determines which indices in the array of STIS associations to
       plot in the preview plots.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from __future__ import division
from builtins import range
#--------------------
# Package Imports
#--------------------
from spec_plots import __version__

#--------------------
def get_association_indices(associations):
    """
    Given a list of associations, determines the indices that will be plotted.
    If n_associations > 3, this is the first, middle, and last associations,
    otherwise it's all of them.

    :param associations: All of the associations for this STIS spectrum.

    :type associations: list

    :returns: list -- The indices of the list that will be plotted.
    """

    n_associations = len(associations)
    if n_associations <= 3:
        subplot_indices = list(range(n_associations))
    else:
        midindex = int(round(n_associations/2.))
        subplot_indices = [0, midindex, n_associations-1]

    return subplot_indices
#--------------------
