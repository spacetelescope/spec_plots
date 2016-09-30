"""
.. module:: __init__

   :synopsis: Used to treat "specutils" as a directory and package.

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
from spec_plots.utils.specutils.specutilserror import SpecUtilsError
from spec_plots.utils.specutils.avoidregion import (AvoidRegion,
                                                    generate_avoid_regions)
from spec_plots.utils.specutils.is_bad_dq import is_bad_dq
from spec_plots.utils.specutils.stitch_components import stitch_components
from spec_plots.utils.specutils.calc_covering_fraction import (
    calc_covering_fraction)
from spec_plots.utils.specutils.calc_plot_metrics import calc_plot_metrics
from spec_plots.utils.specutils.debug_oplot import debug_oplot
from spec_plots.utils.specutils.dq_has_flag import dq_has_flag
from spec_plots.utils.specutils.edge_trim import (edge_trim,
                                                  _set_plot_xrange_test)
from spec_plots.utils.specutils.get_flux_stats import get_flux_stats
from spec_plots.utils.specutils.set_plot_xrange import set_plot_xrange
from spec_plots.utils.specutils.set_plot_yrange import set_plot_yrange
