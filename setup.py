"""
.. module:: setup

   :synopsis: This script is used to setup the pip packages.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

from __future__ import absolute_import
from setuptools import setup
from spec_plots import __version__

setup(name="spec_plots",
      version=__version__,
      description="Create preview plots of HST or JWST spectra.",
      classifiers=["Programming Language :: Python :: 3"],
      url="https://github.com/spacetelescope/spec_plots",
      author="Scott W. Fleming",
      author_email="fleming@stsci.edu",
      license="MIT",
      packages=["spec_plots", "spec_plots.utils", "spec_plots.utils.specutils",
                "spec_plots.utils.specutils_cos",
                "spec_plots.utils.specutils_hasp",
                "spec_plots.utils.specutils_jwst",
                "spec_plots.utils.specutils_stis"],
      install_requires=["astropy>=5.2.2", "matplotlib>=3.7.1", "numpy>=1.24.3",
                        "future>=0.18.3"],
      entry_points={"console_scripts" :
                    ["make_hst_spec_previews = spec_plots.__main__:main",
                     "make_jwst_spec_previews = spec_plots.__main_jwst__:main"
                    ]},
      zip_safe=False)
