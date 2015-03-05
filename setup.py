from setuptools import setup

setup(name="spec_plots",
      version="1.32.0",
      description="Create preview plots of HST spectra.",
      classifiers=["Programming Language :: Python :: 2.7"],
      url="https://github.com/openSAIL/spec_plots",
      author="Scott W. Fleming",
      author_email="fleming@stsci.edu",
      license="MIT",
      packages=["spec_plots", "spec_plots.utils", "spec_plots.utils.specutils", "spec_plots.utils.specutils_cos", "spec_plots.utils.specutils_stis"],
      install_requires=["astropy>=0.4.1", "matplotlib>=1.4.1", "numpy>=1.9.1"],
      entry_points={"console_scripts" : ["make_hst_spec_previews = spec_plots.__main__:main"]},
      zip_safe=False)
