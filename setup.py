from setuptools import setup

setup(name="spec_plots",
      version="1.31",
      description="Create preview plots of HST spectra.",
      classifiers=["Programming Language :: Python :: 2.7"],
      url="https://github.com/openSAIL/spec_plots",
      author="Scott W. Fleming",
      author_email="fleming@stsci.edu",
      license="MIT",
      packages=["spec_plots"],
      install_requires=["astropy", "matplotlib", "numpy"],
      zip_safe=False)
