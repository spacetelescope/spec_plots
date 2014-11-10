from setuptools import setup

setup(name="spec-plots",
      version="1.25",
      description="Create preview plots of HST spectra.",
      classifiers=["Programming Language :: Python :: 2.7"],
      url="https://github.com/openSAIL/spec-plots",
      author="Scott W. Fleming",
      author_email="fleming@stsci.edu",
      license="MIT",
      packages=["spec-plots"],
      install_requires=["astropy", "matplotlib", "numpy"],
      zip_safe=False)
