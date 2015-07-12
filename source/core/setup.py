#!/usr/bin/python

"""Run the following command in a terminal to compile module:

  $ python setup.py build_ext --inplace

"""

from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'Path metrics',
  ext_modules = cythonize("path_metrics.pyx"),
  extra_compile_args=["-O3", "-march=native", "-ffast-math"],
)