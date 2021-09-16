# DEPRECATED
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os
import numpy

CCDC_source_dir = os.path.abspath(os.path.dirname('__file__'))

sccd_extension = Extension(
    name = "pycold",
    sources = ["pycold.pyx"],
    libraries = ["sccd"],
    library_dirs = [CCDC_source_dir]
)

setup(
    name = "pycold",
    description="python implementation of COntinuous monitoring of Land disturbances algorithm",
    ext_modules = cythonize([sccd_extension]),
    include_dirs=[numpy.get_include()],
    author = "Su Ye",
    author_email = "remotesensingsuy@gmail.com"
)

