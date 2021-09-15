# from os.path import exists
# from os.path import join
# from os.path import dirname
# from setuptools import find_packages
from skbuild import setup


setup(
    package_dir={'': 'tool/python'},
    name="pycold",
    description="python implementation of COntinuous monitoring of Land disturbances algorithm",
    # ext_modules = cythonize([sccd_extension]),
    # include_dirs=[numpy.get_include()],
    author="Su Ye",
    author_email="remotesensingsuy@gmail.com"
)
