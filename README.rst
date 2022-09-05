PYCOLD
======

|GithubActions| 

.. .. TODO: uncomment these after docs / pypi / coverage are online
.. .. |ReadTheDocs| |Pypi| |Downloads| |Codecov| 

A Python library for COntinuous monitoring of Land Disturbance (COLD) and its extension algorithms at the speed of C language
=============================================================================================================================

The base algorithms were mostly written using C wrapped in Python, and
have been verified with `MATLAB version <https://github.com/GERSL/COLD>`_.

This library provides: 

1. Original COntinuous monitoring of Land Disturbance (COLD): a upgraded CCDC algorithm proposed by Dr.Zhe Zhu for offline satellite-based time-series analysis 
       
2. Stochastic Continuous Change Detection (S-CCD, a near real-time and short-memory implementation of COLD) 
       
3. Object-based COLD (OB-COLD), integrating spatial information into COLD by using a ‘change object’ view 
           
    
1. Pre-work
-----------
   
Clone github repo to your local code directory for the first use:

.. code:: bash

   git clone https://github.com/GERSL/pycold.git

Or you call pull the recent repo if you want to update the existing
pycold repo:

.. code:: bash

   git pull origin devel:devel

2. Installation
---------------

The steps to install this library in development mode are consolidated
into a single script: ``run_developer_setup.sh``.  On debian-based systems,
this will install all of the developer requirements and ensure you are setup
with a working opencv-python-headless and gdal Python modules, as well as other
requirements and then it will compile and install pycold in editable
development mode.


The following is an overview of these details and alternative choices that
could be made.

2.1 Install Required Libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ZLIB, GSL, and OpenMP libraries are required.

For Ubuntu/Debian systems, they can be installed via:

.. code:: bash

   sudo apt-get update
   sudo apt-get install build-essential  -y
   sudo apt-get install zlib1g-dev -y
   sudo apt-get install gfortran -y
   sudo apt-get install libgsl-dev -y

On CentOS systems run:

.. code:: bash

   sudo apt-get install gcc gcc-c++ make  -y
   sudo apt-get install zlib-devel -y
   sudo apt-get install gcc-gfortran -y
   # Yum provides an gsl 1.5, but we need 2.7
   # sudo apt-get install gsl-devel -y
   curl https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz  > gsl.tar.gz && tar xfv gsl.tar.gz && cd gsl-2.7.1 && ./configure --prefix=/usr --disable-static && make && make install

2.2 Compile and Install PYCOLD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Option 1: Install in development mode**

The following instructure assume you are inside a Python virtual environment
(e.g. via conda or pyenv).

The first step is to install GDAL, which does not currently have a binary
distribution on pypi. The following step will install GDAL from a 
`custom pypi server <https://girder.github.io/large_image_wheels>`_ 
containing precompiled wheels. 


.. code:: bash

    # Install GDAL (note-this must be done manually)
    pip install -r requirements/gdal.txt

We must also manually install the Python build dependencies. This is done via:

.. code:: bash

    # In the future the pyproject.toml should take care of this, but 
    # for now, we must do it manually.
    pip install -r requirements/build.txt
   
While the next command will take care of installing all other non-gdal Python
runtime dependencies, if desired these runtime dependencies can be installed manually via:
``pip install -r requirements.txt`` or via invoking
``pip install -r`` on a specific files in the ``requirements`` subdirectory.

Once GDAL and the Pyton build dependences are installed, pycold can be compiled
and installed in development mode via:

.. code:: bash

    export SETUPTOOLS_ENABLE_FEATURES="legacy-editable"
    pip install --no-build-isolation --verbose -e .[headless]

This will install the Python dependencies, then compile the C dependencies
in-place, and finally install the pycold Python package in development mode.
Note specifying "headless" tells pycold which version of opencv to use. The
alternative is "graphics", which is not recommended for fresh builds. Note: you
may need to remove any existing ``_skbuild`` directly if you encounter a build
error.

NOTE: Due to a bug in scikit-build, the editable install link does not point to
the correct path. This can be corrected via:


.. code:: bash

    # Workaround for a scikit-build editable install bug
    REPO_DPATH=$(pwd)
    SITE_DPATH=$(python -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib())")
    PYCOLD_EGG_LINK_FPATH="$SITE_DPATH"/pycold.egg-link
    EASY_INSTALL_FPATH=$SITE_DPATH/easy-install.pth
    echo "$REPO_DPATH/src/python" > "$PYCOLD_EGG_LINK_FPATH"
    echo "../../" >> "$PYCOLD_EGG_LINK_FPATH"
    mv pycold.egg-info ./src/python/
    echo "$REPO_DPATH/src/python" >> "$EASY_INSTALL_FPATH"

Again we note that the above steps and other minor details are consolidated in
the ``run_developer_setup.sh`` script.

**Option 2: Build and install a wheel**

Scikit-build will invoke CMake and build everything. (you may need to
remove any existing ``_skbuild`` directory).

.. code:: bash

   python -m build --wheel .

Then you can pip install the wheel (the exact path will depend on your system
and version of python).

.. code:: bash

   pip install dist/pycold-0.1.0-cp38-cp38-linux_x86_64.whl


You can also use the ``build_wheels.sh`` script to invoke cibuildwheel to
produce portable wheels that can be installed on different than they were built
on. You must have docker and cibuildwheel installed to use this.


**Option 3: build standalone binaries with CMake by itself (recommended
for C development)**

.. code:: bash

   mkdir -p build
   cd build
   cmake ..
   make 

**Option 4: Use a docker image.**

See `dockerfiles/README.rst <dockerfiles/README.rst>`__ for details.

3. Using pycold for pixel-based processing
------------------------------------------

COLD:

.. code:: python

   from pycold import cold_detect
   cold_result = cold_detect(dates, blues, greens, reds, nirs, swir1s, swir2s, thermals, qas)

S-CCD:

.. code:: python

   # require offline processing for the first time 
   from pycold import sccd_detect, sccd_update
   sccd_pack = sccd_detect(dates, blues, greens, reds, nirs, swir1s, swir2s, thermals, qas)

   # then use sccd_pack to do recursive and short-memory NRT update
   sccd_pack_new = sccd_update(sccd_pack, dates, blues, greens, reds, nirs, swir1s, swir2s, thermals, qas)

Q&A
~~~

Q1: Has pycold been verified with original Matlab codes?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Re: yes, multiple rounds of verification have been done. Comparison
based on two testing tiles shows that pycold and Matlab version have
smaller than <2% differences for breakpoint detection and <2%
differences for harmonic coefficients; the accuracy of pycold was also
tested against the same reference dataset used in the original COLD
paper (Zhu et al., 2020), and pycold reached the same accuracy (27%
omission and 28% commission) showing that the discrepancy doesn’t hurt
accuracy. The primary source for the discrepancy is mainly from the
rounding: MATLAB uses float64 precision, while pycold chose float32 to
save the run-time computing memory and boost efficiency.

Q2: how much time for production of a tile-based disturbance map (5000*5000 pixels) using pycold?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Re: I tested it in UCONN HPC environment (200 EPYC7452 cores): for
processing a 40-year Landsat ARD tile (1982-2021), the stacking
typically takes 15 mins; per-pixel COLD processing costs averagely 1
hour; exporting maps needs 7 mins.

4. Citations
------------

If you make use of the algorithms in this repo (or to read more about them),
please cite (/see) the relevant publications from the following list:

`[COLD] <https://www.sciencedirect.com/science/article/am/pii/S0034425719301002>`_ 
Zhu, Z., Zhang, J., Yang, Z., Aljaddani, A. H., Cohen, W. B., Qiu, S., &
Zhou, C. (2020). Continuous monitoring of land disturbance based on
Landsat time series. *Remote Sensing of Environment*, *238*, 111116.

`[S-CCD] <https://www.sciencedirect.com/science/article/pii/S003442572030540X>`_
Ye, S., Rogan, J., Zhu, Z., & Eastman, J. R. (2021). A near-real-time
approach for monitoring forest disturbance using Landsat time series:
Stochastic continuous change detection. *Remote Sensing of Environment*,
*252*, 112167.

[OB-COLD] Ye, S., Zhu, Z., & Cao, G., (2022). Object-based continuous monitoring
of land disturbance. Submitted to Remote Sensing of Environment


.. |Codecov| image:: https://codecov.io/github/GERSL/pycold/badge.svg?branch=devel&service=github
   :target: https://codecov.io/github/GERSL/pycold?branch=devel
.. |Pypi| image:: https://img.shields.io/pypi/v/pycold.svg
   :target: https://pypi.python.org/pypi/pycold
.. |Downloads| image:: https://img.shields.io/pypi/dm/pycold.svg
   :target: https://pypistats.org/packages/pycold
.. |ReadTheDocs| image:: https://readthedocs.org/projects/pycold/badge/?version=latest
    :target: http://pycold.readthedocs.io/en/latest/
.. |GithubActions| image:: https://github.com/GERSL/pycold/actions/workflows/tests.yml/badge.svg?branch=devel
    :target: https://github.com/GERSL/pycold/actions?query=branch%3Amain