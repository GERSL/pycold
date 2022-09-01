#!/bin/bash
__doc__="
Script that does its best to ensure a robust installation of pycold in
development mode.
"

if [[ "$VIRTUAL_ENV" == "" ]]; then
    echo "WARNING: NOT INSIDE OF A Python VIRTUAL_ENV. This script may not run correctly"
fi 

apt_ensure(){
    __doc__="
    Checks to see if the packages are installed and installs them if needed.

    The main reason to use this over normal apt install is that it avoids sudo
    if we already have all requested packages.

    Args:
        *ARGS : one or more requested packages 

    Example:
        apt_ensure git curl htop 

    Ignore:
        REQUESTED_PKGS=(git curl htop) 
    "
    # Note the $@ is not actually an array, but we can convert it to one
    # https://linuxize.com/post/bash-functions/#passing-arguments-to-bash-functions
    ARGS=("$@")
    MISS_PKGS=()
    HIT_PKGS=()
    for PKG_NAME in "${ARGS[@]}"
    do
        #apt_ensure_single $EXE_NAME
        RESULT=$(dpkg -l "$PKG_NAME" | grep "^ii *$PKG_NAME")
        if [ "$RESULT" == "" ]; then 
            echo "Do not have PKG_NAME='$PKG_NAME'"
            # shellcheck disable=SC2268,SC2206
            MISS_PKGS=(${MISS_PKGS[@]} "$PKG_NAME")
        else
            echo "Already have PKG_NAME='$PKG_NAME'"
            # shellcheck disable=SC2268,SC2206
            HIT_PKGS=(${HIT_PKGS[@]} "$PKG_NAME")
        fi
    done

    if [ "${#MISS_PKGS}" -gt 0 ]; then
        sudo apt install -y "${MISS_PKGS[@]}"
    else
        echo "No missing packages"
    fi
}

###  ENSURE DEPENDENCIES ###

# If on debian/ubuntu ensure the dependencies are installed 
if [[ "$(command -v apt)" != "" ]]; then
    apt_ensure build-essential zlib1g-dev libgsl-dev gfortran
else
    echo "
    WARNING: Check and install of system packages is currently only supported
    on Debian Linux. You will need to verify that ZLIB, GSL, OpenMP are
    installed before running this script.
    "
fi



###  CLEANUP PREVIOUS INSTALLS ###


# Reference: https://stackoverflow.com/questions/59895/bash-script-dir
REPO_DPATH=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Check that specific files exist to indicate we aren't in the wrong place.
# Exit if we are.
echo "REPO_DPATH = $REPO_DPATH"
if [[ ! -f "$REPO_DPATH/setup.py" ]] || [[ ! -f "$REPO_DPATH/src/python/pycold/__init__.py" ]]; then
    echo "NOT RUNNING FROM THE CORRECT DIRECTORY. CANNOT PREFORM POST-INSTALL HACKS. EXITING"
    set -e
    false
fi

if [ -d "./_skbuild" ]; then
    # Cleanup the scikit build directory
    echo "Removing old _skbuild directory"
    rm -rf _skbuild
else
    echo "Detected clean state (no _skbuild dir)"
fi

# Also clean up any shared libraries
rm src/python/pycold/_pycold_cython.*.so
# Clean up old egg links and such
rm -rf src/python/pycold.egg-info
rm -rf pycold.egg-info


# There editable install seems ot have bugs. This is an attempt to fix them.
SITE_DPATH=$(python -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib())")
PYCOLD_EGG_LINK_FPATH="$SITE_DPATH"/pycold.egg-link
echo "PYCOLD_EGG_LINK_FPATH = $PYCOLD_EGG_LINK_FPATH"

# Need to get rid of the easy install entry if it exists
EASY_INSTALL_FPATH=$SITE_DPATH/easy-install.pth
if cat "$EASY_INSTALL_FPATH" | grep "$REPO_DPATH" ; then
    echo "Detected pycold in easy install path. Removing it before we reinstall."
    grep -v "$REPO_DPATH" "$EASY_INSTALL_FPATH" > tmpfile && mv tmpfile "$EASY_INSTALL_FPATH"
else
    echo "Easy install pth seems clean"
fi


###  COMPILE STEP ###
#
# Compile and install PyCold in development mode.
# 
# It is important to have no-build-isolation to work
# around an issue with finding numpy headers.
# 
###
pip install --no-build-isolation --verbose -e . 


# This script will analyze the editable install and detect if anything is wrong
# It require some dependencies so it is commented by default
DEBUG_EDITABLE_INSTALL=0
if [[ "$DEBUG_EDITABLE_INSTALL" == "1" ]]; then
    python dev/setup.py mwe analyize --mod_name=pycold --repo_dpath="."
fi


###  HACKY FIXUP STEP ###
#
# Fix the contents of the egg-link in site-packages. See
# [EggFormatDetails]_ for egg format details.
# 
# References:
#     .. [EggFormatDetails] https://svn.python.org/projects/sandbox/trunk/setuptools/doc/formats.txt
###
echo "PYCOLD_EGG_LINK_FPATH = $PYCOLD_EGG_LINK_FPATH"
cat "$PYCOLD_EGG_LINK_FPATH"
echo "SITE_DPATH = $SITE_DPATH"
echo "$REPO_DPATH/src/python
../../" > "$PYCOLD_EGG_LINK_FPATH"


###
# HACK #2 - Move the .egg-info folder into the correct location
###
# Notes about egg-link files:
mv pycold.egg-info ./src/python/


# Now fixup the easy-install.pth
echo "Detected pycold in easy install path. Need to remove it"
grep -v "$REPO_DPATH" "$EASY_INSTALL_FPATH" > tmpfile && mv tmpfile "$EASY_INSTALL_FPATH"
echo "$REPO_DPATH/src/python" >> "$EASY_INSTALL_FPATH"


## Quick tests that the install worked
echo "Quick pycold tests to verify the install:"
echo "Pycold Version: $(python -c 'import pycold; print(pycold.__version__)')"
echo "Python Package Location: $(python -c 'import pycold; print(pycold.__file__)')"
echo "Compiled Cython Module: $(python -c 'import pycold; print(pycold._pycold_cython.__file__)')"
