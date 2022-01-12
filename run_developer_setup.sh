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
            MISS_PKGS=(${MISS_PKGS[@]} "$PKG_NAME")
        else
            echo "Already have PKG_NAME='$PKG_NAME'"
            HIT_PKGS=(${HIT_PKGS[@]} "$PKG_NAME")
        fi
    done

    if [ "${#MISS_PKGS}" -gt 0 ]; then
        sudo apt install -y "${MISS_PKGS[@]}"
    else
        echo "No missing packages"
    fi
}


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


if [ -d "./_skbuild" ]; then
    # Cleanup the scikit build directory
    echo "Removing old _skbuild directory"
    rm -rf _skbuild
else
    echo "Detected clean state (no _skbuild dir)"
fi

# Compile and install PyCold in development mode.
# pip install --verbose -e . 
pip install --no-build-isolation --verbose -e .

# Quick tests that the install worked
echo "Pycold Version: $(python -c 'import pycold; print(pycold.__version__)')"
echo "Python Package Location: $(python -c 'import pycold; print(pycold.__file__)')"
echo "Compiled Cython Module: $(python -c 'import pycold; print(pycold._pycold_cython.__file__)')"
