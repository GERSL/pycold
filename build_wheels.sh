#!/bin/bash
__doc__="""
SeeAlso:
    pyproject.toml
"""
if ! which docker ; then 
    echo "Missing requirement: docker. Please install docker before running build_wheels.sh"
    exit 1
fi
if ! which cibuildwheel ; then 
    echo "The cibuildwheel module is not installed. Please pip install cibuildwheel before running build_wheels.sh"
    exit 1
fi
cibuildwheel --config-file pyproject.toml --platform linux --arch x86_64 
