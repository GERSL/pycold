#!/bin/bash
flake8 --count --select=E9,F63,F7,F82 --show-source --statistics src/python/pycold
flake8 --count --select=E9,F63,F7,F82 --show-source --statistics ./tests