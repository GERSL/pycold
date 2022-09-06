#!/bin/bash
flake8 --count --select=E9,F63,F7,F82 --show-source --statistics ./src/python/pycold
flake8 --count --select=E9,F63,F7,F82 --show-source --statistics ./tests

#flake8 --max-line-length 79 '--ignore=B006,B007,B009,C401,C405,C408,C409,C414,C416,E123,E126,E127,E201,E202,E203,E221,E222,E241,E26 ,E265,E266,E271,E272,E301,E305,E306,E501,EXE001,EXE002,I100,I201,N801,N802,N803,N804,N805,N806,N807,N811,N812,N813,N814,N817,N818,W503,W504,W602' --statistics ./src/python/pycold
