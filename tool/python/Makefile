# DEPRECATED
#
# For cython to generate python library
default: pycold

pycold: setup.py pycold.pyx libsccd.so
	python3 setup.py build_ext --inplace && rm -f pycold.c && rm -Rf build
	export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):${CURDIR}

clean:
	rm pycold*.so
