# For cython to generate python library
default: pycold

pysccd: setup.py pycold.pyx libsccd.so
	python3 setup.py build_ext --inplace && rm -f pysccd.c && rm -Rf build
	# cp pysccd*.so $(PYTHON_SCCDPROJECT_DIR)
	# cp libsccd.so $(PYTHON_LIB_DIR)
	export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):${CURDIR}
	# cp pysccd.pyx $(PYTHON_PACKAGE_DIR)
	# cp pysccd*.so /home/coloury/Dropbox/Documents/PyCharmProjects/Dissertation
	# cp libsccd.so /home/coloury/miniconda3/envs/sccd/lib

clean:
	rm pycold*.so
	# rm /home/coloury/Dropbox/Documents/PyCharmProjects/Dissertation/pysccd*.so
	# rm /home/coloury/miniconda3/envs/sccd/lib/libsccd.so
	# rm $(PYTHON_SCCDPROJECT_DIR)/pysccd*.so
	# rm libsccd.so
	# rm $(PYTHON_SCCDPROJECT_DIR)/libsccd.so
