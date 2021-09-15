# PYCOLD

# PYthon version of COntinuous monitoring of Land Disturbance algorithm (V0.1)
### Authors: Su Ye (remotesensingsuy@gmail.com)
####  History
- Version 0.1: a Python package based on C library of COLD. The verification with MATLAB COLD shows that  C-COLD has smaller than <2% differences for breakpoint detection and <2% differences for harmonic coefficents

## 1. Pre-work: clone github repo to your local directory
### Clone or pull the latest pycold git repo
```
git clone https://gitlab.kitware.com/smart/pycold.git
```

## 2. installation
### 2.1 Install dependencies
(Ubuntu)
```
sudo apt-get install build-essential
sudo apt-get install libgsl-dev
sudo apt-get install gfortran
```
(MAC)
```
brew install gsl
brew install gcc
```
### 2.2 local installation
It's highly recommended to do all your development & testing in a virtual environment.
```
conda create -n pycoldenv python=3.6 # create a new environmnet
conda activate pycoldenv
conda install --file=/YOUR_PYCOLD_DIRECTORY/tool/python/requirements.txt # install required packages
```

Then run a compilation file which combines the make procedures for original C codes and cython codes. 
```
./YOUR_PYCOLD_DIRECTORY/compilation_procedure
```
Just ignore the warning messages, and if successful, you can find two executable files called 'libsccd.so' and 'pycold.cpython-36m-darwin.so' (this name is for mac; it should be different affix under different system) in the directory /YOUR_PYCOLD_DIRECTORY/tool/python

### 2.3 Move executable files into your python project directory
Copy pycold.cython-36m-darwin.so and libsccd.so into your python project directory
Ok, I admit that this step is ''stupied''. The future work should change it to a pip-installable package.


## 3. running pycold
```python
import pycold

here are the codes for reading your dataset..

cold_result = pycold.pycold(dates, blue, green, red, nir, swir1, swir2, thermal, qa)

```
A full example of running pycold for an example dataset can be found in /YOUR_PYCOLD_DIRECTORY/tool/python/COLDexample_fromcsv.py


# New Cmake Build Instructions:


```
# Install required libraries
sudo apt-get install build-essential
sudo apt-get install libgsl-dev
sudo apt-get install gfortran

# Scikit-build will invoke CMake and build everything
python setup.py build

# Or you can build with CMake by itself
mkdir -p build
cd build
cmake ..
```
