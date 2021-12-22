# PYCOLD

# PYthon version of COntinuous monitoring of Land Disturbance algorithm (V0.1)
####  History
- Version 0.2: add pip-installable functionality by using scikit-build
- Version 0.1: a Python package based on C library of COLD. The verification with MATLAB COLD shows that  C-COLD has smaller than <2% differences for breakpoint detection and <2% differences for harmonic coefficents

## 1. Pre-work: clone github repo to your local directory
### Clone or pull the latest pycold git repo
```
git clone https://gitlab.kitware.com/smart/pycold.git
```

## 2. installation
### 2.1 install required libraries
sudo apt-get install build-essential zlib1g-dev
sudo apt-get install libgsl-dev
sudo apt-get install gfortran

### 2.2 install pycold
**Option 1: Build and install a wheel** 
(Scikit-build will invoke CMake and build everythinig)
```
python setup.py bdist_wheel
```
Then you can pip install the wheel
```
pip install dist/pycold-0.1.0-cp38-cp38-linux_x86_64.whl
```

**Option 2: Install in development mode**
```
pip install -r requirements.txt
pip install -e .
```

**Option 3: build standalone binaries with CMake by itself (recommended for C development)**
```
mkdir -p build
cd build
cmake ..
make 
```

### 3. Using pycold
```python
>>> import pycold
>>> cold_result = pycold(dates, blues, greens, reds, nirs, swir1s, swir2s, thermals, qas)
```

## Please refer to the example in /tool/notebook/pycold_example.ipynb
