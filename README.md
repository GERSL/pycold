# PYCOLD

# A Python library for COntinuous monitoring of Land Disturbance (COLD) and its extension algorithmsat the speed of C language
The base algorithms were mostly written using C wrapped in python, and have been verified with MATLAB version (if it has). This library provides:
  1. Original COntinuous monitoring of Land Disturbance (COLD): a upgraded CCDC algorithm proposed by Dr.Zhe Zhu for offline satellite-based time-series analysis
  2. Stochastic Continuous Change Detection (S-CCD, a near real-time implementation of COLD)
  3. Object-based COLD (OB-COLD, integrating spatial information into COLD by using a 'change object' view)
## 1. Pre-work: clone github repo to your local directory
### Clone or pull the latest pycold git repo

## 2. installation

The steps to install this library in development mode are consolidated into a
single script: `run_developer_setup.sh`. The following is an overview of these
details and alternative choices that could be made.

### 2.1 install required libraries

The ZLIB, GSL, and OpenMP libraries are required.

For Ubuntu/Debian systems, they can be installed via:

```bash
sudo apt-get install build-essential zlib1g-dev -y
sudo apt-get install libgsl-dev -y
sudo apt-get install gfortran -y
```

### 2.2 install pycold
**Option 1: Build and install a wheel** 

Scikit-build will invoke CMake and build everything. (you may need to remove
any existing `_skbuild` directory).

```bash
python setup.py bdist_wheel
```

Then you can pip install the wheel

```bash
pip install dist/pycold-0.1.0-cp38-cp38-linux_x86_64.whl
```

**Option 2: Install in development mode**

```bash
pip install -r requirements.txt
bash run_developer_setup.sh
```

**Option 3: build standalone binaries with CMake by itself (recommended for C development)**

```bash
mkdir -p build
cd build
cmake ..
make 
```

### 3. Using pycold
```python
>>> from pycold import cold_detect
>>> cold_result = cold_detect(dates, blues, greens, reds, nirs, swir1s, swir2s, thermals, qas)
```

### 4. Examples
#### 1) Running single pixel-based time series and plot time series:
/tool/notebook/pycold_example.ipynb
#### 2) Converting original ARD to binary uncompressed stack data (preprocessing):
/tool/python/AutoPrepareDataARD.py
#### 3) Running tile-based COLD to produce change records based on the stack data:
/tool/python/pycold_workflow.py
#### 4) Export disturbance maps from change records:
/tool/python/exportChangeMap.py
