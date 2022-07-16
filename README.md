# PYCOLD

# A Python library for COntinuous monitoring of Land Disturbance (COLD) and its extension algorithms at the speed of C language
The base algorithms were mostly written using C wrapped in python, and have been verified with MATLAB version (if it has). This library provides:
  1. Original COntinuous monitoring of Land Disturbance (COLD): a upgraded CCDC algorithm proposed by Dr.Zhe Zhu for offline satellite-based time-series analysis
  2. Stochastic Continuous Change Detection (S-CCD, a near real-time and short-memory implementation of COLD)
  3. Object-based COLD (OB-COLD), integrating spatial information into COLD by using a 'change object' view
## 1. Pre-work: clone github repo to your local directory

Clone repo for the first use:

```
git clone https://github.com/GERSL/pycold.git
```

Or you call pull the recent repo if you want to update the existing pycold repo:

```
git pull origin devel:devel
```



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
**Option 1: Install in development mode**

```bash
pip install -r requirements.txt
conda install gdal  # the easiest way to install gdal 
bash run_developer_setup.sh
```

**Option 2: Build and install a wheel** 

Scikit-build will invoke CMake and build everything. (you may need to remove
any existing `_skbuild` directory).

```bash
python setup.py bdist_wheel
```

Then you can pip install the wheel

```bash
pip install dist/pycold-0.1.0-cp38-cp38-linux_x86_64.whl
```

**Option 3: build standalone binaries with CMake by itself (recommended for C development)**

```bash
mkdir -p build
cd build
cmake ..
make 
```

## 3. Using pycold for pixel-based processing

COLD:

```python
>>> from pycold import cold_detect
>>> cold_result = cold_detect(dates, blues, greens, reds, nirs, swir1s, swir2s, thermals, qas)
```

S-CCD: 

```
# require offline processing for the first time 
>>> from pycold import sccd_detect, sccd_update
>>> sccd_pack = sccd_detect(dates, blues, greens, reds, nirs, swir1s, swir2s, thermals, qas)

# then use sccd_pack to do recursive and short-memory NRT update
>>> sccd_pack_new = sccd_update(sccd_pack, dates, blues, greens, reds, nirs, swir1s, swir2s, thermals, qas)
```




### Q&A
#### Q1: Has pycold been verified with original Matlab codes?
Re: yes, multiple rounds of verification have been done. Comparison based on two testing tiles shows that pycold and Matlab version have smaller than <2% differences for breakpoint detection and <2% differences for harmonic coefficients; the accuracy of pycold was also tested against the same reference dataset used in the original COLD paper (Zhu et al., 2020), and pycold reached the same accuracy (27% omission and 28% commission) showing that the discrepancy doesn't hurt accuracy. The primary source for the discrepancy is mainly from the rounding: MATLAB uses float64 precision, while pycold chose float32 to save the run-time computing memory and boost efficiency. 

#### Q2: how much time for production of a tile-based disturbance map (5000*5000 pixels) using pycold?
Re: I tested it in UCONN HPC environment (200 EPYC7452 cores): for processing a 40-year Landsat ARD tile (1982-2021), the stacking typically takes 15 mins; per-pixel COLD processing costs averagely 1 hour; exporting maps needs 7 mins.  



## 4. Citation

COLD:

Zhu, Z., Zhang, J., Yang, Z., Aljaddani, A. H., Cohen, W. B., Qiu, S., & Zhou, C. (2020). Continuous monitoring of land disturbance based on Landsat time series. *Remote Sensing of Environment*, *238*, 111116.

S-CCD:

Ye, S., Rogan, J., Zhu, Z., & Eastman, J. R. (2021). A near-real-time approach for monitoring forest disturbance using Landsat time series: Stochastic continuous change detection. *Remote Sensing of Environment*, *252*, 112167.

OB-COLD:

Ye, S., Zhu, Z., & Cao, G., (2022). Object-based continuous monitoring of land disturbance. Submitted to Remote Sensing of Environment
