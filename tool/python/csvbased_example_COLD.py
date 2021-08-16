# this script show how to run COLD from a csv file that store time series information
import pysccd
import numpy as np
import os
import pandas as pd
import time
import seaborn as sns
import datetime
from datetime import date

mode = 4
in_path = '/Users/coloury/Dropbox/Documents/QTProjects/S-CCD/test/spectral_336_3980_obs.csv' # please change as needed
out_path = '/Users/coloury/Dropbox/Documents/QTProjects/S-CCD/test' # please change as needed
probability = 0.95
min_days_conse = 80
row = 1 # won't affect the result for csv-based processing
col = 1 # won't affect the result for csv-based processing
method = 1  # 1 - COLD algorithm, 2 - S-CCD algorithm
user_mask_path = ''
output_mode = 1 # 11: both output; 1: output states; 10: output observation file
bandselection_bit = 62 # the bit value for the default spectral inputs, namely green, red, nir, SWIR1 and SWIR2
sccd_dt = np.dtype([('t_start', np.int32), ('t_end', np.int32), ('t_break', np.int32), ('pos', np.int32), ('num_obs', np.int32),
               ('category', np.int32), ('change_prob', np.int32),
               ('coef', np.double, (7, 8)), ('rmse', np.double, (7, 1)), ('magnitude', np.double, (7, 1))])

ret = pysccd.py_sccd(mode, in_path.encode('ascii'), out_path.encode('ascii'), \
                    row, col, method, bandselection_bit, user_mask_path.encode('ascii'),
                     probability, min_days_conse, output_mode=output_mode)

# extract base name
basename = os.path.splitext(os.path.basename(in_path))[0]

# read .dat file
ccd_plot = np.fromfile(os.path.join(out_path, basename + '_ccd.dat'), dtype=sccd_dt)
print(ccd_plot)