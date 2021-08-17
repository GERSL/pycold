# this script show how to run COLD from a csv file that store time series information
import pysccd
import numpy as np
import os
import pandas as pd
import time
import datetime
from datetime import date
import pandas as pd
import numpy as np
import pycold

in_path = '/Users/coloury/Dropbox/Documents/S-CCD/test/spectral_336_3980_obs.csv'  # please change as needed
data = pd.read_csv(in_path, header=None)
# the sensor columnis not useful for now
data.columns = ['date', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'thermal', 'qa', 'sensor']

ret = pycold.pycold(data['date'].to_numpy(dtype=np.int32), data['blue'].to_numpy(dtype=np.short),
                    data['green'].to_numpy(dtype=np.short), data['red'].to_numpy(dtype=np.short),
                    data['nir'].to_numpy(dtype=np.short), data['swir1'].to_numpy(dtype=np.short),
                    data['swir2'].to_numpy(dtype=np.short), data['thermal'].to_numpy(dtype=np.short),
                    data['qa'].to_numpy(dtype=np.short))


# the below is used validate it using previous dataset
# cold_dt = np.dtype([('t_start', np.int32),
#                ('t_end', np.int32),
#                ('t_break', np.int32),
#                ('pos', np.int32),
#                ('num_obs', np.int32),
#                ('category', np.short),
#                ('change_prob', np.short),
#                ('coefs', np.double, (7, 8)),
#                ('rmse', np.double, 7),
#                ('magnitude', np.double, 7)])
# cold_plot = np.fromfile('/Users/coloury/Dropbox/Documents/spectral_336_3980_obs_ccd.dat', dtype=cold_dt)