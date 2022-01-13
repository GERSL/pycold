import numpy as np

dt = np.dtype([('t_start', np.int32),
                   ('t_end', np.int32),
                   ('t_break', np.int32),
                   ('pos', np.int32),
                   ('num_obs', np.int32),
                   ('category', np.short),
                   ('change_prob', np.short),
                   ('coefs', np.float32, (7, 8)),
                   ('rmse', np.float32, 7),
                   ('magnitude', np.float32, 7)])


with open('tests/resources/test_parameters.yaml', 'r') as yaml_obj:
    test_parameter = yaml.safe_load(yaml_obj)
pyclassifier = PyClassifierHPC(test_parameter, 'tests/resources', 'tests/resources', year_lowbound=2015,
                               year_uppbound=2020)

cold_block = np.load(os.path.join('tests/resources', 'record_change_x1_y1_cold_ori.npy'))
block_features = pyclassifier.predict_features(1, cold_block)

curve_x1_y1 = []
curve_x1_y2 = []
curve_x2_y1 = []
curve_x2_y2 = []
for curve in cold_block:
    i_row = int((curve['pos'] - 1) / 100)
    i_col = int((curve['pos'] - 1) % 100)
    if i_row < 50 and i_col < 50:
        curve_x1_y1.append(curve)
    elif i_row < 50 and i_col >= 50:
        curve_x2_y1.append(curve)
    elif i_row >= 50 and i_col < 50:
        curve_x1_y2.append(curve)
    else:
        curve_x2_y2.append(curve)
curve_x1_y1 = np.array(curve_x1_y1, dtype=dt)
curve_x2_y1 = np.array(curve_x2_y1, dtype=dt)
curve_x1_y2 = np.array(curve_x1_y2, dtype=dt)
curve_x2_y2 = np.array(curve_x2_y2, dtype=dt)

np.save(os.path.join('tests/resources', 'record_change_x1_y1_cold.npy'), curve_x1_y1)
np.save(os.path.join('tests/resources', 'record_change_x1_y2_cold.npy'), curve_x1_y2)
np.save(os.path.join('tests/resources', 'record_change_x2_y1_cold.npy'), curve_x2_y1)
np.save(os.path.join('tests/resources', 'record_change_x2_y2_cold.npy'), curve_x2_y2)

# pyclassifier.hpc_preparation()
cold_block = np.load(os.path.join('tests/resources', 'record_change_x1_y1_cold.npy'))
block_features = pyclassifier.predict_features(1, cold_block)
pyclassifier.step1_save_features(1, block_features)
cold_block = np.load(os.path.join('tests/resources', 'record_change_x2_y1_cold.npy'))
block_features = pyclassifier.predict_features(2, cold_block)
pyclassifier.step1_save_features(2, block_features)
cold_block = np.load(os.path.join('tests/resources', 'record_change_x1_y2_cold.npy'))
block_features = pyclassifier.predict_features(3, cold_block)
pyclassifier.step1_save_features(3, block_features)
cold_block = np.load(os.path.join('tests/resources', 'record_change_x2_y2_cold.npy'))
block_features = pyclassifier.predict_features(4, cold_block)
pyclassifier.step1_save_features(4, block_features)

import matplotlib.pyplot as plt
plt.imshow(b1)

from osgeo import gdal_array
import gdal
reference_path = '/Users/coloury/Dropbox/UCONN/spatial/test_results/h027v009/LCMAP_CU_2001_V01_LCPRI_027009.tif'
ref_image = gdal.Open(reference_path, gdal.GA_ReadOnly)
trans = ref_image.GetGeoTransform()
proj = ref_image.GetProjection()
cols = ref_image.RasterXSize
ref_image_array = gdal_array.LoadFile(reference_path)
rows = ref_image.RasterYSize

outdriver1 = gdal.GetDriverByName("GTiff")
outdata = outdriver1.Create('/Users/coloury/Dropbox/UCONN/spatial/test_results/h027v009/LCMAP_CU_2001_V01_LCPRI_027009_partial.tif',
                            100, 100, 1, gdal.GDT_Byte)
outdata.GetRasterBand(1).WriteArray(ref_image_array[0:100,0:100])
outdata.FlushCache()
outdata.SetGeoTransform(trans)
outdata.FlushCache()
outdata.SetProjection(proj)
outdata.FlushCache()
outdata = None

np.save(full_feature_array)

def getdate_from_obiaresult_fn(filename):
    return int()

import pandas as pd
import os
from os.path import join
obiaresult_path = '/Users/coloury/Dropbox/Documents/OBIAresults'
obia_files = [f for f in os.listdir(obiaresult_path) if f.startswith('obiaresult')]
# sort image files by dates
cm_dates = [int(f[f.find('obiaresult_') + len('obiaresult_'):
                  f.find('obiaresult_') + len('obiaresult_')+6])
            for f in obia_files]


files_date_zip = sorted(zip(cm_dates, obia_files))
obia_files_sorted = [x[1] for x in files_date_zip]
obia_dates_sorted = np.asarray([x[0] for x in files_date_zip])
obia_tstack = [np.load(join(obiaresult_path, f[1])).reshape(53 * 53) * f[0]
               for f in files_date_zip]
obia_tstack = np.dstack(obia_tstack)
