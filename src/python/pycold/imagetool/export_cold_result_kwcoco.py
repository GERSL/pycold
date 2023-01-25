"""
This script is for exporting COLD algorithm results (change vector, coefficients, RMSEs)
to geotiff raster with kwcoco dataset.
See original code: ~/code/pycold/src/python/pycold/imagetool/export_change_map.py
"""

import os
import numpy as np
import pandas as pd
from osgeo import gdal
from osgeo import gdal_array
# import click
import datetime as datetime
from os.path import join
import json
# from mpi4py import MPI
import pickle
import yaml
from collections import namedtuple
import kwcoco, kwimage

# Now write out the videospace image with correct updated transforms.
image_fpath = "~/kwcoco/US_C000_ref.tif" # image_fpath is reference image in image space
ref_image = gdal.Open(image_fpath, gdal.GA_ReadOnly)
trans = ref_image.GetGeoTransform()
proj = ref_image.GetProjection()
cols = ref_image.RasterXSize
rows = ref_image.RasterYSize

# Get original transform from projection to image space
c, a, b, f, d, e = trans
original = kwimage.Affine(np.array([
    [a, b, c],
    [d, e, f],
    [0, 0, 1],
]))

# load kwcoco image
coco_fpath = '~/data/dvc-repos/smart_data_dvc/Aligned-Drop4-2022-08-08-TA1-S2-L8-ACC/data.kwcoco.json'
dset = kwcoco.CocoDataset.coerce(coco_fpath)
video_id = 17
coco_image = dset.images(video_id=video_id).coco_images[0]

# Get the modifier transform to move from image space to video space
warp_vid_from_img = kwimage.Affine.coerce(coco_image.img['warp_img_to_vid']).inv()

# Combine transforms to get a new transform that goes
# from the projection to video space
new_geotrans =  original @ warp_vid_from_img
a, b, c, d, e, f, g, h, i = np.array(new_geotrans).ravel().tolist()
new_gdal_transform = (c, a, b, f, d, e)

vid_w = coco_image.img['width'] 
vid_h = coco_image.img['height'] 

new_gdal_transform = (275011.57468469924, 29.78988341490997, 3.256665984888778e-29, 4332467.113652797, -6.011541802002886e-29, -29.804444582314744)

# FIXME: stack_path should be "dpath / 'stacked' / video_name" defined in prepare_kwcoco.py
# from pycold.imagetool.prepare_kwcoco import grab_demo_kwcoco_dataset
# coco_fpath = grab_demo_kwcoco_dataset()
# coco_dset = kwcoco.CocoDataset.coerce(coco_fpath)
# dset_hashid = coco_dset._cached_hashid()[0:8]
# dpath = (ub.Path.appdir('pycold/kwcoco_prep') / dset_hashid).ensuredir()
# out_dir = (dpath / 'stacked').ensuredir()
#
# For now, I manually define stack_path
stack_path = '~/kwcoco/stacked/US_C000'


coef_names = ['a0', 'c1', 'a1', 'b1', 'a2', 'b2', 'a3', 'b3', 'cv', 'rmse']
band_names = [0, 1, 2, 3, 4, 5, 6]
SLOPE_SCALE = 10000

BAND_INFO = {0: 'blue',
             1: 'green',
             2: 'red',
             3: 'nir',
             4: 'swir1',
             5: 'swir2'}

# copy from /pycold/src/python/pycold/pyclassifier.py because MPI has conflicts with the pycold package in UCONN HPC.
# Dirty approach!
def extract_features(cold_plot, band, ordinal_day_list, nan_val, timestamp, feature_outputs=['a0', 'a1', 'b1']):
    """
    generate features for classification based on a plot-based rec_cg and a list of days to be predicted
    Parameters
    ----------
    cold_plot: nested array
        plot-based rec_cg
    band: integer
        the predicted band number range from 0 to 6
    ordinal_day_list: list
        a list of days that this function will predict every days as a list as output
    nan_val: integer
        NA value assigned to the output
    feature_outputs: a list of outputted feature name
        it must be within [a0, c1, a1, b1, a2, b2, a3, b3, cv, rmse]
    Returns
    -------
        feature: a list (length = n_feature) of 1-array [len(ordinal_day_list)]
    """
    features = [np.full(len(ordinal_day_list), nan_val, dtype=np.double) for x in range(len(feature_outputs))]
    for index, ordinal_day in enumerate(ordinal_day_list):
        for idx, cold_curve in enumerate(cold_plot):
            if idx == len(cold_plot) - 1:
                last_year = pd.Timestamp.fromordinal(cold_plot[idx]['t_end']).year
                max_days = datetime.date(last_year, 12, 31).toordinal()
            else:
                max_days = cold_plot[idx + 1]['t_start']
            break_year = pd.Timestamp.fromordinal(cold_curve['t_break']).year if(cold_curve['t_break'] > 0 and cold_curve['change_prob'] == 100) else -9999

            if cold_curve['t_start'] <= ordinal_day < max_days:
                for n, feature in enumerate(feature_outputs):
                    if feature not in feature_outputs:
                        raise Exception('the outputted feature must be in [a0, c1, a1, b1,a2, b2, a3, b3, cv, rmse]')
                    if feature == 'a0':
                        features[n][index] = cold_curve['coefs'][band][0] + cold_curve['coefs'][band][1] * \
                                             ordinal_day / SLOPE_SCALE
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'c1':
                        features[n][index] = cold_curve['coefs'][band][1] / SLOPE_SCALE
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'a1':
                        features[n][index] = cold_curve['coefs'][band][2]
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'b1':
                        features[n][index] = cold_curve['coefs'][band][3]
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'a2':
                        features[n][index] = cold_curve['coefs'][band][4]
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'b2':
                        features[n][index] = cold_curve['coefs'][band][5]
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'a3':
                        features[n][index] = cold_curve['coefs'][band][6]
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'b3':
                        features[n][index] = cold_curve['coefs'][band][7]
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'rmse':
                        features[n][index] = cold_curve['rmse'][band]
                        if np.isnan(features[n][index]):
                            features[n][index] = 0                    
                break

        if 'cv' in feature_outputs:
            # ordinal_day_years = [pd.Timestamp.fromordinal(day).year for day in ordinal_day_list]
            for index, ordinal_day in enumerate(ordinal_day_list):
                ordinal_year = pd.Timestamp.fromordinal(ordinal_day).year
                for cold_curve in cold_plot:
                    if (cold_curve['t_break'] == 0) or (cold_curve['change_prob'] != 100):
                        continue
                    break_year = pd.Timestamp.fromordinal(cold_curve['t_break']).year
                    if timestamp == True:
                        if ordinal_day == cold_curve['t_break']:
                            features[feature_outputs.index('cv')][index] = cold_curve['magnitude'][band]
                            continue
                    else:
                        if break_year == ordinal_year:
                            features[feature_outputs.index('cv')][index] = cold_curve['magnitude'][band]
                            continue

    return features


def getcategory_cold(cold_plot, i_curve):
    t_c = -200
    if cold_plot[i_curve]['magnitude'][3] > t_c and cold_plot[i_curve]['magnitude'][2] < -t_c and \
            cold_plot[i_curve]['magnitude'][4] < -t_c:
        if cold_plot[i_curve + 1]['coefs'][3, 1] > np.abs(cold_plot[i_curve]['coefs'][3, 1]) and \
                cold_plot[i_curve + 1]['coefs'][2, 1] < -np.abs(cold_plot[i_curve]['coefs'][2, 1]) and \
                cold_plot[i_curve + 1]['coefs'][4, 1] < -np.abs(cold_plot[i_curve]['coefs'][4, 1]):
            return 3  # aforestation
        else:
            return 2  # regrowth
    else:
        return 1  # land disturbance

def getcategory_obcold(cold_plot, i_curve, last_dist_type):
    t_c = -250
    if cold_plot[i_curve]['magnitude'][3] > t_c and cold_plot[i_curve]['magnitude'][2] < -t_c and \
            cold_plot[i_curve]['magnitude'][4] < -t_c:
        if cold_plot[i_curve + 1]['coefs'][3, 1] > np.abs(cold_plot[i_curve]['coefs'][3, 1]) and \
                cold_plot[i_curve + 1]['coefs'][2, 1] < -np.abs(cold_plot[i_curve]['coefs'][2, 1]) and \
                cold_plot[i_curve + 1]['coefs'][4, 1] < -np.abs(cold_plot[i_curve]['coefs'][4, 1]):
            return 3  # aforestation
        else:
            return 2  # regrowth
    else:
        if i_curve > 0:
            if (cold_plot[i_curve]['t_break'] - cold_plot[i_curve - 1]['t_break'] > 365.25 * 5) or (
                    last_dist_type != 1):
                return 1
            flip_count = 0
            for b in range(5):
                if cold_plot[i_curve]['magnitude'][b + 1] * cold_plot[i_curve - 1]['magnitude'][b + 1] < 0:
                    flip_count = flip_count + 1
            if flip_count >= 4:
                return 4
            else:
                return 1
        else:
            return 1  # land disturbance

# @click.command()
# @click.option('--reccg_path', type=str, help='rec_cg folder')
# @click.option('--reference_path', type=str, help='image path used to provide georeference for output images')
# @click.option('--out_path', type=str, help='output folder for saving image')
# @click.option('--method', type=click.Choice(['COLD', 'OBCOLD', 'SCCDOFFLINE']), default='COLD', help='the algorithm used for processing')
# @click.option('--region', type=str, help='kwcoco region, e.g., US_C000')
# @click.option('--probability', type=float, help='COLD probability threshold: e.g., 0.99')
# @click.option('--yaml_path', type=str, help='path for yaml file')
# @click.option('--year_lowbound', type=int, default=1982, help='the starting year for exporting')
# @click.option('--year_highbound', type=int, default=2020, help='the ending year for exporting')
# @click.option('--coefs', type=str, default=None, help='if output coefs layers')
# @click.option('--coefs_bands', type=str, default='0, 1, 2, 3, 4, 5, 6', help='indicate the ba_nds for output coefs_bands,'
#                                                                     'only works when coefs is True; note that the band '
#                                                                     'order is b,g,r,n,s1,s2,t')
# @click.option('--stack_path', type=str, help='stack folder')
# @click.option('--timestamp', is_flag=True, show_default=True, default=False, help='True: exporting cold result by timestamp,'
#                                                                                   'False: exporting cold result by year'
#                                                                                   'Default is False')

def main(stack_path, reccg_path, reference_path, out_path, method, region, probability, conse, year_lowbound,
         year_highbound, coefs, coefs_bands, timestamp):
# def main():

    ## MPI mode
    # comm = MPI.COMM_WORLD
    # rank = comm.Get_rank()
    # n_process = comm.Get_size()
    if method == 'OBCOLD':
        reccg_path = os.path.join(reccg_path, 'obcold')
        if timestamp == True:
            out_path = os.path.join(out_path, 'obcold_maps', 'by_timestamp')
        else:
            out_path = os.path.join(out_path, 'obcold_maps', 'by_year')

    elif method == 'COLD' or method == 'HybridCOLD':
        if timestamp == True:
            out_path = os.path.join(out_path, 'cold_maps', 'by_timestamp')
        else:
            out_path = os.path.join(out_path, 'cold_maps', 'by_year')
    if coefs is not None:
        try:
            coefs = list(coefs.split(","))
            coefs = [str(coef) for coef in coefs]
        except:
            print("Illegal coefs inputs: example, --coefs='a0, c1, a1, b1, a2, b2, a3, b3, cv, rmse'")

        try:
            coefs_bands = list(coefs_bands.split(","))
            coefs_bands = [int(coefs_band) for coefs_band in coefs_bands]
        except:
            print("Illegal coefs_bands inputs: example, --coefs_bands='0, 1, 2, 3, 4, 5, 6'")

    dt = np.dtype([('t_start', np.int32),
                   ('t_end', np.int32),
                   ('t_break', np.int32),
                   ('pos', np.int32),
                   ('num_obs', np.int32),
                   ('category', np.short),
                   ('change_prob', np.short),
                   ('coefs', np.float32, (7, 8)),   # note that the slope coefficient was scaled up by 10000
                   ('rmse', np.float32, 7),
                   ('magnitude', np.float32, 7)])

    if coefs is not None:
        assert all(elem in coef_names for elem in coefs)
        assert all(elem in band_names for elem in coefs_bands)

    if rank == 0:
        if not os.path.exists(out_path):
            os.makedirs(out_path)

        # with open(yaml_path, 'r') as yaml_obj:
        #     config = yaml.safe_load(yaml_obj)

        config = {'n_block_x': 20,
                  'n_block_y': 20,
                  'padded_n_cols': 660,
                  'padded_n_rows': 780
        }

        config['block_width'] = int(config['n_cols'] / config['n_block_x'])  # width of a block
        config['block_height'] = int(config['n_rows'] / config['n_block_y'])  # height of a block
        config['n_blocks'] = config['n_block_x'] * config['n_block_y']  # total number of blocks
        ref_image = gdal.Open(reference_path, gdal.GA_ReadOnly)
        trans = ref_image.GetGeoTransform()
        proj = ref_image.GetProjection()
        cols = ref_image.RasterXSize
        rows = ref_image.RasterYSize
    else:
        trans = None
        proj = None
        cols = None
        rows = None
        config = None

    ## MPI mode
    # trans = comm.bcast(trans, root=0)
    # proj = comm.bcast(proj, root=0)
    # cols = comm.bcast(cols, root=0)
    # rows = comm.bcast(rows, root=0)
    # config = comm.bcast(config, root=0)

    # stack_folder = '/gpfs/scratchfs1/zhz18039/jws18003/kwcoco/stacked_KR_R002_drop4_2paths/KR_R002'

    ranks_percore = int(np.ceil(config['n_blocks'] / n_process))
    for i in range(ranks_percore):
        iblock = n_process * i + rank
        if iblock >= config['n_blocks']:
            break
        current_block_y = int(np.floor(iblock / config['n_block_x'])) + 1
        current_block_x = iblock % config['n_block_x'] + 1
        if method == 'OBCOLD':
            filename = 'record_change_x{}_y{}_obcold.npy'.format(current_block_x, current_block_y)
        elif method == 'COLD':
            filename = 'record_change_x{}_y{}_cold.npy'.format(current_block_x, current_block_y)
        elif method == 'HybridCOLD':
            filename = 'record_change_x{}_y{}_hybridcold.npy'.format(current_block_x, current_block_y)

        block_folder = os.path.join(stack_path, 'block_x{}_y{}'.format(current_block_x, current_block_y))

        if timestamp == False:
            year_list_to_predict = list(range(year_lowbound, year_highbound + 1))
            ordinal_day_list = [pd.Timestamp.toordinal(datetime.date(year, 7, 1)) for year
                                in year_list_to_predict]
            results_block = [np.full((config['block_height'], config['block_width']), -9999, dtype=np.int16)
                             for t in range(year_highbound - year_lowbound + 1)]
            if coefs is not None:
                results_block_coefs = np.full(
                    (config['block_height'], config['block_width'], len(coefs) * len(coefs_bands),
                     year_highbound - year_lowbound + 1), -9999, dtype=np.float32)

            print('processing the rec_cg file {}'.format(os.path.join(reccg_path, filename)))
            if not os.path.exists(os.path.join(reccg_path, filename)):
                print('the rec_cg file {} is missing'.format(os.path.join(reccg_path, filename)))
                for year in range(year_lowbound, year_highbound + 1):
                    outfile = os.path.join(out_path, 'tmp_map_block{}_{}.npy'.format(iblock + 1, year))
                    np.save(outfile, results_block[year - year_lowbound])
                continue
        else:
            meta_files = [m for m in os.listdir(block_folder) if m.endswith('.json')]

            # sort image files by ordinal dates
            img_dates = []
            img_files = []

            # read metadata and
            for meta in meta_files:
                metadata = open(join(block_folder, meta))
                meta_config = json.load(metadata)
                ordinal_date = meta_config['ordinal_date']
                img_name = meta_config['image_name'] + '.npy'
                img_dates.append(ordinal_date)
                img_files.append(img_name)

            if year_lowbound > 0:
                year_low_ordinal = pd.Timestamp.toordinal(datetime.datetime(int(year_lowbound), 1, 1))
                img_dates, img_files = zip(*filter(lambda x: x[0] >= year_low_ordinal,
                                                   zip(img_dates, img_files)))
            if year_highbound > 0:
                year_high_ordinal = pd.Timestamp.toordinal(datetime.datetime(int(year_highbound + 1), 1, 1))
                img_dates, img_files = zip(*filter(lambda x: x[0] < year_high_ordinal,
                                                   zip(img_dates, img_files)))
            img_dates = sorted(img_dates)
            ordinal_day_list = img_dates
            results_block = [np.full((config['block_height'], config['block_width']), -9999, dtype=np.int16)
                             for t in range(len(ordinal_day_list))]

            if coefs is not None:
                results_block_coefs = np.full(
                    (config['block_height'], config['block_width'], len(coefs) * len(coefs_bands),
                     len(ordinal_day_list)), -9999, dtype=np.float32)

            print('processing the rec_cg file {}'.format(os.path.join(reccg_path, filename)))
            if not os.path.exists(os.path.join(reccg_path, filename)):
                print('the rec_cg file {} is missing'.format(os.path.join(reccg_path, filename)))

                for day in range(len(ordinal_day_list)):
                    outfile = os.path.join(out_path, 'tmp_map_block{}_{}.npy'.format(iblock + 1, ordinal_day_list[day]))
                    # if not os.path.exists(outfile):
                    np.save(outfile, results_block[day])
                continue

          cold_block = np.array(np.load(os.path.join(reccg_path, filename)), dtype=dt)

          cold_block.sort(order='pos')
          current_processing_pos = cold_block[0]['pos']
          current_dist_type = 0

          for count, curve in enumerate(cold_block):
              if curve['pos'] != current_processing_pos:
                  current_processing_pos = curve['pos']
                  current_dist_type = 0

              if curve['change_prob'] < 100 or curve['t_break'] == 0 or count == (len(cold_block) - 1):  # last segment
                  continue

              i_col = int((curve["pos"] - 1) % config['n_cols']) - \
                      (current_block_x - 1) * config['block_width']
              i_row = int((curve["pos"] - 1) / config['n_cols']) - \
                      (current_block_y - 1) * config['block_height']
              if i_col < 0:
                  dat_pth = '?'
                  print('Processing {} failed: i_row={}; i_col={} for {}'.format(filename, i_row, i_col, dat_pth))
                  return

              if method == 'OBCOLD':
                  current_dist_type = getcategory_obcold(cold_block, count, current_dist_type)
              else:
                  current_dist_type = getcategory_cold(cold_block, count)
              break_year = pd.Timestamp.fromordinal(curve['t_break']).year
              if break_year < year_lowbound or break_year > year_highbound:
                  continue
              results_block[break_year - year_lowbound][i_row][i_col] = current_dist_type * 1000 + curve['t_break'] - \
                  (pd.Timestamp.toordinal(datetime.date(break_year, 1, 1))) + 1

          if coefs is not None:
              cold_block_split = np.split(cold_block, np.argwhere(np.diff(cold_block['pos']) != 0)[:, 0] + 1)
              for element in cold_block_split:
                  # the relative column number in the block
                  i_col = int((element[0]["pos"] - 1) % config['n_cols']) - \
                          (current_block_x - 1) * config['block_width']
                  i_row = int((element[0]["pos"] - 1) / config['n_cols']) - \
                          (current_block_y - 1) * config['block_height']

                  for band_idx, band in enumerate(coefs_bands):
                      feature_row = extract_features(element, band, ordinal_day_list, -9999, timestamp,
                                                     feature_outputs=coefs)
                      for index, coef in enumerate(coefs):
                          results_block_coefs[i_row][i_col][index + band_idx * len(coefs)][:] = \
                              feature_row[index]
          # e.g., 1315 means that disturbance happens at doy of 315
          # save the temp dataset out
          if timestamp == False:
              for year in range(year_lowbound, year_highbound + 1):
                  outfile = os.path.join(out_path,
                                         'tmp_map_block{}_{}.npy'.format(iblock + 1, year))
                  np.save(outfile, results_block[year - year_lowbound])
                  if coefs is not None:
                      outfile = os.path.join(out_path,
                                             'tmp_coefmap_block{}_{}.npy'.format(iblock + 1, year))
                      np.save(outfile, results_block_coefs[:, :, :, year - year_lowbound])
          else:
              for day in range(len(ordinal_day_list)):
                  if coefs is not None:
                      outfile = os.path.join(out_path,
                                             'tmp_coefmap_block{}_{}.npy'.format(iblock + 1,
                                                                                 ordinal_day_list[day]))
                      np.save(outfile, results_block_coefs[:, :, :, day])

    # MPI mode (wait for all processes)
    # comm.Barrier()

    if rank == 0:
        # assemble
        if timestamp is False:
            for year in range(year_lowbound, year_highbound + 1):
                tmp_map_blocks = [np.load(os.path.join(out_path, 'tmp_map_block{}_{}.npy'.format(x + 1, year)))
                                  for x in range(config['n_blocks'])]

                results = np.hstack(tmp_map_blocks)
                results = np.vstack(np.hsplit(results, config['n_block_x']))

                for x in range(config['n_blocks']):
                    os.remove(os.path.join(out_path, 'tmp_map_block{}_{}.npy'.format(x + 1, year)))
                outname = '%s_%s_prob_%s_conse_%s_%s_break_map.tif' % (region, method, probability, conse, year)
                outfile = os.path.join(out_path, outname)
                outdriver1 = gdal.GetDriverByName("GTiff")
                outdata = outdriver1.Create(outfile, vid_w, vid_h, 1, gdal.GDT_Int16)
                outdata.GetRasterBand(1).WriteArray(results[:vid_h, :vid_w])
                outdata.FlushCache()
                outdata.SetGeoTransform(new_gdal_transform)
                outdata.FlushCache()
                outdata.SetProjection(proj)
                outdata.FlushCache()

            # output recent disturbance year
            recent_dist = np.full((vid_h, vid_w), 0, dtype=np.int16)
            for year in range(year_lowbound, year_highbound + 1):
                outname = '%s_%s_prob_%s_conse_%s_%s_break_map.tif' % (region, method, probability, conse, year)
                breakmap = gdal_array.LoadFile(os.path.join(out_path, outname))
                recent_dist[
                    (breakmap / 1000).astype(np.byte) == 1] = year
            outname = "%s_%s_prob_%s_conse_%s_recent_disturbance_map.tif" % (region, method, probability, conse)
            outfile = os.path.join(out_path, outname)
            outdriver1 = gdal.GetDriverByName("GTiff")
            outdata = outdriver1.Create(outfile, vid_w, vid_h, 1, gdal.GDT_Int16)
            outdata.GetRasterBand(1).WriteArray(recent_dist[:vid_h, :vid_w])
            outdata.FlushCache()
            outdata.SetGeoTransform(new_gdal_transform)
            outdata.FlushCache()
            outdata.SetProjection(proj)
            outdata.FlushCache()

            first_dist = np.full((vid_h, vid_w), 0, dtype=np.int16)
            for year in range(year_highbound, year_lowbound - 1, -1):
                outname = '%s_%s_prob_%s_conse_%s_%s_break_map.tif' % (region, method, probability, conse, year)
                breakmap = gdal_array.LoadFile(os.path.join(out_path, outname))
                first_dist[(breakmap / 1000).astype(np.byte) == 1] = year
            outname = "%s_%s_prob_%s_conse_%s_first_disturbance_map.tif" % (region, method, probability, conse)
            outfile = os.path.join(out_path, outname)
            outdriver1 = gdal.GetDriverByName("GTiff")
            outdata = outdriver1.Create(outfile, vid_w, vid_h, 1, gdal.GDT_Int16)
            outdata.GetRasterBand(1).WriteArray(first_dist[:vid_h, :vid_w])
            outdata.FlushCache()
            outdata.SetGeoTransform(new_gdal_transform)
            outdata.FlushCache()
            outdata.SetProjection(proj)
            outdata.FlushCache()

            if coefs is not None:
                for year in range(year_lowbound, year_highbound + 1):
                    tmp_map_blocks = [np.load(os.path.join(out_path, 'tmp_coefmap_block{}_{}.npy'.format(x + 1, year)))
                                      for x in range(config['n_blocks'])]
                    results = np.hstack(tmp_map_blocks)
                    results = np.vstack(np.hsplit(results, config['n_block_x']))
                    ninput = 0
                    for band_idx, band_name in enumerate(coefs_bands):
                        for coef_index, coef in enumerate(coefs):
                            # FIXME: Best name for outputs...?
                            band = BAND_INFO[band_name]
                            if coef == 'cv':
                                results[results == -9999.0] = 0
                            outname = '%s_%s_prob_%s_conse_%s_%s_%s_%s.tif' % (
                            region, method, probability, conse, year, band,
                            coef)
                            outfile = os.path.join(out_path, outname)
                            outdriver1 = gdal.GetDriverByName("GTiff")
                            outdata = outdriver1.Create(outfile, vid_w, vid_h, 1, gdal.GDT_Float32)
                            outdata.GetRasterBand(1).WriteArray(results[:vid_h, :vid_w, ninput])
                            outdata.FlushCache()
                            outdata.SetGeoTransform(new_gdal_transform)
                            outdata.FlushCache()
                            outdata.SetProjection(proj)
                            outdata.FlushCache()
                            ninput = ninput + 1
                    for x in range(config['n_blocks']):
                        os.remove(
                            os.path.join(out_path, 'tmp_coefmap_block{}_{}.npy'.format(x + 1, year)))

        else:
            if coefs is not None:
                for day in range(len(ordinal_day_list)):
                    tmp_map_blocks = [np.load(
                        os.path.join(out_path, 'tmp_coefmap_block{}_{}.npy'.format(x + 1, ordinal_day_list[day])))
                        for x in range(config['n_blocks'])]

                    results = np.hstack(tmp_map_blocks)
                    results = np.vstack(np.hsplit(results, config['n_block_x']))
                    ninput = 0
                    for band_idx, band_name in enumerate(coefs_bands):
                        for coef_index, coef in enumerate(coefs):
                            # FIXME: Best name for outputs...?
                            date = str(datetime.date.fromordinal(ordinal_day_list[day]))
                            band = BAND_INFO[band_name]
                            outname = '%s_%s_prob_%s_conse_%s_%s_%s_%s.tif' % (
                            region, method, probability, conse, date, band, coef)
                            outfile = os.path.join(out_path, outname)
                            outdriver1 = gdal.GetDriverByName("GTiff")
                            outdata = outdriver1.Create(outfile, vid_w, vid_h, 1, gdal.GDT_Float32)
                            outdata.GetRasterBand(1).WriteArray(results[:vid_h, :vid_w, ninput])
                            outdata.FlushCache()
                            outdata.SetGeoTransform(new_gdal_transform)
                            outdata.FlushCache()
                            outdata.SetProjection(proj)
                            outdata.FlushCache()
                            ninput = ninput + 1

                    for x in range(config['n_blocks']):
                        os.remove(
                            os.path.join(out_path, 'tmp_coefmap_block{}_{}.npy'.format(x + 1, ordinal_day_list[day])))

if __name__ == '__main__':
    # main()
    rank = 0
    n_process = 1
    method = "COLD"
    region = "US_C000"
    probability = 0.99
    stack_path = '~/kwcoco/stacked/US_C000_metadata'
    out_path = "~/kwcoco/COLD_US_C000"
    reccg_path = "~/kwcoco/COLD_US_C000"
    reference_path = "~/kwcoco/US_C000_ref.tif"
    yaml_path = "~/pycold-uconnhpc/config_watch.yaml"
    coefs = ['cv']
    band_names = [0, 1, 2, 3, 4, 5]
    year_lowbound = 2017
    year_highbound = 2021
    coefs_bands = [5]
    timestamp = True
    main(stack_path, reccg_path, reference_path, out_path, method, region, probability, year_lowbound, year_highbound, yaml_path, coefs, coefs_bands, timestamp)