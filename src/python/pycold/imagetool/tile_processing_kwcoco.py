"""
This script is for running COLD algorithm with kwcoco dataset.
See original code: ~/code/pycold/src/python/pycold/imagetool/tile_processing.py
"""

import os
import json
from os.path import join
import pandas as pd
import datetime
import numpy as np
import pycold
from datetime import datetime
from pytz import timezone
# import click
import time
from pycold import cold_detect
from scipy.stats import chi2
import ubelt as ub
from pycold.utils import get_rowcol_intile
import kwcoco

def tileprocessing_report(result_log_path, stack_path, version, algorithm, config, startpoint, cold_timepoint, tz,
                          n_cores, starting_date=0, n_cm_maps=0, year_lowbound=0, year_highbound=0):
    """
    output tile-based processing report
    Parameters
    ----------
    result_log_path: string
        outputted log path
    stack_path: string
        stack path
    version: string
    algorithm: string
    config: dictionary structure
    startpoint: a time point, when the program starts
    tz: string, time zone
    n_cores: the core number used
    starting_date: the first date of the total dataset
    n_cm_maps: the number of snapshots
    year_lowbound: the low bound of year range
    year_highbound: the upper bound of year range
    Returns
    -------
    """
    endpoint = datetime.now(tz)
    file = open(result_log_path, "w")
    file.write("PYCOLD V{} \n".format(version))
    file.write("Author: Su Ye(remoteseningsuy@gmail.com)\n")
    file.write("Algorithm: {} \n".format(algorithm))
    file.write("Starting_time: {}\n".format(startpoint.strftime('%Y-%m-%d %H:%M:%S')))
    file.write("Change probability threshold: {}\n".format(config['probability_threshold']))
    file.write("Conse: {}\n".format(config['conse']))
    file.write("stack_path: {}\n".format(stack_path))
    file.write("The number of requested cores: {}\n".format(n_cores))
    file.write("The program starts at {}\n".format(startpoint.strftime('%Y-%m-%d %H:%M:%S')))
    file.write("The COLD ends at {}\n".format(cold_timepoint.strftime('%Y-%m-%d %H:%M:%S')))
    file.write("The program ends at {}\n".format(endpoint.strftime('%Y-%m-%d %H:%M:%S')))
    file.write("The program lasts for {:.2f}mins\n".format((endpoint - startpoint) / datetime.timedelta(minutes=1)))
    file.close()

def is_finished_cold_blockfinished(result_path, nblocks):
    """
    check if the COLD algorithm finishes all blocks
    Parameters
    ----------
    result_path: the path that save COLD results
    nblocks: the block number

    Returns: boolean
    -------
        True -> all block finished
    """
    for n in range(nblocks):
        if not os.path.exists(os.path.join(result_path, 'COLD_block{}_finished.txt'.format(n + 1))):
            return False
    return True

def get_stack_date(block_x, block_y, stack_path, year_lowbound=0, year_highbound=0, nband=8):
    """
    :param block_x: block id at x axis
    :param block_y: block id at y axis
    :param stack_path: stack path
    :param year_lowbound: ordinal data of low bounds of selection date range
    :param year_highbound: ordinal data of upper bounds of selection date range
    :return:
        img_tstack, img_dates_sorted
        img_tstack - 3-d array (block_width * block_height, nband, nimage)
    """
    block_folder = join(stack_path, 'block_x{}_y{}'.format(block_x, block_y))
    meta_files = [m for m in os.listdir(block_folder) if m.endswith('.json')]

    # sort image files by ordinal dates
    img_dates = []
    img_files = []

    # read metadata and
    for meta in meta_files:
        metadata = open(join(block_folder, meta))
        config = json.load(metadata)
        ordinal_date = config['ordinal_date']
        img_name = config['image_name'] + '.npy'
        img_dates.append(ordinal_date)
        img_files.append(img_name)

    if len(img_files) == 0:
        return None, None

    sample_np = np.load(join(block_folder, img_files[0]))
    block_width = sample_np.shape[1]
    block_height = sample_np.shape[0]

    if year_lowbound > 0:
        year_low_ordinal = pd.Timestamp.toordinal(datetime(int(year_lowbound), 1, 1))
        img_dates, img_files = zip(*filter(lambda x: x[0] >= year_low_ordinal,
                                           zip(img_dates, img_files)))
    if year_highbound > 0:
        year_high_ordinal = pd.Timestamp.toordinal(datetime(int(year_highbound + 1), 1, 1))
        img_dates, img_files = zip(*filter(lambda x: x[0] < year_high_ordinal,
                                           zip(img_dates, img_files)))

    files_date_zip = sorted(zip(img_dates, img_files))
    img_files_sorted = [x[1] for x in files_date_zip]
    img_dates_sorted = np.asarray([x[0] for x in files_date_zip])
    img_tstack = np.dstack([np.load(join(block_folder, f)).reshape(block_width * block_height, nband)
                            for f in img_files_sorted])
    return img_tstack, img_dates_sorted

# @click.command()
# @click.option('--rank', type=int, default=0, help='the rank id')
# @click.option('--n_cores', type=int, default=0, help='the total cores assigned')
# @click.option('--stack_path', type=str, default=None, help='the path for stack data')
# @click.option('--result_path', type=str, default=None, help='the path for storing results')
# @click.option('--method', type=click.Choice(['COLD', 'OBCOLD', 'SCCDOFFLINE']), help='COLD, OBCOLD, SCCD-OFFLINE')
# @click.option('--year_lowbound', type=str, default=None, help='min year bound of image selection for processing. '
#                                                               'Example - 2015-01-01')
# @click.option('--year_highbound', type=str, default=None, help='max year bound of image selection for processing.'
#                                                                 'Example - 2021-12-31')
# @click.option('--b_c2', is_flag=True, show_default=True, default=False, help='indicate if it is c2 or hls; if so, '
#                                                                              'thermal bands will not be considered for observation selection ')

# def main(rank, n_cores, stack_path, result_path, method, year_lowbound, year_highbound, b_c2):

# debugging mode
def main():
    rank = 1
    n_cores = 1
    stack_path = '~/kwcoco/stacked/US_C000'
    result_path = '~/kwcoco/COLD/US_C000'
    yaml_path = '/home/jws18003/Document/pycold-uconnhpc/config_watch.yaml'
    method = 'COLD'
    year_lowbound = 2017
    year_highbound = 2022
    b_c2 = True

    # FIXME: stack_path should be "dpath / 'stacked' / video_name" defined in prepare_kwcoco.py
    # from pycold.imagetool.prepare_kwcoco import grab_demo_kwcoco_dataset
    # coco_fpath = grab_demo_kwcoco_dataset()
    # coco_dset = kwcoco.CocoDataset.coerce(coco_fpath)
    # dset_hashid = coco_dset._cached_hashid()[0:8]
    # dpath = (ub.Path.appdir('pycold/kwcoco_prep') / dset_hashid).ensuredir()
    # out_dir = (dpath / 'stacked').ensuredir()
    #
    # # For now, I manually define stack_path
    # stack_path = '~/kwcoco/stacked/US_C000'

    tz = timezone('US/Eastern')
    start_time = datetime.now(tz)

    # Define year_low_ordinal and year_high_ordinal to filter year for COLD processing
    if year_lowbound is None:
        year_lowbound = 0
    else:
        year_low_ordinal = pd.Timestamp.toordinal(datetime(int(year_lowbound), 1, 1))
    if year_highbound is None:
        year_highbound = 0
    else:
        year_high_ordinal = pd.Timestamp.toordinal(datetime(int(year_highbound + 1), 1, 1))

    # Reading/Defining config
    config = {'n_cols': 660,
              'n_rows': 780,
              'n_block_x': 20,
              'n_block_y': 20,
              'probability_threshold': 0.99,
              'conse': 6
              }
    # set up some additional config
    block_width = int(config['n_cols'] / config['n_block_x'])  # width of a block
    block_height = int(config['n_rows'] / config['n_block_y'])  # height of a block
    if (config['n_cols'] % block_width != 0) or (config['n_rows'] % block_height != 0):
        print('n_cols, n_rows must be divisible respectively by block_width, block_height! Please double '
              'check your config yaml')
        exit()

    # logging and folder preparation
    if rank == 1:
        if not os.path.exists(result_path):
            os.makedirs(result_path)

        print("The per-pixel time series processing begins: {}".format(start_time.strftime('%Y-%m-%d %H:%M:%S')))

        if not os.path.exists(stack_path):
            print("Failed to locate stack folders. The program ends: {}".format(
                datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))
            return

    #########################################################################
    #                        per-pixel COLD procedure                       #
    #########################################################################
    threshold = chi2.ppf(config['probability_threshold'], 5)
    nblock_eachcore = int(np.ceil(config['n_block_x'] * config['n_block_y'] * 1.0 / n_cores))
    for i in range(nblock_eachcore):
        block_id = n_cores * i + rank  # started from 1, i.e., rank, rank + n_cores, rank + 2 * n_cores
        if block_id > config['n_block_x'] * config['n_block_y']:
            break
        block_y = int((block_id - 1) / config['n_block_x']) + 1  # note that block_x and block_y start from 1
        block_x = int((block_id - 1) % config['n_block_x']) + 1
        if os.path.exists(join(result_path, 'COLD_block{}_finished.txt'.format(block_id))):
            print("Per-pixel COLD processing is finished for block_x{}_y{} ({})".format(block_x, block_y,
                                                                                        datetime.now(tz).strftime(
                                                                                            '%Y-%m-%d %H:%M:%S')))
            continue
        img_tstack, img_dates_sorted = get_stack_date(block_x, block_y, stack_path, year_lowbound,
                                                      year_highbound)

        # Define empty list
        result_collect = []
        date_collect = []

        if method == "COLD":
            # start looping every pixel in the block
            for pos in range(block_width * block_height):
                original_row, original_col = get_rowcol_intile(pos, block_width,
                                                               block_height, block_x, block_y)
                try:
                    # if method == 'COLD':
                    cold_result = cold_detect(img_dates_sorted,
                                              img_tstack[pos, 0, :].astype(np.int64),
                                              img_tstack[pos, 1, :].astype(np.int64),
                                              img_tstack[pos, 2, :].astype(np.int64),
                                              img_tstack[pos, 3, :].astype(np.int64),
                                              img_tstack[pos, 4, :].astype(np.int64),
                                              img_tstack[pos, 5, :].astype(np.int64),
                                              img_tstack[pos, 6, :].astype(np.int64),
                                              img_tstack[pos, 7, :].astype(np.int64),
                                              t_cg=threshold,
                                              conse=config['conse'], b_c2=b_c2,
                                              pos=config['n_cols'] * (original_row - 1) + original_col)
                # except RuntimeError:
                except:
                    print("COLD fails at original_row {}, original_col {} ({})".format(original_row, original_col,
                                                                                       datetime.now(tz)
                                                                                       .strftime(
                                                                                           '%Y-%m-%d %H:%M:%S')))
                else:
                    result_collect.append(cold_result)

            # save the dataset
            if len(result_collect) > 0:
                np.save(join(result_path, 'record_change_x{}_y{}_cold.npy'.format(block_x, block_y)),
                        np.hstack(result_collect))

        with open(join(result_path, 'COLD_block{}_finished.txt'.format(block_id)), 'w'):
            pass

        print("Per-pixel COLD processing is finished for block_x{}_y{} ({})"
              .format(block_x, block_y, datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

    # wait for all cores to be finished
    while not is_finished_cold_blockfinished(result_path, config['n_block_x'] * config['n_block_y']):
        time.sleep(30)

    if rank == 1:
        cold_timepoint = datetime.now(tz)

    if rank == 1:
        # tile_based report
        tileprocessing_report(join(result_path, 'tile_processing_report.log'), stack_path, pycold.__version__,
                              method, config, start_time, cold_timepoint, tz, n_cores)
        print("The whole procedure finished: {}".format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))


if __name__ == '__main__':
    main()
