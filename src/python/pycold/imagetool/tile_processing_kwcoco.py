"""
This script is for running COLD algorithm with kwcoco dataset.
See original code: ~/code/pycold/src/python/pycold/imagetool/tile_processing.py
"""

import os, time
import json
from os.path import join
import pandas as pd
import datetime
import numpy as np
from datetime import datetime, timedelta
from pytz import timezone
# import click
import yaml
from scipy.stats import chi2
import ubelt as ub
import pycold
from pycold import cold_detect, sccd_detect
from pycold.utils import get_rowcol_intile, get_doy, assemble_cmmaps, unindex_sccdpack
from pycold.ob_analyst import ObjectAnalystHPC
from pycold.pyclassifier import PyClassifierHPC
from pycold.app import defaults
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

#########################################################################
#                     function for OB-COLD procedure                    #
#########################################################################

def reading_start_dates_nmaps(stack_path, year_lowbound, year_highbound, cm_interval):
    """
    Parameters
    ----------
    stack_path: string
        stack_path for saving starting_last_dates.txt
    cm_interval: interval
        day interval for outputting change magnitudes
    Returns
    -------
        (starting_date, n_cm_maps)
        starting_date - starting date is the first date of the whole dataset,
        n_cm_maps - the number of change magnitudes to be outputted per pixel per band
    """
    # read starting and ending dates, note that all blocks only has one starting and last date (mainly for obcold)
    try:
        block_folder = join(stack_path, 'block_x1_y1')
        meta_files = [m for m in os.listdir(block_folder) if m.endswith('.json')]

        # sort image files by ordinal dates
        img_dates = []

        # read metadata and
        for meta in meta_files:
            metadata = open(join(block_folder, meta))
            config = json.load(metadata)
            ordinal_date = config['ordinal_date']
            img_dates.append(ordinal_date)
        sorted_img_dates = sorted(img_dates)

        if year_lowbound > 0:
            year_low_ordinal = pd.Timestamp.toordinal(datetime(int(year_lowbound), 1, 1))
            img_dates = (lambda x: x >= year_low_ordinal, img_dates)
        if year_highbound > 0:
            year_high_ordinal = pd.Timestamp.toordinal(datetime(int(year_highbound + 1), 1, 1))
            img_dates = (lambda x: x < year_high_ordinal, img_dates)

    except IOError:
        raise
    else:
        starting_date = sorted_img_dates[0]
        ending_date = sorted_img_dates[-1]
        n_cm_maps = int((ending_date - starting_date + 1) / cm_interval) + 1
        return starting_date, n_cm_maps

def is_finished_assemble_cmmaps(cmmap_path, n_cm, starting_date, cm_interval):
    """
    Parameters
    ----------
    cmmap_path: the path for saving change magnitude maps
    n_cm: the number of change magnitudes outputted per pixel
    starting_date: the starting date of the whole dataset
    cm_interval: the day interval for outputting change magnitudes
    Returns: boolean
    -------
    True -> assemble finished
    """
    for count in range(n_cm):
        ordinal_date = starting_date + count * cm_interval
        if not os.path.exists(join(cmmap_path,
                                   'CM_maps_{}_{}{}.npy'.format(str(ordinal_date),
                                                                pd.Timestamp.fromordinal(ordinal_date).year,
                                                                get_doy(ordinal_date)))):
            return False
        if not os.path.exists(join(cmmap_path,
                                   'CM_date_maps_{}_{}{}.npy'.format(str(ordinal_date),
                                                                     pd.Timestamp.fromordinal(ordinal_date).year,
                                                                     get_doy(ordinal_date)))):
            return False
    return True

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
    # rank = 1
    # n_cores = 1
    # stack_path = '~/kwcoco/stacked/US_C000'
    # result_path = '~/kwcoco/COLD/US_C000'
    # yaml_path = '/home/jws18003/Document/pycold-uconnhpc/config_watch.yaml'
    # method = 'COLD'
    # year_lowbound = 2017
    # year_highbound = 2022
    # b_c2 = True

    rank = 1
    n_cores = 1
    # stack_path = '/gpfs/scratchfs1/zhz18039/jws18003/kwcoco/stacked/ASI/US_C000'
    stack_path = '/gpfs/scratchfs1/zhz18039/jws18003/kwcoco/stacked/US_C000_metadata'
    # result_path = '/gpfs/scratchfs1/zhz18039/jws18003/kwcoco/COLD_result/ASI/US_C000_prob099_conse6' ##"/gpfs/scratchfs1/zhz18039/jws18003/kwcoco/COLD_KR_R002_drop4_115034_no_boundary"
    result_path = '/gpfs/scratchfs1/zhz18039/jws18003/kwcoco/COLD_result/COLD_parameter_test/OBCOLD_US_C000_prob099_conse4'
    yaml_path = '/home/jws18003/Document/pycold-uconnhpc/config_watch.yaml'
    method = 'OBCOLD'
    seedmap_path = None
    year_lowbound = None
    year_highbound = None
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
              'conse': 6,
              'CM_OUTPUT_INTERVAL': 60
              }

    with open(yaml_path, 'r') as yaml_obj:
        config = yaml.safe_load(yaml_obj)

    conse = int(config['conse'])
    cm_output_interval = int(config['CM_OUTPUT_INTERVAL'])
    probability_threshold = config['probability_threshold']
    block_width = int(config['n_cols'] / config['n_block_x'])  # width of a block
    block_height = int(config['n_rows'] / config['n_block_y'])  # height of a block
    if (config['n_cols'] % block_width != 0) or (config['n_rows'] % block_height != 0):
        print('n_cols, n_rows must be divisible respectively by block_width, block_height! Please double '
              'check your config yaml')
        exit()

    # set up some additional config
    block_width = int(config['n_cols'] / config['n_block_x'])  # width of a block
    block_height = int(config['n_rows'] / config['n_block_y'])  # height of a block
    if (config['n_cols'] % block_width != 0) or (config['n_rows'] % block_height != 0):
        print('n_cols, n_rows must be divisible respectively by block_width, block_height! Please double '
              'check your config yaml')
        exit()

    # set up additional parameters for obcold
    if method == 'OBCOLD':
        # we need read 'global starting date' to save CM which will be only used for ob-cold
        try:
            starting_date, n_cm_maps = reading_start_dates_nmaps(stack_path, year_lowbound, year_highbound,
                                                                 config['CM_OUTPUT_INTERVAL'])
            year_lowbound = pd.Timestamp.fromordinal(starting_date).year
            year_highbound = pd.Timestamp.fromordinal(
                starting_date + (n_cm_maps - 1) * config['CM_OUTPUT_INTERVAL']).year
        except IOError:
            print("reading start dates errors: {}".format(start_time.strftime('%Y-%m-%d %H:%M:%S')))
            exit()

    # logging and folder preparation
    if rank == 1:
        if not os.path.exists(result_path):
            os.makedirs(result_path)
        if method == 'OBCOLD':
            if not os.path.exists(join(result_path, 'cm_maps')):
                os.makedirs(join(result_path, 'cm_maps'))
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

        if img_tstack is None:  # empty block
            if method == 'OBCOLD':
                for pos in range(block_width * block_height):
                    CM_collect.append(np.full(n_cm_maps, -9999, dtype=np.short))
                    date_collect.append(np.full(n_cm_maps, -9999, dtype=np.short))
                np.save(join(result_path, 'CM_date_x{}_y{}.npy'.format(block_x, block_y)), np.hstack(date_collect))
                np.save(join(result_path, 'CM_x{}_y{}.npy'.format(block_x, block_y)), np.hstack(CM_collect))

        else:
            # start looping every pixel in the block
            if method == "COLD" or method == "HybridCOLD" or method == "OBCOLD":
                for pos in range(block_width * block_height):
                    original_row, original_col = get_rowcol_intile(pos, block_width,
                                                                   block_height, block_x, block_y)
                    try:
                        if method == 'OBCOLD':
                            [cold_result, CM, CM_date] = cold_detect(img_dates_sorted,
                                                                     img_tstack[pos, 0, :].astype(np.int64),
                                                                     img_tstack[pos, 1, :].astype(np.int64),
                                                                     img_tstack[pos, 2, :].astype(np.int64),
                                                                     img_tstack[pos, 3, :].astype(np.int64),
                                                                     img_tstack[pos, 4, :].astype(np.int64),
                                                                     img_tstack[pos, 5, :].astype(np.int64),
                                                                     img_tstack[pos, 6, :].astype(np.int64),
                                                                     img_tstack[pos, 7, :].astype(np.int64),
                                                                     pos=config['n_cols'] * (original_row - 1) +
                                                                         original_col,
                                                                     conse=config['conse'],
                                                                     starting_date=starting_date,
                                                                     n_cm=n_cm_maps, b_c2=b_c2,
                                                                     cm_output_interval=config['CM_OUTPUT_INTERVAL'],
                                                                     b_output_cm=True)
                        else:
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

                    except RuntimeError:
                        print("COLD fails at original_row {}, original_col {} ({})".format(original_row, original_col,
                                                                                           datetime.now(tz)
                                                                                           .strftime('%Y-%m-%d %H:%M:%S')))
                    except Exception:
                        if method == 'OBCOLD':
                            CM = np.full(n_cm_maps, -9999, dtype=np.short)
                            CM_date = np.full(n_cm_maps, -9999, dtype=np.short)
                    else:
                        result_collect.append(cold_result)
                    finally:
                        if method == 'OBCOLD':
                            CM_collect.append(CM)
                            date_collect.append(CM_date)

                # save the dataset
                if len(result_collect) > 0:
                    if method == 'HybridCOLD':
                        np.save(join(result_path, 'record_change_x{}_y{}_hybridcold.npy'.format(block_x, block_y)),
                                np.hstack(result_collect))
                    else:
                        np.save(join(result_path, 'record_change_x{}_y{}_cold.npy'.format(block_x, block_y)),
                                np.hstack(result_collect))
                if method == 'OBCOLD':
                    np.save(join(result_path, 'CM_date_x{}_y{}.npy'.format(block_x, block_y)), np.hstack(date_collect))
                    np.save(join(result_path, 'CM_x{}_y{}.npy'.format(block_x, block_y)), np.hstack(CM_collect))

            with open(join(result_path, 'COLD_block{}_finished.txt'.format(block_id)), 'w'):
                pass

            print("Per-pixel COLD processing is finished for block_x{}_y{} ({})"
                  .format(block_x, block_y, datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

    # wait for all cores to be finished
    if method == 'OBCOLD':
        while not is_finished_cold_blockfinished(result_path, config['n_block_x'] * config['n_block_y']):
            time.sleep(30)

    if rank == 1:
        cold_timepoint = datetime.now(tz)

    #################################################################################
    #                        the below is object-based process                      #
    #################################################################################
    if method == 'OBCOLD':
        if seedmap_path is None:
            ob_analyst = ObjectAnalystHPC(config, starting_date=starting_date, stack_path=stack_path,
                                          result_path=result_path)
        else:
            pyclassifier = PyClassifierHPC(config, record_path=result_path, year_lowbound=year_lowbound,
                                           year_uppbound=year_highbound, seedmap_path=seedmap_path)
            ob_analyst = ObjectAnalystHPC(config, starting_date=starting_date, stack_path=stack_path,
                                          result_path=result_path, thematic_path=join(result_path, 'feature_maps'))
        if rank == 1:
            # need to create folders first
            if seedmap_path is not None:
                pyclassifier.hpc_preparation()
            ob_analyst.hpc_preparation()

        #########################################################################
        #                        reorganize cm snapshots                        #
        #########################################################################

        if not is_finished_assemble_cmmaps(join(result_path, 'cm_maps'), n_cm_maps,
                                           starting_date, config['CM_OUTPUT_INTERVAL']):
            if rank == 1:
                assemble_cmmaps(config, result_path, join(result_path, 'cm_maps'), starting_date, n_cm_maps, 'CM',
                                clean=False)
            elif rank == 2:
                assemble_cmmaps(config, result_path, join(result_path, 'cm_maps'), starting_date, n_cm_maps,
                                'CM_date',
                                clean=False)

            while not is_finished_assemble_cmmaps(join(result_path, 'cm_maps'), n_cm_maps,
                                                  starting_date, config['CM_OUTPUT_INTERVAL']):
                time.sleep(15)

        #########################################################################
        #                      producing classification maps                    #
        #########################################################################
        if seedmap_path is not None:  # we used thematic info
            if not pyclassifier.is_finished_step4_assemble():
                if rank == 1:
                    print("Starts predicting features: {}".format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))
                for i in range(nblock_eachcore):
                    if n_cores * i + rank > config['n_block_x'] * config['n_block_y']:
                        break
                    pyclassifier.step1_feature_generation(block_id=n_cores * i + rank)

                if rank == 1:  # serial mode for producing rf
                    pyclassifier.step2_train_rf()
                    print("Training rf ends: {}".format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

                for i in range(nblock_eachcore):
                    if n_cores * i + rank > config['n_block_x'] * config['n_block_y']:
                        break
                    pyclassifier.step3_classification(block_id=n_cores * i + rank)

                if rank == 1:  # serial mode for assemble
                    pyclassifier.step4_assemble()
            while not pyclassifier.is_finished_step4_assemble():
                time.sleep(15)
            if rank == 1:
                print("Assemble classification map ends: {}".format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))
        #########################################################################
        #                      object-based image analysis                      #
        #########################################################################
        if not ob_analyst.is_finished_object_analysis(np.arange(starting_date,
                                                                starting_date + config[
                                                                    'CM_OUTPUT_INTERVAL'] * n_cm_maps,
                                                                config['CM_OUTPUT_INTERVAL'])):
            n_map_percore = int(np.ceil(n_cm_maps / n_cores))
            max_date = starting_date + (n_cm_maps - 1) * config['CM_OUTPUT_INTERVAL']
            for i in range(n_map_percore):
                if starting_date + (rank - 1 + i * n_cores) * config['CM_OUTPUT_INTERVAL'] > max_date:
                    break
                date = starting_date + (rank - 1 + i * n_cores) * config['CM_OUTPUT_INTERVAL']
                ob_analyst.obia_execute(date)

            while not ob_analyst.is_finished_object_analysis(np.arange(starting_date,
                                                                       starting_date + config[
                                                                           'CM_OUTPUT_INTERVAL'] * n_cm_maps,
                                                                       config['CM_OUTPUT_INTERVAL'])):
                time.sleep(15)
        if rank == 1:
            print("OBIA ends: {}".format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

        #########################################################################
        #                        reconstruct change records                     #
        #########################################################################
        for i in range(nblock_eachcore):
            block_id = n_cores * i + rank  # started from 1, i.e., rank, rank + n_cores, rank + 2 * n_cores
            if block_id > config['n_block_x'] * config['n_block_y']:
                break
            block_y = int((block_id - 1) / config['n_block_y']) + 1  # note that block_x and block_y start from 1
            block_x = int((block_id - 1) % config['n_block_x']) + 1
            img_tstack, img_dates_sorted = get_stack_date(block_x, block_y, stack_path)
            result_collect = ob_analyst.reconstruct_reccg(block_id=block_id,
                                                          img_stack=img_tstack,
                                                          img_dates_sorted=img_dates_sorted)
            ob_analyst.save_obcoldrecords(block_id=block_id, result_collect=result_collect)

    if rank == 1:
        # tile_based report
        if method == 'OBCOLD':
            tileprocessing_report(join(result_path, 'tile_processing_report.log'),
                                  stack_path, pycold.__version__, method, config, start_time, cold_timepoint, tz,
                                  n_cores, starting_date, n_cm_maps, year_lowbound, year_highbound)
        else:
            tileprocessing_report(join(result_path, 'tile_processing_report.log'), stack_path, pycold.__version__,
                                  method, config, start_time, cold_timepoint, tz, n_cores)
        print("The whole procedure finished: {}".format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

if __name__ == '__main__':
    main()
