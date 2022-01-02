# This script is an example for running pycold in a HPC environment
# This script requires 'AutoPrepareDataARD.py' to produce intermediate block-based stack. This scripts assign each
# core to independently process each block, and produce python-native change record formats (.npy) for each block.

# Three types of logging files are used: 1) print() to save block-based  processing info into individual slurm file;
# 2) system.log is mainly used to log system information; 3) tile_processing_report.log records parameters for each tile
# 2) and 3) are only called when rank == 1

# For a 42-year Landsat ARD C1 tile (~3000 images), this script averagely produces ~8G change record
# files, and takes ~1.5 hours to finish if 200 EPYC 7452 cores are used.
import yaml
import os
from os.path import join
import pandas as pd
import datetime as dt
import numpy as np
import pycold
import logging
from datetime import datetime
from pytz import timezone
from datetime import timedelta
import click
import time
from pycold import cold_detect


def get_rowcol_intile(pos, block_width, block_height, block_x, block_y):
    """
    calculate row and col in original images based on pos index and block location
    :param pos: position id of the pixel (i.e., i_row * n_cols + i_col)
    :param block_width: the width of each block
    :param block_height: the height of each block
    :param block_x: block location at x direction
    :param block_y: block location at y direction
    :return: row and col number (starting from 1) in original image (e.g., Landsat ARD 5000*5000)
    """
    original_row = int(pos / block_width + (block_y - 1) * block_height + 1)
    original_col = int(pos % block_width + (block_x - 1) * block_width + 1)
    return original_row, original_col


def tileprocessing_report(result_log_path, stack_path, version, algorithm, params, startpoint, tz):
    """
    output tile-based processing report
    :param result_log_path: outputted log path
    :param stack_path: stack data folder
    :param version: version of current pyscold software
    :param algorithm: COLD, OBCOLD, or S-CCD
    :param params: a structure of inputted parameters
    :param startpoint: timepoint that the program starts
    :param tz: time zone
    :return:
    """
    endpoint = datetime.now(tz)
    file = open(result_log_path, "w")
    file.write("S-COLD V{} \n".format(version))
    file.write("Author: Su Ye(remoteseningsuy@gmail.com)\n")
    file.write("Algorithm: {} \n".format(algorithm))
    file.write("Starting_time: {}\n".format(startpoint.strftime('%Y-%m-%d %H:%M:%S')))
    file.write("Change probability threshold: {}\n".format(params['probability_threshold']))
    file.write("Conse: {}\n".format(params['conse']))
    file.write("First date: {}\n".format(params['starting_date']))
    file.write("n_cm_maps: {}\n".format(params['n_cm_maps']))
    file.write("stack_path: {}\n".format(stack_path))
    file.write("The program starts at {}\n".format(startpoint.strftime('%Y-%m-%d %H:%M:%S')))
    file.write("The program ends at {}\n".format(endpoint.strftime('%Y-%m-%d %H:%M:%S')))
    file.write("The program lasts for {:.2f}mins\n".format((endpoint - startpoint) / timedelta(minutes=1)))
    if algorithm == 'OBCOLD':
        file.write("Land-cover-specific parameters:\n")
        file.write("  C1_threshold: {}\n".format(params['C1_threshold']))
        file.write("  C1_sizeslope: {}\n".format(params['C1_sizeslope']))
        file.write("  C2_threshold: {}\n".format(params['C2_threshold']))
        file.write("  C2_sizeslope: {}\n".format(params['C2_sizeslope']))
        file.write("  C3_threshold: {}\n".format(params['C3_threshold']))
        file.write("  C3_sizeslope: {}\n".format(params['C3_sizeslope']))
        file.write("  C4_threshold: {}\n".format(params['C4_threshold']))
        file.write("  C4_sizeslope: {}\n".format(params['C4_sizeslope']))
        file.write("  C5_threshold: {}\n".format(params['C5_threshold']))
        file.write("  C5_sizeslope: {}\n".format(params['C5_sizeslope']))
        file.write("  C6_threshold: {}\n".format(params['C6_threshold']))
        file.write("  C6_sizeslope: {}\n".format(params['C6_sizeslope']))
        file.write("  C7_threshold: {}\n".format(params['C7_threshold']))
        file.write("  C7_sizeslope: {}\n".format(params['C7_sizeslope']))
        file.write("  C8_threshold: {}\n".format(params['C8_threshold']))
        file.write("  C8_sizeslope: {}\n".format(params['C8_sizeslope']))
    file.close()


def get_time_now(tz):
    """
    :param tz: time zone
    :return: return readable format of current time
    """
    return datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')


def reading_start_end_dates(stack_path):
    """
    :param stack_path: stack_path for saving starting_last_dates.txt
    :return:
    """
    # read starting and ending dates, note that all blocks only has one starting and last date (mainly for obcold)
    try:
        f = open(join(stack_path, "starting_last_dates.txt"),
                 "r")  # read starting and ending date info from stack files
    except IOError as e:
        raise
    else:
        starting_date = int(f.readline().rstrip('\n'))
        ending_date = int(f.readline().rstrip('\n'))
        f.close()
        return starting_date, ending_date


def check_cold_blockfinished(result_path, nblocks):
    """
    check if the COLD algorithm finishes all blocks
    :param result_path: the path to save signal text
    :param nblocks: the block number
    :return:
    """
    for n in range(nblocks):
        if not os.path.exists(os.path.join(result_path, 'COLD_block{}_finished.txt'.format(n+1))):
            return False
    return True


@click.command()
@click.option('--rank', type=int, default=0, help='the rank id')
@click.option('--n_cores', type=int, default=0, help='the total cores assigned')
@click.option('--stack_path', type=str, default=None, help='the path for stack data')
@click.option('--result_path', type=str, default=None, help='the path for storing results')
@click.option('--yaml_path', type=str, default=None, help='YAML path')
@click.option('--method', type=str, default=None, help='COLD, OBCOLD-COLD, OBCOLD-SPATIAL, OBCOLD-COMPLETE'
                                                       'SCCD_OFFLINE')
def main(rank, n_cores, stack_path, result_path, yaml_path, method):
    tz = timezone('US/Eastern')
    start_time = datetime.now(tz)
    if rank == 1:
        if not os.path.exists(result_path):
            os.makedirs(result_path)
        # system logging
        logging.basicConfig(filename='{}/system.log'.format(result_path),
                            filemode='w', level=logging.INFO)
        logger = logging.getLogger(__name__)

        logger.info("The per-pixel COLD algorithm begins: {}".format(start_time.strftime('%Y-%m-%d %H:%M:%S')))

        if not os.path.exists(stack_path):
            logger.error("Failed to locate stack folders. The program ends: {}"
                         .format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))
            exit()

    b_outputCM = False if method == 'COLD' or method == 'SCCD_OFFLINE' else True

    # Reading parameters
    with open(yaml_path, 'r') as yaml_obj:
        parameters = yaml.safe_load(yaml_obj)
    params = parameters['obcold']
    params.update(parameters['common'])

    # set up some additional parameters
    try:
        starting_date, ending_date = reading_start_end_dates(stack_path)
    except IOError as e:
        exit()

    params['block_width'] = int(params['n_cols'] / params['n_block_x'])  # width of a block
    params['block_height'] = int(params['n_rows'] / params['n_block_y'])  # height of a block
    params['starting_date'] = starting_date
    params['n_cm_maps'] = int((ending_date - starting_date + 1) / params['CM_OUTPUT_INTERVAL'])

    nblock_eachcore = int(np.ceil(params['n_block_x'] * params['n_block_y'] / n_cores))
    for i in range(nblock_eachcore):
        block_id = n_cores * i + rank
        if block_id > params['n_block_x'] * params['n_block_y']:
            break

        block_y = int((block_id - 1) / params['n_block_x']) + 1  # note that block_x and block_y start from 1
        block_x = int((block_id - 1) % params['n_block_x']) + 1

        # skip the block if the change record block has been created
        if os.path.exists(join(result_path, 'record_change_x{}_y{}_cold.npy'.format(block_x, block_y))):
            continue

        block_folder = join(stack_path, 'block_x{}_y{}'.format(block_x, block_y))

        img_files = [f for f in os.listdir(block_folder) if f.startswith('L')]

        # sort image files by dates
        img_dates = [pd.Timestamp.toordinal(dt.datetime(int(folder_name[9:13]), 1, 1) +
                                            dt.timedelta(int(folder_name[13:16]) - 1)) + 366
                     for folder_name in img_files]
        files_date_zip = sorted(zip(img_dates, img_files))
        img_files_sorted = [x[1] for x in files_date_zip]
        img_dates_sorted = np.asarray([x[0] for x in files_date_zip])
        img_tstack = [np.load(join(block_folder, f)).reshape(params['block_width'] * params['block_height'],
                                                             params['TOTAL_IMAGE_BANDS'] + 1)
                      for f in img_files_sorted]
        img_tstack = np.dstack(img_tstack)

        # initialize a list (may need to changing it to generator in future)
        result_collect = []
        CM_collect = []
        direction_collect = []
        date_collect = []

        # start looping every pixel in the block
        for pos in range(params['block_width'] * params['block_height']):
            original_row, original_col = get_rowcol_intile(pos, params['block_width'],
                                                           params['block_height'], block_x, block_y)
            try:
                if b_outputCM:
                    [cold_result, CM, CM_direction, CM_date] = cold_detect(img_dates_sorted,
                                                                           img_tstack[pos, 0, :].astype(np.int64),
                                                                           img_tstack[pos, 1, :].astype(np.int64),
                                                                           img_tstack[pos, 2, :].astype(np.int64),
                                                                           img_tstack[pos, 3, :].astype(np.int64),
                                                                           img_tstack[pos, 4, :].astype(np.int64),
                                                                           img_tstack[pos, 5, :].astype(np.int64),
                                                                           img_tstack[pos, 6, :].astype(np.int64),
                                                                           img_tstack[pos, 7, :].astype(np.int64),
                                                                           col_pos=original_col,
                                                                           row_pos=original_row,
                                                                           t_cg=params['default_threshold'],
                                                                           conse=params['conse'],
                                                                           num_samples=params['n_cols'],
                                                                           starting_date=params['starting_date'],
                                                                           CM_OUTPUT_INTERVAL=params[
                                                                               'CM_OUTPUT_INTERVAL'],
                                                                           b_outputCM=b_outputCM)
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
                                              col_pos=original_col,
                                              row_pos=original_row,
                                              t_cg=params['default_threshold'],
                                              conse=params['conse'],
                                              num_samples=params['n_cols'],
                                              starting_date=params['starting_date'],
                                              CM_OUTPUT_INTERVAL=params['CM_OUTPUT_INTERVAL'],
                                              b_outputCM=b_outputCM)
            except RuntimeError:
                print("COLD fails at original_row {}, original_col {} ({})".format(original_row, original_col,
                                                                                   datetime.now(tz)
                                                                                   .strftime('%Y-%m-%d %H:%M:%S')))
            except Exception as e:
                print(e)
            else:
                result_collect.append(cold_result)
                if b_outputCM:
                    CM_collect.append(CM)
                    direction_collect.append(CM_direction)
                    date_collect.append(CM_date)

        # save the dataset
        np.save(join(result_path, 'CM_x{}_y{}.npy'.format(block_x, block_y)), np.hstack(CM_collect))
        if b_outputCM:
            np.save(join(result_path, 'CM_date_x{}_y{}.npy'.format(block_x, block_y)), np.hstack(date_collect))
            np.save(join(result_path, 'CM_direction_x{}_y{}.npy'.format(block_x, block_y)), np.hstack(direction_collect))
            np.save(join(result_path, 'record_change_x{}_y{}_cold.npy'.format(block_x, block_y)), np.hstack(result_collect))

        with open(os.path.join(result_path, 'COLD_block{}_finished.txt'.format(block_id)), 'w') as fp:
            pass

        print("Per-pixel COLD processing is finished for block_x{}_y{} ({})"
              .format(block_x, block_y, datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

    while not check_cold_blockfinished(result_path, params['n_block_x'] * params['n_block_y']):
        time.sleep(15)

    if rank == 1:
        # tile_based report
        tileprocessing_report(join(result_path, 'tile_processing_report.log'),
                              stack_path, pycold.__version__, method, params, start_time, tz)
        # if b_outputCM:
        #     deleted_files = [file for file in os.listdir(result_path) if file.startswith('CM_')]
        #     for file in deleted_files:
        #         os.remove(os.path.join(result_path, file))
        # system log
        logger.info("The per-pixel COLD algorithm ends: {}".format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))


if __name__ == '__main__':
    main()

