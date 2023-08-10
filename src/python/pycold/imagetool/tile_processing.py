# Author: Su Ye
# This script is an example for running pyscold in UCONN job array environment
# Due to the independence features of job array, we use writing disk files for communication
# Two types of logging files are used: 1) print() to save block-based  processing info into individual slurm file;
# 2) tile_processing_report.log records config for each tile
# 2) is only called when rank == 1
import yaml
import os
from os.path import join
import pandas as pd
import datetime as dt
import numpy as np
import pycold
from datetime import datetime
from pytz import timezone
import click
import time
from pycold import cold_detect, sccd_detect
from scipy.stats import chi2
import pickle
from dateutil.parser import parse

from pycold.utils import (
    assemble_cmmaps,
    get_rowcol_intile,
    get_doy,
    unindex_sccdpack,
    class_from_dict,
)
from pycold.ob_analyst import ObjectAnalystHPC
from pycold.pyclassifier import PyClassifierHPC
from pycold.app import defaults
from pycold.common import DatasetInfo

TZ = timezone("US/Eastern")


def tileprocessing_report(
    result_log_path,
    stack_path,
    version,
    algorithm,
    dataset_info,
    startpoint,
    cold_timepoint,
    tz,
    n_cores,
    probability_threshold,
    conse,
    starting_date=0,
    n_cm_maps=0,
    year_lowbound=0,
    year_uppbound=0,
):
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
    dataset_info: dictionary structure
    startpoint: a time point, when the program starts
    tz: string, time zone
    n_cores: the core number used
    starting_date: the first date of the total dataset
    n_cm_maps: the number of snapshots
    year_lowbound: the low bound of year range
    year_uppbound: the upper bound of year range
    Returns
    -------

    """
    endpoint = datetime.now(TZ)
    file = open(result_log_path, "w")
    file.write("PYCOLD V{} \n".format(version))
    file.write("Author: Su Ye(remoteseningsuy@gmail.com)\n")
    file.write("Algorithm: {} \n".format(algorithm))
    file.write("Starting_time: {}\n".format(startpoint.strftime("%Y-%m-%d %H:%M:%S")))
    file.write("Change probability threshold: {}\n".format(probability_threshold))
    file.write("Conse: {}\n".format(conse))
    file.write("stack_path: {}\n".format(stack_path))
    file.write("The number of requested cores: {}\n".format(n_cores))
    file.write("The program starts at {}\n".format(startpoint.strftime("%Y-%m-%d %H:%M:%S")))
    file.write("The COLD ends at {}\n".format(cold_timepoint.strftime("%Y-%m-%d %H:%M:%S")))
    file.write("The program ends at {}\n".format(endpoint.strftime("%Y-%m-%d %H:%M:%S")))
    file.write(
        "The program lasts for {:.2f}mins\n".format(
            (endpoint - startpoint) / dt.timedelta(minutes=1)
        )
    )
    if algorithm == "OBCOLD":
        file.write("Starting date of the dataset: {}\n".format(starting_date))
        file.write("Number of snapshots being processed: {}\n".format(n_cm_maps))
        file.write("Lower bound of year range: {}\n".format(year_lowbound))
        file.write("Upper bound of year range: {}\n".format(year_uppbound))
        file.write("Land-cover-specific config:\n")
        file.write("  C1_threshold: {}\n".format(defaults["OBCOLD"]["C1_threshold"]))
        file.write("  C1_sizeslope: {}\n".format(defaults["OBCOLD"]["C1_sizeslope"]))
        file.write("  C2_threshold: {}\n".format(defaults["OBCOLD"]["C2_threshold"]))
        file.write("  C2_sizeslope: {}\n".format(defaults["OBCOLD"]["C2_sizeslope"]))
        file.write("  C3_threshold: {}\n".format(defaults["OBCOLD"]["C3_threshold"]))
        file.write("  C3_sizeslope: {}\n".format(defaults["OBCOLD"]["C3_sizeslope"]))
        file.write("  C4_threshold: {}\n".format(defaults["OBCOLD"]["C4_threshold"]))
        file.write("  C4_sizeslope: {}\n".format(defaults["OBCOLD"]["C4_sizeslope"]))
        file.write("  C5_threshold: {}\n".format(defaults["OBCOLD"]["C5_threshold"]))
        file.write("  C5_sizeslope: {}\n".format(defaults["OBCOLD"]["C5_sizeslope"]))
        file.write("  C6_threshold: {}\n".format(defaults["OBCOLD"]["C6_threshold"]))
        file.write("  C6_sizeslope: {}\n".format(defaults["OBCOLD"]["C6_sizeslope"]))
        file.write("  C7_threshold: {}\n".format(defaults["OBCOLD"]["C7_threshold"]))
        file.write("  C7_sizeslope: {}\n".format(defaults["OBCOLD"]["C7_sizeslope"]))
        file.write("  C8_threshold: {}\n".format(defaults["OBCOLD"]["C8_threshold"]))
        file.write("  C8_sizeslope: {}\n".format(defaults["OBCOLD"]["C8_sizeslope"]))
    file.close()


def reading_start_dates_nmaps(stack_path, cm_interval):
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
        # read starting and ending date info from stack files
        f = open(join(stack_path, "starting_last_dates.txt"), "r")
    except IOError:
        raise
    else:
        starting_date = int(f.readline().rstrip("\n"))
        ending_date = int(f.readline().rstrip("\n"))
        n_cm_maps = int((ending_date - starting_date + 1) / cm_interval) + 1
        f.close()
        return starting_date, n_cm_maps


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
        if not os.path.exists(os.path.join(result_path, "COLD_block{}_finished.txt".format(n + 1))):
            return False
    return True


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
        if not os.path.exists(
            join(
                cmmap_path,
                "CM_maps_{}_{}{}.npy".format(
                    str(ordinal_date),
                    pd.Timestamp.fromordinal(ordinal_date).year,
                    get_doy(ordinal_date),
                ),
            )
        ):
            return False
        if not os.path.exists(
            join(
                cmmap_path,
                "CM_date_maps_{}_{}{}.npy".format(
                    str(ordinal_date),
                    pd.Timestamp.fromordinal(ordinal_date).year,
                    get_doy(ordinal_date),
                ),
            )
        ):
            return False
    return True


def get_stack_date(
    dataset_info: DatasetInfo,
    block_x: int,
    block_y: int,
    stack_path: str,
    low_datebound: int = 0,
    high_datebound: int = 0,
    nband: int = 8,
):
    """
    :param dataset_info:
    :param block_x: block id at x axis
    :param block_y: block id at y axis
    :param stack_path: stack path
    :param low_datebound: ordinal data of low bounds of selection date range
    :param high_datebound: ordinal data of upper bounds of selection date range
    :return:
        img_tstack, img_dates_sorted
        img_tstack - 3-d array (dataset_info.block_width * dataset_info.block_height, nband, nimage)
    """
    dataset_info.block_width = int(dataset_info.n_cols / dataset_info.n_block_x)  # width of a block
    dataset_info.block_height = int(
        dataset_info.n_rows / dataset_info.n_block_y
    )  # height of a block
    block_folder = join(stack_path, "block_x{}_y{}".format(block_x, block_y))
    img_files = [f for f in os.listdir(block_folder) if f.startswith("L") or f.startswith("S")]
    if len(img_files) == 0:
        return None, None
    # sort image files by dates
    img_dates = [
        pd.Timestamp.toordinal(
            dt.datetime(int(folder_name[9:13]), 1, 1) + dt.timedelta(int(folder_name[13:16]) - 1)
        )
        for folder_name in img_files
    ]
    if low_datebound > 0:
        img_dates, img_files = zip(
            *filter(lambda x: x[0] >= low_datebound, zip(img_dates, img_files))
        )
    if high_datebound > 0:
        img_dates, img_files = zip(
            *filter(lambda x: x[0] <= high_datebound, zip(img_dates, img_files))
        )
    files_date_zip = sorted(zip(img_dates, img_files))
    img_files_sorted = [x[1] for x in files_date_zip]
    img_dates_sorted = np.asarray([x[0] for x in files_date_zip])
    img_tstack = np.dstack(
        [
            np.load(join(block_folder, f)).reshape(
                dataset_info.block_width * dataset_info.block_height, nband
            )
            for f in img_files_sorted
        ]
    )
    return img_tstack, img_dates_sorted


@click.command()
@click.option("--rank", type=int, default=0, help="the rank id")
@click.option("--n_cores", type=int, default=1, help="the total cores assigned")
@click.option("--stack_path", type=str, help="the path for stack data")
@click.option("--result_path", type=str, help="the path for storing results")
@click.option("--yaml_path", type=str, help="YAML path")
@click.option(
    "--seedmap_path",
    type=str,
    default=None,
    help="an existing label map path; " "none means not using thematic info",
)
@click.option(
    "--method",
    type=click.Choice(["COLD", "OBCOLD", "SCCDOFFLINE"]),
    help="COLD, OBCOLD, SCCD-OFFLINE",
)
@click.option(
    "--low_datebound",
    type=str,
    default=None,
    help="low date bound of image selection for processing. " "Example - 2015-01-01",
)
@click.option(
    "--upper_datebound",
    type=str,
    default=None,
    help="upper date bound of image selection for processing." "Example - 2021-12-31",
)
@click.option(
    "--b_c2",
    is_flag=True,
    show_default=True,
    default=False,
    help="indicate if it is c2 or hls; if so, "
    "thermal bands will not be considered for observation selection ",
)
@click.option(
    "--seedmap_year",
    type=int,
    show_default=True,
    default=None,
    help="the year of coefficients to generate classification model. If none, it will use the default from pycold/src/python/pycold/ob_parameters.yaml",
)
def main(
    rank,
    n_cores,
    stack_path,
    result_path,
    yaml_path,
    method,
    seedmap_path,
    low_datebound,
    upper_datebound,
    b_c2,
    seedmap_year,
):
    # def main():
    #     rank = 29
    #     n_cores = 200
    #     stack_path = '/gpfs/sharedfs1/zhulab/Jiwon/HLS/Stack/10SEG/10SEG_stack_0.2'
    #     result_path = '/scratch/suy20004/suy20004/test'
    #     yaml_path = '/home/suy20004/Document/pycold-uconnhpc/config_hls.yaml'
    #     method = 'COLD'
    #     seedmap_path = None
    #     low_datebound = None
    #     upper_datebound = None
    #     b_c2 =True
    start_time = datetime.now(TZ)

    if low_datebound is None:
        low_datebound = 0
    else:
        low_datebound = parse(low_datebound).toordinal()

    if upper_datebound is None:
        upper_datebound = 0
    else:
        upper_datebound = parse(upper_datebound).toordinal()

    # Reading config
    with open(yaml_path, "r") as yaml_obj:
        config_general = yaml.safe_load(yaml_obj)
    dataset_info = class_from_dict(DatasetInfo, config_general["DATASETINFO"])
    conse = int(config_general["ALGORITHMINFO"]["conse"])
    cm_output_interval = int(config_general["ALGORITHMINFO"]["CM_OUTPUT_INTERVAL"])
    probability_threshold = config_general["ALGORITHMINFO"]["probability_threshold"]

    if (dataset_info.n_cols % dataset_info.block_width != 0) or (
        dataset_info.n_rows % dataset_info.block_height != 0
    ):
        print(
            "n_cols, n_rows must be divisible respectively by dataset_info.block_width, dataset_info.block_height! Please double "
            "check your config yaml"
        )
        exit()

    # set up additional parameters for obcold
    if method == "OBCOLD":
        # we need read 'global starting date' to save CM which will be only used for ob-cold
        try:
            starting_date, n_cm_maps = reading_start_dates_nmaps(stack_path, cm_output_interval)
            year_lowbound = pd.Timestamp.fromordinal(starting_date).year
            year_uppbound = pd.Timestamp.fromordinal(
                starting_date + (n_cm_maps - 1) * cm_output_interval
            ).year
        except IOError:
            print("reading start dates errors: {}".format(start_time.strftime("%Y-%m-%d %H:%M:%S")))
            exit()

    # logging and folder preparation
    if rank == 1:
        if not os.path.exists(result_path):
            os.makedirs(result_path)
        if method == "OBCOLD":
            if not os.path.exists(join(result_path, "cm_maps")):
                os.makedirs(join(result_path, "cm_maps"))
        print(
            "The per-pixel time series processing begins: {}".format(
                start_time.strftime("%Y-%m-%d %H:%M:%S")
            )
        )

        if not os.path.exists(stack_path):
            print(
                "Failed to locate stack folders. The program ends: {}".format(
                    datetime.now(TZ).strftime("%Y-%m-%d %H:%M:%S")
                )
            )
            return

    #########################################################################
    #                        per-pixel COLD procedure                       #
    #########################################################################
    threshold = chi2.ppf(probability_threshold, 5)
    nblock_eachcore = int(np.ceil(dataset_info.n_block_x * dataset_info.n_block_y * 1.0 / n_cores))
    for i in range(nblock_eachcore):
        block_id = (
            n_cores * i + rank
        )  # started from 1, i.e., rank, rank + n_cores, rank + 2 * n_cores
        if block_id > dataset_info.n_block_x * dataset_info.n_block_y:
            break
        # note that block_x and block_y start from 1
        block_y = int((block_id - 1) / dataset_info.n_block_x) + 1
        block_x = int((block_id - 1) % dataset_info.n_block_x) + 1
        if os.path.exists(join(result_path, "COLD_block{}_finished.txt".format(block_id))):
            print(
                "Per-pixel COLD processing is finished for block_x{}_y{} ({})".format(
                    block_x, block_y, datetime.now(TZ).strftime("%Y-%m-%d %H:%M:%S")
                )
            )
            continue
        img_tstack, img_dates_sorted = get_stack_date(
            dataset_info, block_x, block_y, stack_path, low_datebound, upper_datebound
        )

        # initialize a list (may better change it to generator for future)
        result_collect = []
        CM_collect = []
        date_collect = []

        if img_tstack is None:  # empty block
            if method == "OBCOLD":
                for pos in range(dataset_info.block_width * dataset_info.block_height):
                    CM_collect.append(np.full(n_cm_maps, -9999, dtype=np.short))
                    date_collect.append(np.full(n_cm_maps, -9999, dtype=np.short))
                np.save(
                    join(result_path, "CM_date_x{}_y{}.npy".format(block_x, block_y)),
                    np.hstack(date_collect),
                )
                np.save(
                    join(result_path, "CM_x{}_y{}.npy".format(block_x, block_y)),
                    np.hstack(CM_collect),
                )
        else:
            # for sccd, as the output is heterogeneous, we continuously save the sccd pack for each pixel
            if method == "SCCDOFFLINE":
                # block_status = np.full((dataset_info.block_width, dataset_info.block_height), 0, dtype=np.int8)
                # block_last_change_date = np.full((dataset_info.block_width, dataset_info.block_height), 0, dtype=np.int32)
                f = open(
                    join(result_path, "record_change_x{}_y{}_sccd.npy".format(block_x, block_y)),
                    "wb+",
                )
                # start looping every pixel in the block
                for pos in range(dataset_info.block_width * dataset_info.block_height):
                    original_row, original_col = get_rowcol_intile(
                        pos, dataset_info.block_width, dataset_info.block_height, block_x, block_y
                    )
                    try:
                        sccd_result = sccd_detect(
                            img_dates_sorted,
                            img_tstack[pos, 0, :].astype(np.int64),
                            img_tstack[pos, 1, :].astype(np.int64),
                            img_tstack[pos, 2, :].astype(np.int64),
                            img_tstack[pos, 3, :].astype(np.int64),
                            img_tstack[pos, 4, :].astype(np.int64),
                            img_tstack[pos, 5, :].astype(np.int64),
                            img_tstack[pos, 6, :].astype(np.int64),
                            img_tstack[pos, 7, :].astype(np.int64),
                            t_cg=threshold,
                            conse=conse,
                            b_c2=b_c2,
                            pos=dataset_info.n_cols * (original_row - 1) + original_col,
                        )
                    except RuntimeError:
                        print(
                            "S-CCD fails at original_row {}, original_col {} ({})".format(
                                original_row,
                                original_col,
                                datetime.now(TZ).strftime("%Y-%m-%d %H:%M:%S"),
                            )
                        )
                    else:
                        # replace structural array to list for saving storage space
                        pickle.dump(unindex_sccdpack(sccd_result), f)
                f.close()

                # from pycold.utils import save_obs2csv
                # save_obs2csv('/home/coloury/Dropbox/UCONN/NRT/test_src/sccd_packid{}.csv'.format(pos),
                #              pd.DataFrame(np.vstack((img_dates_sorted, img_tstack[pos, :, :])).T))

                # free memory
                del img_tstack
                del img_dates_sorted
                # del block_status
                # del block_last_change_date
            else:
                # start looping every pixel in the block
                # for pos in [17682]:
                for pos in range(dataset_info.block_width * dataset_info.block_height):
                    original_row, original_col = get_rowcol_intile(
                        pos, dataset_info.block_width, dataset_info.block_height, block_x, block_y
                    )
                    try:
                        if method == "OBCOLD":
                            [cold_result, CM, CM_date] = cold_detect(
                                img_dates_sorted,
                                img_tstack[pos, 0, :].astype(np.int64),
                                img_tstack[pos, 1, :].astype(np.int64),
                                img_tstack[pos, 2, :].astype(np.int64),
                                img_tstack[pos, 3, :].astype(np.int64),
                                img_tstack[pos, 4, :].astype(np.int64),
                                img_tstack[pos, 5, :].astype(np.int64),
                                img_tstack[pos, 6, :].astype(np.int64),
                                img_tstack[pos, 7, :].astype(np.int64),
                                pos=dataset_info.n_cols * (original_row - 1) + original_col,
                                conse=conse,
                                starting_date=starting_date,
                                n_cm=n_cm_maps,
                                b_c2=b_c2,
                                cm_output_interval=cm_output_interval,
                                b_output_cm=True,
                            )
                        else:
                            cold_result = cold_detect(
                                img_dates_sorted,
                                img_tstack[pos, 0, :].astype(np.int64),
                                img_tstack[pos, 1, :].astype(np.int64),
                                img_tstack[pos, 2, :].astype(np.int64),
                                img_tstack[pos, 3, :].astype(np.int64),
                                img_tstack[pos, 4, :].astype(np.int64),
                                img_tstack[pos, 5, :].astype(np.int64),
                                img_tstack[pos, 6, :].astype(np.int64),
                                img_tstack[pos, 7, :].astype(np.int64),
                                t_cg=threshold,
                                conse=conse,
                                b_c2=b_c2,
                                pos=dataset_info.n_cols * (original_row - 1) + original_col,
                            )
                    except RuntimeError:
                        print(
                            "COLD fails at original_row {}, original_col {} ({})".format(
                                original_row,
                                original_col,
                                datetime.now.strftime("%Y-%m-%d %H:%M:%S"),
                            )
                        )
                    except Exception:
                        if method == "OBCOLD":
                            CM = np.full(n_cm_maps, -9999, dtype=np.short)
                            CM_date = np.full(n_cm_maps, -9999, dtype=np.short)
                    else:
                        result_collect.append(cold_result)
                    finally:
                        if method == "OBCOLD":
                            CM_collect.append(CM)
                            date_collect.append(CM_date)

                # save the dataset
                if len(result_collect) > 0:
                    np.save(
                        join(
                            result_path, "record_change_x{}_y{}_cold.npy".format(block_x, block_y)
                        ),
                        np.hstack(result_collect),
                    )

                if method == "OBCOLD":
                    np.save(
                        join(result_path, "CM_date_x{}_y{}.npy".format(block_x, block_y)),
                        np.hstack(date_collect),
                    )
                    np.save(
                        join(result_path, "CM_x{}_y{}.npy".format(block_x, block_y)),
                        np.hstack(CM_collect),
                    )

        with open(join(result_path, "COLD_block{}_finished.txt".format(block_id)), "w"):
            pass

        print(
            "Per-pixel COLD processing is finished for block_x{}_y{} ({})".format(
                block_x, block_y, datetime.now(TZ).strftime("%Y-%m-%d %H:%M:%S")
            )
        )
    # wait for all cores to be finished
    if method == "OBCOLD":
        while not is_finished_cold_blockfinished(
            result_path, dataset_info.n_block_x * dataset_info.n_block_y
        ):
            time.sleep(30)

    if rank == 1:
        cold_timepoint = datetime.now(TZ)

    #################################################################################
    #                        the below is object-based process                      #
    #################################################################################
    if method == "OBCOLD":
        if seedmap_path is None:
            ob_analyst = ObjectAnalystHPC(
                dataset_info,
                starting_date=starting_date,
                stack_path=stack_path,
                result_path=result_path,
            )
        else:
            pyclassifier = PyClassifierHPC(
                dataset_info,
                record_path=result_path,
                year_list_to_predict=list(range(year_lowbound, year_uppbound)),
                seedmap_path=seedmap_path,
            )
            ob_analyst = ObjectAnalystHPC(
                dataset_info,
                starting_date=starting_date,
                stack_path=stack_path,
                result_path=result_path,
                thematic_path=join(result_path, "feature_maps"),
            )
        if rank == 1:
            # need to create folders first
            if seedmap_path is not None:
                pyclassifier.hpc_preparation()
            ob_analyst.hpc_preparation()
        #########################################################################
        #                        reorganize cm snapshots                        #
        #########################################################################

        if not is_finished_assemble_cmmaps(
            join(result_path, "cm_maps"), n_cm_maps, starting_date, cm_output_interval
        ):
            if rank == 1:
                assemble_cmmaps(
                    dataset_info,
                    result_path,
                    join(result_path, "cm_maps"),
                    starting_date,
                    n_cm_maps,
                    "CM",
                    cm_output_interval,
                    clean=False,
                )
            if rank == 16:
                assemble_cmmaps(
                    dataset_info,
                    result_path,
                    join(result_path, "cm_maps"),
                    starting_date,
                    n_cm_maps,
                    "CM_date",
                    cm_output_interval,
                    clean=False,
                )

            while not is_finished_assemble_cmmaps(
                join(result_path, "cm_maps"), n_cm_maps, starting_date, cm_output_interval
            ):
                time.sleep(15)

        #########################################################################
        #                      producing classification maps                    #
        #########################################################################
        if seedmap_path is not None:  # we used thematic info
            if not pyclassifier.is_finished_step4_assemble():
                if rank == 1:
                    print(
                        "Starts predicting features: {}".format(
                            datetime.now(TZ).strftime("%Y-%m-%d %H:%M:%S")
                        )
                    )
                if not pyclassifier.is_finished_step1_predict_features():
                    for i in range(nblock_eachcore):
                        if n_cores * i + rank > dataset_info.n_block_x * dataset_info.n_block_y:
                            break
                        pyclassifier.step1_feature_generation(block_id=n_cores * i + rank)

                if rank == 1:  # serial mode for producing rf
                    pyclassifier.step2_train_rf(ref_year=seedmap_year)
                    print(
                        "Training rf ends: {}".format(
                            datetime.now(TZ).strftime("%Y-%m-%d %H:%M:%S")
                        )
                    )

                for i in range(nblock_eachcore):
                    if n_cores * i + rank > dataset_info.n_block_x * dataset_info.n_block_y:
                        break
                    pyclassifier.step3_classification(block_id=n_cores * i + rank)

                if rank == 1:  # serial mode for assemble
                    pyclassifier.step4_assemble()
            while not pyclassifier.is_finished_step4_assemble():
                time.sleep(15)
            if rank == 1:
                print(
                    "Assemble classification map ends: {}".format(
                        datetime.now(TZ).strftime("%Y-%m-%d %H:%M:%S")
                    )
                )
        #########################################################################
        #                      object-based image analysis                      #
        #########################################################################
        if not ob_analyst.is_finished_object_analysis(
            np.arange(
                starting_date, starting_date + cm_output_interval * n_cm_maps, cm_output_interval
            )
        ):
            n_map_percore = int(np.ceil(n_cm_maps / n_cores))
            max_date = starting_date + (n_cm_maps - 1) * cm_output_interval
            for i in range(n_map_percore):
                if starting_date + (rank - 1 + i * n_cores) * cm_output_interval > max_date:
                    break
                date = starting_date + (rank - 1 + i * n_cores) * cm_output_interval
                ob_analyst.obia_execute(date, cm_output_interval=cm_output_interval)

            while not ob_analyst.is_finished_object_analysis(
                np.arange(
                    starting_date,
                    starting_date + cm_output_interval * n_cm_maps,
                    cm_output_interval,
                )
            ):
                time.sleep(15)
        if rank == 1:
            print("OBIA ends: {}".format(datetime.now(TZ).strftime("%Y-%m-%d %H:%M:%S")))

        #########################################################################
        #                        reconstruct change records                     #
        #########################################################################
        for i in range(nblock_eachcore):
            block_id = (
                n_cores * i + rank
            )  # started from 1, i.e., rank, rank + n_cores, rank + 2 * n_cores
            if block_id > dataset_info.n_block_x * dataset_info.n_block_y:
                break
            # note that block_x and block_y start from 1
            block_y = int((block_id - 1) / dataset_info.n_block_x) + 1
            block_x = int((block_id - 1) % dataset_info.n_block_x) + 1
            img_tstack, img_dates_sorted = get_stack_date(
                dataset_info, block_x, block_y, stack_path
            )
            result_collect = ob_analyst.reconstruct_reccg(
                block_id=block_id,
                img_stack=img_tstack,
                img_dates_sorted=img_dates_sorted,
                conse=conse,
            )
            ob_analyst.save_obcoldrecords(block_id=block_id, result_collect=result_collect)

    if rank == 1:
        # tile_based report
        if method == "OBCOLD":
            tileprocessing_report(
                join(result_path, "tile_processing_report.log"),
                stack_path,
                pycold.__version__,
                method,
                dataset_info,
                start_time,
                cold_timepoint,
                TZ,
                n_cores,
                probability_threshold,
                conse,
                starting_date,
                n_cm_maps,
                year_lowbound,
                year_uppbound,
            )
        else:
            tileprocessing_report(
                join(result_path, "tile_processing_report.log"),
                stack_path,
                pycold.__version__,
                method,
                dataset_info,
                start_time,
                cold_timepoint,
                TZ,
                n_cores,
                probability_threshold,
                conse,
            )
        print(
            "The whole procedure finished: {}".format(
                datetime.now(TZ).strftime("%Y-%m-%d %H:%M:%S")
            )
        )


if __name__ == "__main__":
    main()
