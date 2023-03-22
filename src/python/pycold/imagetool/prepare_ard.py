# Author: Su Ye
# This script is an example for generating block-based stack files from original ARD zip as
# intermediate inputs to the COLD algorithm in a HPC environment. As preparation, you need to
# download '_BT' and '_SR' for all Landsat ARD collection 1.
# This script has 4 steps: 1) warp single-path array to limit the observation inputs for each pixel
# if single_path is set True; 2) unzip all images and unpack bit-based QA bands; 3)
# partition each 5000*5000 temporal images to blocks and eliminate those image blocks if no clear,
# water or snow pixel were in it (so to save disk space and enable independent IO for individual
# block in later time-series analysis); 4) save each image block to python-native binary format
# (.npy) into its block folders

# For a 42-year Landsat ARD C1 tile (~3000 images), this script averagely produces ~350 G intermediate disk
# files, and takes ~12 mins to finish if 200 EPYC 7452 cores are used.
import os
import sys
import shutil
import tarfile
import logging
import time
from typing import Optional
import datetime as dt
import xml.etree.ElementTree as ET
from logging import Logger
import yaml
import click
from pytz import timezone
from datetime import datetime
from os import listdir
from os.path import isfile, join, isdir
from pathlib import Path
from importlib import resources as importlib_resources

import pandas as pd
import numpy as geek
import numpy as np
from glob import glob
from dateutil.parser import parse
from osgeo import gdal_array
from osgeo import gdal
import fiona
from pycold.common import DatasetInfo
from pycold.utils import class_from_dict


# define constant here
QA_CLEAR = 0
QA_WATER = 1
QA_SHADOW = 2
QA_SNOW = 3
QA_CLOUD = 4
QA_FILL = 255

QA_CIRRUS_HLS = 0
QA_WATER_HLS = 5
QA_SHADOW_HLS = 3
QA_SNOW_HLS = 4
QA_CLOUDADJACENT_HLS = 2
QA_CLOUD_HLS = 1
res = 30


s2_stack_bands = ["B02", "B03", "B04", "B8A", "B11", "B12", "Fmask"]
l8_stack_bands = ["B02", "B03", "B04", "B05", "B06", "B07", "Fmask"]


def mask_value(vector, val):
    """
    Build a boolean mask around a certain value in the vector.

    Args:
        vector: 1-d ndarray of values
        val: values to mask on
    Returns:
        1-d boolean ndarray
    """
    return vector == val


def qabitval_array_HLS(packedint_array: np.ndarray):
    """
    Institute a hierarchy of qa values that may be flagged in the bitpacked
    value.

    fill > cloud > shadow > snow > water > clear

    Args:
        packedint: int value to bit check
    Returns:
        offset value to use
    """
    unpacked = np.full(packedint_array.shape, 0)
    QA_CLOUD_unpacked = geek.bitwise_and(packedint_array, 1 << 1)
    QA_CLOUD_ADJ = geek.bitwise_and(packedint_array, 1 << 2)
    QA_SHADOW_unpacked = geek.bitwise_and(packedint_array, 1 << 3)
    QA_SNOW_unpacked = geek.bitwise_and(packedint_array, 1 << 4)
    QA_WATER_unpacked = geek.bitwise_and(packedint_array, 1 << 5)

    unpacked[QA_WATER_unpacked > 0] = QA_WATER
    unpacked[QA_SNOW_unpacked > 0] = QA_SNOW
    unpacked[QA_SHADOW_unpacked > 0] = QA_SHADOW
    unpacked[QA_CLOUD_ADJ > 0] = QA_CLOUD
    unpacked[QA_CLOUD_unpacked > 0] = QA_CLOUD
    unpacked[packedint_array == QA_FILL] = QA_FILL

    return unpacked


def qabitval_array(packedint_array: np.ndarray):
    """
    Institute a hierarchy of qa values that may be flagged in the bitpacked
    value.
    fill > cloud > shadow > snow > water > clear
    Args:
        packedint: int value to bit check
    Returns:
        offset value to use
    """
    unpacked = np.full(packedint_array.shape, QA_FILL)
    QA_CLOUD_unpacked = geek.bitwise_and(packedint_array, 1 << (QA_CLOUD + 1))
    QA_SHADOW_unpacked = geek.bitwise_and(packedint_array, 1 << (QA_SHADOW + 1))
    QA_SNOW_unpacked = geek.bitwise_and(packedint_array, 1 << (QA_SNOW + 1))
    QA_WATER_unpacked = geek.bitwise_and(packedint_array, 1 << (QA_WATER + 1))
    QA_CLEAR_unpacked = geek.bitwise_and(packedint_array, 1 << (QA_CLEAR + 1))

    unpacked[QA_CLEAR_unpacked > 0] = QA_CLEAR
    unpacked[QA_WATER_unpacked > 0] = QA_WATER
    unpacked[QA_SNOW_unpacked > 0] = QA_SNOW
    unpacked[QA_SHADOW_unpacked > 0] = QA_SHADOW
    unpacked[QA_CLOUD_unpacked > 0] = QA_CLOUD
    return unpacked


def qabitval_array_c2(packedint_array: np.ndarray):
    """
    Institute a hierarchy of qa values that may be flagged in the bitpacked
    value for c2

    fill > cloud > shadow > snow > water > clear

    Args:
        packedint: int value to bit check
    Returns:
        offset value to use
    """
    unpacked = np.full(packedint_array.shape, QA_FILL)
    QA_CLEAR_unpacked = geek.bitwise_and(packedint_array, 1 << 6)
    QA_SHADOW_unpacked = geek.bitwise_and(packedint_array, 1 << 4)
    QA_CLOUD_unpacked = geek.bitwise_and(packedint_array, 1 << 3)
    QA_DILATED_unpacked = geek.bitwise_and(packedint_array, 1 << 1)
    QA_SNOW_unpacked = geek.bitwise_and(packedint_array, 1 << 5)
    QA_WATER_unpacked = geek.bitwise_and(packedint_array, 1 << 7)

    unpacked[QA_CLEAR_unpacked > 0] = QA_CLEAR
    unpacked[QA_WATER_unpacked > 0] = QA_WATER
    unpacked[QA_SNOW_unpacked > 0] = QA_SNOW
    unpacked[QA_SHADOW_unpacked > 0] = QA_SHADOW
    unpacked[QA_CLOUD_unpacked > 0] = QA_CLOUD
    unpacked[QA_DILATED_unpacked > 0] = QA_CLOUD

    return unpacked


# def load_data(file_name, gdal_driver='GTiff'):
#     '''
#     Converts a GDAL compatable file into a numpy array and associated geodata.
#     The rray is provided so you can run with your processing - the geodata consists of the geotransform and gdal dataset object
#     if you're using an ENVI binary as input, this willr equire an associated .hdr file otherwise this will fail.
#     This needs modifying if you're dealing with multiple bands.

#     VARIABLES
#     file_name : file name and path of your file

#     RETURNS
#     image array
#     (geotransform, inDs)
#     '''
#     driver_t = gdal.GetDriverByName(gdal_driver)  # http://www.gdal.org/formats_list.html
#     driver_t.Register()

#     inDs = gdal.Open(file_name, gdal.GA_ReadOnly)
#     # print(inDs)
#     if inDs is None:
#         print('Couldnt open this file {}'.format(file_name))
#         sys.exit("Try again!")

#     # Extract some info form the inDs
#     geotransform = inDs.GetGeoTransform()

#     # Get the data as a numpy array
#     band = inDs.GetRasterBand(1)
#     cols = inDs.RasterXSize
#     rows = inDs.RasterYSize
#     image_array = band.ReadAsArray(0, 0, cols, rows)

#     return image_array, (geotransform, inDs)


def single_image_stacking_hls(
    source_dir: str,
    out_dir: str,
    logger: Logger,
    dataset_info: DatasetInfo,
    folder: str,
    is_partition: bool = True,
    clear_threshold: float = 0.0,
    low_date_bound: Optional[str] = None,
    upp_date_bound: Optional[str] = None,
):
    """
    unzip single image, convert bit-pack qa to byte value, and save as numpy
    :param source_dir: the parent folder to save image 'folder'
    :param out_dir: the folder to save result
    :param folder: the folder name of image
    :param logger: the handler of logger file
    :param data_info: data info data class
    :param is_partition: True, partition each image into blocks; False, save original size of image
    :param clear_threshold: threshold of clear pixel percentage, if lower than threshold, won't be processed
    :param low_date_bound: the lower date of user interested date range
    :param upp_date_bound: the upper date of user interested date range
    :return:
    """
    try:
        QA_band = gdal_array.LoadFile(join(join(source_dir, folder), "{}.Fmask.tif".format(folder)))
    except ValueError as e:
        # logger.error('Cannot open QA band for {}: {}'.format(folder, e))
        logger.error("Cannot open QA band for {}: {}".format(folder, e))
        return False

    # convertQA = np.vectorize(qabitval)
    QA_band_unpacked = qabitval_array_HLS(QA_band).astype(np.short)
    if clear_threshold > 0:
        clear_ratio = np.sum(
            np.logical_or(QA_band_unpacked == QA_CLEAR, QA_band_unpacked == QA_WATER)
        ) / np.sum(QA_band_unpacked != QA_FILL)
    else:
        clear_ratio = 1

    if clear_ratio > clear_threshold:
        [collection, sensor, tile_id, imagetime, version1, version2] = folder.rsplit(".")
        year = imagetime[0:4]
        doy = imagetime[4:7]
        file_name = sensor + tile_id + year + doy + collection + version1
        if low_date_bound is not None:
            if (dt.datetime(int(year), 1, 1) + dt.timedelta(int(doy) - 1)) < parse(low_date_bound):
                return True
        if upp_date_bound is not None:
            if (dt.datetime(int(year), 1, 1) + dt.timedelta(int(doy) - 1)) > parse(upp_date_bound):
                return True

        if sensor == "L30":
            try:
                B1 = gdal_array.LoadFile(
                    join(join(source_dir, folder), "{}.B02.tif".format(folder))
                )
                B2 = gdal_array.LoadFile(
                    join(join(source_dir, folder), "{}.B03.tif".format(folder))
                )
                B3 = gdal_array.LoadFile(
                    join(join(source_dir, folder), "{}.B04.tif".format(folder))
                )
                B4 = gdal_array.LoadFile(
                    join(join(source_dir, folder), "{}.B05.tif".format(folder))
                )
                B5 = gdal_array.LoadFile(
                    join(join(source_dir, folder), "{}.B06.tif".format(folder))
                )
                B6 = gdal_array.LoadFile(
                    join(join(source_dir, folder), "{}.B07.tif".format(folder))
                )
                B7 = np.full(B6.shape, 0)  # assign zero

            except Exception as e:
                # logger.error('Cannot open spectral bands for {}: {}'.format(folder, e))
                logger.error("Cannot open Landsat bands for {}: {}".format(folder, e))
                return False
        elif sensor == "S30":
            try:
                B1 = gdal_array.LoadFile(
                    join(join(source_dir, folder), "{}.B02.tif".format(folder))
                )
                B2 = gdal_array.LoadFile(
                    join(join(source_dir, folder), "{}.B03.tif".format(folder))
                )
                B3 = gdal_array.LoadFile(
                    join(join(source_dir, folder), "{}.B04.tif".format(folder))
                )
                B4 = gdal_array.LoadFile(
                    join(join(source_dir, folder), "{}.B8A.tif".format(folder))
                )
                B5 = gdal_array.LoadFile(
                    join(join(source_dir, folder), "{}.B11.tif".format(folder))
                )
                B6 = gdal_array.LoadFile(
                    join(join(source_dir, folder), "{}.B12.tif".format(folder))
                )
                B7 = np.full(B6.shape, 0)

            except Exception as e:
                # logger.error('Cannot open spectral bands for {}: {}'.format(folder, e))
                logger.error("Cannot open Landsat bands for {}: {}".format(folder, e))
                return False

        if (
            (B1 is None)
            or (B2 is None)
            or (B3 is None)
            or (B4 is None)
            or (B5 is None)
            or (B6 is None)
        ):
            logger.error("Reading Landsat band fails for {}".format(folder))
            return False

        if is_partition is True:
            # width of a block
            bytesize = 2  # short16 = 2 * byte
            # source: https://towardsdatascience.com/efficiently-splitting-an-image-into-tiles-in-python-using-numpy-d1bf0dd7b6f7
            B1_blocks = np.lib.stride_tricks.as_strided(
                B1,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B2_blocks = np.lib.stride_tricks.as_strided(
                B2,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B3_blocks = np.lib.stride_tricks.as_strided(
                B3,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B4_blocks = np.lib.stride_tricks.as_strided(
                B4,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B5_blocks = np.lib.stride_tricks.as_strided(
                B5,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B6_blocks = np.lib.stride_tricks.as_strided(
                B6,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B7_blocks = np.lib.stride_tricks.as_strided(
                B7,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            QA_blocks = np.lib.stride_tricks.as_strided(
                QA_band_unpacked,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            for i in range(dataset_info.n_block_y):
                for j in range(dataset_info.n_block_x):
                    # check if no valid pixels in the chip, then eliminate
                    qa_unique = np.unique(QA_blocks[i][j])

                    # skip blocks are all cloud, shadow or filled values
                    # in DHTC, we also don't need to save pixel that has qa value of 'QA_CLOUD',
                    # 'QA_SHADOW', or FILLED value (255)
                    if (
                        QA_CLEAR not in qa_unique
                        and QA_WATER not in qa_unique
                        and QA_SNOW not in qa_unique
                    ):
                        continue

                    block_folder = "block_x{}_y{}".format(j + 1, i + 1)
                    np.save(
                        join(join(out_dir, block_folder), file_name),
                        np.dstack(
                            [
                                B1_blocks[i][j],
                                B2_blocks[i][j],
                                B3_blocks[i][j],
                                B4_blocks[i][j],
                                B5_blocks[i][j],
                                B6_blocks[i][j],
                                B7_blocks[i][j],
                                QA_blocks[i][j],
                            ]
                        ).astype(np.int16),
                    )

        else:
            np.save(
                join(out_dir, file_name),
                np.dstack([B1, B2, B3, B4, B5, B6, B7, QA_band_unpacked]).astype(np.int16),
            )
        # scene_list.append(folder_name)
    else:
        # logger.info('Not enough clear observations for {}'.format(folder[0:len(folder) - 3]))
        logger.warn("Not enough clear observations for {}".format(folder))

    return True


def single_image_stacking_hls14(
    out_dir: str,
    logger: Logger,
    dataset_info: DatasetInfo,
    folder: str,
    is_partition: bool = True,
    clear_threshold: int = 0,
    low_date_bound: Optional[str] = None,
    upp_date_bound: Optional[str] = None,
):
    """
    unzip single image, convert bit-pack qa to byte value, and save as numpy
    :param source_dir: the parent folder to save image 'folder'
    :param out_dir: the folder to save result
    :param folder: the folder name of image
    :param logger: the handler of logger file
    :param dataset_info: dataset information dataclass
    :param is_partition: True, partition each image into blocks; False, save original size of image
    :param clear_threshold: threshold of clear pixel percentage, if lower than threshold, won't be processed
    :param low_date_bound: the lower bound of user interested date range
    :param upp_date_bound: the upper bound of user interested date range
    :return:
    """
    try:
        hdf_ds = gdal.Open(folder, gdal.GA_ReadOnly).GetSubDatasets()
        qa_ds = gdal.Open(hdf_ds[-1][0])
        QA_band = qa_ds.ReadAsArray()
    except (ValueError, AttributeError) as e:
        # logger.error('Cannot open QA band for {}: {}'.format(folder, e))
        logger.error("Cannot open QA band for {}: {}".format(folder, e))
        return False

    # convertQA = np.vectorize(qabitval)
    QA_band_unpacked = qabitval_array_HLS(QA_band).astype(np.short)
    if clear_threshold > 0:
        clear_ratio = np.sum(
            np.logical_or(QA_band_unpacked == QA_CLEAR, QA_band_unpacked == QA_WATER)
        ) / np.sum(QA_band_unpacked != QA_FILL)
    else:
        clear_ratio = 1

    if clear_ratio > clear_threshold:
        [collection, sensor, tile_id, imagetime, version1, version2, ext] = folder.split("/")[
            -1
        ].rsplit(".")
        year = imagetime[0:4]
        doy = imagetime[4:7]
        file_name = sensor + tile_id + year + doy + collection + version1
        if low_date_bound is not None:
            if (dt.datetime(int(year), 1, 1) + dt.timedelta(int(doy) - 1)) < parse(low_date_bound):
                return True

        if upp_date_bound is not None:
            if (dt.datetime(int(year), 1, 1) + dt.timedelta(int(doy) - 1)) > parse(upp_date_bound):
                return True

        if sensor == "L30":
            try:
                B1 = gdal.Open(hdf_ds[1][0]).ReadAsArray()
                B2 = gdal.Open(hdf_ds[2][0]).ReadAsArray()
                B3 = gdal.Open(hdf_ds[3][0]).ReadAsArray()
                B4 = gdal.Open(hdf_ds[4][0]).ReadAsArray()
                B5 = gdal.Open(hdf_ds[5][0]).ReadAsArray()
                B6 = gdal.Open(hdf_ds[6][0]).ReadAsArray()
                B7 = np.full(B6.shape, 0)  # assign zero

            except Exception as e:
                # logger.error('Cannot open spectral bands for {}: {}'.format(folder, e))
                logger.error("Cannot open Landsat bands for {}: {}".format(folder, e))
                return False
        elif sensor == "S30":
            try:
                B1 = gdal.Open(hdf_ds[1][0]).ReadAsArray()
                B2 = gdal.Open(hdf_ds[2][0]).ReadAsArray()
                B3 = gdal.Open(hdf_ds[3][0]).ReadAsArray()
                B4 = gdal.Open(hdf_ds[8][0]).ReadAsArray()
                B5 = gdal.Open(hdf_ds[11][0]).ReadAsArray()
                B6 = gdal.Open(hdf_ds[12][0]).ReadAsArray()
                B7 = np.full(B6.shape, 0)

            except Exception as e:
                # logger.error('Cannot open spectral bands for {}: {}'.format(folder, e))
                logger.error("Cannot open Landsat bands for {}: {}".format(folder, e))
                return False

        if (
            (B1 is None)
            or (B2 is None)
            or (B3 is None)
            or (B4 is None)
            or (B5 is None)
            or (B6 is None)
        ):
            logger.error("Reading Landsat band fails for {}".format(folder))
            return False

        if is_partition is True:
            # width of a block
            bytesize = 2  # short16 = 2 * byte
            # source: https://towardsdatascience.com/efficiently-splitting-an-image-into-tiles-in-python-using-numpy-d1bf0dd7b6f7
            B1_blocks = np.lib.stride_tricks.as_strided(
                B1,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B2_blocks = np.lib.stride_tricks.as_strided(
                B2,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B3_blocks = np.lib.stride_tricks.as_strided(
                B3,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B4_blocks = np.lib.stride_tricks.as_strided(
                B4,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B5_blocks = np.lib.stride_tricks.as_strided(
                B5,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B6_blocks = np.lib.stride_tricks.as_strided(
                B6,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B7_blocks = np.lib.stride_tricks.as_strided(
                B7,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            QA_blocks = np.lib.stride_tricks.as_strided(
                QA_band_unpacked,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            for i in range(dataset_info.n_block_y):
                for j in range(dataset_info.n_block_x):
                    # check if no valid pixels in the chip, then eliminate
                    qa_unique = np.unique(QA_blocks[i][j])

                    # skip blocks are all cloud, shadow or filled values
                    # in DHTC, we also don't need to save pixel that has qa value of 'QA_CLOUD',
                    # 'QA_SHADOW', or FILLED value (255)
                    if (
                        QA_CLEAR not in qa_unique
                        and QA_WATER not in qa_unique
                        and QA_SNOW not in qa_unique
                    ):
                        continue

                    block_folder = "block_x{}_y{}".format(j + 1, i + 1)
                    np.save(
                        join(join(out_dir, block_folder), file_name),
                        np.dstack(
                            [
                                B1_blocks[i][j],
                                B2_blocks[i][j],
                                B3_blocks[i][j],
                                B4_blocks[i][j],
                                B5_blocks[i][j],
                                B6_blocks[i][j],
                                B7_blocks[i][j],
                                QA_blocks[i][j],
                            ]
                        ).astype(np.int16),
                    )

        else:
            np.save(
                join(out_dir, file_name),
                np.dstack([B1, B2, B3, B4, B5, B6, B7, QA_band_unpacked]).astype(np.int16),
            )
        # scene_list.append(folder_name)
    else:
        # logger.info('Not enough clear observations for {}'.format(folder[0:len(folder) - 3]))
        logger.warn("Not enough clear observations for {}".format(folder))

    return True


def single_image_stacking(
    tmp_path: str,
    source_dir: str,
    out_dir: str,
    folder: str,
    clear_threshold: float,
    path_array: np.ndarray,
    logger: Logger,
    dataset_info: DatasetInfo,
    is_partition: bool = True,
    low_date_bound: Optional[str] = None,
    upp_date_bound: Optional[str] = None,
    b_c2: bool = False,
):
    """
    unzip single image, convert bit-pack qa to byte value, and save as numpy
    :param tmp_path: tmp folder to save unzip image
    :param source_dir: image folder save source zipped files
    :param out_dir: the folder to save result
    :param folder: the folder name of image
    :param clear_threshold: threshold of clear pixel percentage, if lower than threshold, won't be processed
    :param path_array: path array has the same dimension of inputted image, and the pixel value indicates
                      the path which the pixel belongs to; if path_array == none, we will use all path
    :param logger: the handler of logger file
    :param dataset_info: dataset information
    :param is_partition: True, partition each image into blocks; False, save original size of image
    :param low_date_bound: the lower bound of user interested year range
    :param upp_date_bound: the upper bound of user interested year range
    :param b_c2: False
    :return:
    """
    # unzip SR
    if os.path.exists(join(tmp_path, folder)):
        shutil.rmtree(join(tmp_path, folder), ignore_errors=True)
    if os.path.exists(join(tmp_path, folder.replace("SR", "BT"))):
        shutil.rmtree(join(tmp_path, folder.replace("SR", "BT")), ignore_errors=True)

    try:
        with tarfile.open(join(source_dir, folder + ".tar")) as tar_ref:
            try:
                tar_ref.extractall(join(tmp_path, folder))
            except Exception:
                # logger.warning('Unzip fails for {}'.format(folder))
                logger.error("Unzip fails for {}".format(folder))
                return
    except IOError as e:
        logger.error("Unzip fails for {}: {}".format(folder, e))
        # return

    # unzip BT
    try:
        with tarfile.open(join(source_dir, folder.replace("SR", "BT") + ".tar")) as tar_ref:
            try:
                tar_ref.extractall(join(tmp_path, folder.replace("SR", "BT")))
            except Exception:
                # logger.warning('Unzip fails for {}'.format(folder.replace("SR", "BT")))
                logger.error("Unzip fails for {}".format(folder.replace("SR", "BT")))
                return
    except IOError as e:
        logger.error("Unzip fails for {}: {}".format(folder.replace("SR", "BT"), e))
        return

    if not isdir(join(tmp_path, folder.replace("SR", "BT"))):
        logger.error("Fail to locate BT folder for {}".format(folder))
        return

    try:
        if b_c2 is False:
            QA_band = gdal_array.LoadFile(
                join(join(tmp_path, folder), "{}_PIXELQA.tif".format(folder[0 : len(folder) - 3]))
            )
        else:
            QA_band = gdal_array.LoadFile(
                join(join(tmp_path, folder), "{}_QA_PIXEL.TIF".format(folder[0 : len(folder) - 3]))
            )
    except ValueError as e:
        # logger.error('Cannot open QA band for {}: {}'.format(folder, e))
        logger.error("Cannot open QA band for {}: {}".format(folder, e))
        return

    # convertQA = np.vectorize(qabitval)
    if b_c2 is False:
        QA_band_unpacked = qabitval_array(QA_band).astype(np.short)
    else:
        QA_band_unpacked = qabitval_array_c2(QA_band).astype(np.short)
    if clear_threshold > 0:
        clear_ratio = np.sum(
            np.logical_or(QA_band_unpacked == QA_CLEAR, QA_band_unpacked == QA_WATER)
        ) / np.sum(QA_band_unpacked != QA_FILL)
    else:
        clear_ratio = 1

    if clear_ratio > clear_threshold:
        if folder[3] == "5":
            sensor = "LT5"
        elif folder[3] == "7":
            sensor = "LE7"
        elif folder[3] == "8":
            sensor = "LC8"
        elif folder[3] == "4":
            sensor = "LT4"
        elif folder[3] == "9":
            sensor = "LC9"
        else:
            logger.error("Sensor is not correctly formatted for the scene {}".format(folder))

        col = folder[8:11]
        row = folder[11:14]
        year = folder[15:19]
        doy = datetime(int(year), int(folder[19:21]), int(folder[21:23])).strftime("%j")
        collection = "C{}".format(folder[35:36])
        if b_c2 is False:
            version = folder[37:40]
        else:
            version = folder[34:35]
        file_name = sensor + col + row + year + doy + collection + version
        if low_date_bound is not None:
            if (dt.datetime(int(year), 1, 1) + dt.timedelta(int(doy) - 1)) < parse(low_date_bound):
                return
        if upp_date_bound is not None:
            if (dt.datetime(int(year), 1, 1) + dt.timedelta(int(doy) - 1)) > parse(upp_date_bound):
                return
        if b_c2 is False:
            if sensor == "LT5" or sensor == "LE7" or sensor == "LT4":
                try:
                    B1 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}B1.tif".format(folder))
                    )
                    B2 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}B2.tif".format(folder))
                    )
                    B3 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}B3.tif".format(folder))
                    )
                    B4 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}B4.tif".format(folder))
                    )
                    B5 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}B5.tif".format(folder))
                    )
                    B6 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}B7.tif".format(folder))
                    )
                    B7 = gdal_array.LoadFile(
                        join(
                            join(tmp_path, "{}_BT".format(folder[0 : len(folder) - 3])),
                            "{}_BTB6.tif".format(folder[0 : len(folder) - 3]),
                        )
                    )
                except ValueError as e:
                    # logger.error('Cannot open spectral bands for {}: {}'.format(folder, e))
                    logger.error("Cannot open Landsat bands for {}: {}".format(folder, e))
                    return
            elif sensor == "LC8" or "LC9":
                try:
                    B1 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}B2.tif".format(folder))
                    )
                    B2 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}B3.tif".format(folder))
                    )
                    B3 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}B4.tif".format(folder))
                    )
                    B4 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}B5.tif".format(folder))
                    )
                    B5 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}B6.tif".format(folder))
                    )
                    B6 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}B7.tif".format(folder))
                    )
                    B7 = gdal_array.LoadFile(
                        join(
                            join(tmp_path, "{}_BT".format(folder[0 : len(folder) - 3])),
                            "{}_BTB10.tif".format(folder[0 : len(folder) - 3]),
                        )
                    )
                except ValueError as e:
                    # logger.error('Cannot open spectral bands for {}: {}'.format(folder, e))
                    logger.error("Cannot open Landsat bands for {}: {}".format(folder, e))
                    return
        else:
            if sensor == "LT5" or sensor == "LE7" or sensor == "LT4":
                try:
                    B1 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}_B1.TIF".format(folder))
                    )
                    B2 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}_B2.TIF".format(folder))
                    )
                    B3 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}_B3.TIF".format(folder))
                    )
                    B4 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}_B4.TIF".format(folder))
                    )
                    B5 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}_B5.TIF".format(folder))
                    )
                    B6 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}_B7.TIF".format(folder))
                    )
                    B7 = gdal_array.LoadFile(
                        join(
                            join(tmp_path, "{}_BT".format(folder[0 : len(folder) - 3])),
                            "{}_BT_B6.TIF".format(folder[0 : len(folder) - 3]),
                        )
                    )
                except ValueError as e:
                    # logger.error('Cannot open spectral bands for {}: {}'.format(folder, e))
                    logger.error("Cannot open Landsat bands for {}: {}".format(folder, e))
                    return
            elif sensor == "LC8" or "LC9":
                try:
                    B1 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}_B2.TIF".format(folder))
                    )
                    B2 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}_B3.TIF".format(folder))
                    )
                    B3 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}_B4.TIF".format(folder))
                    )
                    B4 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}_B5.TIF".format(folder))
                    )
                    B5 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}_B6.TIF".format(folder))
                    )
                    B6 = gdal_array.LoadFile(
                        join(join(tmp_path, folder), "{}_B7.TIF".format(folder))
                    )
                    B7 = gdal_array.LoadFile(
                        join(
                            tmp_path,
                            "{}_BT".format(folder[0 : len(folder) - 3]),
                            "{}_BT_B10.TIF".format(folder[0 : len(folder) - 3]),
                        )
                    )
                except ValueError as e:
                    # logger.error('Cannot open spectral bands for {}: {}'.format(folder, e))
                    logger.error("Cannot open Landsat bands for {}: {}".format(folder, e))
                    return

        if (
            (B1 is None)
            or (B2 is None)
            or (B3 is None)
            or (B4 is None)
            or (B5 is None)
            or (B6 is None)
            or (B7 is None)
        ):
            logger.error("Reading Landsat band fails for {}".format(folder))
            return

        if b_c2 is True:
            B1 = (10000 * (B1 * 2.75e-05 - 0.2)).astype(np.int16)
            B2 = (10000 * (B2 * 2.75e-05 - 0.2)).astype(np.int16)
            B3 = (10000 * (B3 * 2.75e-05 - 0.2)).astype(np.int16)
            B4 = (10000 * (B4 * 2.75e-05 - 0.2)).astype(np.int16)
            B5 = (10000 * (B5 * 2.75e-05 - 0.2)).astype(np.int16)
            B7 = (10000 * (B7 * 2.75e-05 - 0.2)).astype(np.int16)
            B6 = (10 * (B6 * 0.00341802 + 149)).astype(np.int16)

        # if path_array is not None, we will eliminate those observation that has different path with its assigned path
        if path_array is not None:  # meaning that single-path processing
            if not os.path.exists(join(join(tmp_path, folder), folder.replace("_SR", ".xml"))):
                logger.error(
                    "Cannot find xml file for {}".format(
                        join(join(tmp_path, folder), folder.replace("_SR", ".xml"))
                    )
                )
                return
            tree = ET.parse(join(join(tmp_path, folder), folder.replace("_SR", ".xml")))

            # get root element
            root = tree.getroot()
            if b_c2 is False:
                elements = root.findall(
                    "./{https://landsat.usgs.gov/ard/v1}scene_metadata/{https://landsat.usgs.gov/"
                    "ard/v1}global_metadata/{https://landsat.usgs.gov/ard/v1}wrs"
                )
            else:
                elements = root.findall("SCENE_METADATA/IMAGE_ATTRIBUTES")
            if len(elements) == 0:
                logger.error("Parsing xml fails for {}".format(folder))
                return
            if b_c2 is False:
                pathid = int(elements[0].attrib["path"])
            else:
                if elements[0][3].text is not None:
                    pathid = int(elements[0][3].text)
                else:
                    raise ValueError("path id could not be located.")
            # print(pathid)

            # assign filled value to the pixels has different path id so won't be processed
            QA_band_unpacked[path_array != pathid] = QA_FILL

        if is_partition is True:
            # width of a block
            bytesize = 2  # short16 = 2 * byte
            # source: https://towardsdatascience.com/efficiently-splitting-an-image-into-tiles-in-python-using-numpy-d1bf0dd7b6f7
            B1_blocks = np.lib.stride_tricks.as_strided(
                B1,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B2_blocks = np.lib.stride_tricks.as_strided(
                B2,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B3_blocks = np.lib.stride_tricks.as_strided(
                B3,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B4_blocks = np.lib.stride_tricks.as_strided(
                B4,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B5_blocks = np.lib.stride_tricks.as_strided(
                B5,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B6_blocks = np.lib.stride_tricks.as_strided(
                B6,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B7_blocks = np.lib.stride_tricks.as_strided(
                B7,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            QA_blocks = np.lib.stride_tricks.as_strided(
                QA_band_unpacked,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            for i in range(dataset_info.n_block_y):
                for j in range(dataset_info.n_block_x):
                    # check if no valid pixels in the chip, then eliminate
                    qa_unique = np.unique(QA_blocks[i][j])

                    # skip blocks are all cloud, shadow or filled values
                    # in DHTC, we also don't need to save pixel that has qa value of 'QA_CLOUD',
                    # 'QA_SHADOW', or FILLED value (255)
                    if (
                        QA_CLEAR not in qa_unique
                        and QA_WATER not in qa_unique
                        and QA_SNOW not in qa_unique
                    ):
                        continue

                    block_folder = "block_x{}_y{}".format(j + 1, i + 1)
                    np.save(
                        join(join(out_dir, block_folder), file_name),
                        np.dstack(
                            [
                                B1_blocks[i][j],
                                B2_blocks[i][j],
                                B3_blocks[i][j],
                                B4_blocks[i][j],
                                B5_blocks[i][j],
                                B6_blocks[i][j],
                                B7_blocks[i][j],
                                QA_blocks[i][j],
                            ]
                        ).astype(np.int16),
                    )

        else:
            np.save(
                join(out_dir, file_name),
                np.dstack([B1, B2, B3, B4, B5, B6, B7, QA_band_unpacked]).astype(np.int16),
            )
        # scene_list.append(folder_name)
    else:
        # logger.info('Not enough clear observations for {}'.format(folder[0:len(folder) - 3]))
        logger.warn("Not enough clear observations for {}".format(folder[0 : len(folder) - 3]))

    # delete unzip folder
    shutil.rmtree(join(tmp_path, folder), ignore_errors=True)
    shutil.rmtree(join(tmp_path, folder.replace("SR", "BT")), ignore_errors=True)


def single_image_stacking_collection2(
    tmp_path: str,
    source_dir: str,
    out_dir: str,
    folder: str,
    clear_threshold: float,
    logger: Logger,
    dataset_info: DatasetInfo,
    bounds: list,
    is_partition: bool = True,
    low_date_bound: Optional[str] = None,
    upp_date_bound: Optional[str] = None,
):
    """
    for collection 2
    :param tmp_path: tmp folder to save unzip image
    :param source_dir: image folder save source zipped files
    :param out_dir: the folder to save result
    :param folder: the folder name of image
    :param clear_threshold: threshold of clear pixel percentage, if lower than threshold, won't be processed
    :param logger: the handler of logger file
    :param dataset_info:DatasetInfo
    :param is_partition: True, partition each image into blocks; False, save original size of image
    :param low_date_bound: the lower bound of user interested date range
    :param upp_date_bound: the upper bound of user interested date range
    :param bounds
    :return:
    """
    # unzip SR
    if os.path.exists(join(tmp_path, folder)):
        shutil.rmtree(join(tmp_path, folder), ignore_errors=True)

    try:
        with tarfile.open(join(source_dir, folder + ".tar")) as tar_ref:
            try:
                tar_ref.extractall(join(tmp_path, folder))
            except Exception:
                # logger.warning('Unzip fails for {}'.format(folder))
                logger.error("Unzip fails for {}".format(folder))
                return
    except IOError as e:
        logger.error("Unzip fails for {}: {}".format(folder, e))
        # return

    try:
        QA_band = gdal.Warp(
            os.path.join(tmp_path, "_tmp_img"),
            gdal.Open(join(join(tmp_path, folder), "{}_QA_PIXEL.TIF".format(folder))),
            outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
            xRes=res,
            yRes=res,
            dstNodata=1,
            srcNodata=1,
            outputType=gdal.GDT_UInt16,
        ).ReadAsArray()
    except ValueError as e:
        # logger.error('Cannot open QA band for {}: {}'.format(folder, e))
        logger.error("Cannot open QA band for {}: {}".format(folder, e))
        return

    # convertQA = np.vectorize(qabitval)
    QA_band_unpacked = qabitval_array_c2(QA_band).astype(np.short)
    if clear_threshold > 0:
        clear_ratio = np.sum(
            np.logical_or(QA_band_unpacked == QA_CLEAR, QA_band_unpacked == QA_WATER)
        ) / np.sum(QA_band_unpacked != QA_FILL)
    else:
        clear_ratio = 1

    if clear_ratio > clear_threshold:
        if folder[3] == "5":
            sensor = "LT5"
        elif folder[3] == "7":
            sensor = "LE7"
        elif folder[3] == "8":
            sensor = "LC8"
        elif folder[3] == "9":
            sensor = "LC9"
        elif folder[3] == "4":
            sensor = "LT4"
        else:
            logger.error("Sensor is not correctly formatted for the scene {}".format(folder))

        path = folder[10:13]
        row = folder[13:16]
        year = folder[17:21]
        doy = datetime(int(year), int(folder[21:23]), int(folder[23:25])).strftime("%j")
        collection = "C2"
        version = folder[len(folder) - 2 : len(folder)]
        file_name = sensor + path + row + year + doy + collection + version
        if low_date_bound is not None:
            if (dt.datetime(int(year), 1, 1) + dt.timedelta(int(doy) - 1)) < parse(low_date_bound):
                return True

        if upp_date_bound is not None:
            if (dt.datetime(int(year), 1, 1) + dt.timedelta(int(doy) - 1)) > parse(upp_date_bound):
                return True

        if sensor == "LT5" or sensor == "LE7" or sensor == "LT4":
            try:
                # B1 = gdal_array.LoadFile(join(join(tmp_path, folder),
                #                               "{}_SR_B1.TIF".format(folder)))
                B1 = gdal.Warp(
                    os.path.join(tmp_path, "_tmp_img"),
                    gdal.Open(join(join(tmp_path, folder), "{}_SR_B1.TIF".format(folder))),
                    outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
                    xRes=res,
                    yRes=res,
                    dstNodata=0,
                    srcNodata=0,
                    outputType=gdal.GDT_UInt16,
                ).ReadAsArray()
                B2 = gdal.Warp(
                    os.path.join(tmp_path, "_tmp_img"),
                    gdal.Open(join(join(tmp_path, folder), "{}_SR_B2.TIF".format(folder))),
                    outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
                    xRes=res,
                    yRes=res,
                    dstNodata=0,
                    srcNodata=0,
                    outputType=gdal.GDT_UInt16,
                ).ReadAsArray()
                B3 = gdal.Warp(
                    os.path.join(tmp_path, "_tmp_img"),
                    gdal.Open(join(join(tmp_path, folder), "{}_SR_B3.TIF".format(folder))),
                    outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
                    xRes=res,
                    yRes=res,
                    dstNodata=0,
                    srcNodata=0,
                    outputType=gdal.GDT_UInt16,
                ).ReadAsArray()
                B4 = gdal.Warp(
                    os.path.join(tmp_path, "_tmp_img"),
                    gdal.Open(join(join(tmp_path, folder), "{}_SR_B4.TIF".format(folder))),
                    outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
                    xRes=res,
                    yRes=res,
                    dstNodata=0,
                    srcNodata=0,
                    outputType=gdal.GDT_UInt16,
                ).ReadAsArray()
                B5 = gdal.Warp(
                    os.path.join(tmp_path, "_tmp_img"),
                    gdal.Open(join(join(tmp_path, folder), "{}_SR_B5.TIF".format(folder))),
                    outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
                    xRes=res,
                    yRes=res,
                    dstNodata=0,
                    srcNodata=0,
                    outputType=gdal.GDT_UInt16,
                ).ReadAsArray()

                B6 = gdal.Warp(
                    os.path.join(tmp_path, "_tmp_img"),
                    gdal.Open(join(join(tmp_path, folder), "{}_SR_B7.TIF".format(folder))),
                    outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
                    xRes=res,
                    yRes=res,
                    dstNodata=0,
                    srcNodata=0,
                    outputType=gdal.GDT_UInt16,
                ).ReadAsArray()
                B7 = gdal.Warp(
                    os.path.join(tmp_path, "_tmp_img"),
                    gdal.Open(join(join(tmp_path, folder), "{}_ST_B6.TIF".format(folder))),
                    outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
                    xRes=res,
                    yRes=res,
                    dstNodata=0,
                    srcNodata=0,
                    outputType=gdal.GDT_UInt16,
                ).ReadAsArray()
            except ValueError as e:
                # logger.error('Cannot open spectral bands for {}: {}'.format(folder, e))
                logger.error("Cannot open Landsat bands for {}: {}".format(folder, e))
                return
        elif sensor == "LC8" or "LC9":
            try:
                B1 = gdal.Warp(
                    os.path.join(tmp_path, "_tmp_img"),
                    gdal.Open(join(join(tmp_path, folder), "{}_SR_B2.TIF".format(folder))),
                    outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
                    xRes=res,
                    yRes=res,
                    dstNodata=0,
                    srcNodata=0,
                    outputType=gdal.GDT_UInt16,
                ).ReadAsArray()
                B2 = gdal.Warp(
                    os.path.join(tmp_path, "_tmp_img"),
                    gdal.Open(join(join(tmp_path, folder), "{}_SR_B3.TIF".format(folder))),
                    outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
                    xRes=res,
                    yRes=res,
                    dstNodata=0,
                    srcNodata=0,
                    outputType=gdal.GDT_UInt16,
                ).ReadAsArray()
                B3 = gdal.Warp(
                    os.path.join(tmp_path, "_tmp_img"),
                    gdal.Open(join(join(tmp_path, folder), "{}_SR_B4.TIF".format(folder))),
                    outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
                    xRes=res,
                    yRes=res,
                    dstNodata=0,
                    srcNodata=0,
                    outputType=gdal.GDT_UInt16,
                ).ReadAsArray()
                B4 = gdal.Warp(
                    os.path.join(tmp_path, "_tmp_img"),
                    gdal.Open(join(join(tmp_path, folder), "{}_SR_B5.TIF".format(folder))),
                    outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
                    xRes=res,
                    yRes=res,
                    dstNodata=0,
                    srcNodata=0,
                    outputType=gdal.GDT_UInt16,
                ).ReadAsArray()
                B5 = gdal.Warp(
                    os.path.join(tmp_path, "_tmp_img"),
                    gdal.Open(join(join(tmp_path, folder), "{}_SR_B6.TIF".format(folder))),
                    outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
                    xRes=res,
                    yRes=res,
                    dstNodata=0,
                    srcNodata=0,
                    outputType=gdal.GDT_UInt16,
                ).ReadAsArray()
                B6 = gdal.Warp(
                    os.path.join(tmp_path, "_tmp_img"),
                    gdal.Open(join(join(tmp_path, folder), "{}_SR_B7.TIF".format(folder))),
                    outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
                    xRes=res,
                    yRes=res,
                    dstNodata=0,
                    srcNodata=0,
                    outputType=gdal.GDT_UInt16,
                ).ReadAsArray()
                B7 = gdal.Warp(
                    os.path.join(tmp_path, "_tmp_img"),
                    gdal.Open(join(join(tmp_path, folder), "{}_ST_B10.TIF".format(folder))),
                    outputBounds=[bounds[0], bounds[1], bounds[2], bounds[3]],
                    xRes=res,
                    yRes=res,
                    dstNodata=0,
                    srcNodata=0,
                    outputType=gdal.GDT_UInt16,
                ).ReadAsArray()
            except ValueError as e:
                # logger.error('Cannot open spectral bands for {}: {}'.format(folder, e))
                logger.error("Cannot open Landsat bands for {}: {}".format(folder, e))
                return

        if (
            (B1 is None)
            or (B2 is None)
            or (B3 is None)
            or (B4 is None)
            or (B5 is None)
            or (B6 is None)
            or (B7 is None)
        ):
            logger.error("Reading Landsat band fails for {}".format(folder))
            return

        # source: https://www.usgs.gov/faqs/how-do-i-use-scale-factor-landsat-level-2-science-products?qt-
        # news_science_products=0#qt-news_science_products recommended by yongquan
        B1 = (10000 * (B1 * 2.75e-05 - 0.2)).astype(np.int16)
        B2 = (10000 * (B2 * 2.75e-05 - 0.2)).astype(np.int16)
        B3 = (10000 * (B3 * 2.75e-05 - 0.2)).astype(np.int16)
        B4 = (10000 * (B4 * 2.75e-05 - 0.2)).astype(np.int16)
        B5 = (10000 * (B5 * 2.75e-05 - 0.2)).astype(np.int16)
        B7 = (10000 * (B7 * 2.75e-05 - 0.2)).astype(np.int16)
        B6 = (10 * (B6 * 0.00341802 + 149)).astype(np.int16)

        if is_partition is True:
            # width of a block
            bytesize = 2  # short16 = 2 * byte
            # source: https://towardsdatascience.com/efficiently-splitting-an-image-into-tiles-
            # in-python-using-numpy-d1bf0dd7b6f7
            B1_blocks = np.lib.stride_tricks.as_strided(
                B1,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B2_blocks = np.lib.stride_tricks.as_strided(
                B2,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B3_blocks = np.lib.stride_tricks.as_strided(
                B3,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B4_blocks = np.lib.stride_tricks.as_strided(
                B4,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B5_blocks = np.lib.stride_tricks.as_strided(
                B5,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B6_blocks = np.lib.stride_tricks.as_strided(
                B6,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            B7_blocks = np.lib.stride_tricks.as_strided(
                B7,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            QA_blocks = np.lib.stride_tricks.as_strided(
                QA_band_unpacked,
                shape=(
                    dataset_info.n_block_y,
                    dataset_info.n_block_x,
                    dataset_info.block_height,
                    dataset_info.block_width,
                ),
                strides=(
                    dataset_info.n_cols * dataset_info.block_height * bytesize,
                    dataset_info.block_width * bytesize,
                    dataset_info.n_cols * bytesize,
                    bytesize,
                ),
            )
            for i in range(dataset_info.n_block_y):
                for j in range(dataset_info.n_block_x):
                    # check if no valid pixels in the chip, then eliminate
                    qa_unique = np.unique(QA_blocks[i][j])

                    # skip blocks are all cloud, shadow or filled values
                    # in DHTC, we also don't need to save pixel that has qa value of 'QA_CLOUD',
                    # 'QA_SHADOW', or FILLED value (255)
                    if (
                        QA_CLEAR not in qa_unique
                        and QA_WATER not in qa_unique
                        and QA_SNOW not in qa_unique
                    ):
                        continue

                    block_folder = "block_x{}_y{}".format(j + 1, i + 1)
                    np.save(
                        join(join(out_dir, block_folder), file_name),
                        np.dstack(
                            [
                                B1_blocks[i][j],
                                B2_blocks[i][j],
                                B3_blocks[i][j],
                                B4_blocks[i][j],
                                B5_blocks[i][j],
                                B6_blocks[i][j],
                                B7_blocks[i][j],
                                QA_blocks[i][j],
                            ]
                        ).astype(np.int16),
                    )

        else:
            np.save(
                join(out_dir, file_name),
                np.dstack([B1, B2, B3, B4, B5, B6, B7, QA_band_unpacked]).astype(np.int16),
            )
        # scene_list.append(folder_name)
    else:
        # logger.info('Not enough clear observations for {}'.format(folder[0:len(folder) - 3]))
        logger.warn("Not enough clear observations for {}".format(folder[0 : len(folder) - 3]))

    # delete unzip folder
    shutil.rmtree(join(tmp_path, folder), ignore_errors=True)
    shutil.rmtree(join(tmp_path, folder.replace("SR", "BT")), ignore_errors=True)


def checkfinished_step1(tmp_path):
    """
    :param tmp_path:
    :return:
    """
    if not os.path.exists(tmp_path):
        return False
    return True


def checkfinished_step2(out_dir, n_cores):
    """
    :param out_dir:
    :param n_cores:
    :return:
    """
    for i in range(n_cores):
        if not os.path.exists(join(out_dir, "rank{}_finished.txt".format(i + 1))):
            return False
    return True


def checkfinished_step3_partition(out_dir):
    """
    :param out_dir:
    :return:
    """
    if not os.path.exists(join(out_dir, "starting_last_dates.txt")):
        return False
    else:
        return True


def checkfinished_step3_nopartition(out_dir):
    """
    :param out_dir:
    :return:
    """
    if not os.path.exists(join(out_dir, "scene_list.txt")):
        return False
    return True


def get_extent(extent_geojson, res, buf=0):
    """
    read shapefile of a tile from an S3 bucket, and convert projection to be aligned with sample image.
    arg:
        'extent_geojson': sharply geojson object
        res: planet resolution
    return:
        (float, float, float, float), (int, int)) tuple
    """
    # txmin = min([row[0] for row in extent_geojson['coordinates'][0]]) - res / 2.0
    # txmax = max([row[0] for row in extent_geojson['coordinates'][0]]) + res / 2.0
    # tymin = min([row[1] for row in extent_geojson['coordinates'][0]]) - res / 2.0
    # tymax = max([row[1] for row in extent_geojson['coordinates'][0]]) + res / 2.0
    txmin = extent_geojson["bbox"][0] - res * buf
    txmax = extent_geojson["bbox"][2] + res * buf
    tymin = extent_geojson["bbox"][1] - res * buf
    tymax = extent_geojson["bbox"][3] + res * buf
    n_row = np.ceil((tymax - tymin) / res)
    n_col = np.ceil((txmax - txmin) / res)
    txmin_new = (txmin + txmax) / 2 - n_col / 2 * res
    txmax_new = (txmin + txmax) / 2 + n_col / 2 * res
    tymin_new = (tymin + tymax) / 2 - n_row / 2 * res
    tymax_new = (tymin + tymax) / 2 + n_row / 2 * res
    return (txmin_new, txmax_new, tymin_new, tymax_new), (n_row, n_col)


def get_feature(shp, id):
    for feature in shp:
        if feature["properties"]["id"] == id:
            return feature


def explode(coords):
    """Explode a GeoJSON geometry's coordinates object and yield coordinate tuples.
    As long as the input is conforming, the type of the geometry doesn't matter."""
    for e in coords:
        if isinstance(e, (float, int)):
            yield coords
            break
        else:
            for f in explode(e):
                yield f


def bbox(f):
    x, y = zip(*list(explode(f["geometry"]["coordinates"])))
    return min(x), min(y), max(x), max(y)


@click.command()
@click.option(
    "--source_dir",
    type=str,
    default=None,
    help="the folder directory of Landsat tar files " "downloaded from USGS website",
)
@click.option("--out_dir", type=str, default=None, help="the folder directory for storing stacks")
@click.option(
    "--clear_threshold", type=float, default=0, help="user-defined clear pixel proportion"
)
@click.option(
    "--single_path", type=bool, default=True, help="indicate if using single_path or sidelap"
)
@click.option("--rank", type=int, default=1, help="the rank id")
@click.option("--n_cores", type=int, default=1, help="the total cores assigned")
@click.option("--is_partition/--continuous", default=True, help="partition the output to blocks")
@click.option("--yaml_path", type=str, help="yaml file path")
@click.option("--hpc/--dhtc", default=True, help="if it is set for HPC or DHTC environment")
@click.option(
    "--low_date_bound",
    type=str,
    default=None,
    help="the lower bound of the date range of user interest",
)
@click.option(
    "--upp_date_bound",
    type=str,
    default=None,
    help="the upper bound of the date range of user interest",
)
@click.option(
    "--collection",
    type=click.Choice(["ARD", "C2", "HLS", "HLS14", "ARD-C2"]),
    default="ARD",
    help="image source category",
)
@click.option("--shapefile_path", type=str, default=None)
@click.option(
    "--id",
    type=int,
    default=0,
    help=" the id of shapefile to crop image, only for C2 image category",
)
def main(
    source_dir,
    out_dir,
    clear_threshold,
    single_path,
    rank,
    n_cores,
    is_partition,
    yaml_path,
    hpc,
    low_date_bound,
    upp_date_bound,
    collection,
    shapefile_path,
    id,
):
    if not os.path.exists(source_dir):
        print("Source directory not exists!")

    # select only _SR
    if collection == "ARD" or collection == "ARD-C2":
        folder_list = [
            f[0 : len(f) - 4]
            for f in listdir(source_dir)
            if (
                isfile(join(source_dir, f))
                and f.endswith(".tar")
                and f[len(f) - 6 : len(f) - 4] == "SR"
            )
        ]
    elif collection == "C2":
        folder_list = [
            f[0 : len(f) - 4]
            for f in listdir(source_dir)
            if (isfile(join(source_dir, f)) and f.endswith(".tar"))
        ]
    elif collection == "HLS":
        folder_list = [f for f in listdir(source_dir) if f.startswith("HLS")]
    elif collection == "HLS14":
        folder_list = [y for x in os.walk(source_dir) for y in glob(os.path.join(x[0], "*.hdf"))]

    tmp_path = join(out_dir, "tmp{}".format(rank))

    tz = timezone("US/Eastern")
    with open(yaml_path) as yaml_obj:
        config_general = yaml.safe_load(yaml_obj)

    dataset_info = class_from_dict(DatasetInfo, config_general["DATASETINFO"])

    logger = None
    # create needed folders
    if rank == 1:
        # mode = w enables the log file to be overwritten
        logging.basicConfig(
            filename=join(os.getcwd(), "prepare_ard.log"), filemode="w", level=logging.INFO
        )
        logger = logging.getLogger(__name__)

        logger.info(
            "AutoPrepareDataARD starts: {}".format(datetime.now(tz).strftime("%Y-%m-%d %H:%M:%S"))
        )

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        if is_partition is True:
            for i in range(dataset_info.n_block_y):
                for j in range(dataset_info.n_block_x):
                    block_folder = "block_x{}_y{}".format(j + 1, i + 1)
                    if not os.path.exists(join(out_dir, block_folder)):
                        os.mkdir(join(out_dir, block_folder))

        if not os.path.exists(tmp_path):
            os.mkdir(tmp_path)

        if hpc is True:
            if collection == "ARD" or collection == "ARD-C2":
                # warp a tile-based single path tif
                with importlib_resources.path(
                    "pycold.imagetool", "singlepath_landsat_conus.tif"
                ) as conus_image_fpath:
                    # conus_image_fpath = (Path(__file__).parent / 'singlepath_landsat_conus.tif').resolve()
                    conus_image = gdal.Open(os.fspath(conus_image_fpath))

                if os.path.exists(join(tmp_path, folder_list[0])):
                    shutil.rmtree(join(tmp_path, folder_list[0]), ignore_errors=True)
                with tarfile.open(join(source_dir, folder_list[0] + ".tar")) as tar_ref:
                    try:
                        tar_ref.extractall(join(tmp_path, folder_list[0]))
                    except Exception:
                        logger.error("Unzip fails for {}".format(folder_list[0]))
                if collection == "ARD":
                    ref_image = gdal.Open(
                        join(tmp_path, folder_list[0], "{}B1.tif".format(folder_list[0]))
                    )
                else:
                    ref_image = gdal.Open(
                        join(tmp_path, folder_list[0], "{}_B1.TIF".format(folder_list[0]))
                    )
                trans = ref_image.GetGeoTransform()
                proj = ref_image.GetProjection()
                xmin = trans[0]
                ymax = trans[3]
                xmax = xmin + trans[1] * ref_image.RasterXSize
                ymin = ymax + trans[5] * ref_image.RasterYSize
                params = gdal.WarpOptions(
                    dstSRS=proj,
                    outputBounds=[xmin, ymin, xmax, ymax],
                    width=ref_image.RasterXSize,
                    height=ref_image.RasterYSize,
                )

                out_fpath = join(out_dir, "singlepath_landsat_tile.tif")
                dst = gdal.Warp(out_fpath, conus_image, options=params)
                # must close the dst
                dst = None  # NOQA
                out_img = None  # NOQA
                shutil.rmtree(join(tmp_path, folder_list[0]), ignore_errors=True)

        if collection == "ARD" or collection == "ARD-C2":
            ordinal_dates = [
                pd.Timestamp.toordinal(
                    dt.datetime(int(folder[15:19]), int(folder[19:21]), int(folder[21:23]))
                )
                for folder in folder_list
            ]
        elif collection == "C2":
            ordinal_dates = [
                pd.Timestamp.toordinal(
                    dt.datetime(int(folder[17:21]), int(folder[21:23]), int(folder[23:25]))
                )
                for folder in folder_list
            ]
        elif collection == "HLS":
            ordinal_dates = [
                pd.Timestamp.toordinal(dt.datetime(int(folder[15:19]), 1, 1))
                + int(folder[19:22])
                - 1
                for folder in folder_list
            ]
        elif collection == "HLS14":
            ordinal_dates = [
                pd.Timestamp.toordinal(dt.datetime(int(folder.split("/")[-1][15:19]), 1, 1))
                + int(folder.split("/")[-1][19:22])
                - 1
                for folder in folder_list
            ]
        ordinal_dates.sort()
        file = open(join(out_dir, "starting_last_dates.txt"), "w+")  # need to save out starting and
        # lasting date for this tile
        if low_date_bound is not None:
            tmp = pd.Timestamp.toordinal(parse(low_date_bound))
        else:
            tmp = 0
        file.writelines("{}\n".format(str(np.max([ordinal_dates[0], tmp]))))
        if upp_date_bound is not None:
            tmp = pd.Timestamp.toordinal(parse(upp_date_bound))
        else:
            tmp = 999999
        file.writelines("{}\n".format(str(np.min([ordinal_dates[-1], tmp]))))
        file.close()
    else:
        # mode = w enables the log file to be overwritten
        logging.basicConfig(
            filename=join(os.getcwd(), "prepare_ard.log"), filemode="a", level=logging.INFO
        )
        logger = logging.getLogger(__name__)

    # use starting_last_dates.txt to indicate all folder have been created. Very stupid way. Need improve it in future.
    while not checkfinished_step1(join(out_dir, "starting_last_dates.txt")):
        time.sleep(5)

    # read a general path file which can indicate which pixel is assigned to which path
    path_array = None
    if single_path is True:
        if collection == "ARD" or collection == "ARD-C2":
            # parse tile h and v from folder name
            folder_name = os.path.basename(source_dir)
            if hpc is True:
                path_array = gdal_array.LoadFile(join(out_dir, "singlepath_landsat_tile.tif"))
            else:
                try:
                    tile_h = int(folder_name[folder_name.find("h") + 1 : folder_name.find("h") + 4])
                    tile_v = int(folder_name[folder_name.find("v") + 1 : folder_name.find("v") + 4])
                    path_array = gdal_array.LoadFile(
                        join(
                            Path(os.path.realpath(__file__)).parent,
                            "singlepath_landsat_tile_crop_compress.tif",
                        ),
                        xoff=tile_h * dataset_info.n_cols,
                        yoff=tile_v * dataset_info.n_rows,
                        xsize=dataset_info.n_cols,
                        ysize=dataset_info.n_rows,
                    )  # read a partial array
                    # from a large file
                except IOError as e:
                    logger.error(e)
                    exit()

    if collection == "C2":
        # gpd_tile = gpd.read_file(shapefile_path)
        bounds = bbox(get_feature(fiona.open(shapefile_path), id))

        for i in range(int(np.ceil(len(folder_list) / n_cores))):
            new_rank = rank - 1 + i * n_cores
            # means that all folder has been processed
            if new_rank > (len(folder_list) - 1):
                break
            folder = folder_list[new_rank]
            single_image_stacking_collection2(
                tmp_path,
                source_dir,
                out_dir,
                folder,
                clear_threshold,
                logger,
                dataset_info,
                bounds,
                is_partition=is_partition,
                low_date_bound=low_date_bound,
                upp_date_bound=upp_date_bound,
            )
    elif collection == "ARD":
        # assign files to each core
        for i in range(int(np.ceil(len(folder_list) / n_cores))):
            new_rank = rank - 1 + i * n_cores
            # means that all folder has been processed
            if new_rank > (len(folder_list) - 1):
                break
            folder = folder_list[new_rank]
            single_image_stacking(
                tmp_path,
                source_dir,
                out_dir,
                folder,
                clear_threshold,
                path_array,
                logger,
                dataset_info,
                is_partition=is_partition,
                low_date_bound=low_date_bound,
                upp_date_bound=upp_date_bound,
            )
    elif collection == "ARD-C2":
        # assign files to each core
        for i in range(int(np.ceil(len(folder_list) / n_cores))):
            new_rank = rank - 1 + i * n_cores
            # means that all folder has been processed
            if new_rank > (len(folder_list) - 1):
                break
            folder = folder_list[new_rank]
            single_image_stacking(
                tmp_path,
                source_dir,
                out_dir,
                folder,
                clear_threshold,
                path_array,
                logger,
                dataset_info,
                is_partition=is_partition,
                low_date_bound=low_date_bound,
                upp_date_bound=upp_date_bound,
                b_c2=True,
            )
    elif collection == "HLS":
        # assign files to each core
        for i in range(int(np.ceil(len(folder_list) / n_cores))):
            new_rank = rank - 1 + i * n_cores
            # means that all folder has been processed
            if new_rank > (len(folder_list) - 1):
                break
            folder = folder_list[new_rank]
            single_image_stacking_hls(
                source_dir,
                out_dir,
                logger,
                dataset_info,
                folder,
                clear_threshold=clear_threshold,
                is_partition=is_partition,
                low_date_bound=low_date_bound,
                upp_date_bound=upp_date_bound,
            )
    elif collection == "HLS14":
        # assign files to each core
        for i in range(int(np.ceil(len(folder_list) / n_cores))):
            new_rank = rank - 1 + i * n_cores
            # means that all folder has been processed
            if new_rank > (len(folder_list) - 1):
                break
            folder = folder_list[new_rank]
            single_image_stacking_hls14(
                out_dir,
                logger,
                dataset_info,
                folder,
                clear_threshold=clear_threshold,
                is_partition=is_partition,
                low_date_bound=low_date_bound,
                upp_date_bound=upp_date_bound,
            )
    # create an empty file for signaling the core that has been finished
    with open(os.path.join(out_dir, "rank{}_finished.txt".format(rank)), "w"):
        pass

    # wait for other cores assigned
    # while not checkfinished_step2(out_dir, n_cores):
    #     time.sleep(5)
    shutil.rmtree(tmp_path, ignore_errors=True)
    # create scene list after stacking is finished and remove folders. Not need it in DHTC
    if rank == 1:
        # remove tmp folder
        logger.info(
            "Stacking procedure finished: {}".format(datetime.now(tz).strftime("%Y-%m-%d %H:%M:%S"))
        )


if __name__ == "__main__":
    main()
