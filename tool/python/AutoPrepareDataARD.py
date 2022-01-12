# This script is an example for generating block-based stack files from original ARD zip as intermediate inputs to the
# COLD algorithm in a HPC environment. As preparation, you need to download '_BT' and '_SR' for all Landsat
# ARD collection 1.
# This script has 4 steps: 1) warp single-path array to limit the observation inputs for each pixel
# if single_path is set True; 2) unzip all images and unpack bit-based QA bands; 3) partition each 5000*5000 temporal
# images to blocks and eliminate those image blocks if no clear, water or snow pixel were in it (so to save disk space
# and  enable independent IO for individual block in later time-series analysis); 4) save each image block to
# python-native binary format (.npy) into its block folders

# For a 42-year Landsat ARD C1 tile (~3000 images), this script averagely produces ~350 G intermediate disk
# files, and takes ~12 mins to finish if 200 EPYC 7452 cores are used.

import warnings
warnings.filterwarnings("ignore")
import os
from osgeo import gdal_array
import numpy as np
import gdal
import tarfile
from os import listdir
import logging
import numpy as geek
from datetime import datetime
import datetime as dt
import click
import shutil
from pytz import timezone
import time
import xml.etree.ElementTree as ET
import yaml
import pandas as pd
from os.path import isfile, join, isdir

# define constant here
QA_CLEAR = 0
QA_WATER = 1
QA_SHADOW = 2
QA_SNOW = 3
QA_CLOUD = 4
QA_FILL = 255

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


def qabitval_array(packedint_array):
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
    QA_CLOUD_unpacked = geek.bitwise_and(packedint_array, 1 << (QA_CLOUD+1))
    QA_SHADOW_unpacked = geek.bitwise_and(packedint_array, 1 << (QA_SHADOW+1))
    QA_SNOW_unpacked = geek.bitwise_and(packedint_array, 1 << (QA_SNOW+1))
    QA_WATER_unpacked = geek.bitwise_and(packedint_array, 1 << (QA_WATER+1))
    QA_CLEAR_unpacked = geek.bitwise_and(packedint_array, 1 << (QA_CLEAR+1))

    unpacked[QA_CLEAR_unpacked > 0] = QA_CLEAR
    unpacked[QA_WATER_unpacked > 0] = QA_WATER
    unpacked[QA_SNOW_unpacked > 0] = QA_SNOW
    unpacked[QA_SHADOW_unpacked > 0] = QA_SHADOW
    unpacked[QA_CLOUD_unpacked > 0] = QA_CLOUD
    return unpacked


def load_data(file_name, gdal_driver='GTiff'):
    '''
    Converts a GDAL compatable file into a numpy array and associated geodata.
    The rray is provided so you can run with your processing - the geodata consists of the geotransform and gdal dataset object
    if you're using an ENVI binary as input, this willr equire an associated .hdr file otherwise this will fail.
    This needs modifying if you're dealing with multiple bands.

    VARIABLES
    file_name : file name and path of your file

    RETURNS
    image array
    (geotransform, inDs)
    '''
    driver_t = gdal.GetDriverByName(gdal_driver) ## http://www.gdal.org/formats_list.html
    driver_t.Register()

    inDs = gdal.Open(file_name, gdal.GA_ReadOnly)
    # print(inDs)
    if inDs is None:
        print('Couldnt open this file {}'.format(file_name))
        sys.exit("Try again!")

    # Extract some info form the inDs
    geotransform = inDs.GetGeoTransform()

    # Get the data as a numpy array
    band = inDs.GetRasterBand(1)
    cols = inDs.RasterXSize
    rows = inDs.RasterYSize
    image_array = band.ReadAsArray(0, 0, cols, rows)

    return image_array, (geotransform, inDs)


def single_image_processing(tmp_path, source_dir, out_dir, folder, clear_threshold, path_array, logger, parameters,
                            is_partition=True, low_year_bound=1, upp_year_bound=9999):
    """
    unzip single image, convert bit-pack qa to byte value, and save as numpy as
    :param tmp_path: tmp folder to save unzip image
    :param source_dir: image folder save source zipped files
    :param out_dir: the folder to save result
    :param folder: the folder name of image
    :param clear_threshold: threshold of clear pixel percentage, if lower than threshold, won't be processed
    :param path_array: path array has the same dimension of inputted image, and the pixel value indicates
                      the path which the pixel belongs to; if path_array == none, we will use all path
    :param logger: the handler of logger file
    :param parameters
    :param is_partition: True, partition each image into blocks; False, save original size of image
    :param low_year_bound: the lower bound of user interested year range
    :param upp_year_bound: the upper bound of user interested year range
    :return:
    """
    # unzip SR
    if os.path.exists(join(tmp_path, folder)):
        shutil.rmtree(join(tmp_path, folder), ignore_errors=True)
    if os.path.exists(join(tmp_path, folder.replace("SR", "BT"))):
        shutil.rmtree(join(tmp_path, folder.replace("SR", "BT")), ignore_errors=True)

    try:
        with tarfile.open(join(source_dir, folder+'.tar')) as tar_ref:
            try:
                tar_ref.extractall(join(tmp_path, folder))
            except:
                # logger.warning('Unzip fails for {}'.format(folder))
                logger.warn('Unzip fails for {}'.format(folder))
                return
    except IOError as e:
        logger.warn('Unzip fails for {}: {}'.format(folder, e))
        # return

    # unzip BT
    try:
        with tarfile.open(join(source_dir, folder.replace("SR", "BT")+'.tar')) as tar_ref:
            try:
                tar_ref.extractall(join(tmp_path, folder.replace("SR", "BT")))
            except:
                # logger.warning('Unzip fails for {}'.format(folder.replace("SR", "BT")))
                logger.warn('Unzip fails for {}'.format(folder.replace("SR", "BT")))
                return
    except IOError as e:
        logger.warn('Unzip fails for {}: {}'.format(folder.replace("SR", "BT"), e))
        return

    if not isdir(join(tmp_path, folder.replace("SR", "BT"))):
        logger.warn('Fail to locate BT folder for {}'.format(folder))
        return

    try:
        QA_band = gdal_array.LoadFile(join(join(tmp_path, folder),
                                                   "{}_PIXELQA.tif".format(folder[0:len(folder) - 3])))
    except ValueError as e:
        # logger.error('Cannot open QA band for {}: {}'.format(folder, e))
        logger.warn('Cannot open QA band for {}: {}'.format(folder, e))
        return

    # convertQA = np.vectorize(qabitval)
    QA_band_unpacked = qabitval_array(QA_band).astype(np.short)
    if clear_threshold > 0:
        clear_ratio = np.sum(np.logical_or(QA_band_unpacked == QA_CLEAR,
                                           QA_band_unpacked == QA_WATER)) \
                      / np.sum(QA_band_unpacked != QA_FILL)
    else:
        clear_ratio = 1

    if clear_ratio > clear_threshold:
        if folder[3] == '5':
            sensor = 'LT5'
        elif folder[3] == '7':
            sensor = 'LE7'
        elif folder[3] == '8':
            sensor = 'LC8'
        elif folder[3] == '4':
            sensor = 'LT4'
        else:
            logger.warn('Sensor is not correctly formatted for the scene {}'.format(folder))

        col = folder[8:11]
        row = folder[11:14]
        year = folder[15:19]
        doy = datetime(int(year), int(folder[19:21]), int(folder[21:23])).strftime('%j')
        collection = "C{}".format(folder[35:36])
        version = folder[37:40]
        file_name = sensor + col + row + year + doy + collection + version
        if low_year_bound != 0:
            if int(year) < low_year_bound:
                return
        if upp_year_bound != 9999:
            if int(year) > upp_year_bound:
                return


        if sensor == 'LT5' or sensor == 'LE7' or sensor == 'LT4':
            try:
                B1 = gdal_array.LoadFile(join(join(tmp_path, folder),
                                         "{}B1.tif".format(folder)))
                B2 = gdal_array.LoadFile(join(join(tmp_path, folder),
                                         "{}B2.tif".format(folder)))
                B3 = gdal_array.LoadFile(join(join(tmp_path, folder),
                                         "{}B3.tif".format(folder)))
                B4 = gdal_array.LoadFile(join(join(tmp_path, folder),
                                         "{}B4.tif".format(folder)))
                B5 = gdal_array.LoadFile(join(join(tmp_path, folder),
                                         "{}B5.tif".format(folder)))
                B6 = gdal_array.LoadFile(join(join(tmp_path, folder),
                                         "{}B7.tif".format(folder)))
                B7 = gdal_array.LoadFile(
                    join(join(tmp_path, "{}_BT".format(folder[0:len(folder) - 3])),
                         "{}_BTB6.tif".format(folder[0:len(folder) - 3])))
            except ValueError as e:
                # logger.error('Cannot open spectral bands for {}: {}'.format(folder, e))
                logger.warn('Cannot open Landsat bands for {}: {}'.format(folder, e))
                return
        elif sensor == 'LC8':
            try:
                B1 = gdal_array.LoadFile(join(join(tmp_path, folder),
                                              "{}B2.tif".format(folder)))
                B2 = gdal_array.LoadFile(join(join(tmp_path, folder),
                                              "{}B3.tif".format(folder)))
                B3 = gdal_array.LoadFile(join(join(tmp_path, folder),
                                              "{}B4.tif".format(folder)))
                B4 = gdal_array.LoadFile(join(join(tmp_path, folder),
                                              "{}B5.tif".format(folder)))
                B5 = gdal_array.LoadFile(join(join(tmp_path, folder),
                                              "{}B6.tif".format(folder)))
                B6 = gdal_array.LoadFile(join(join(tmp_path, folder),
                                              "{}B7.tif".format(folder)))
                B7 = gdal_array.LoadFile(
                    join(join(tmp_path, "{}_BT".format(folder[0:len(folder) - 3])),
                         "{}_BTB10.tif".format(folder[0:len(folder) - 3])))
            except ValueError as e:
                # logger.error('Cannot open spectral bands for {}: {}'.format(folder, e))
                logger.warn('Cannot open Landsat bands for {}: {}'.format(folder, e))
                return

        if (B1 is None) or (B2 is None) or (B3 is None) or (B4 is None) or (B5 is None) or (B6 is None) or \
                (B7 is None):
            logger.warn('Reading Landsat band fails for {}'.format(folder))
            return

        # if path_array is not None, we will eliminate those observation that has different path with its assigned path
        if path_array is not None:  # meaning that single-path processing
            if not os.path.exists(join(join(tmp_path, folder), folder.replace("_SR", ".xml"))):
                logger.warn('Cannot find xml file for {}'.format(join(join(tmp_path, folder),
                                                                      folder.replace("_SR", ".xml"))))
                return
            tree = ET.parse(join(join(tmp_path, folder), folder.replace("_SR", ".xml")))

            # get root element
            root = tree.getroot()
            elements = root.findall(
                './{https://landsat.usgs.gov/ard/v1}scene_metadata/{https://landsat.usgs.gov/'
                'ard/v1}global_metadata/{https://landsat.usgs.gov/ard/v1}wrs')
            if len(elements) == 0:
                logger.warn('Parsing xml fails for {}'.format(folder))
                return
            pathid = int(elements[0].attrib['path'])

            # assign filled value to the pixels has different path id so won't be processed
            QA_band_unpacked[path_array != pathid] = QA_FILL

        if is_partition is True:
            b_width = int(parameters['n_cols'] / parameters['n_block_x'])  # width of a block
            b_height = int(parameters['n_rows'] / parameters['n_block_y'])
            bytesize = 2  # short16 = 2 * byte
            # source: https://towardsdatascience.com/efficiently-splitting-an-image-into-tiles-in-python-using-numpy-d1bf0dd7b6f7
            B1_blocks = np.lib.stride_tricks.as_strided(B1, shape=(parameters['n_block_y'],
                                                        parameters['n_block_x'], b_height, b_width),
                                                        strides=(parameters['n_cols'] * b_height * bytesize,
                                                                 b_width * bytesize,
                                                                 parameters['n_cols'] * bytesize, bytesize))
            B2_blocks = np.lib.stride_tricks.as_strided(B2, shape=(parameters['n_block_y'],
                                                        parameters['n_block_x'], b_height, b_width),
                                                        strides=(parameters['n_cols'] * b_height * bytesize,
                                                                 b_width * bytesize,
                                                                 parameters['n_cols'] * bytesize, bytesize))
            B3_blocks = np.lib.stride_tricks.as_strided(B3, shape=(parameters['n_block_y'],
                                                        parameters['n_block_x'], b_height, b_width),
                                                        strides=(parameters['n_cols'] * b_height * bytesize,
                                                                 b_width * bytesize,
                                                                 parameters['n_cols'] * bytesize, bytesize))
            B4_blocks = np.lib.stride_tricks.as_strided(B4, shape=(parameters['n_block_y'],
                                                        parameters['n_block_x'], b_height, b_width),
                                                        strides=(parameters['n_cols'] * b_height * bytesize,
                                                                 b_width * bytesize,
                                                                 parameters['n_cols'] * bytesize, bytesize))
            B5_blocks = np.lib.stride_tricks.as_strided(B5, shape=(parameters['n_block_y'],
                                                        parameters['n_block_x'], b_height, b_width),
                                                        strides=(parameters['n_cols'] * b_height * bytesize,
                                                                 b_width * bytesize,
                                                                 parameters['n_cols'] * bytesize, bytesize))
            B6_blocks = np.lib.stride_tricks.as_strided(B6, shape=(parameters['n_block_y'],
                                                        parameters['n_block_x'], b_height, b_width),
                                                        strides=(parameters['n_cols'] * b_height * bytesize,
                                                                 b_width * bytesize,
                                                                 parameters['n_cols'] * bytesize, bytesize))
            B7_blocks = np.lib.stride_tricks.as_strided(B7, shape=(parameters['n_block_y'],
                                                        parameters['n_block_x'], b_height, b_width),
                                                        strides=(parameters['n_cols'] * b_height * bytesize,
                                                                 b_width * bytesize,
                                                                 parameters['n_cols'] * bytesize, bytesize))
            QA_blocks = np.lib.stride_tricks.as_strided(QA_band_unpacked,
                                                        shape=(parameters['n_block_y'],
                                                               parameters['n_block_x'], b_height,
                                                               b_width),
                                                        strides=(parameters['n_cols']*b_height*bytesize,
                                                                 b_width * bytesize,
                                                                 parameters['n_cols']*bytesize,
                                                                 bytesize))
            for i in range(parameters['n_block_y']):
                for j in range(parameters['n_block_x']):
                    # check if no valid pixels in the chip, then eliminate
                    qa_unique = np.unique(QA_blocks[i][j])

                    # skip blocks are all cloud, shadow or filled values
                    # in DHTC, we also don't need to save pixel that has qa value of 'QA_CLOUD',
                    # 'QA_SHADOW', or FILLED value (255)
                    if QA_CLEAR not in qa_unique and \
                            QA_WATER not in qa_unique and \
                            QA_SNOW not in qa_unique:
                        continue

                    block_folder = 'block_x{}_y{}'.format(j + 1, i + 1)
                    np.save(join(join(out_dir, block_folder), file_name),
                            np.dstack([B1_blocks[i][j], B2_blocks[i][j], B3_blocks[i][j], B4_blocks[i][j],
                                       B5_blocks[i][j], B6_blocks[i][j], B7_blocks[i][j], QA_blocks[i][j]]))

        else:
            np.save(join(join(out_dir, block_folder), file_name), np.dstack([B1, B2, B3, B4, B5,
                                                                             B6, B7, QA_band_unpacked]))
        # scene_list.append(folder_name)
    else:
        # logger.info('Not enough clear observations for {}'.format(folder[0:len(folder) - 3]))
        logger.warn('Not enough clear observations for {}'.format(folder[0:len(folder) - 3]))


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
        if not os.path.exists(join(out_dir, 'rank{}_finished.txt'.format(i+1))):
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


@click.command()
@click.option('--source_dir', type=str, default=None, help='the folder directory of Landsat tar files '
                                                           'downloaded from USGS website')
@click.option('--out_dir', type=str, default=None, help='the folder directory for storing stacks')
@click.option('--clear_threshold', type=float, default=0, help='user-defined clear pixel proportion')
@click.option('--single_path', type=bool, default=True, help='indicate if using single_path or sidelap')
@click.option('--rank', type=int, default=1, help='the rank id')
@click.option('--n_cores', type=int, default=1, help='the total cores assigned')
@click.option('--is_partition', type=bool, default=True, help='partition the output to blocks')
@click.option('--yaml_path', type=str, help='yaml file path')
@click.option('--hpc', type=bool, default=True, help='if it is set for HPC or DHTC environment')
@click.option('--low_year_bound', type=int, default=1, help='the lower bound of the year range of user interest')
@click.option('--upp_year_bound', type=int, default=9999, help='the upper bound of the year range of user interest')
def main(source_dir, out_dir, clear_threshold, single_path, rank, n_cores, is_partition, yaml_path, hpc, low_year_bound,
         upp_year_bound):
    if not os.path.exists(source_dir):
        print('Source directory not exists!')

    # select only _SR
    folder_list = [f[0:len(f) - 4] for f in listdir(source_dir) if
                   (isfile(join(source_dir, f)) and f.endswith('.tar')
                    and f[len(f) - 6:len(f) - 4] == 'SR')]
    tmp_path = join(out_dir, 'tmp')

    tz = timezone('US/Eastern')
    with open(yaml_path) as yaml_obj:
        parameters = yaml.safe_load(yaml_obj)

    logger = None
    # create needed folders
    if rank == 1:
        logging.basicConfig(filename=join(os.getcwd(), 'AutoPrepareDataARD.log'),
                            filemode='w', level=logging.INFO)  # mode = w enables the log file to be overwritten
        logger = logging.getLogger(__name__)

        logger.info('AutoPrepareDataARD starts: {}'.format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        if is_partition is True:
            for i in range(parameters['n_block_y']):
                for j in range(parameters['n_block_x']):
                    block_folder = 'block_x{}_y{}'.format(j + 1, i + 1)
                    if not os.path.exists(join(out_dir, block_folder)):
                        os.mkdir(join(out_dir, block_folder))
        if not os.path.exists(tmp_path):
            os.mkdir(tmp_path)

        if hpc is True:
            # warp a tile-based single path tif
            conus_path_image = gdal.Open(join(os.getcwd(), 'singlepath_landsat_conus.tif'))
            if os.path.exists(join(tmp_path, folder_list[0])):
                shutil.rmtree(join(tmp_path, folder_list[0]), ignore_errors=True)
            with tarfile.open(join(source_dir, folder_list[0] + '.tar')) as tar_ref:
                try:
                    tar_ref.extractall(join(tmp_path, folder_list[0]))
                except:
                    logger.warning('Unzip fails for {}'.format(folder))
            ref_image = gdal.Open(join(join(tmp_path, folder_list[0]), "{}B1.tif".format(folder_list[0])))
            trans = ref_image.GetGeoTransform()
            proj = ref_image.GetProjection()
            xmin = trans[0]
            ymax = trans[3]
            xmax = xmin + trans[1] * ref_image.RasterXSize
            ymin = ymax + trans[5] * ref_image.RasterYSize
            params = gdal.WarpOptions(dstSRS=proj, outputBounds=[xmin, ymin, xmax, ymax],
                                      width=ref_image.RasterXSize, height=ref_image.RasterYSize)
            dst = gdal.Warp(join(out_dir, 'singlepath_landsat_tile.tif'), conus_path_image,
                            options=params)
            # must close the dst
            dst = None
            out_img = None
            shutil.rmtree(join(tmp_path, folder_list[0]), ignore_errors=True)
        
        # saving starting and ending date for this dataset, note that python ordinal dates has 366 days with Matlab
        # ordinal dates, so we add 366 to facilitate verification with MATLAB results
        ordinal_Dates = [pd.Timestamp.toordinal(dt.datetime(int(folder[15:19]), int(folder[19:21]),
                                                            int(folder[21:23]))) + 366
                         for folder in folder_list]
        ordinal_Dates.sort()
        file = open(join(out_dir, "starting_last_dates.txt"), "w+")  # need to save out starting and
        # lasting date for this tile
        file.writelines("{}\n".format(str(np.max([ordinal_Dates[0],
                                                 pd.Timestamp.toordinal(dt.datetime(low_year_bound, 1, 1)) + 366]))))
        file.writelines("{}\n".format(str(np.min([ordinal_Dates[-1],
                                                 pd.Timestamp.toordinal(dt.datetime(upp_year_bound, 12, 31)) + 366]))))
        file.close()
    else:
        logging.basicConfig(filename=join(os.getcwd(), 'AutoPrepareDataARD.log'),
                            filemode='a', level=logging.INFO)  # mode = w enables the log file to be overwritten
        logger = logging.getLogger(__name__)

    # use starting_last_dates.txt to indicate all folder have been created. Very stupid way. Need improve it in future.
    while not checkfinished_step1(join(out_dir, "starting_last_dates.txt")):
        time.sleep(5)

    # read a general path file which can indicate which pixel is assigned to which path
    path_array = None
    if single_path is True:
        # parse tile h and v from folder name
        folder_name = os.path.basename(source_dir)
        tile_h = int(folder_name[folder_name.find('h') + 1: folder_name.find('h') + 4])
        tile_v = int(folder_name[folder_name.find('v') + 1: folder_name.find('v') + 4])
        if hpc is True:
            path_array = gdal_array.LoadFile(join(out_dir, 'singlepath_landsat_tile.tif'))
        else:
            try:
                path_array = gdal_array.LoadFile(join(os.getcwd(), 'singlepath_landsat_tile_crop_compress.tif'),
                                                 xoff=tile_h * parameters['n_cols'],
                                                 yoff=tile_v * parameters['n_rows'],
                                                 xsize=parameters['n_cols'],
                                                 ysize=parameters['n_rows'])  # read a partial array
                # from a large file
            except IOError as e:
                logger.error(e)
                exit()
   
    # assign files to each core
    for i in range(int(np.ceil(len(folder_list) / n_cores))):
        new_rank = rank - 1 + i * n_cores
        if new_rank > (len(folder_list) - 1):  # means that all folder has been processed
            break
        folder = folder_list[new_rank]
        single_image_processing(tmp_path, source_dir, out_dir, folder, clear_threshold, path_array, logger,
                                parameters, is_partition=is_partition, low_year_bound=low_year_bound,
                                upp_year_bound=upp_year_bound)

    # create an empty file for signaling the core that has been finished
    with open(os.path.join(out_dir, 'rank{}_finished.txt'.format(rank)), 'w') as fp:
        pass

    # wait for other cores assigned
    while not checkfinished_step2(out_dir, n_cores):
        time.sleep(5)

    # create scene list after stacking is finished and remove folders. Not need it in DHTC
    if rank == 1:
        # remove tmp folder
        shutil.rmtree(tmp_path, ignore_errors=True)
        logger.info('Stacking procedure finished: {}'.format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))


if __name__ == '__main__':
    main()
