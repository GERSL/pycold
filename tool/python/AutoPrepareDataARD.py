# import warnings
# warnings.filterwarnings("ignore")
import os
import sys
import pandas  # NOQA  # must import pandas ahead of gdal; uconn hpc specific issue
from osgeo import gdal_array
import numpy as np
import gdal
import tarfile
from os import listdir
from os.path import isfile, join, isdir
import logging
import numpy as geek
from datetime import datetime
import click
import shutil
from pytz import timezone
# from mpi4py import MPI
from fixed_thread_pool_executor import FixedThreadPoolExecutor
import multiprocessing
import time
import xml.etree.ElementTree as ET

from pycold import const


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


class Parameters(dict):
    def __init__(self, params):

        super(Parameters, self).__init__(params)

    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError('No such attribute: ' + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError('No such attribute: ' + name)


def qabitval_array(packedint_array, proc_params):
    """
    Institute a hierarchy of qa values that may be flagged in the bitpacked
    value.

    fill > cloud > shadow > snow > water > clear

    Args:
        packedint: int value to bit check
        proc_params: dictionary of processing parameters
    Returns:
        offset value to use
    """
    unpacked = np.full(packedint_array.shape, 255)
    QA_FILL_unpacked = geek.bitwise_and(packedint_array, 1 << proc_params.QA_FILL)
    QA_CLOUD_unpacked = geek.bitwise_and(packedint_array, 1 << proc_params.QA_CLOUD)
    QA_SHADOW_unpacked = geek.bitwise_and(packedint_array, 1 << proc_params.QA_SHADOW)
    QA_SNOW_unpacked = geek.bitwise_and(packedint_array, 1 << proc_params.QA_SNOW)
    QA_WATER_unpacked = geek.bitwise_and(packedint_array, 1 << proc_params.QA_WATER)
    QA_CLEAR_unpacked = geek.bitwise_and(packedint_array, 1 << proc_params.QA_CLEAR)
    # QA_CIRRUS1 = geek.bitwise_and(packedint_array, 1 << proc_params.QA_CIRRUS1)
    # QA_CIRRUS2 = geek.bitwise_and(packedint_array, 1 << proc_params.QA_CIRRUS2)
    # QA_OCCLUSION = geek.bitwise_and(packedint_array, 1 << proc_params.QA_OCCLUSION)

    # unpacked[QA_OCCLUSION > 0] = proc_params.QA_CLEAR - 1
    # unpacked[np.logical_and(QA_CIRRUS1 > 0, QA_CIRRUS2 > 0)] = proc_params.QA_CLEAR - 1
    unpacked[QA_CLEAR_unpacked > 0] = proc_params.QA_CLEAR - 1
    unpacked[QA_WATER_unpacked > 0] = proc_params.QA_WATER - 1
    unpacked[QA_SNOW_unpacked > 0] = proc_params.QA_SNOW - 1
    unpacked[QA_SHADOW_unpacked > 0] = proc_params.QA_SHADOW - 1
    unpacked[QA_CLOUD_unpacked > 0] = proc_params.QA_CLOUD - 1
    unpacked[QA_FILL_unpacked > 0] = 255
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
    driver_t = gdal.GetDriverByName(gdal_driver)  # http://www.gdal.org/formats_list.html
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


def single_image_processing(tmp_path, source_dir, out_dir, folder, clear_threshold, width, height, band_count,
                            image_count, total_image_count, single_path, logger):

    # unzip SR
    if os.path.exists(join(tmp_path, folder)):
        shutil.rmtree(join(tmp_path, folder), ignore_errors=True)
    if os.path.exists(join(tmp_path, folder.replace("SR", "BT"))):
        shutil.rmtree(join(tmp_path, folder.replace("SR", "BT")), ignore_errors=True)

    try:
        with tarfile.open(join(source_dir, folder + '.tar')) as tar_ref:
            try:
                tar_ref.extractall(join(tmp_path, folder))
            except Exception:
                # logger.warning('Unzip fails for {}'.format(folder))
                logger.warn('Unzip fails for {}'.format(folder))
                return
    except IOError as e:
        logger.warn('Unzip fails for {}: {}'.format(folder, e))
        return

    # unzip BT
    try:
        with tarfile.open(join(source_dir, folder.replace("SR", "BT") + '.tar')) as tar_ref:
            try:
                tar_ref.extractall(join(tmp_path, folder.replace("SR", "BT")))
            except Exception:
                # logger.warning('Unzip fails for {}'.format(folder.replace("SR", "BT")))
                logger.warn('Unzip fails for {}'.format(folder.replace("SR", "BT")))
                return
    except IOError as e:
        logger.warn('Unzip fails for {}: {}'.format(folder.replace("SR", "BT"), e))
        return

    driver = gdal.GetDriverByName('ENVI')
    if isdir(join(tmp_path, folder.replace("SR", "BT"))):
        try:
            QA_band = gdal_array.LoadFile(join(join(tmp_path, folder),
                                                       "{}_PIXELQA.tif".format(folder[0:len(folder) - 3])))
        except ValueError as e:
            # logger.error('Cannot open QA band for {}: {}'.format(folder, e))
            logger.warn('Cannot open QA band for {}: {}'.format(folder, e))
            return

        # convertQA = np.vectorize(qabitval)
        proc_params = Parameters(const.defaults)
        QA_band_unpacked = qabitval_array(QA_band, proc_params)
        if clear_threshold > 0:
            clear_ratio = np.sum(np.logical_or(QA_band_unpacked == proc_params.QA_CLEAR - 1,
                                               QA_band_unpacked == proc_params.QA_WATER - 1)) / np.sum(QA_band_unpacked != 255)
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
                logger.warn('Sensor is not correctly formated for the scene {}'.format(folder))

            col = folder[8:11]
            row = folder[11:14]
            year = folder[15:19]
            doy = datetime(int(year), int(folder[19:21]), int(folder[21:23])).strftime('%j')
            collection = "C{}".format(folder[35:36])
            version = folder[37:40]
            folder_name = sensor + col + row + year + doy + collection + version

            file_name = folder_name + '_MTLstack'

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
                return

            if not os.path.exists(join(out_dir, folder_name)):
                os.makedirs(join(out_dir, folder_name))

            outDs = driver.Create(join(join(out_dir, folder_name), file_name), width, height,
                                  band_count,
                                  gdal.GDT_Int16, options=["INTERLEAVE=BIP"])
            outDs.GetRasterBand(1).WriteArray(B1)
            outDs.GetRasterBand(2).WriteArray(B2)
            outDs.GetRasterBand(3).WriteArray(B3)
            outDs.GetRasterBand(4).WriteArray(B4)
            outDs.GetRasterBand(5).WriteArray(B5)
            outDs.GetRasterBand(6).WriteArray(B6)
            outDs.GetRasterBand(7).WriteArray(B7)
            # the last step assign filled value for sidelap region if single path is true
            if single_path is True:
                singlepath_tile = gdal_array.LoadFile(join(out_dir, 'singlepath_landsat_tile.tif'))
                if not os.path.exists(join(join(tmp_path, folder), folder.replace("_SR", ".xml"))):
                    logger.warn('Cannot find xml file for {}'.format(join(join(tmp_path, folder),
                                                                            folder.replace("_SR", ".xml"))))
                    return
                tree = ET.parse(join(join(tmp_path, folder), folder.replace("_SR", ".xml")))
                # get root element
                root = tree.getroot()
                elements = root.findall('./{https://landsat.usgs.gov/ard/v1}scene_metadata/{https://landsat.usgs.gov/'
                                        'ard/v1}global_metadata/{https://landsat.usgs.gov/ard/v1}wrs')
                if len(elements) == 0:
                    logger.warn('Parsing xml fails for {}'.format(folder))
                    return
                pathid = int(elements[0].attrib['path'])
                # assign the region has different pathid to filled value so won't be processed
                QA_band_unpacked[singlepath_tile != pathid] = 255

            outDs.GetRasterBand(8).WriteArray(QA_band_unpacked)
            # print(join(join(tmp_path, folder), "{}B1.tif".format(folder)))
            srsdata, srsgeodata = load_data(join(join(tmp_path, folder), "{}B1.tif".format(folder)))
            original_geotransform, inDs = srsgeodata

            # srsgeodata[0][1] is resolution
            outDs.SetGeoTransform([original_geotransform[0], srsgeodata[0][1], 0.0, original_geotransform[3],
                                   0.0, -srsgeodata[0][1]])
            outDs.SetProjection(inDs.GetProjection())

            outDs.FlushCache()
            outDs = None
            # scene_list.append(folder_name)
        else:
            # logger.info('Not enough clear observations for {}'.format(folder[0:len(folder) - 3]))
            logger.warn('Not enough clear observations for {}'.format(folder[0:len(folder) - 3]))
    else:
        # logger.warning('Fail to locate BT folder for {}'.format(folder))
        logger.warn('Fail to locate BT folder for {}'.format(folder))
        return

    # delete unzip folder
    shutil.rmtree(join(tmp_path, folder), ignore_errors=True)
    shutil.rmtree(join(tmp_path, folder.replace("SR", "BT")), ignore_errors=True)

    # logger.info("Finished processing {} th scene in total {} scene ".format(image_count, total_image_count))
    print("Finished processing {} th scene in total {} scene ".format(image_count, total_image_count))


def checkfinished_step1(out_dir):
    """
    :param out_dir:
    :return:
    """
    if not os.path.exists(join(out_dir, 'singlepath_landsat_tile.tif')):
        return False
    return True


def checkfinished_step2(out_dir, n_cores):
    """
    :param out_dir:
    :param n_cores:
    :return:
    """
    for i in range(n_cores):
        if not os.path.exists(join(out_dir, 'rank{}_finished.txt'.format(i + 1))):
            return False
    return True


def checkfinished_step3(out_dir):
    """
    :param out_dir:
    :return:
    """
    if not os.path.exists(join(out_dir, "scene_list.txt")):
        return False
    return True


@click.command()
@click.option('--source_dir', type=str, default=None, help='the folder directory of Landsat tar files downloaded from USGS website')
@click.option('--out_dir', type=str, default=None, help='the folder directory for ENVI outputs')
@click.option('--threads_number', type=int, default=0, help='user-defined thread number')
@click.option('--parallel_mode', type=str, default='desktop', help='desktop or HPC')
@click.option('--clear_threshold', type=float, default=0, help='user-defined clear pixel proportion')
@click.option('--single_path', type=bool, default=True, help='indicate if using single_path or sidelap')
@click.option('--rank', type=int, default=0, help='the rank id')
@click.option('--n_cores', type=int, default=0, help='the total cores assigned')
@click.option('--n_cores_step0', type=int, default=0, help='the cores only used for this step')
def main(source_dir, out_dir, threads_number, parallel_mode, clear_threshold, single_path, rank, n_cores, n_cores_step0):
    # source_dir = '/Users/coloury/Dropbox/transfer_landsat'
    # out_dir = '/Users/coloury/sccd_test'
    # clear_threshold = 0.2
    # single_path = True
    if not os.path.exists(source_dir):
        print('Source directory not exists!')

    if parallel_mode == 'desktop':
        tz = timezone('US/Eastern')
        logging.basicConfig(filename=join(os.getcwd(), 'AutoPrepareDataARD_{}.log'.format(datetime.now(tz).
                                                                                          strftime('%Y-%m-%d %H:%M:%S'))),
                            filemode='w+', level=logging.INFO)

        logger = logging.getLogger(__name__)

        tmp_path = join(out_dir, 'tmp')

        if os.path.exists(tmp_path) is False:
            os.mkdir(tmp_path)

        if threads_number == 0:
            threads_number = multiprocessing.cpu_count()
        else:
            threads_number = int(threads_number)

        print('The thread number to be paralleled is {}'.format(threads_number))

        folder_list = [f[0:len(f) - 4] for f in listdir(source_dir) if
                       (isfile(join(source_dir, f)) and f.endswith('.tar')
                        and f[len(f) - 6:len(f) - 4] == 'SR')]
        if single_path is True:
            conus_pathfile_path = join(os.getcwd(), 'singlepath_landsat_conus.tif')
            if os.path.exists(join(tmp_path, folder_list[0])) is False:  # open first folder as reference
                with tarfile.open(join(source_dir, folder_list[0] + '.tar')) as tar_ref:
                    try:
                        tar_ref.extractall(join(tmp_path, folder_list[0]))
                    except Exception:
                        # logger.warning('Unzip fails for {}'.format(folder))
                        print('Unzip fails for {} and gdal-warp single-path fail'.format(folder_list[0]))
            ref_image = gdal.Open(
                join(join(tmp_path, folder_list[0]), "{}B1.tif".format(folder_list[0])))
            trans = ref_image.GetGeoTransform()
            proj = ref_image.GetProjection()
            conus_path_image = gdal.Open(conus_pathfile_path)
            xmin = trans[0]
            ymax = trans[3]
            xmax = xmin + trans[1] * ref_image.RasterXSize
            ymin = ymax + trans[5] * ref_image.RasterYSize
            out_img = gdal.Warp(join(out_dir, 'singlepath_landsat_tile.tif'), conus_path_image,
                                outputBounds=[xmin, ymin, xmax, ymax],
                                width=ref_image.RasterXSize, height=ref_image.RasterYSize, dstNodata=0,
                                outputType=gdal.GDT_Byte, dstSRS=proj)
            shutil.rmtree(join(tmp_path, folder_list[0]), ignore_errors=True)
            print('gdal-warp for single-path array succeed: {}'.format(datetime.now(tz)
                                                                       .strftime('%Y-%m-%d %H:%M:%S')))

        width = 5000
        height = 5000
        band_count = 8

        prepare_executor = FixedThreadPoolExecutor(size=threads_number)

        for count, folder in enumerate(folder_list):
            print("it is processing {} th scene in total {} scene ".format(count + 1, len(folder_list)))
            prepare_executor.submit(single_image_processing, tmp_path, source_dir, out_dir, folder, clear_threshold, width,
                                    height, band_count, count + 1, len(folder_list), single_path)

        # await all tile finished
        prepare_executor.drain()

        # await thread pool to stop
        prepare_executor.close()

        logger.info("Final report: finished preparation task ({})"
                    .format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

        # count_valid = len(scene_list)
        # logger.warning("Total processing scene number is {}; valid scene number is {}".format(count, count_valid))

        # remove tmp folder
        shutil.rmtree(tmp_path, ignore_errors=True)

    else:  # for HPC mode
        width = 5000
        height = 5000
        band_count = 8
        tz = timezone('US/Eastern')

        # select only _SR
        folder_list = [f[0:len(f) - 4] for f in listdir(source_dir) if
                       (isfile(join(source_dir, f)) and f.endswith('.tar')
                        and f[len(f) - 6:len(f) - 4] == 'SR')]
        tmp_path = join(out_dir, 'tmp')

        # print('AutoPrepareDataARD starts for {}: {}'.format(source_dir, datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))
        # step 1: create folders and single-tile path
        if rank == 1:
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            # if tmp path exists, delete path
            if not os.path.exists(tmp_path):
                os.mkdir(tmp_path)

            logging.basicConfig(filename=join(os.getcwd(), 'LOG_AutoPrepareDataARD.log'),
                                filemode='w', level=logging.INFO)   # mode = w enables the log file to be overwritten
            logger = logging.getLogger(__name__)
            logger.info('AutoPrepareDataARD starts: {}'.format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

            # if single path, gdalwarp an image from conus single path image
            conus_pathfile_path = join(os.getcwd(), 'singlepath_landsat_conus.tif')
            conus_path_image = gdal.Open(conus_pathfile_path)
            if os.path.exists(join(tmp_path, folder_list[0])):
                shutil.rmtree(join(tmp_path, folder_list[0]), ignore_errors=True)
            with tarfile.open(join(source_dir, folder_list[0] + '.tar')) as tar_ref:
                try:
                    tar_ref.extractall(join(tmp_path, folder_list[0]))
                except Exception:
                    logger.warning('Unzip fails for {}'.format(folder))
                    # print('Unzip fails for {} and gdal-warp single-path fail'.format(folder_list[0]))
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
            dst = None  # NOQA
            out_img = None  # NOQA
            shutil.rmtree(join(tmp_path, folder_list[0]), ignore_errors=True)
            logger.info('gdal-warp for single-path array succeed: {}'.format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))
            logger = None

        logging.basicConfig(filename=join(os.getcwd(), 'LOG_AutoPrepareDataARD.log'),
                            filemode='a', level=logging.INFO)
        logger = logging.getLogger(__name__)

        while not checkfinished_step1(out_dir):
            time.sleep(5)

        n_map_percore = int(np.ceil(len(folder_list) / n_cores_step0))

        # only the core in the early assignment (n_cores_step0) will be used
        if rank <= n_cores_step0:
            for i in range(n_map_percore):
                new_rank = rank - 1 + i * n_cores_step0
                if new_rank > (len(folder_list) - 1):
                    break
                folder = folder_list[new_rank]
                single_image_processing(tmp_path, source_dir, out_dir, folder, clear_threshold, width, height,
                                        band_count, new_rank + 1, len(folder_list), single_path, logger)

            # create an empty file for signaling
            with open(os.path.join(out_dir, 'rank{}_finished.txt'.format(rank)), 'w') as fp:  # NOQA
                pass

        # wait for all cores assigned
        while not checkfinished_step2(out_dir, n_cores_step0):
            time.sleep(5)

        # create scene list after stacking is finished and remove folders
        if rank == 1:
            # remove tmp folder
            shutil.rmtree(tmp_path, ignore_errors=True)

            # delete signal files
            # signal_filenames = [file for file in os.listdir(out_dir) if file.startswith('rank')]
            # for file in signal_filenames:
            #     os.remove(os.path.join(out_dir, file))

            # out_dir = '/scratch/suy20004/h016v008_stack'
            scene_list = [f for f in os.listdir(out_dir) if (os.path.isdir(join(out_dir, f))) and (f.startswith('L'))]
            scene_file = open(join(out_dir, "scene_list.txt"), "w+")
            for L in scene_list:
                scene_file.writelines("{}\n".format(L))
            scene_file.close()
            logger.info('Stacking procedure finished: {}'.format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

        while not checkfinished_step3(out_dir):
            time.sleep(5)

        time.sleep(3)  # sleep 3 seconds to guarantee scene_list writting finished

        # os.remove(join(out_dir, 'singlepath_landsat_tile.tif'))

        # logger.info("Final report: finished preparation task ({})"
        #             .format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

    # count_valid = len(scene_list)
    # logger.warning("Total processing scene number is {}; valid scene number is {}".format(count, count_valid))

if __name__ == '__main__':
    main()
