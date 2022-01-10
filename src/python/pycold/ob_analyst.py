import warnings
warnings.filterwarnings("ignore") # mainly to filter out lzma warnings
import numpy as np
import pandas as pd
import os
import gdal
import datetime as dt
# import datetime
from osgeo import gdal_array
from cv2 import floodFill
import cv2 as cv2
from skimage.feature import peak_local_max
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
import time
from scipy.io import loadmat
from sklearn.ensemble import RandomForestClassifier
import joblib
from scipy.stats import chi2
from os.path import join, exists
from utils import gdal_save_file_1band
from pycold.app import defaults, logging
from pycold.utils import get_block_y, get_block_x, get_col_index, get_row_index

logger = logging.getLogger(__name__)


def filename_fromordinaldate(ordinal_date):
    """
    get file name from ordinate date
    Args:
        ordinal_date: the inputted ordinate date
    Returns:
        CM file name
    """
    return 'CM_maps_{}_{}{}'.format(str(ordinal_date),
                             pd.Timestamp.fromordinal(ordinal_date - 366).year,
                             str(pd.Timestamp.fromordinal(ordinal_date - 366).timetuple().tm_yday).zfill(3))


def directionfilename_fromordinaldate(ordinal_date):
    """
    get direction file name from ordinate date
    Args:
        ordinal_date: the inputted ordinate date
    Returns:
        CM direction file name
    """
    return 'CM_direction_maps_{}_{}{}'.format(str(ordinal_date), pd.Timestamp.fromordinal(ordinal_date - 366).year,
                                              str(pd.Timestamp.fromordinal(ordinal_date
                                                                           - 366).timetuple().tm_yday).zfill(3))


def datefilename_fromordinaldate(ordinal_date):
    """
    get date file name from ordinate date
    Args:
        ordinal_date: the inputted ordinate date
    Returns:
        CM date file name
    """
    return 'CM_date_maps_{}_{}{}'.format(str(ordinal_date), pd.Timestamp.fromordinal(ordinal_date - 366).year,
                                              str(pd.Timestamp.fromordinal(ordinal_date
                                                                           - 366).timetuple().tm_yday).zfill(3))


def is_change_object(config, b_thematic, stats_lut_row, uniform_threshold, uniform_sizeslope, keyword):
    """
    config
    ----------
    config:
        obcold parameter structure
    b_thematic: True -> that has thematic inputs; False-> default config will be applied
    stats_lut_row: a table of ['id', 'change magnitude', 'mode of lc category', 'pixel number']
    uniform_threshold: grid searching usage, overwriting default config
    uniform_sizeslope: grid searching usage, overwriting default config
    keyword: string, keyword to calculate change magnitude, 'cm_average' or 'cm_median'

    Returns
    -------

    """
    # LCMAP category id
    # developed_id: 1
    # cropland_id: 2
    # grass_id: 3
    # tree_id: 4
    # water_id: 5
    # wetland_id: 6
    # ice: 7
    # barren: 8

    log10_size = np.log10(int(stats_lut_row['npixels']))
    intercept = 1
    if b_thematic == False:
        # intercept = 1 - config['default_sizeslope']
        if log10_size * config['default_sizeslope'] + intercept < 2:
            scale = log10_size * config['default_sizeslope'] + intercept
        else:
            scale = 2
        if np.double(stats_lut_row[keyword]) * scale > config['default_threshold']:
            return True
        else:
            return False
    else:
        if uniform_threshold != config['NAN_VAL'] and uniform_sizeslope != config['NAN_VAL']:
            # intercept = 1 - uniform_sizeslope
            if log10_size * uniform_sizeslope + intercept < 2:
                scale = log10_size * uniform_sizeslope + intercept
            else:
                scale = 2
            if np.double(stats_lut_row[keyword]) * scale > uniform_threshold:
                return True
            else:
                return False
        else:
            if int(stats_lut_row['mode']) == 255:
                # intercept = 1 - config['default_sizeslope']
                if log10_size * config['default_sizeslope'] + intercept < 2:
                    scale = log10_size * config['default_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > config['default_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 1:  # develop in LCMAP
                # intercept = 1 - config['C1_sizeslope']
                if log10_size * config['C1_sizeslope'] + intercept < 2:
                    scale = log10_size * config['C1_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > config['C1_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 2:  # cropland in LCMAP
                # intercept = 1 - config['C2_sizeslope']
                if log10_size * config['C2_sizeslope'] + intercept < 2:
                    scale = log10_size * config['C2_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > config['C2_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 3:  # grassland in LCMAP
                # intercept = 1 - config['C3_sizeslope']
                if log10_size * config['C3_sizeslope'] + intercept < 2:
                    scale = log10_size * config['C3_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > config['C3_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 4:  # forest in LCMAP
                # intercept = 1 - config['C4_sizeslope']
                if log10_size * config['C4_sizeslope'] + intercept < 2:
                    scale = log10_size * config['C4_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > config['C4_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 5:  # water in LCMAP
                # intercept = 1 - config['C5_sizeslope']
                if log10_size * config['C5_sizeslope'] + intercept < 2:
                    scale = log10_size * config['C5_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > config['C5_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 6:  # wetland in LCMAP
                # intercept = 1 - config['C6_sizeslope']
                if log10_size * config['C6_sizeslope'] + intercept < 2:
                    scale = log10_size * config['C6_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > config['C6_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 7:  # ice in LCMAP
                # intercept = 1 - config['C7_sizeslope']
                if log10_size * config['C7_sizeslope'] + intercept < 2:
                    scale = log10_size * config['C7_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > config['C7_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 8:  # barren in LCMAP
                # intercept = 1 - config['C8_sizeslope']
                if log10_size * config['C8_sizeslope'] + intercept < 2:
                    scale = log10_size * config['C8_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > config['C8_threshold']:
                    return True
                else:
                    return False


def modeby(input_array, index_array):
    """
    calculate modes of input_array groupped by index_array.
    Args:
        input_array: input array to group
        index_array: index array of input array

    Returns:
        a list of modes for splitted subarray following ascending order of unique id
        modified from: https://stackoverflow.com/questions/49372918/group-numpy-into-multiple-sub-arrays-using-an-array-of-values
    """
    # Get argsort indices, to be used to sort a and b in the next steps
    # input_array = classification_map
    # index_array = object_map
    sidx = index_array.argsort(kind='mergesort')
    a_sorted = input_array[sidx]
    b_sorted = index_array[sidx]

    # Get the group limit indices (start, stop of groups)
    cut_idx = np.flatnonzero(np.r_[True, b_sorted[1:] != b_sorted[:-1], True])

    # Split input array with those start, stop ones
    split = [a_sorted[i:j] for i, j in zip(cut_idx[:-1], cut_idx[1:])]
    mode_list = [stats.mode(x)[0][0] for x in split]
    return mode_list


def mode_median_by(input_array_mode, input_array_median, index_array):
    """
    Args:
        input_array_mode: input array for calculating mode
        input_array_median: input array for calculating mode
        index_array: index array of input array

    Returns:
        a list of modes for splitted subarray following ascending order of unique id
    """
    sidx = index_array.argsort(kind='mergesort')
    a1_sorted = input_array_mode[sidx]
    a2_sorted = input_array_median[sidx]
    b_sorted = index_array[sidx]

    # Get the group limit indices (start, stop of groups)
    cut_idx = np.flatnonzero(np.r_[True, b_sorted[1:] != b_sorted[:-1], True])

    # Split input array with those start, stop ones
    split_mode = [a1_sorted[i:j] for i, j in zip(cut_idx[:-1], cut_idx[1:])]
    split_median = [a2_sorted[i:j] for i, j in zip(cut_idx[:-1], cut_idx[1:])]
    mode_list = [stats.mode(x)[0][0] for x in split_mode]
    median_list = [np.median(x[~np.isnan(x)]) for x in split_median]
    return mode_list, median_list


class ObjectAnalyst:
    def __init__(self, config, results_path, model_ready=False, thematic_src=None):
        """
        parameters
        ----------
        config
        results_path
        year_lowbound
        year_uppbound
        thematic_mode
        thematic_src
        """
        self.config = config
        self.results_path = results_path
        self.model_ready = model_ready
        self.thematic_src = thematic_src
        self.obiaresults_path = join(results_path, 'OBIAresults')

    @staticmethod
    def _check_inputs(parameters, stack_path, results_path, model_ready,thematic_src):
        if type(parameters['n_rows']) != int or parameters['n_rows'] < 0:
            raise ValueError('n_rows must be positive integer')
        if type(parameters['n_cols']) != int or parameters['n_cols'] < 0:
            raise ValueError('n_cols must be positive integer')
        if type(parameters['n_block_x']) != int or parameters['n_block_x'] < 0:
            raise ValueError('n_block_x must be positive integer')
        if type(parameters['n_block_y']) != int or parameters['n_block_y'] < 0:
            raise ValueError('n_block_y must be positive integer')
        if type(parameters['n_block_y']) != int or parameters['n_block_y'] < 0:
            raise ValueError('n_block_y must be positive integer')

        if os.path.isdir(stack_path) is False:
            raise FileExistsError('No such directory: {}'.format(stack_path))
        if os.path.isdir(results_path) is False:
            raise FileExistsError('No such directory: {}'.format(results_path))
        if model_ready:
            if thematic_src is None:
                raise ValueError('thematic_src must be assigned as value if model_ready is True')


    def hpc_preparation(self):
        """
        Returns
        -------

        """
        if not exists(self.obiaresults_path):
            os.makedirs(self.obiaresults_path)

    def get_classification(self, date):
        if pd.Timestamp.fromordinal(date - 366).year - 1 < self.year_lowbound:  # we used the the year before the date
            classified_year = self.year_lowbound
        else:
            classified_year = pd.Timestamp.fromordinal(date - 366).year - 1
        if self.thematic_mode == 'single_map' or 'rf_model':
            classification_map = np.load(join(self.out_thematic_path,
                                              'yearlyclassification_{}.npy'.format(classified_year)))
        elif self.thematic_mode == 'multi_maps':
            classification_map = gdal_array.LoadFile(
                join(self.thematic_src, 'yearlyclassification_{}.tif'.format(classified_year)))
        else:
            classification_map = np.full((self.config['n_rows'], self.config['n_cols']),
                                         self.config['NAN_VAL'],
                                         dtype=np.int32)
        return classification_map

    def segmentation(self, date, cm_array, cm_direction_array, cm_date_array, cm_array_l1, cm_array_l1_direction,
                     cm_array_l1_date, classification_map, uniform_threshold=-9999, uniform_sizeslope=-9999,
                     devel_mode=False):
        filenm = filename_fromordinaldate(date)
        peak_threshold = chi2.ppf(0.95, 5)

        # assign valid CM values in the stack into current cm array where NA values are
        cm_array[cm_array == self.config['NAN_VAL']] = cm_array_l1[cm_array == self.config['NAN_VAL']]

        # config['NAN_VAL_UINT8'], i.e., 255, is also the NA value of direction.
        cm_direction_array[cm_direction_array == self.config['NAN_VAL_UINT8']] = \
            cm_array_l1_direction[cm_direction_array == self.config['NAN_VAL_UINT8']]

        cm_date_array[cm_date_array == self.config['NAN_VAL']] = \
            cm_array_l1_date[cm_date_array == self.config['NAN_VAL']]

        if devel_mode is True:
            gdal_save_file_1band(join(self.obiaresults_path, filenm + '_cm_array.tif'), cm_array,
                                 gdal.GDT_Float32, self.trans, self.proj, self.config['n_cols'],
                                 self.config['n_rows'])
            gdal_save_file_1band(join(self.obiaresults_path, filenm + '_cm_direction_array.tif'), cm_direction_array,
                                 gdal.GDT_Byte, self.trans, self.proj, self.config['n_cols'],
                                 self.config['n_rows'])

        #######################################################################################
        #                               Scale 1: change superpixel                            #
        #######################################################################################
        cm_array[cm_array < peak_threshold * (1 - self.config['floodfill_thres'])] = np.nan
        bandwidth = 1
        
        # using gaussian kernel ( with 1 sigma value) to smooth images in preparation for floodfill
        kernel = Gaussian2DKernel(x_stddev=bandwidth, y_stddev=bandwidth)
        cm_array_gaussian_s1 = convolve(cm_array, kernel, boundary='extend', preserve_nan=True)
        cm_array_gaussian_s1[np.isnan(cm_array)] = self.config['NAN_VAL']

        # cm_array_gaussian_s1 = cm_array  # try not using gaussian
        if devel_mode is True:
            gdal_save_file_1band(join(self.obiaresults_path, filenm + '_cm_array_gaussian_s1.tif'),
                                 cm_array_gaussian_s1,
                                 gdal.GDT_Float32, self.trans, self.proj, self.config['n_cols'],
                                 self.config['n_rows'])

        seed_index = peak_local_max(cm_array_gaussian_s1, threshold_abs=peak_threshold,
                                    exclude_border=False, min_distance=0)
        # seed_index = np.flip(seed_index, axis=0)

        seed_labels = np.zeros((self.config['n_rows'], self.config['n_cols']))
        seed_labels[tuple(np.transpose(seed_index))] = 1
        if devel_mode:
            gdal_save_file_1band(join(self.obiaresults_path, filenm + '_seed.tif'), seed_labels,
                                 gdal.GDT_Byte, self.trans, self.proj, self.config['n_cols'],
                                 self.config['n_rows'])

        mask_s1 = np.zeros((self.config['n_rows'] + 2, self.config['n_cols'] + 2)).astype(np.uint8)
        mask_label_s1 = np.zeros((self.config['n_rows'] + 2, self.config['n_cols'] + 2))
        floodflags_base = 8
        floodflags_base |= cv2.FLOODFILL_MASK_ONLY
        floodflags_base |= cv2.FLOODFILL_FIXED_RANGE
        no = 0
        i = 0
        cm_stack = np.dstack([cm_array_gaussian_s1, cm_direction_array, cm_date_array]).astype(np.float32)
        # config['floodfill_thres'] = 0.5
        for i in range(len(seed_index)):
            # print(i)
            remainder = i % 255
            floodflags = floodflags_base | ((remainder + 1) << 8)
            seedcm = cm_array_gaussian_s1[tuple(seed_index[i])]
            num, im, mask_s1, rect = floodFill(cm_stack, mask_s1, tuple(reversed(seed_index[i])), 0,
                                               loDiff=[seedcm * self.config['floodfill_thres'], 0,
                                                       self.config['CM_OUTPUT_INTERVAL']],
                                               upDiff=[seedcm * self.config['floodfill_thres'], 0,
                                                       self.config['CM_OUTPUT_INTERVAL']],
                                               flags=floodflags)
            # the opencv mask only supports 8-bit, so every 255 seed needs to update values in mask_label
            if remainder == 254:
                no = int(i / 255)
                mask_label_s1[(mask_label_s1 == 0) & (mask_s1 > 0)] = mask_s1[
                                                                          (mask_label_s1 == 0) & (mask_s1 > 0)].astype(
                    int) \
                                                                      + no * 255
            # flood_fill(cm_array_log, tuple(reversed(seed_index[i])), new_value=-np.ceil(floodfill_threshold)*2-i,
            #            tolerance=floodfill_threshold, connectivity=2, in_place=True)
        if devel_mode:
            print("segmentation finished for {}".format(date))

        if len(seed_index) > 0:
            # processing the rest that hasn't be recorded into mask_label
            no = int(i / 255)
            mask_label_s1[(mask_label_s1 == 0) & (mask_s1 > 0)] = mask_s1[(mask_label_s1 == 0) & (mask_s1 > 0)].astype(
                int) + no * 255

        object_map_s1 = mask_label_s1[1:self.config['n_rows'] + 1, 1:self.config['n_cols'] + 1]

        if devel_mode:
            gdal_save_file_1band(
                os.path.join(self.obiaresults_path, filenm + '_floodfill_gaussian_{}_s1.tif'.format(bandwidth)),
                object_map_s1, gdal.GDT_Int32, self.trans, self.proj, self.config['n_rows'],
                self.config['n_cols'])

        #######################################################################################
        #                                 Scale 2: change patch                              #
        #######################################################################################

        # create a object-based change
        unq_s1, ids_s1, count_s1 = np.unique(object_map_s1, return_inverse=True, return_counts=True)
        mean_list = np.bincount(ids_s1.astype(int), weights=cm_array_gaussian_s1.reshape(ids_s1.shape)) / count_s1
        mean_list[unq_s1 == 0] = self.config['NAN_VAL']  # force mean of unchanged objects to be -9999
        lut_dict_s1 = dict(zip(unq_s1, mean_list))
        cm_array_s2 = np.vectorize(lut_dict_s1.get)(object_map_s1)

        bandwidth = 1
        # using gaussian kernel ( with 0.5 sigma value) to smooth images in preparation for floodfill
        # cm_array_s2[cm_array_s2 == config['NAN_VAL']] = np.nan
        # kernel = Gaussian2DKernel(x_stddev=bandwidth, y_stddev=bandwidth)
        # cm_array_gaussian_s2 = convolve(cm_array_s2, kernel, boundary='extend', preserve_nan=True)
        # cm_array_gaussian_s2[np.isnan(cm_array_s2)] = config['NAN_VAL']

        cm_array_gaussian_s2 = cm_array_s2  # not using gaussian
        if devel_mode:
            gdal_save_file_1band(join(self.obiaresults_path,  filenm + '_cm_array_gaussian_s2.tif'),
                                 cm_array_gaussian_s2,
                                 gdal.GDT_Float32, self.trans, self.proj, self.config['n_cols'],
                                 self.config['n_rows'])

        mask_s2 = np.zeros((self.config['n_rows'] + 2, self.config['n_cols'] + 2)).astype(np.uint8)
        mask_label_s2 = np.zeros((self.config['n_rows'] + 2, self.config['n_cols'] + 2))
        floodflags_base = 8
        floodflags_base |= cv2.FLOODFILL_MASK_ONLY
        i = 0
        cm_array_gaussian_s2 = cm_array_gaussian_s2.astype(np.float32)
        # cm_stack = np.dstack([cm_array_gaussian_s2, cm_direction_array, cm_direction_array]).astype(np.float32)
        # config['floodfill_thres'] = 0.5
        for i in range(len(seed_index)):
            # print(i)
            remainder = i % 255
            floodflags = floodflags_base | ((remainder + 1) << 8)
            num, im, mask_s2, rect = floodFill(cm_array_gaussian_s2, mask_s2, tuple(reversed(seed_index[i])), 0,
                                               loDiff=[1000],
                                               upDiff=[1000],
                                               flags=floodflags)  # assign an extremely large value
            # to connect any polygon adjacent to each other
            # the opencv mask only supports 8-bit, so every 255 seed needs to update values in mask_label
            if remainder == 254:
                no = int(i / 255)
                mask_label_s2[(mask_label_s2 == 0) & (mask_s2 > 0)] = mask_s2[
                                                                          (mask_label_s2 == 0) & (mask_s2 > 0)].astype(
                    int) + \
                                                                      no * 255
            # flood_fill(cm_array_log, tuple(reversed(seed_index[i])), new_value=-np.ceil(floodfill_threshold)*2-i,
            #            tolerance=floodfill_threshold, connectivity=2, in_place=True)
        if devel_mode:
            print("segmentation finished ")

        if len(seed_index) > 0:
            # processing the rest that hasn't be recorded into mask_label
            no = int(i / 255)
            mask_label_s2[(mask_label_s2 == 0) & (mask_s2 > 0)] = mask_s2[(mask_label_s2 == 0) & (mask_s2 > 0)].astype(
                int) + no * 255

        object_map_s2 = mask_label_s2[1:self.config['n_rows'] + 1, 1:self.config['n_cols'] + 1]
        if devel_mode:
            gdal_save_file_1band(
                join(self.obiaresults_path,  'test_floodfill_gaussian_{}_s2.tif'.format(bandwidth)),
                object_map_s2, gdal.GDT_Int32, self.trans, self.proj, self.config['n_rows'],
                self.config['n_cols'])

        return object_map_s1, object_map_s2, ids_s1, unq_s1, mean_list

    def object_analysis(self, object_map_s1, object_map_s2, classification_map, ids_s1, unq_s1, mean_list,
                        devel_mode=False):
        #######################################################################################
        #                                 object-based decision                               #
        #######################################################################################
        change_map = np.zeros((self.config['n_rows'], self.config['n_cols'])).astype(np.uint8)
        change_group = []

        if devel_mode:
            start_time = time.time()

        unq_s2, ids_s2, count_s2 = np.unique(object_map_s2, return_inverse=True, return_counts=True)
        lut_dict_s2 = dict(zip(unq_s2, count_s2))
        size_map = np.vectorize(lut_dict_s2.get)(object_map_s2)  # we use size of scale 2 to represent object szie
        size_map[object_map_s2 == 0] = 0

        mode_list = modeby(classification_map.reshape(ids_s1.shape), object_map_s1.reshape(ids_s1.shape))
        size_list = modeby(size_map.reshape(ids_s1.shape),
                           object_map_s1.reshape(ids_s1.shape))  # the mode of s1 objects in sizemap
        # mode_list, median_list = mode_median_by(classification_map.reshape(ids_s1.shape),
        #                                         cm_array_gaussian_s1.reshape(ids_s1.shape),
        #                                         object_map_s1.reshape(ids_s1.shape))

        # Use the tags/IDs to perform ID based summation of array2 elems and
        # thus divide by the ID counts to get ID based average values
        stats_lut = pd.DataFrame(np.column_stack((unq_s1, mean_list,
                                                  mode_list, size_list)),
                                 columns=['id', 'cm_average', 'mode', 'npixels'])

        class_labels = unq_s1[unq_s1 > 0]
        if len(class_labels) > 0:  # if has valid change objects to output
            for index in class_labels:
                stats_lut_row = stats_lut.loc[stats_lut['id'] == index]
                if self.thematic_src is 'no_input':
                    b_thematic = False
                else:
                    b_thematic = True
                # print(count)
                if is_change_object(self.config, b_thematic, stats_lut_row, uniform_threshold, uniform_sizeslope,
                                    'cm_average'):
                    change_group.append(index)

            if devel_mode:
                print("--- %s seconds ---" % (time.time() - start_time))
            change_map[np.isin(object_map_s1, change_group)] = 1

        # if devel_mode:
        #     gdal_save_file_1band(join(self.obiaresults_path, filenm + '_OBIAresult.tif'), change_map,
        #                          gdal.GDT_Byte, self.trans, self.proj, self.config['n_cols'],
        #                          self.config['n_rows'])
        # else:
        #     np.save(join(self.obiaresults_path, filenm + '_OBIAresult.tif'), change_map)

        if devel_mode:
            print("change object analysis finished for {}".format(date))

        return change_map

    def hpc_save_obia(self, date, change_map, outputformat):
        filenm = filename_fromordinaldate(date)
        if outputformat == 'tiff':
            gdal_save_file_1band(join(self.obiaresults_path, filenm + '_OBIAresult.tif'), change_map,
                                 gdal.GDT_Byte, self.trans, self.proj, defaults['n_cols'], defaults['n_rows'])
        else:
            np.save(join(self.obiaresults_path, filenm + '_OBIAresult.npy'), change_map)
        
    def is_finished_object_analysis(self):
        """
        :return:
        """
        for n in range(self.config['n_cm_maps']):
            if not exists(join(self.obiaresults_path, filename_fromordinaldate(self.config['starting_date'] +
                                                                        n * self.config['CM_OUTPUT_INTERVAL']) +
                                               '_OBIAresult.tif')):
                return False
        return True




