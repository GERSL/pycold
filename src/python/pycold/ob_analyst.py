import warnings
warnings.filterwarnings("ignore") # mainly to filter out lzma warnings
import numpy as np
import pandas as pd
import os
import gdal
# import datetime
from cv2 import floodFill
import cv2 as cv2
from skimage.feature import peak_local_max
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from scipy import stats
from scipy.stats import chi2
from os.path import join, exists
from pycold.utils import gdal_save_file_1band, get_block_x, get_block_y, read_blockdata, get_rowcol_intile
from pycold.app import defaults, logging
from pycold import obcold_reconstruct

logger = logging.getLogger(__name__)


def cmname_fromdate(ordinal_date):
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


def cmdirectionname_fromdate(ordinal_date):
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

def cmdatename_fromdate(ordinal_date):
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


def obiaresname_fromdate(ordinal_date):
    """
    get file name from ordinate date
    Args:
        ordinal_date: the inputted ordinate date
    Returns:
        CM file name
    """
    return 'obiaresult_{}_{}{}'.format(str(ordinal_date),
                                    pd.Timestamp.fromordinal(ordinal_date - 366).year,
                                    str(pd.Timestamp.fromordinal(ordinal_date - 366).timetuple().tm_yday).zfill(3))




def is_change_object(stats_lut_row, uniform_threshold, uniform_sizeslope, keyword, classification_map):
    """
    parameters
    ----------
    is_thematic: True -> that has thematic inputs; False-> default config will be applied
    stats_lut_row: a table of ['id', 'change magnitude', 'mode of lc category', 'pixel number']
    uniform_threshold: grid searching usage, overwriting default config
    uniform_sizeslope: grid searching usage, overwriting default config
    keyword: string, keyword to calculate change magnitude, 'cm_average' or 'cm_median'

    Returns: boolean
    -------
        True -> is change object
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
    if classification_map is None:
        if log10_size * ['default_sizeslope'] + intercept < 2:
            scale = log10_size * defaults['default_sizeslope'] + intercept
        else:
            scale = 2
        if np.double(stats_lut_row[keyword]) * scale > defaults['default_threshold']:
            return True
        else:
            return False
    else:
        if uniform_threshold is not None and uniform_sizeslope is not None:
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
                if log10_size * defaults['default_sizeslope'] + intercept < 2:
                    scale = log10_size * defaults['default_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > defaults['default_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 1:  # develop in LCMAP
                # intercept = 1 - config['C1_sizeslope']
                if log10_size * defaults['C1_sizeslope'] + intercept < 2:
                    scale = log10_size * defaults['C1_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > defaults['C1_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 2:  # cropland in LCMAP
                # intercept = 1 - config['C2_sizeslope']
                if log10_size * defaults['C2_sizeslope'] + intercept < 2:
                    scale = log10_size * defaults['C2_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > defaults['C2_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 3:  # grassland in LCMAP
                # intercept = 1 - config['C3_sizeslope']
                if log10_size * defaults['C3_sizeslope'] + intercept < 2:
                    scale = log10_size * defaults['C3_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > defaults['C3_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 4:  # forest in LCMAP
                # intercept = 1 - config['C4_sizeslope']
                if log10_size * defaults['C4_sizeslope'] + intercept < 2:
                    scale = log10_size * defaults['C4_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > defaults['C4_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 5:  # water in LCMAP
                # intercept = 1 - config['C5_sizeslope']
                if log10_size * defaults['C5_sizeslope'] + intercept < 2:
                    scale = log10_size * defaults['C5_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > defaults['C5_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 6:  # wetland in LCMAP
                # intercept = 1 - config['C6_sizeslope']
                if log10_size * defaults['C6_sizeslope'] + intercept < 2:
                    scale = log10_size * defaults['C6_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > defaults['C6_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 7:  # ice in LCMAP
                # intercept = 1 - config['C7_sizeslope']
                if log10_size * defaults['C7_sizeslope'] + intercept < 2:
                    scale = log10_size * defaults['C7_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > defaults['C7_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 8:  # barren in LCMAP
                # intercept = 1 - config['C8_sizeslope']
                if log10_size * defaults['C8_sizeslope'] + intercept < 2:
                    scale = log10_size * defaults['C8_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > defaults['C8_threshold']:
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


def segmentation(cm_array, cm_direction_array, cm_date_array, cm_array_l1=None, cm_array_l1_direction=None,
                 cm_array_l1_date=None, cm_interval=32, floodfill_ratio=None, devel_mode=False,
                 filenm=None, out_path=None, trans=None, proj=None):
    cm_array = cm_array.astype(float) / defaults['cm_scale']
    cm_array_l1_date = cm_array_l1_date.astype(float) / defaults['cm_scale']

    peak_threshold = chi2.ppf(0.95, 5)
    [n_rows, n_cols] = cm_array.shape
    if floodfill_ratio is None:
        floodfill_ratio = defaults['floodfill_thres']
    if filenm is None:
        filenm = 'test'
    if cm_array_l1 is None:
        cm_array_l1 = np.full((n_rows, n_cols), defaults['NAN_VAL'], dtype=np.int16)
    if cm_array_l1_direction is None:
        cm_array_l1_direction = np.full((n_rows, n_cols), defaults['NAN_VAL_UINT8'], dtype=np.uint8)
    if cm_array_l1_date is None:
        cm_array_l1_date = np.full((n_rows, n_cols), defaults['NAN_VAL'], dtype=np.int32)

    # assign valid CM values in the stack into current cm array where NA values are
    cm_array[cm_array == defaults['NAN_VAL']] = cm_array_l1[cm_array == defaults['NAN_VAL']]

    # defaults['NAN_VAL_UINT8'], i.e., 255, is also the NA value of direction.
    cm_direction_array[cm_direction_array == defaults['NAN_VAL_UINT8']] = \
        cm_array_l1_direction[cm_direction_array == defaults['NAN_VAL_UINT8']]

    cm_date_array[cm_date_array == defaults['NAN_VAL']] = \
        cm_array_l1_date[cm_date_array == defaults['NAN_VAL']]

    if devel_mode is True:
        gdal_save_file_1band(join(out_path, filenm + '_cm_array.tif'), cm_array,
                             gdal.GDT_Float32, trans, proj, n_cols, n_rows)
        gdal_save_file_1band(join(out_path, filenm + '_cm_direction_array.tif'), cm_direction_array,
                             gdal.GDT_Byte, trans, proj, n_cols, n_rows)

    #######################################################################################
    #                               Scale 1: change superpixel                            #
    #######################################################################################
    cm_array[cm_array < peak_threshold * (1 - floodfill_ratio)] = np.nan
    bandwidth = 1

    # using gaussian kernel ( with 1 sigma value) to smooth images in hpc_preparation for floodfill
    kernel = Gaussian2DKernel(x_stddev=bandwidth, y_stddev=bandwidth)
    cm_array_gaussian_s1 = convolve(cm_array, kernel, boundary='extend', preserve_nan=True)
    cm_array_gaussian_s1[np.isnan(cm_array)] = defaults['NAN_VAL']

    # cm_array_gaussian_s1 = cm_array  # try not using gaussian
    if devel_mode is True:
        gdal_save_file_1band(join(out_path, filenm + '_cm_array_gaussian_s1.tif'),
                             cm_array_gaussian_s1,
                             gdal.GDT_Float32, trans, proj, n_cols, n_rows)

    seed_index = peak_local_max(cm_array_gaussian_s1, threshold_abs=peak_threshold,
                                exclude_border=False, min_distance=0)
    # seed_index = np.flip(seed_index, axis=0)

    seed_labels = np.zeros((n_rows, n_cols))
    seed_labels[tuple(np.transpose(seed_index))] = 1
    if devel_mode:
        gdal_save_file_1band(join(out_path, filenm + '_seed.tif'), seed_labels,
                             gdal.GDT_Byte, trans, proj, n_cols, n_rows)

    mask_s1 = np.zeros((n_rows + 2, n_cols + 2)).astype(np.uint8)
    mask_label_s1 = np.zeros((n_rows + 2, n_cols + 2))
    floodflags_base = 8
    floodflags_base |= cv2.FLOODFILL_MASK_ONLY
    floodflags_base |= cv2.FLOODFILL_FIXED_RANGE
    no = 0
    i = 0
    # cm_stack = np.dstack([cm_array_gaussian_s1, cm_direction_array, cm_date_array]).astype(np.float32)
    cm_stack = np.dstack([cm_array_gaussian_s1, cm_direction_array, cm_direction_array]).astype(np.float32)
    for i in range(len(seed_index)):
        # print(i)
        remainder = i % 255
        floodflags = floodflags_base | ((remainder + 1) << 8)
        seedcm = cm_array_gaussian_s1[tuple(seed_index[i])]
        num, im, mask_s1, rect = floodFill(cm_stack, mask_s1, tuple(reversed(seed_index[i])), 0,
                                           loDiff=[seedcm * floodfill_ratio, 0, 0],
                                           upDiff=[seedcm * floodfill_ratio, 0, 0],
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
        print("segmentation finished for {}".format(filenm))

    if len(seed_index) > 0:
        # processing the rest that hasn't be recorded into mask_label
        no = int(i / 255)
        mask_label_s1[(mask_label_s1 == 0) & (mask_s1 > 0)] = mask_s1[(mask_label_s1 == 0) & (mask_s1 > 0)].astype(
            int) + no * 255

    object_map_s1 = mask_label_s1[1:n_rows + 1, 1:n_cols + 1]

    if devel_mode:
        gdal_save_file_1band(
            os.path.join(out_path, filenm + '_floodfill_gaussian_{}_s1.tif'.format(bandwidth)),
            object_map_s1, gdal.GDT_Int32, trans, proj, n_rows, n_cols)

    #######################################################################################
    #                                 Scale 2: change patch                              #
    #######################################################################################

    # create an object-based change
    unq_s1, ids_s1, count_s1 = np.unique(object_map_s1, return_inverse=True, return_counts=True)
    mean_list = np.bincount(ids_s1.astype(int), weights=cm_array_gaussian_s1.reshape(ids_s1.shape)) / count_s1
    mean_list[unq_s1 == 0] = defaults['NAN_VAL']  # force mean of unchanged objects to be -9999
    lut_dict_s1 = dict(zip(unq_s1, mean_list))
    cm_array_s2 = np.vectorize(lut_dict_s1.get)(object_map_s1)

    bandwidth = 1
    # using gaussian kernel ( with 0.5 sigma value) to smooth images in hpc_preparation for floodfill
    # cm_array_s2[cm_array_s2 == config['NAN_VAL']] = np.nan
    # kernel = Gaussian2DKernel(x_stddev=bandwidth, y_stddev=bandwidth)
    # cm_array_gaussian_s2 = convolve(cm_array_s2, kernel, boundary='extend', preserve_nan=True)
    # cm_array_gaussian_s2[np.isnan(cm_array_s2)] = config['NAN_VAL']

    cm_array_gaussian_s2 = cm_array_s2  # not using gaussian
    if devel_mode:
        gdal_save_file_1band(join(out_path,  filenm + '_cm_array_gaussian_s2.tif'),
                             cm_array_gaussian_s2,
                             gdal.GDT_Float32, trans, proj, n_cols,
                             n_rows)

    mask_s2 = np.zeros((n_rows + 2, n_cols + 2)).astype(np.uint8)
    mask_label_s2 = np.zeros((n_rows + 2, n_cols + 2))
    floodflags_base = 8
    floodflags_base |= cv2.FLOODFILL_MASK_ONLY
    i = 0
    cm_array_gaussian_s2 = cm_array_gaussian_s2.astype(np.float32)
    # cm_stack = np.dstack([cm_array_gaussian_s2, cm_direction_array, cm_direction_array]).astype(np.float32)
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
                                                                      (mask_label_s2 == 0) & (mask_s2 > 0)].\
                                                                      astype(int) + no * 255

    if devel_mode:
        print("segmentation finished ")

    if len(seed_index) > 0:
        # processing the rest that hasn't be recorded into mask_label
        no = int(i / 255)
        mask_label_s2[(mask_label_s2 == 0) & (mask_s2 > 0)] = mask_s2[(mask_label_s2 == 0) & (mask_s2 > 0)].astype(
            int) + no * 255

    object_map_s2 = mask_label_s2[1:n_rows + 1, 1:n_cols + 1]
    if devel_mode:
        gdal_save_file_1band(
            join(out_path,  '{}_floodfill_gaussian_{}_s2.tif'.format(filenm, bandwidth)),
            object_map_s2, gdal.GDT_Int32, trans, proj, n_rows, n_cols)
    return object_map_s1, object_map_s2, unq_s1, mean_list


def object_analysis(object_map_s1, object_map_s2, unq_s1, mean_list, classification_map=None,
                    uniform_threshold=None, uniform_sizeslope=None):
    # if classification_map is not None:
    #     assert classification_map.shape == object_map_s1.shape
    [n_rows, n_cols] = object_map_s1.shape
    change_map = np.zeros((n_rows, n_cols)).astype(np.uint8)

    unq_s2, ids_s2, count_s2 = np.unique(object_map_s2, return_inverse=True, return_counts=True)
    lut_dict_s2 = dict(zip(unq_s2, count_s2))
    size_map = np.vectorize(lut_dict_s2.get)(object_map_s2)  # we use size of scale 2 to represent object szie
    size_map[object_map_s2 == 0] = 0

    if classification_map is None:
        mode_list = [0] * len(unq_s1)
    else:
        mode_list = modeby(classification_map.reshape(ids_s2.shape), object_map_s1.reshape(ids_s2.shape))
    size_list = modeby(size_map.reshape(ids_s2.shape),
                       object_map_s1.reshape(ids_s2.shape))  # the mode of s1 objects in sizemap
    # mode_list, median_list = mode_median_by(classification_map.reshape(ids_s2.shape),
    #                                         cm_array_gaussian_s1.reshape(ids_s2.shape),
    #                                         object_map_s1.reshape(ids_s2.shape))

    # Use the tags/IDs to perform ID based summation of array2 elems and
    # thus divide by the ID counts to get ID based average values
    stats_lut = pd.DataFrame(np.column_stack((unq_s1, mean_list,
                                              mode_list, size_list)),
                             columns=['id', 'cm_average', 'mode', 'npixels'])

    class_labels = unq_s1[unq_s1 > 0]
    change_group = []
    if len(class_labels) > 0:  # if has valid change objects to output
        for index in class_labels:
            stats_lut_row = stats_lut.loc[stats_lut['id'] == index]
            if is_change_object(stats_lut_row, uniform_threshold, uniform_sizeslope,
                                'cm_average', classification_map):
                change_group.append(index)
        change_map[np.isin(object_map_s1, change_group)] = 1

    return change_map


class ObjectAnalystHPC:
    def __init__(self, config, stack_path, result_path, starting_date, cmmap_path=None, obia_path=None,
                 obcold_recg_path=None, thematic_path=None):
        self.config = config
        self.config['block_width'] = int(self.config['n_cols'] / self.config['n_block_x'])
        self.config['block_height'] = int(self.config['n_rows'] / self.config['n_block_y'])
        self.config['n_blocks'] = self.config['n_block_x'] * self.config['n_block_y']
        self.thematic_path = thematic_path
        if cmmap_path is None:  # default path
            self.cmmap_path = join(result_path, 'cm_maps')
        else:
            self.cmmap_path = cmmap_path
        if obia_path is None:  # default path
            self.obia_path = join(result_path, 'OBIAresults')
        else:
            self.obia_path = obia_path
        if obcold_recg_path is None:  # default path
            self.obcold_recg_path = join(result_path, 'obcold')
        else:
            self.obcold_recg_path = obcold_recg_path
        self.stack_path = stack_path
        self.starting_date = starting_date
        self.year_lowbound = pd.Timestamp.fromordinal(starting_date - 366).year

    @staticmethod
    def _check_inputs(config, stack_path, result_path, cmmap_path, obia_path, obcold_recg_path,
                      thematic_path):
        if type(config['n_rows']) != int or config['n_rows'] < 0:
            raise ValueError('n_rows must be positive integer')
        if type(config['n_cols']) != int or config['n_cols'] < 0:
            raise ValueError('n_cols must be positive integer')
        if type(config['n_block_x']) != int or config['n_block_x'] < 0:
            raise ValueError('n_block_x must be positive integer')
        if type(config['n_block_y']) != int or config['n_block_y'] < 0:
            raise ValueError('n_block_y must be positive integer')
        if type(config['n_block_y']) != int or config['n_block_y'] < 0:
            raise ValueError('n_block_y must be positive integer')

        if os.path.isdir(stack_path) is False:
            raise FileExistsError('No such directory: {}'.format(stack_path))

        if os.path.isdir(result_path) is False:
            raise FileExistsError('No such directory: {}'.format(result_path))

    def hpc_preparation(self):
        if exists(self.cmmap_path) is False:
            try:
                os.mkdir(self.cmmap_path)
            except IOError as e:
                raise e

        if exists(self.obia_path) is False:
            try:
                os.mkdir(self.obia_path)
            except IOError as e:
                raise e

        if exists(self.obcold_recg_path) is False:
            try:
                os.mkdir(self.obcold_recg_path)
            except IOError as e:
                raise e

    def get_lastyear_cmap_fromdate(self, date):
        if pd.Timestamp.fromordinal(date - 366).year - 1 < self.year_lowbound:  # we used the the year before the date
            classified_year = self.year_lowbound
        else:
            classified_year = pd.Timestamp.fromordinal(date - 366).year - 1
        try:
            classification_map = np.load(join(self.thematic_path, 'yearlyclassification_{}.npy'.format(classified_year)))
        except IOError as e:
            raise e
        return classification_map

    def obia_execute(self, date, floodfill_ratio=None, uniform_threshold=None, uniform_sizeslope=None):
        if date - self.config['CM_OUTPUT_INTERVAL'] < self.starting_date:
            change_map = np.full((self.config['n_rows'], self.config['n_cols']), 0,
                                  dtype=np.byte)
        else:
            [object_map_s1, object_map_s2, unq_s1, mean_list] \
                = segmentation(np.load(join(self.cmmap_path, cmname_fromdate(date)+'.npy')),
                               np.load(join(self.cmmap_path, cmdirectionname_fromdate(date)+'.npy')),
                               np.load(join(self.cmmap_path, cmdatename_fromdate(date)+'.npy')),
                               np.load(join(self.cmmap_path, cmname_fromdate(date -
                                                                             self.config['CM_OUTPUT_INTERVAL'])+'.npy')),
                               np.load(join(self.cmmap_path, cmdirectionname_fromdate(date -
                                                                                      self.config['CM_OUTPUT_INTERVAL'])+'.npy')),
                               np.load(join(self.cmmap_path, cmdatename_fromdate(date -
                                                                                 self.config['CM_OUTPUT_INTERVAL'])+'.npy')))

            if self.thematic_path is not None:
                classification_map = self.get_lastyear_cmap_fromdate(date)
                change_map = object_analysis(object_map_s1, object_map_s2, unq_s1, mean_list,
                                             classification_map=classification_map,
                                             uniform_threshold=uniform_threshold, uniform_sizeslope=uniform_sizeslope)
            else:
                change_map = object_analysis(object_map_s1, object_map_s2, unq_s1, mean_list,
                                             uniform_threshold=uniform_threshold, uniform_sizeslope=uniform_sizeslope)
        self.save_obiaresult(date, change_map)

    def save_obiaresult(self, date, change_map):
        # if outputformat == 'tiff':
        #     gdal_save_file_1band(join(self.obia_path, filenm + '_OBIAresult.tif'), change_map,
        #                          gdal.GDT_Byte, trans, proj, defaults['n_cols'], defaults['n_rows'])
        # else:
        filenm = obiaresname_fromdate(date)
        bytesize = 1
        result_blocks = np.lib.stride_tricks.as_strided(change_map,
                                                        shape=(self.config['n_block_y'], self.config['n_block_x'],
                                                               self.config['block_height'], self.config['block_width']),
                                                        strides=(self.config['n_cols'] * self.config['block_height'] * bytesize,
                                                                 self.config['block_width'] * bytesize,
                                                                 self.config['n_cols'] * bytesize, bytesize))
        for i in range(self.config['n_block_y']):
            for j in range(self.config['n_block_x']):
                np.save(join(self.obia_path, filenm + '_x{}_y{}.npy'.format(j + 1, i + 1)),
                        result_blocks[i][j])

    def is_finished_object_analysis(self, date_list):
        """
        :return:
        """
        for date in date_list:
            for i in range(self.config['n_block_y']):
                for j in range(self.config['n_block_x']):
                    if not exists(join(self.obia_path, '{}_{}'.format(cmname_fromdate(date),
                                                                      'OBIAresult_x{}_y{}.npy'.format(j + 1, i + 1)))):
                        return False
        return True

    def get_allobiaresult_asarray(self, block_x, block_y):
        obia_files = [f for f in os.listdir(self.obia_path) if f.startswith('obiaresult')
                      and f.endswith('_x{}_y{}.npy'.format(block_x, block_y))]
        # sort image files by dates
        cm_dates = [int(f[f.find('obiaresult_') + len('obiaresult_'):
                          f.find('obiaresult_') + len('obiaresult_')+6])
                    for f in obia_files]
        files_date_zip = sorted(zip(cm_dates, obia_files))
        cm_dates = [x[0] for x in files_date_zip]

        # load date block
        cm_date_block = np.dstack([np.load(join(self.cmmap_path, cmdatename_fromdate(date)+'.npy'))
                                   for date in cm_dates])
        cm_date_block = np.reshape(cm_date_block, (self.config['block_width'] * self.config['block_height'],
                                                      cm_date_block.shape[2]))

        assert cm_date_block.shape[1] == len(obia_files)
        obia_tstack = [np.multiply(np.load(join(self.obia_path, f[1])).reshape(self.config['block_width']
                                                                               * self.config['block_height']),
                                   cm_date_block[:, count])
                       for count, f in enumerate(files_date_zip)]
        return np.vstack(obia_tstack)

    def reconstruct_reccg(self, block_id, img_stack=None, img_dates_sorted=None):
        block_y = get_block_y(block_id, self.config['n_block_x'])
        block_x = get_block_x(block_id, self.config['n_block_x'])
        if img_stack is None and img_dates_sorted is None:
            block_folder = join(self.stack_path, 'block_x{}_y{}'.format(block_x, block_y))
            img_stack, img_dates_sorted = read_blockdata(block_folder, self.config['n_cols']*self.config['n_rows'],
                                                         defaults['TOTAL_IMAGE_BANDS']+1)
        obia_breaks = self.get_allobiaresult_asarray(block_x, block_y)
        result_collect = []
        for pos in range(self.config['block_width'] * self.config['block_height']):
            original_row, original_col = get_rowcol_intile(pos, self.config['block_width'],
                                                           self.config['block_height'], block_x, block_y)

            try:
                obcold_result = obcold_reconstruct(img_dates_sorted,
                                                   img_stack[pos, 0, :].astype(np.int64),
                                                   img_stack[pos, 1, :].astype(np.int64),
                                                   img_stack[pos, 2, :].astype(np.int64),
                                                   img_stack[pos, 3, :].astype(np.int64),
                                                   img_stack[pos, 4, :].astype(np.int64),
                                                   img_stack[pos, 5, :].astype(np.int64),
                                                   img_stack[pos, 6, :].astype(np.int64),
                                                   img_stack[pos, 7, :].astype(np.int64),
                                                   obia_breaks[:, pos].astype(np.int64),
                                                   conse=self.config['conse'],
                                                   pos=self.config['n_cols'] * (original_row - 1) + original_col)
            except RuntimeError:
                logger.error("COLD fails at original_row {}, original_col {}".format(original_row, original_col))
            else:
                result_collect.append(obcold_result)
        return result_collect

    def save_obcoldresults(self, block_id, result_collect):
        block_x = get_block_x(block_id, self.config['n_block_x'])
        block_y = get_block_y(block_id, self.config['n_block_x'])
        np.save(join(self.obcold_recg_path, 'record_change_x{}_y{}_obcold.npy'.format(block_x, block_y)),
                np.hstack(result_collect))
