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
from shared import gdal_save_file_1band
from pycold.app import defaults, logging
from pycold.utils import get_block_y, get_block_x

logger = logging.getLogger(__name__)


thematic_mode_types = ['no_input', 'single_map', 'rf_model', 'multi_maps']


def extract_features(cold_plot, band, ordinal_day_list, nan_val, slope_scale, n_features, ismat=False):
    """
    Parameters
    ----------
    cold_plot: nested array, plot-based rec_cg
    band: integer, the band number range from 0 to 6
    ordinal_day_list: a list of ordinal days that the output will predict
    nan_val:
    slope_scale:
    n_features: integer the number of features for classification, from an order
    ismat

    Returns
    -------

    """
    features = [np.full(len(ordinal_day_list), nan_val) for x in range(n_features)]
    for index, ordinal_day in enumerate(ordinal_day_list):
        # print(index)
        for idx, cold_curve in enumerate(cold_plot):
            if idx == len(cold_plot) - 1:
                max_days = cold_plot[idx]['t_end']
            else:
                max_days = cold_plot[idx + 1]['t_start']
            if cold_curve['t_start'] <= ordinal_day < max_days:
                if ismat:
                    for n in range(n_features):
                        if n == 0:
                            # if cold_curve['t_start'] <= ordinal_day < cold_curve['t_end']:
                            features[n][index] = cold_curve['coefs'][0][band] + cold_curve['coefs'][1][band] * ordinal_day
                        else:
                            features[n][index] = cold_curve['coefs'][n+1][band]
                    break
                else:
                    for n in range(n_features):
                        if n == 0:
                            # if cold_curve['t_start'] <= ordinal_day < cold_curve['t_end']:
                            features[n][index] = cold_curve['coefs'][band][0] + cold_curve['coefs'][band][1] * ordinal_day / \
                                          slope_scale
                        else:
                            features[n][index] = cold_curve['coefs'][band][n+1]  # n + 1 is because won't need slope as output
                    break
    return features


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


def is_change_object(parameters, b_thematic, stats_lut_row, uniform_threshold, uniform_sizeslope, keyword):
    """
    Parameters
    ----------
    parameters:
        obcold parameter structure
    b_thematic: True -> that has thematic inputs; False-> default parameters will be applied
    stats_lut_row: a table of ['id', 'change magnitude', 'mode of lc category', 'pixel number']
    uniform_threshold: grid searching usage, overwriting default parameters
    uniform_sizeslope: grid searching usage, overwriting default parameters
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
        # intercept = 1 - parameters['default_sizeslope']
        if log10_size * parameters['default_sizeslope'] + intercept < 2:
            scale = log10_size * parameters['default_sizeslope'] + intercept
        else:
            scale = 2
        if np.double(stats_lut_row[keyword]) * scale > parameters['default_threshold']:
            return True
        else:
            return False
    else:
        if uniform_threshold != parameters['NAN_VAL'] and uniform_sizeslope != parameters['NAN_VAL']:
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
                # intercept = 1 - parameters['default_sizeslope']
                if log10_size * parameters['default_sizeslope'] + intercept < 2:
                    scale = log10_size * parameters['default_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > parameters['default_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 1:  # develop in LCMAP
                # intercept = 1 - parameters['C1_sizeslope']
                if log10_size * parameters['C1_sizeslope'] + intercept < 2:
                    scale = log10_size * parameters['C1_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > parameters['C1_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 2:  # cropland in LCMAP
                # intercept = 1 - parameters['C2_sizeslope']
                if log10_size * parameters['C2_sizeslope'] + intercept < 2:
                    scale = log10_size * parameters['C2_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > parameters['C2_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 3:  # grassland in LCMAP
                # intercept = 1 - parameters['C3_sizeslope']
                if log10_size * parameters['C3_sizeslope'] + intercept < 2:
                    scale = log10_size * parameters['C3_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > parameters['C3_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 4:  # forest in LCMAP
                # intercept = 1 - parameters['C4_sizeslope']
                if log10_size * parameters['C4_sizeslope'] + intercept < 2:
                    scale = log10_size * parameters['C4_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > parameters['C4_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 5:  # water in LCMAP
                # intercept = 1 - parameters['C5_sizeslope']
                if log10_size * parameters['C5_sizeslope'] + intercept < 2:
                    scale = log10_size * parameters['C5_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > parameters['C5_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 6:  # wetland in LCMAP
                # intercept = 1 - parameters['C6_sizeslope']
                if log10_size * parameters['C6_sizeslope'] + intercept < 2:
                    scale = log10_size * parameters['C6_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > parameters['C6_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 7:  # ice in LCMAP
                # intercept = 1 - parameters['C7_sizeslope']
                if log10_size * parameters['C7_sizeslope'] + intercept < 2:
                    scale = log10_size * parameters['C7_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > parameters['C7_threshold']:
                    return True
                else:
                    return False
            elif int(stats_lut_row['mode']) == 8:  # barren in LCMAP
                # intercept = 1 - parameters['C8_sizeslope']
                if log10_size * parameters['C8_sizeslope'] + intercept < 2:
                    scale = log10_size * parameters['C8_sizeslope'] + intercept
                else:
                    scale = 2
                if float(stats_lut_row[keyword]) * scale > parameters['C8_threshold']:
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


def generate_sample_num(label, parameters):
    """
    generate sample number for each land cover category using the method from 'Optimizing selection of training and
    auxiliary data for operational land cover classification for the LCMAP initiative'
    Args:
        label: an array of land-cover categories, i.e., map
        parameters: obcold parameter structure
    Returns:

    """
    unique_category, unique_counts_category = np.unique(label[label <= parameters['total_landcover_category']], return_counts=True)
    counts_category = np.array([0] * parameters['total_landcover_category'])
    for x in range(len(unique_counts_category)):
        counts_category[unique_category[x]-1] = unique_counts_category[x]
    percate_samples = np.array([round(x * parameters['total_samples'] / sum(counts_category)) for x in counts_category])
    percate_samples[percate_samples > parameters['max_category_samples']] = parameters['max_category_samples']
    percate_samples[percate_samples < parameters['min_category_samples']] = parameters['min_category_samples']
    percate_samples = np.minimum(percate_samples, counts_category) # needs to check not exceed the total category pixels
    return percate_samples



def _check_inputs(parameters, stack_path, results_path, year_lowbound, year_uppbound, thematic_mode, thematic_src,
                  devel_mode):
    """
    Args:
        parameters:
        stack_path:
        results_path:
        year_lowbound:
        year_uppbound:
        logger:
        thematic_mode:
        thematic_src:
        devel_mode:

    Returns:

    """
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

    if year_lowbound is None and year_uppbound is None:
        if os.path.exists(join(stack_path, 'starting_last_dates.txt')):
            raise FileExistsError('No starting_last_dates.txt, please specify year_lowbound and year_uppbound')

    if thematic_mode in thematic_mode_types is False:
        raise ValueError('thematic mode must be in no_input, single_map, rf_model, multi_maps')

    if thematic_mode != 'no_input':
        if thematic_src is None:
            raise ValueError('thematic_src must be assigned as value if thematic_mode is not no_input')


class ObjectAnalyst:
    def __init__(self, parameters, stack_path, results_path, year_lowbound=None, year_uppbound=None,
                 thematic_mode='no_input', thematic_src=None, devel_mode=False):
        """
        Parameters
        ----------
        parameters
        stack_path
        results_path
        year_lowbound
        year_uppbound
        thematic_mode
        thematic_src
        devel_mode
        """

        try:
            _check_inputs(parameters, stack_path, results_path, year_lowbound, year_uppbound,
                          thematic_mode, thematic_src, devel_mode)
        except ValueError or FileExistsError as e:
            raise e

        self.parameters = parameters
        self.parameters['block_width'] = int(self.parameters['n_cols'] / self.parameters['n_block_x'])
        self.parameters['block_height'] = int(self.parameters['n_rows'] / self.parameters['n_block_y'])
        self.parameters['n_blocks'] = self.parameters['n_block_x'] * self.parameters['n_block_y']
        self.results_path = results_path
        self.stack_path = stack_path
        self.out_cmmap_path = join(results_path, 'cm_maps')
        self.out_yearly_path = join(results_path, 'feature_maps')
        self.obiaresults_path = join(results_path, 'OBIAresults')
        self.breakdatemaps_path = join(results_path, 'breakdate_maps')
        self.obcoldreccg_path = join(results_path, 'obcold')
        if year_uppbound is None and year_lowbound is None:
            with open(join(stack_path, 'starting_last_dates.txt')) as f:
                date_list = f.readlines()
            self.year_lowbound = pd.Timestamp.fromordinal(int(date_list[0]) - 366).year
            self.year_uppbound = pd.Timestamp.fromordinal(int(date_list[1]) - 366).year
        else:
            self.year_lowbound = year_lowbound
            self.year_uppbound = year_uppbound
        self.thematic_mode = thematic_mode
        self.thematic_src = thematic_src
        self.devel_mode = devel_mode

    def hpc_preparation(self):
        """
        Returns
        -------

        """
        if not exists(self.out_cmmap_path):
            os.makedirs(self.out_cmmap_path)
        if not exists(self.out_yearly_path):
            os.makedirs(self.out_yearly_path)
        if not exists(self.obiaresults_path):
            os.makedirs(self.obiaresults_path)
        if not exists(self.breakdatemaps_path):
            os.makedirs(self.breakdatemaps_path)
        if not exists(self.obcoldreccg_path):
            os.makedirs(self.obcoldreccg_path)

    def step1_predict_features(self, block_id, ismat=False, row=None):
        """
        Parameters
        ----------
        block_id
        ismat
        row

        Returns
        -------

        """
        feature_list = [np.full(((self.year_uppbound - self.year_lowbound + 1) * self.parameters['block_height'],
                                 self.parameters['block_width']), defaults['NAN_VAL'], dtype=np.float32)
                        for x in range(defaults['TOTAL_IMAGE_BANDS'] * defaults['N_FEATURES'])]

        ordinal_day_list = [pd.Timestamp.toordinal(dt.date(year, 7, 1)) + 366 for year
                            in range(self.year_lowbound, self.year_uppbound + 1)]
        current_block_x = get_block_x(block_id, self.parameters['n_block_x'])
        current_block_y = get_block_y(block_id, self.parameters['n_block_x'])
        if ismat:
            rec_cg_path = os.path.join(self.results_path, 'record_change_r{}.mat'.format(str(row).zfill(5)))
            ccd_scanline = loadmat(rec_cg_path)['rec_cg'][0]
        else:
            rec_cg_path = join(self.results_path, 'record_change_x{}_y{}_cold.npy'
                               .format(current_block_x, current_block_y))
            ccd_scanline = np.array(np.load(rec_cg_path))

        if len(ccd_scanline) == 0:
            logger.warning('the rec_cg file {} is missing'.format(rec_cg_path))
            return feature_list

        ccd_scanline.sort(order='pos')
        processing_pos = ccd_scanline[0]['pos']
        identical_pos_curve = []
        for count, curve in enumerate(ccd_scanline):
            if curve['pos'] != processing_pos:  # means all curve have been collected for processing_pos
                if len(identical_pos_curve) == 0:
                    continue
                identical_pos_curve_np = np.array(identical_pos_curve)

                if identical_pos_curve_np[0]["pos"].size == 0:  # MATLAB empty array size == 0
                    continue

                # i_col and i_row is the local location in each block
                i_col = int((identical_pos_curve_np[0]["pos"] - 1) % self.parameters['n_cols']) - \
                            (current_block_x - 1) * self.parameters['block_width']
                i_row = int((identical_pos_curve_np[0]["pos"] - 1) / self.parameters['n_cols']) - \
                            (current_block_y - 1) * self.parameters['block_height']
                if i_col < 0:
                    self.logger.warning('Processing failed: i_row={}; i_col={} for {}'.format(i_row, i_col, rec_cg_path))
                for band in range(defaults['TOTAL_IMAGE_BANDS']):
                    feature_row = extract_features(identical_pos_curve_np, band, ordinal_day_list,
                                                   defaults['NAN_VAL'], defaults['SLOPE_SCALE'],
                                                   defaults['N_FEATURES'], ismat)
                    for index in range(defaults['N_FEATURES']):
                        feature_list[band * defaults['N_FEATURES'] + index][i_row * (self.year_uppbound -
                                                                                               self.year_lowbound + 1):
                                                                                 (i_row + 1) * (self.year_uppbound -
                                                                                                self.year_lowbound + 1),
                        i_col] = \
                            feature_row[index]
                # clear each list
                identical_pos_curve.clear()
                identical_pos_curve.append(
                    curve)  # need to append current curve as the first curve for the next position
                # go to the next pos
                processing_pos = curve['pos']
            else:
                identical_pos_curve.append(curve)

            # need to make the codes more elegant in future
            if count == (len(ccd_scanline) - 1):  # processing the last element,

                if len(identical_pos_curve) == 0:
                    continue

                if ismat:
                    identical_pos_curve_np = np.array(identical_pos_curve)
                    if identical_pos_curve_np[0]["pos"].size == 0:  # MATLAB empty array size == 0
                        continue

                i_col = int((identical_pos_curve_np[0]["pos"] - 1) % self.parameters['n_cols']) - \
                    (current_block_x - 1) * self.parameters['block_width']
                i_row = int((identical_pos_curve_np[0]["pos"] - 1) / self.parameters['n_cols']) - \
                    (current_block_y - 1) * self.parameters['block_height']
                if i_col < 0:
                    self.logger.warning('Processing failed: i_row={}; i_col={} for {}'.format(i_row, i_col, rec_cg_path))

                for band in range(defaults['TOTAL_IMAGE_BANDS']):
                    feature_row = extract_features(identical_pos_curve_np, band, ordinal_day_list,
                                                   defaults['NAN_VAL'], defaults['SLOPE_SCALE'],
                                                   defaults['N_FEATURES'], ismat)
                    for index in range(defaults['N_FEATURES']):
                        feature_list[band * defaults['N_FEATURES'] + index][i_row * (self.year_uppbound -
                                     self.year_lowbound + 1): (i_row + 1) * (self.year_uppbound - self.year_lowbound + 1),
                                     i_col] = feature_row[index]

        return feature_list

    def checkfinished_step1_predict_features(self):
        """
        :return: True or false
        """

        for iblock in range(self.parameters['n_blocks']):
            for x in range(self.parameters['TOTAL_IMAGE_BANDS'] * self.parameters['N_FEATURES']):
                if not os.path.exists(join(self.out_yearly_path, 'tmp_feature_B{}_block{}.npy').format(x + 1,
                                                                                                          iblock + 1)):
                    return False
        return True

    def step2_train_rfmodel(self, year):
        """
        :return:
        """
        # copy from def step2_train_rfmodel(parameters, trans, proj, year_lowbound, year_uppbound, out_yearly_path,
        # thematic_src, logger):
        while not self.checkfinished_step1_predict_features():
            time.sleep(5)

        tmp_feature_filenames = [file for file in os.listdir(self.out_yearly_path) if file.startswith('tmp_feature_B')]
        if len(tmp_feature_filenames) < self.parameters['n_blocks'] * self.parameters['TOTAL_IMAGE_BANDS'] * \
                self.parameters['N_FEATURES']:
            self.logger.warning('tmp features are incomplete! should have {}; but actually have {} feature images'.
                                format(self.parameters['n_blocks'] * self.parameters['TOTAL_IMAGE_BANDS'] *
                                       self.parameters['N_FEATURES'],
                len(tmp_feature_filenames)))
        for x in range(self.parameters['TOTAL_IMAGE_BANDS'] * self.parameters['N_FEATURES']):
            tmp_feature_filenames_sub = [file for file in tmp_feature_filenames if
                                         int(file[13:file.find('_', 13, 16)]) ==
                                         x + 1]
            tmp_feature_filenames_sub.sort(
                key=lambda t: int(t[t.find('block') + 5: t.find('.npy')]))  # sorted by row number

            # tmp_feature_array_sub is an array of (block_height * (year_uppbound - year_lowbound + 1), block_width)
            tmp_feature_array_sub = [np.load(os.path.join(self.out_yearly_path, file)) for file in
                                     tmp_feature_filenames_sub]
            # full_feature_array = full_feature_array.reshape(parameters['n_rows'],
            #                                                 parameters['n_cols'],
            #                                                 (year_uppbound - year_lowbound + 1))
            full_feature_array = np.hstack(tmp_feature_array_sub)  # (10000, 100000)
            full_feature_array = np.vstack(np.hsplit(full_feature_array, self.parameters['n_block_y']))  # (200000, 5000)
            full_feature_array = full_feature_array.reshape(self.parameters['n_cols'],
                                                            (self.year_uppbound - self.year_lowbound + 1),
                                                            self.parameters['n_rows'])  # (5000, 40, 5000)
            full_feature_array_list = np.split(full_feature_array, (self.year_uppbound - self.year_lowbound + 1),
                                               axis=1)

            outfile = os.path.join(self.out_yearly_path, 'feature_B{}_{}.npy'.format(x + 1,
                                                                                     self.parameters['classification_year']))
            np.save(outfile, full_feature_array_list[year - self.year_lowbound].reshape(self.parameters['n_rows'],
                    self.parameters['n_cols']))

        label = gdal_array.LoadFile(self.thematic_src)
        samplecount = generate_sample_num(label, self.parameters)
        index_list = []
        label_list = []
        for i in range(self.parameters['total_landcover_category']):
            index = np.argwhere(label == i + 1)
            # print(samplecount[i])
            np.random.seed(42)  # set random seed to reproduce the same result
            index_sample = index[np.random.choice(len(index), int(samplecount[i]), replace=False)]
            index_list.append(index_sample)
            label_list.append(np.array([i + 1] * len(index_sample)))
        index_list = np.vstack(index_list)
        label_list = np.hstack(label_list)
        feature_extraction = np.zeros(
            (len(index_list), self.parameters['TOTAL_IMAGE_BANDS'] * self.parameters['N_FEATURES']))
        for b in range(self.parameters['TOTAL_IMAGE_BANDS'] * self.parameters['N_FEATURES']):
            band = np.load(join(self.out_yearly_path, 'feature_B{}_{}.tif'.format(b + 1,
                                                                                  self.parameters['classification_year'])))
            feature_extraction[:, b] = [band[tuple(x)] for x in index_list]
        rf_model = RandomForestClassifier(random_state=42)
        rf_model.fit(feature_extraction, label_list)
        return rf_model

    def checkfinished_step2_train_rfmodel(self):
        """
        :param rfmodel_path: the path of rf model
        :return: True of false
        """
        if not os.path.exists(self.rfmodel_path):
            return False
        return True

    def step2_classification(self, block_id):
        """
        classifying each tmp feature patches
        :param block_id: the id of blocks. e.g., 1-500 if 500 cores
        """
        # thematic_mode: None, 'single_map', 'rf_model', 'multi_maps'
        if self.thematic_mode != 'no_input':
            if self.thematic_mode == 'rf_model':
                rf_model = joblib.load(self.thematic_src)
            elif self.thematic_src == 'single_map':
                rf_model = joblib.load(os.path.join(self.results_path, 'rf.model'))

        # tmp_feature_{}_block{}.npy is a 2-d array of (block_height * (year_uppbound - year_lowbound + 1), block_width)
        tmp_feature_sub = [
            np.load(os.path.join(self.out_yearly_path, 'tmp_feature_{}_block{}.npy'.format('B{}'.format(x + 1),
                                                                                      block_id))).flatten()
            for x in range(self.parameters['TOTAL_IMAGE_BANDS'] * self.parameters['N_FEATURES'])]

        # tmp_feature_sub_stack is a 2-d array of (block_width * (year_uppbound - year_lowbound + 1) * block_height, nbands)
        tmp_feature_sub_stack = np.vstack(tmp_feature_sub).transpose()
        # image_flat = np.zeros((parameters['n_rows'] * parameters['n_cols'], parameters['TOTAL_IMAGE_BANDS']))
        label = rf_model.predict(tmp_feature_sub_stack).reshape(self.parameters['block_height'],
                                                                (self.year_uppbound - self.year_lowbound + 1),
                                                                self.parameters['block_width'])
        mask = np.all(tmp_feature_sub_stack == self.parameters['NAN_VAL'], axis=1).reshape(self.parameters['block_height'],
                                                                                         (self.year_uppbound -
                                                                                          self.year_lowbound + 1),
                                                                                          self.parameters['block_width'])
        label[mask] = 255
        annualmap_list = np.split(label, (self.year_uppbound - self.year_lowbound + 1), axis=1)
        return annualmap_list

    def checkfinished_step2_classification(self):
        """
        :return: True or false
        """
        for iblock in range(self.parameters['n_blocks']):
            for year in range(self.year_lowbound, self.year_uppbound + 1):
                if not os.path.exists(os.path.join(self.out_yearly_path,
                                                   'tmp_yearlyclassification{}_block{}.npy'.format(year, iblock + 1))):
                    return False
        return True

    def step2_assemble_landcovermap(self, year):
        """
        assemble land cover maps from tmp yearly classification tmp files to a targeted year
        :return:
        """
        while not self.checkfinished_step2_classification():
            time.sleep(5)

        tmp_yearlyclass_filenames = [file for file in os.listdir(self.out_yearly_path)
                                     if file.startswith('tmp_yearlyclassification{}'.format(year))]

        # sort to guarantee order follows low to high rows

        tmp_yearlyclass_filenames.sort(key=lambda x: int(x[x.find('block') + 5: x.find('.npy')]))
        full_yearlyclass_array = np.hstack([np.load(join(self.out_yearly_path, file))
                                            for file in tmp_yearlyclass_filenames])
        full_yearlyclass_array = np.vstack(np.hsplit(full_yearlyclass_array, self.parameters['n_block_v']))

        if (full_yearlyclass_array.shape[1] != self.parameters['n_cols']) or (full_yearlyclass_array.shape[0] !=
                                                                            self.parameters['n_rows']):
            self.logger.warning('The assemble category map is incomplete for {}: nrow = {}; ncol = {}'.
                                format(year, full_yearlyclass_array.shape[0], full_yearlyclass_array.shape[1]))
            # outfile = os.path.join(self.out_yearly_path, 'yearlyclassification_{}.tif'.format(year))
            # np.save(outfile, full_yearlyclass_array)
        return full_yearlyclass_array
        # remove tmp files
        # tmp_yearlyclass_filenames = [file for file in os.listdir(self.out_yearly_path)
        #                              if file.startswith('tmp_')]
        # for file in tmp_yearlyclass_filenames:
        #     os.remove(os.path.join(self.out_yearly_path, file))

    def step3_objectanalysis(self, date, cm_array, cm_direction_array, cm_date_array, cm_array_l1, cm_array_l1_direction,
                             cm_array_l1_date, uniform_threshold=-9999, uniform_sizeslope=-9999,
                             candidate_gatethreshold=-9999):
        """
        segmentation and object-based decision for disturbance object
        :param date
        :param cm_array
        :param cm_direction_array
        :param cm_date_array
        :param cm_array_l1
        :param cm_direction_array_l1
        :param cm_date_array_l1
        :param uniform_threshold: only used for parameter test
        :param uniform_sizeslope: only used for parameter test
        :param candidate_gatethreshold: only used for parameter test
        :return:
        """
        peak_threshold = chi2.ppf(0.95, 5)

        # assign valid CM values in the stack into current cm array where NA values are
        cm_array[cm_array == self.parameters['NAN_VAL']] = cm_array_l1[cm_array == self.parameters['NAN_VAL']]

        # parameters['NAN_VAL_UINT8'], i.e., 255, is also the NA value of direction.
        cm_direction_array[cm_direction_array == self.parameters['NAN_VAL_UINT8']] = \
            cm_array_l1_direction[cm_direction_array == self.parameters['NAN_VAL_UINT8']]

        cm_date_array[cm_date_array == self.parameters['NAN_VAL']] = \
            cm_array_l1_date[cm_date_array == self.parameters['NAN_VAL']]

        if self.devel_mode is True:
            gdal_save_file_1band(join(self.obiaresults_path, filenm + '_cm_array.tif'), cm_array,
                                 gdal.GDT_Float32, self.trans, self.proj, self.parameters['n_cols'],
                                 parameters['n_rows'])
            gdal_save_file_1band(join(self.obiaresults_path, filenm + '_cm_direction_array.tif'), cm_direction_array,
                                 gdal.GDT_Byte, self.trans, self.proj, self.parameters['n_cols'],
                                 self.parameters['n_rows'])

        #######################################################################################
        #                               Scale 1: change superpixel                            #
        #######################################################################################
        cm_array[cm_array < peak_threshold * (1 - self.parameters['floodfill_thres'])] = np.nan
        bandwidth = 1
        
        # using gaussian kernel ( with 1 sigma value) to smooth images in preparation for floodfill
        kernel = Gaussian2DKernel(x_stddev=bandwidth, y_stddev=bandwidth)
        cm_array_gaussian_s1 = convolve(cm_array, kernel, boundary='extend', preserve_nan=True)
        cm_array_gaussian_s1[np.isnan(cm_array)] = self.parameters['NAN_VAL']

        # cm_array_gaussian_s1 = cm_array  # try not using gaussian
        if self.devel_mode is True:
            gdal_save_file_1band(join(self.obiaresults_path, filenm + '_cm_array_gaussian_s1.tif'),
                                 cm_array_gaussian_s1,
                                 gdal.GDT_Float32, self.trans, self.proj, self.parameters['n_cols'],
                                 self.parameters['n_rows'])

        seed_index = peak_local_max(cm_array_gaussian_s1, threshold_abs=peak_threshold,
                                    exclude_border=False, min_distance=0)
        # seed_index = np.flip(seed_index, axis=0)

        seed_labels = np.zeros((self.parameters['n_rows'], self.parameters['n_cols']))
        seed_labels[tuple(np.transpose(seed_index))] = 1
        if self.devel_mode:
            gdal_save_file_1band(join(self.obiaresults_path, filenm + '_seed.tif'), seed_labels,
                                 gdal.GDT_Byte, self.trans, self.proj, self.parameters['n_cols'],
                                 self.parameters['n_rows'])

        mask_s1 = np.zeros((self.parameters['n_rows'] + 2, self.parameters['n_cols'] + 2)).astype(np.uint8)
        mask_label_s1 = np.zeros((self.parameters['n_rows'] + 2, self.parameters['n_cols'] + 2))
        floodflags_base = 8
        floodflags_base |= cv2.FLOODFILL_MASK_ONLY
        floodflags_base |= cv2.FLOODFILL_FIXED_RANGE
        no = 0
        i = 0
        cm_stack = np.dstack([cm_array_gaussian_s1, cm_direction_array, cm_date_array]).astype(np.float32)
        # parameters['floodfill_thres'] = 0.5
        for i in range(len(seed_index)):
            # print(i)
            remainder = i % 255
            floodflags = floodflags_base | ((remainder + 1) << 8)
            seedcm = cm_array_gaussian_s1[tuple(seed_index[i])]
            num, im, mask_s1, rect = floodFill(cm_stack, mask_s1, tuple(reversed(seed_index[i])), 0,
                                               loDiff=[seedcm * self.parameters['floodfill_thres'], 0,
                                                       self.parameters['CM_OUTPUT_INTERVAL']],
                                               upDiff=[seedcm * self.parameters['floodfill_thres'], 0,
                                                       self.parameters['CM_OUTPUT_INTERVAL']],
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
        if self.devel_mode:
            print("segmentation finished for {}".format(filenm))

        if len(seed_index) > 0:
            # processing the rest that hasn't be recorded into mask_label
            no = int(i / 255)
            mask_label_s1[(mask_label_s1 == 0) & (mask_s1 > 0)] = mask_s1[(mask_label_s1 == 0) & (mask_s1 > 0)].astype(
                int) + no * 255

        object_map_s1 = mask_label_s1[1:self.parameters['n_rows'] + 1, 1:self.parameters['n_cols'] + 1]

        if self.devel_mode:
            gdal_save_file_1band(
                os.path.join(self.obiaresults_path, filenm + '_floodfill_gaussian_{}_s1.tif'.format(bandwidth)),
                object_map_s1, gdal.GDT_Int32, self.trans, self.proj, self.parameters['n_rows'],
                self.parameters['n_cols'])

        del mask_label_s1
        del mask_s1

        #######################################################################################
        #                                 Scale 2: change patch                              #
        #######################################################################################

        # create a object-based change
        unq_s1, ids_s1, count_s1 = np.unique(object_map_s1, return_inverse=True, return_counts=True)
        mean_list = np.bincount(ids_s1.astype(int), weights=cm_array_gaussian_s1.reshape(ids_s1.shape)) / count_s1
        mean_list[unq_s1 == 0] = self.parameters['NAN_VAL']  # force mean of unchanged objects to be -9999
        lut_dict_s1 = dict(zip(unq_s1, mean_list))
        cm_array_s2 = np.vectorize(lut_dict_s1.get)(object_map_s1)

        bandwidth = 1
        # using gaussian kernel ( with 0.5 sigma value) to smooth images in preparation for floodfill
        # cm_array_s2[cm_array_s2 == parameters['NAN_VAL']] = np.nan
        # kernel = Gaussian2DKernel(x_stddev=bandwidth, y_stddev=bandwidth)
        # cm_array_gaussian_s2 = convolve(cm_array_s2, kernel, boundary='extend', preserve_nan=True)
        # cm_array_gaussian_s2[np.isnan(cm_array_s2)] = parameters['NAN_VAL']

        cm_array_gaussian_s2 = cm_array_s2  # not using gaussian
        if self.devel_mode:
            gdal_save_file_1band(join(self.obiaresults_path,  'test_cm_array_gaussian_s2.tif'),
                                 cm_array_gaussian_s2,
                                 gdal.GDT_Float32, self.trans, self.proj, self.parameters['n_cols'],
                                 self.parameters['n_rows'])

        mask_s2 = np.zeros((self.parameters['n_rows'] + 2, self.parameters['n_cols'] + 2)).astype(np.uint8)
        mask_label_s2 = np.zeros((self.parameters['n_rows'] + 2, self.parameters['n_cols'] + 2))
        floodflags_base = 8
        floodflags_base |= cv2.FLOODFILL_MASK_ONLY
        i = 0
        cm_array_gaussian_s2 = cm_array_gaussian_s2.astype(np.float32)
        # cm_stack = np.dstack([cm_array_gaussian_s2, cm_direction_array, cm_direction_array]).astype(np.float32)
        # parameters['floodfill_thres'] = 0.5
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
        if self.devel_mode:
            print("segmentation finished ")

        if len(seed_index) > 0:
            # processing the rest that hasn't be recorded into mask_label
            no = int(i / 255)
            mask_label_s2[(mask_label_s2 == 0) & (mask_s2 > 0)] = mask_s2[(mask_label_s2 == 0) & (mask_s2 > 0)].astype(
                int) + no * 255

        object_map_s2 = mask_label_s2[1:self.parameters['n_rows'] + 1, 1:self.parameters['n_cols'] + 1]
        if self.devel_mode:
            gdal_save_file_1band(
                join(self.obiaresults_path,  'test_floodfill_gaussian_{}_s2.tif'.format(bandwidth)),
                object_map_s2, gdal.GDT_Int32, self.trans, self.proj, self.parameters['n_rows'],
                self.parameters['n_cols'])

        del mask_label_s2
        del mask_s2

        #######################################################################################
        #                                 object-based decision                               #
        #######################################################################################
        if pd.Timestamp.fromordinal(date - 366).year - 1 < self.year_lowbound:  # we used the the year before the date
            # processed as the year to be classified
            classified_year = self.year_lowbound
        else:
            classified_year = pd.Timestamp.fromordinal(date - 366).year - 1
        if self.thematic_mode == 'single_map' or 'rf_model':
            classification_map = gdal_array.LoadFile(
                join(self.out_yearly_path, 'yearlyclassification_{}.npy'.format(classified_year)))
        elif self.thematic_mode == 'multi_maps':
            classification_map = gdal_array.LoadFile(
                join(self.thematic_src, 'yearlyclassification_{}.tif'.format(classified_year)))
        else:
            classification_map = np.full((self.parameters['n_rows'], self.parameters['n_cols']),
                                         self.parameters['NAN_VAL'],
                                         dtype=np.int32)

        change_map = np.zeros((self.parameters['n_rows'], self.parameters['n_cols'])).astype(np.uint8)
        change_group = []

        if self.devel_mode:
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
                if is_change_object(self.parameters, b_thematic, stats_lut_row, uniform_threshold, uniform_sizeslope,
                                    'cm_average'):
                    change_group.append(index)

            if self.devel_mode:
                print("--- %s seconds ---" % (time.time() - start_time))
            change_map[np.isin(object_map_s1, change_group)] = 1

        # if self.devel_mode:
        #     gdal_save_file_1band(join(self.obiaresults_path, filenm + '_OBIAresult.tif'), change_map,
        #                          gdal.GDT_Byte, self.trans, self.proj, self.parameters['n_cols'],
        #                          self.parameters['n_rows'])
        # else:
        #     np.save(join(self.obiaresults_path, filenm + '_OBIAresult.tif'), change_map)

        if self.devel_mode:
            print("change object analysis finished for {}".format(filenm))

        return change_map

    def checkfinished_step3_objectanalysis(self):
        """
        :return:
        """
        for n in range(self.parameters['n_cm_maps']):
            if not exists(join(self.obiaresults_path, filename_fromordinaldate(self.parameters['starting_date'] +
                                                                        n * self.parameters['CM_OUTPUT_INTERVAL']) +
                                               '_OBIAresult.tif')):
                return False
        return True




