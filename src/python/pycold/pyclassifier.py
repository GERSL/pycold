import warnings
warnings.filterwarnings("ignore") # mainly to filter out lzma warnings
import numpy as np
import pandas as pd
import os
import datetime as dt
warnings.filterwarnings("ignore", category=DeprecationWarning)
from sklearn.ensemble import RandomForestClassifier
from os.path import join, exists
from pycold.app import defaults, logging
from pycold.utils import get_block_y, get_block_x, get_col_index, get_row_index
import joblib
import time
from osgeo import gdal_array

logger = logging.getLogger(__name__)


def extract_features(cold_plot, band, ordinal_day_list, nan_val, n_features_perband, ismat=False):
    """
    generate features for classification based on a plot-based rec_cg and a list of days to be predicted
    Parameters
    ----------
    cold_plot: nested array
        plot-based rec_cg
    band: integer
        the predicted band number range from 0 to 6
    ordinal_day_list: list
        a list of days that this function will predict every days as a list as output
    nan_val: integer
        NA value assigned to the output
    n_features_perband: integer
        the number of features per band
    ismat: bool
        True -> the input is MATLAB rec_cg
    Returns
    -------
        feature: a list (length = n_feature) of 1-array [len(ordinal_day_list)]
    """
    features = [np.full(len(ordinal_day_list), nan_val, dtype=np.double) for x in range(n_features_perband)]
    for index, ordinal_day in enumerate(ordinal_day_list):
        # print(index)
        for idx, cold_curve in enumerate(cold_plot):
            if idx == len(cold_plot) - 1:
                max_days = cold_plot[idx]['t_end']
            else:
                max_days = cold_plot[idx + 1]['t_start']
            if cold_curve['t_start'] <= ordinal_day < max_days:
                if ismat:
                    for n in range(n_features_perband):
                        if n == 0:
                            # if cold_curve['t_start'] <= ordinal_day < cold_curve['t_end']:
                            features[n][index] = cold_curve['coefs'][0][band] + cold_curve['coefs'][1][band] * \
                                                 ordinal_day
                        else:
                            features[n][index] = cold_curve['coefs'][n+1][band]
                    break
                else:
                    for n in range(n_features_perband):
                        if n == 0:
                            # if cold_curve['t_start'] <= ordinal_day < cold_curve['t_end']:
                            features[n][index] = cold_curve['coefs'][band][0] + cold_curve['coefs'][band][1] * \
                                                 ordinal_day / defaults['SLOPE_SCALE']
                        else:
                            features[n][index] = cold_curve['coefs'][band][n+1]  # n + 1 is because won't need slope as output
                    break
    return features


def generate_sample_num(label, sample_parameters):
    """
    generate sample number for each land cover category using the method from 'Optimizing selection of training and
    auxiliary data for operational land cover classification for the LCMAP initiative'
    Args:
        label: an array of land-cover categories, i.e., map
        sample_parameters: obcold parameter structure
    Returns:

    """
    unique_category, unique_counts_category = np.unique(label[label <= sample_parameters['total_landcover_category']],
                                                        return_counts=True)
    counts_category = np.array([0] * sample_parameters['total_landcover_category'])
    for x in range(len(unique_counts_category)):
        counts_category[unique_category[x]-1] = unique_counts_category[x]
    percate_samples = np.array([round(x * sample_parameters['total_samples'] / sum(counts_category)) for x in
                                counts_category])
    percate_samples[percate_samples > sample_parameters['max_category_samples']] = \
        sample_parameters['max_category_samples']
    percate_samples[percate_samples < sample_parameters['min_category_samples']] = \
        sample_parameters['min_category_samples']
    percate_samples = np.minimum(percate_samples, counts_category)  # needs to check not exceed the total category pixels
    return percate_samples


class PyClassifier:
    def __init__(self, config, n_features=None):
        """
        Parameters
        ----------
        config: from config.yaml
        n_features: the total feature number for a single pixel, it may be overwritten when accessory bands are added
        """
        self.config = parameters
        self.config['block_width'] = int(self.config['n_cols'] / self.config['n_block_x'])
        self.config['block_height'] = int(self.config['n_rows'] / self.config['n_block_y'])
        self.config['n_blocks'] = self.config['n_block_x'] * self.config['n_block_y']
        if n_features is None:
            self.n_features = defaults['TOTAL_IMAGE_BANDS'] * defaults['N_FEATURES']
        else:
            self.n_features = n_features

    def predict_features(self, block_id, cold_block, year_lowbound, year_uppbound, ismat=False):
        """
        Parameters
        ----------
        block_id: integer
            the block id
        cold_block: nested array of colddt datatype
            the block-based change records produced by cold algorithms
        year_lowbound: integer
            the lower bound of the analysis year range
        year_uppbound: integer
            the upper bound of the analysis year range
            Note that the reason for not parsing cold_block to get year bounds is that the year ranges of blocks
            may vary from each other, so the year bounds are required to be defined from the tile level, not block level
            such as from 'starting_end_date.txt'
        ismat:boolean
            True -> the input is MATLAB-outputted change records

        Returns
        -------
        an array [year_uppbound-year_lowbound+1, block_width*block_height, n_features]
        """
        block_features = np.full(((year_uppbound - year_lowbound + 1),
                                  self.config['block_width'] * self.config['block_height'],
                                  self.n_features),
                                 defaults['NAN_VAL'], dtype=np.float32)
        ordinal_day_list = [pd.Timestamp.toordinal(dt.date(year, 7, 1)) + 366 for year
                            in range(year_lowbound, year_uppbound + 1)]
        if len(cold_block) == 0:
            logger.warning('the rec_cg file for block_id'.format(block_id))
            return block_features

        cold_block_split = np.split(cold_block, np.argwhere(np.diff(cold_block['pos']) != 0)[:, 0] + 1)
        for element in cold_block_split:
            # the relative column number in the block
            i_col = get_col_index(element[0]["pos"], self.config['n_cols'],
                                  get_block_x(block_id, self.config['n_block_x']),
                                  self.config['block_width'])
            i_row = get_row_index(element[0]["pos"], self.config['n_cols'],
                                  get_block_y(block_id, self.config['n_block_x']),
                                  self.config['block_height'])

            for band in range(defaults['TOTAL_IMAGE_BANDS']):
                feature_row = extract_features(element, band, ordinal_day_list, defaults['NAN_VAL'],
                                               defaults['N_FEATURES'], ismat)  # N_FEATURE *[year_number]
                for index in range(defaults['N_FEATURES']):
                    block_features[:, i_row * self.config['block_width'] + i_col, band * defaults['N_FEATURES'] +
                                                                                    index] = feature_row[index]

        return block_features

    def assemble_array(self, array_list):
        """
        assemble a list of block-based array to a bigger array that aligns with the dimension of an ARD tile
        Parameters
        ----------
        array_list

        Returns
        -------
        an array [nrows, ncols]
        """
        full_feature_array = np.hstack(array_list)
        full_feature_array = np.vstack(np.hsplit(full_feature_array, self.config['n_block_x']))  # (nrows, ncols, nfeatures)
        return full_feature_array

    def train_rfmodel(self, full_feature_array, label):
        """
        Parameters
        ----------
        full_feature_array: (nrows, ncols, feature_number) array
            feature array for the whole tile
        label: (nrows, ncols) array
            reference array for the whole tile
        Returns
        -------
            a sklearn random forest model
        """
        assert label.shape == (self.config['n_rows'], self.config['n_cols'])
        samplecount = generate_sample_num(label, defaults)
        index_list = []
        label_list = []
        for i in range(defaults['total_landcover_category']):
            index = np.argwhere(label == i + 1)
            np.random.seed(42)  # set random seed to reproduce the same result
            index_sample = index[np.random.choice(len(index), int(samplecount[i]), replace=False)]
            index_list.append(index_sample)
            label_list.append(np.array([i + 1] * len(index_sample)))
        index_list = np.vstack(index_list)
        label_list = np.hstack(label_list)
        feature_extraction = np.array([full_feature_array[tuple(x)] for x in index_list])
        rf_model = RandomForestClassifier(random_state=42)
        rf_model.fit(feature_extraction, label_list)
        return rf_model

    def classification_block(self, rf_model, tmp_feature):
        """
        make classification for temp feature array for a single year, single block
        Parameters
        ----------
        rf_model: sklearn random forest model
        tmp_feature:
        Return: (block_height, block_width) array
        -------
            cmap, the feature temp array file for block_id and year
        """
        # image_flat = np.zeros((parameters['n_rows'] * parameters['n_cols'], parameters['TOTAL_IMAGE_BANDS']))
        cmap = rf_model.predict(tmp_feature).reshape(self.config['block_height'], self.config['block_width'])
        mask = np.all(tmp_feature == defaults['NAN_VAL'], axis=1).reshape(self.config['block_height'],
                                                                          self.config['block_width'])
        cmap[mask] = 255
        return cmap


class PyClassifierHPC(PyClassifier):
    """
    this class adds new IO functions in the HPC environment for the base class
    """
    def __init__(self, config, record_path, year_lowbound=1982, year_uppbound=2021, tmp_path=None, output_path=None,
                 n_features=None, labelmap_path=None, rf_path=None):
        """
        Parameters
        ----------
        config: configuration structure from config.yaml
        record_path
        year_lowbound
        year_uppbound
        labelmap_path
        """
        try:
            self.__check_inputs_thematic(config, record_path, year_lowbound, year_uppbound, tmp_path,
                                         labelmap_path, rf_path)
        except ValueError or FileExistsError as e:
            raise e

        self.config = config
        self.config['block_width'] = int(self.config['n_cols'] / self.config['n_block_x'])
        self.config['block_height'] = int(self.config['n_rows'] / self.config['n_block_y'])
        self.config['n_blocks'] = self.config['n_block_x'] * self.config['n_block_y']
        self.record_path = record_path

        if tmp_path is None:
            self.tmp_path = join(record_path, 'feature_maps')  # default path
        else:
            self.tmp_path = tmp_path

        if output_path is None:
            self.output_path = join(record_path, 'feature_maps')  # default path
        else:
            self.output_path = tmp_path

        if n_features is None:
            self.n_features = defaults['TOTAL_IMAGE_BANDS'] * defaults['N_FEATURES']
        else:
            self.n_features = n_features

        self.year_lowbound = year_lowbound
        self.year_uppbound = year_uppbound
        self.labelmap_path = labelmap_path
        if rf_path is None:
            self.rf_path = join(self.output_path, 'rf.model')  # default path
        else:
            self.rf_path = rf_path

    @staticmethod
    def __check_inputs_thematic(config, record_path, year_lowbound, year_uppbound, tmp_path,  labelmap_path,
                                rf_path):
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

        if os.path.isdir(record_path) is False:
            raise FileExistsError('No such directory: {}'.format(record_path))

        if labelmap_path is not None:
            if os.path.isfile(labelmap_path) is False:
                raise FileExistsError('No such file: {}'.format(labelmap_path))

        if rf_path is not None:
            if os.path.isfile(rf_path) is False:
                raise FileExistsError('No such file: {}'.format(rf_path))

    def hpc_preparation(self):
        if exists(self.tmp_path) is False:
            try:
                os.mkdir(self.tmp_path)
            except IOError as e:
                raise e

        if exists(self.output_path) is False:
            try:
                os.mkdir(self.output_path)
            except IOError as e:
                raise e

    def save_features(self, block_id, block_features):
        '''
        Parameters
        ----------
        block_id: integer
        block_features
        Returns
        -------

        '''
        for i in range(block_features.shape[0]):
            np.save(os.path.join(self.tmp_path, 'tmp_feature_year{}_block{}.npy').format(self.year_lowbound+i,
                                                                                                  block_id),
                    block_features[i, :, :])

    def is_finished_step1_predict_features(self):
        return any(exists(join(self.tmp_path, 'tmp_feature_year{}_block{}.npy').format(year, block + 1))
                   for block in range(self.config['n_blocks']) for year in range(self.year_lowbound,
                                                                                 self.year_uppbound+1))

    def get_fullfeature_forcertainyear(self, year):
        tmp_feature_filenames = [file for file in os.listdir(self.tmp_path)
                                 if file.startswith('tmp_feature_year{}'.format(year))]
        if len(tmp_feature_filenames) < self.config['n_blocks']:
            logger.warning('tmp features are incomplete! should have {}; but actually have {} feature images'.
                           format(self.config['n_blocks'], len(tmp_feature_filenames)))

        tmp_feature_filenames.sort(
            key=lambda t: t[t.find('block'): t.find('.npy')])  # sorted by row number
        full_feature_array = self.assemble_array([np.load(join(self.tmp_path, file))
                                                 .reshape(self.config['block_height'],
                                                          self.config['block_width'],
                                                          self.n_features)
                                                  for file in tmp_feature_filenames])
        if (full_feature_array.shape[1] != self.config['n_cols']) or (full_feature_array.shape[0] !=
                                                                              self.config['n_rows']):
            logger.error('The feature image is incomplete for {}'.format(year))
        return full_feature_array

    @staticmethod
    def save_rf_model(rf_model, rf_path):
        joblib.dump(rf_model, rf_path)

    def is_finished_step2_train_rfmodel(self):
        return exists(self.rf_path)

    @staticmethod
    def get_features(path):
        """
        Parameters
        ----------
        path:
        Returns: (block_width*block_height, total feature) array
        -------
            the feature temp array file for block_id and year
        """
        return np.load(path)

    def get_rf_model(self):
        return joblib.load(self.rf_path)

    def save_yearlyclassification_maps(self, block_id, year, cmap):
        outfile = join(self.tmp_path, 'tmp_yearlyclassification{}_block{}.npy'.format(year, block_id))
        np.save(outfile, cmap)

    def is_finished_step3_classification(self):
        """
        :return: True or false
        """
        return any(exists(join(self.tmp_path, 'tmp_yearlyclassification{}_block{}.npy'.format(year, iblock + 1)))
                   for iblock in range(self.config['n_blocks']) for year in range(self.year_lowbound,
                                                                                      self.year_uppbound + 1))

    def assemble_covermap(self, year):
        tmp_yearlyclass_filenames = [file for file in os.listdir(self.tmp_path)
                                     if file.startswith('tmp_yearlyclassification{}'.format(year))]

        # sort to guarantee order follows low to high rows
        tmp_yearlyclass_filenames.sort(key=lambda t: t[t.find('block'): t.find('.npy')])
        full_yearlyclass_array = self.assemble_array([np.load(join(self.tmp_path, file))
                                                     .reshape(self.config['block_height'],
                                                              self.config['block_width'], 1)
                                                      for file in tmp_yearlyclass_filenames])

        if (full_yearlyclass_array.shape[1] != self.config['n_cols']) or (full_yearlyclass_array.shape[0] !=
                                                                              self.config['n_rows']):
            logger.error('The assemble category map is incomplete for {}'.format(year))
        return full_yearlyclass_array[:, :, 0]  # we only return the first channel

    def save_covermaps(self, full_yearlyclass_array, year):
        np.save(join(self.output_path, 'yearlyclassification_{}.npy'.format(year)), full_yearlyclass_array)

    def clean(self):
        tmp_yearlyclass_filenames = [file for file in os.listdir(self.tmp_path)
                                     if file.startswith('tmp_')]
        for file in tmp_yearlyclass_filenames:
            os.remove(join(self.tmp_path, file))

    def step1_feature_generation(self, block_id):
        cold_block = np.load(join(self.record_path, 'record_change_x{}_y{}_cold.npy').format(
                             get_block_x(block_id, self.config['n_block_x']),
                             get_block_y(block_id, self.config['n_block_x'])))
        block_features = self.predict_features(block_id, cold_block, self.year_lowbound, self.year_uppbound)
        self.save_features(block_id, block_features)

    def step2_train_rf(self, ref_year=None):
        if ref_year is None:
            ref_year = defaults['classification_year']
        while not self.is_finished_step1_predict_features():
            time.sleep(5)
        full_feature_array = self.get_fullfeature_forcertainyear(ref_year)
        rf_model = self.train_rfmodel(full_feature_array, gdal_array.LoadFile(self.labelmap_path))
        self.save_rf_model(rf_model, self.rf_path)

    def step3_classification(self, block_id):
        while not self.is_finished_step2_train_rfmodel():
            time.sleep(5)

        try:
            rf_model = self.get_rf_model()
        except IOError as e:
            raise ("Please double check your rf model file directory or generate random forest model first:"
                   " {}".format(e))

        for year in range(self.year_lowbound, self.year_uppbound + 1):
            tmp_feature_block = self.get_features(join(self.tmp_path, 'tmp_feature_year{}_block{}.npy'.format(year,
                                                                                                        block_id)))
            cmap = self.classification_block(rf_model, tmp_feature_block)
            self.save_yearlyclassification_maps(block_id, year, cmap)

    def step4_assemble(self):
        while not self.is_finished_step3_classification():
            time.sleep(5)
        for year in range(self.year_lowbound, self.year_uppbound + 1):
            full_yearlyclass_array = self.assemble_covermap(year)
            self.save_covermaps(full_yearlyclass_array, year)
        self.clean()  # clean all temp files






