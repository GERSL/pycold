import yaml
import os
from pycold.pyclassifier import PyClassifierHPC
import shutil
import numpy as np
import joblib
from osgeo import gdal_array
from pycold.utils import assemble_array

with open('tests/resources/test_config_pyclassifier.yaml', 'r') as yaml_obj:
    test_config = yaml.safe_load(yaml_obj)

ref_path = 'tests/resources/LCMAP_CU_2001_V01_LCPRI_027009_100_100.tif'
year_lowbound = 2015
year_uppbound = 2020
testing_folder = 'tests/resources/tmp'


def test_step1_4():
    pyclassifier = PyClassifierHPC(test_config, record_path='tests/resources', tmp_path=testing_folder,
                                   output_path=testing_folder,
                                   year_lowbound=year_lowbound,
                                   year_uppbound=year_uppbound,
                                   seedmap_path=ref_path)

    # delete testing folder if it has
    if os.path.exists(pyclassifier.tmp_path):
        shutil.rmtree(pyclassifier.tmp_path)

    pyclassifier.hpc_preparation()

    for iblock in range(pyclassifier.config['n_blocks']):
        pyclassifier.step1_feature_generation(iblock + 1)
        for year in range(year_lowbound, year_uppbound + 1):
            assert os.path.exists(os.path.join(pyclassifier.tmp_path,
                                               'tmp_feature_year{}_block{}.npy'.format(year, iblock + 1)))

    pyclassifier.step2_train_rf(2017)
    assert os.path.exists(os.path.join(pyclassifier.output_path, 'rf.model'))

    for iblock in range(pyclassifier.config['n_blocks']):
        pyclassifier.step3_classification(iblock + 1)
        for year in range(year_lowbound, year_uppbound + 1):
            assert os.path.exists(os.path.join(pyclassifier.tmp_path,
                                               'tmp_yearlyclassification{}_block{}.npy'.format(year, iblock + 1)))
    pyclassifier.step4_assemble()
    for year in range(year_lowbound, year_uppbound + 1):
        assert os.path.exists(os.path.join(pyclassifier.output_path, 'yearlyclassification_{}.npy'.format(year)))

    # import matplotlib.pyplot as plt
    # plt.imshow(np.load(os.path.join(os.getcwd(),'tests/resources/tmp/yearlyclassification_2017.npy')))
    shutil.rmtree(pyclassifier.tmp_path)


def test_predict_features():
    pyclassifier = PyClassifierHPC(test_config, record_path='tests/resources', year_lowbound=year_lowbound,
                                     year_uppbound=year_uppbound,
                                     seedmap_path=ref_path)
    cold_block = np.load(os.path.join('tests/resources', 'record_change_x1_y1_cold.npy'))
    block_features = pyclassifier.predict_features(1, cold_block, 2015, 2020)
    assert np.shape(block_features) == (6, pyclassifier.config['block_width']*pyclassifier.config['block_height'],
                                        21)   # 6 years, 2500 pixels (50, 50), 21 features


def test_train_rf_model():
    pyclassifier = PyClassifierHPC(test_config, record_path='tests/resources', tmp_path='tests/resources/feature_maps',
                                     year_lowbound=year_lowbound,
                                     year_uppbound=year_uppbound,
                                     seedmap_path=ref_path)

    full_feature_2017 = assemble_array(pyclassifier.get_fullfeature_forcertainyear(2017),
                                       pyclassifier.config['n_block_x'])
    rf_model = pyclassifier.train_rfmodel(full_feature_2017, gdal_array.LoadFile(ref_path))
    # pyclassifier.save_rf_model(rf_model)
    assert rf_model is not None


def test_classification_block():
    pyclassifier = PyClassifierHPC(test_config, record_path='tests/resources', tmp_path='tests/resources/feature_maps',
                                   rf_path='tests/resources/feature_maps/rf.model',
                                   year_lowbound=year_lowbound,
                                   year_uppbound=year_uppbound)
    rf_model = joblib.load('tests/resources/feature_maps/rf.model')
    tmp_feature_block = np.load('tests/resources/feature_maps/tmp_feature_year2017_block1.npy')
    cmap = pyclassifier.classification_block(rf_model, tmp_feature_block)
    print(cmap.shape)
    assert cmap.shape == (50, 50)

