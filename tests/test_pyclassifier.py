import yaml
import os
from pycold.pyclassifier import PyClassifierHPC
import numpy as np
from osgeo import gdal_array
import shutil

with open('tests/resources/test_config.yaml', 'r') as yaml_obj:
    test_config = yaml.safe_load(yaml_obj)

ref_path = 'tests/resources/LCMAP_CU_2001_V01_LCPRI_027009_partial.tif'
year_lowbound = 2015
year_uppbound = 2020
testing_folder = 'tests/resources/tmp'

pyclassifier = PyClassifierHPC(test_config, 'tests/resources', 'tests/resources', year_lowbound=year_lowbound,
                               year_uppbound=year_uppbound,
                               thematic_src=ref_path)


def test_workflow():
    # create a test folder
    if not os.path.exists(testing_folder):
        os.mkdir(testing_folder)
    else:
        shutil.rmtree(testing_folder)
        os.mkdir(testing_folder)

    for block_id in range(pyclassifier.config['n_blocks']):
        pyclassifier.step1_feature_generation(block_id, in_path='tests/resources', out_path=testing_folder)
        for year in range(year_lowbound, year_uppbound + 1):
            assert os.path.exists(os.path.join(testing_folder, 'tmp_feature_year{}_block{}.npy'.format(year, block_id)))

    pyclassifier.step2_train_rf(2017, in_path=testing_folder, out_path=testing_folder)
    assert os.path.exists(os.path.join(testing_folder, 'rf.model'))

    for block_id in range(pyclassifier.config['n_blocks']):
        pyclassifier.step3_classification(block_id, rf_path=os.path.join(testing_folder, 'rf.model'),
                                          in_path=testing_folder,
                                          out_path=testing_folder)
        for year in range(year_lowbound, year_uppbound + 1):
            assert os.path.exists(os.path.join(testing_folder, 'tmp_yearlyclassification{}_block{}.npy'.format(year,
                                                                                                               block_id)))
    pyclassifier.step4_assemble(in_path=testing_folder, out_path=testing_folder)
    for year in range(year_lowbound, year_uppbound + 1):
        assert os.path.exists(os.path.join(testing_folder, 'yearlyclassification_{}.npy'.format(year)))

    # import matplotlib.pyplot as plt
    # plt.imshow(np.load(os.path.join(os.getcwd(),'tests/resources/tmp/yearlyclassification_2017.npy')))
    shutil.rmtree(testing_folder)


def test_predict_features():
    cold_block = np.load(os.path.join('tests/resources', 'record_change_x1_y1_cold.npy'))
    block_features = pyclassifier.predict_features(1, cold_block, 2015, 2020)
    assert np.shape(block_features) == (6, pyclassifier.config['block_width']*pyclassifier.config['block_height'],
                                        21)   # 6 years, 2500 pixels (50, 50), 21 features


def test_produce_rf_model():
    full_feature_2017 = pyclassifier.get_fullfeature_forcertainyear('tests/resources', 2017)
    assert full_feature_2017.shape == (pyclassifier.config['n_rows'],
                                        pyclassifier.config['n_cols'], 21)
    rf_model = pyclassifier.train_rfmodel(full_feature_2017, gdal_array.LoadFile(ref_path))
    # pyclassifier.save_rf_model(rf_model)
    assert rf_model is not None


def test_classification():
    rf_model = pyclassifier.get_rf_model(os.path.join('tests/resources', 'rf.model'))
    tmp_feature_block = pyclassifier.get_features(2017, 1, in_path='tests/resources')
    cmap = pyclassifier.classification_block(rf_model, tmp_feature_block)
    print(cmap.shape)
    assert cmap.shape == (50,50)

