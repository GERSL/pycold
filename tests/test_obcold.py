import yaml
import os
from pycold.obcold import get_block_y
import sys
from pycold.obcold import ObjectAnalyst

os.chdir('/Users/coloury/Dropbox/Documents/pycold')


def test_obcold():
    with open('tests/resources/test_parameters.yaml', 'r') as yaml_obj:
        test_parameter = yaml.safe_load(yaml_obj)
    object_analyst = ObjectAnalyst(test_parameter, 'tests/resources', 'tests/resources', year_lowbound=2015,
                                   year_uppbound=2020)
    feature_list = object_analyst.step1_predict_features(1)
    assert feature_list is not None

