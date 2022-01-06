import yaml
import os
from pycold import object_analyst
import sys

os.chdir('/Users/coloury/Dropbox/Documents/pycold')


def main_test_object_analyst():
    with open('tests/resources/parameters.yaml', 'r') as yaml_obj:
        parameters = yaml.safe_load(yaml_obj)
    logger = logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    object_analyst(obcold_params, 'tests/resources', 'tests/resources', 2015, 2020, logger)
