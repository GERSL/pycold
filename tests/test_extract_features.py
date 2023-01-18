import datetime
import numpy as np
import pathlib
from pycold.imagetool.export_change_map import extract_features

try:
    TEST_RESOURCE_DPATH = (pathlib.Path(__file__).parent / 'resources').resolve()
except NameError:
    TEST_RESOURCE_DPATH = pathlib.Path('tests/resources').resolve()

cold_plot = np.load(TEST_RESOURCE_DPATH / 'record_change_pid15043028_cold.npy')
ordinal_day_list = [datetime.date(year, 7, 1).toordinal() for year
                    in range(1985, 2023, 1)]


def test_extract_features():
    features = extract_features(
        cold_plot, 0, ordinal_day_list, -9999,
        feature_outputs=['cv'])
    assert features[0] == -9999
