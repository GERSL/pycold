from pycold import sccd_detect, sccd_update
from pycold.utils import read_data
import pickle
from collections import namedtuple
SccdOutput = namedtuple("SccdOutput", "position rec_cg min_rmse nrt_mode nrt_model nrt_queue")


def test_sccd_detect():
    # running COLD for a Landsat time series provided by a csv
    in_path = 'tests/resources/spectral_336_3980_obs.csv'
    data = read_data(in_path)
    dates, blues, greens, reds, nirs, swir1s, swir2s, thermals, qas, sensor = data.copy()  # exclude header
    sccd_result = sccd_detect(dates, blues, greens, reds, nirs, swir1s, swir2s, thermals, qas)
    assert sccd_result.rec_cg['t_break'] == 735016
    assert sccd_result.rec_cg['num_obs'] == 141
    assert sccd_result.nrt_model['num_obs'] == 63
    assert sccd_result.nrt_mode == 1


def test_sccd_update():
    pickle_file = open('tests/resources/spectral_15043028_sccdplot.pickle', 'rb')
    sccd_plot = SccdOutput(*pickle.load(pickle_file))
    in_path = 'tests/resources/spectral_15043028_ext.csv'
    data = read_data(in_path)
    dates, blues, greens, reds, nirs, swir1s, swir2s, thermals, qas, sensor = data.copy()  # exclude header
    sccd_plot_new = sccd_update(sccd_plot, dates, blues, greens, reds, nirs, swir1s, swir2s, thermals, qas)
    assert sccd_plot_new.nrt_model['num_obs'] == 64
    assert sccd_plot_new.nrt_model['obs_date_since1982'][0][3] == 9059
