from ._colds_cython import _cold_detect, _obcold_reconstruct, _sccd_update, SccdOutput, SccdWrapper
import numpy as np
from ._param_validation import validate_parameter_constraints, Interval, Integral, Real, check_consistent_length, \
    check_1d
import copy


_parameter_constraints: dict = {
    "t_cg": [Interval(Real, 0.0, None, closed="neither")],
    "pos": [Interval(Integral, 0, None, closed="neither")],
    "b_output_cm": ["boolean"],
    "starting_date": [Interval(Integral, 0, None, closed="left")],
    "n_cm": [Interval(Integral, 0, None, closed="left")],
    "conse": [Interval(Integral, 0, 12, closed="neither")],
    "b_output_cm": ["boolean"],
    "gap_days": [Interval(Real, 0.0, None, closed="left")],
    "b_pinpoint": ["boolean"],
    "gate_tcg": [Interval(Real, 0.0, None, closed="neither")],
    "b_monitor_init": ["boolean"]
}


NUM_FC = 40  # define the maximum number of outputted curves
NUM_FC_SCCD = 40
NUM_NRT_QUEUE = 240


def _validate_params(func_name, **kwargs):
    """Validate types and values of constructor parameters
    The expected type and values must be defined in the `_parameter_constraints`
    class attribute, which is a dictionary `param_name: list of constraints`. See
    the docstring of `validate_parameter_constraints` for a description of the
    accepted constraints.
    """
    # params = dict()
    # for key in kwargs:
    #     dict[key] = kwargs[key]

    validate_parameter_constraints(
        _parameter_constraints,
        kwargs,
        caller_name=func_name,
    )


def _validate_data(dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas, break_dates=None):
    """
    validate and forcibly change the data format
    Parameters
    ----------
    dates: 1d array of shape(observation numbers), list of ordinal dates
    ts_b: 1d array of shape(observation numbers), time series of blue band.
    ts_g: 1d array of shape(observation numbers), time series of green band
    ts_r: 1d array of shape(observation numbers), time series of red band
    ts_n: 1d array of shape(observation numbers), time series of nir band
    ts_s1: 1d array of shape(observation numbers), time series of swir1 band
    ts_s2: 1d array of shape(observation numbers), time series of swir2 band
    ts_t: 1d array of shape(observation numbers), time series of thermal band
    qas: 1d array, the QA cfmask bands. '0' - clear; '1' - water; '2' - shadow; '3' - snow; '4' - cloud

    Returns
    ----------
    """
    check_consistent_length(dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas)
    check_1d(dates, 'dates')
    check_1d(ts_b, 'ts_b')
    check_1d(ts_g, 'ts_g')
    check_1d(ts_r, 'ts_r')
    check_1d(ts_n, 'ts_n')
    check_1d(ts_s1, 'ts_s1')
    check_1d(ts_s2, 'ts_s2')
    check_1d(qas, 'qas')
    if break_dates is not None:
        check_1d(break_dates, 'break_dates')

    dates = dates.astype(dtype=np.int64, order="C")
    ts_b = ts_b.astype(dtype=np.int64, order="C")
    ts_g = ts_g.astype(dtype=np.int64, order="C")
    ts_r = ts_r.astype(dtype=np.int64, order="C")
    ts_n = ts_n.astype(dtype=np.int64, order="C")
    ts_s1 = ts_s1.astype(dtype=np.int64, order="C")
    ts_s2 = ts_s2.astype(dtype=np.int64, order="C")
    ts_t = ts_t.astype(dtype=np.int64, order="C")
    qas = qas.astype(dtype=np.int64, order="C")
    if break_dates is not None:
        break_dates = break_dates.astype(dtype=np.int64, order="C")
        return dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas, break_dates
    else:
        return dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas


def cold_detect(dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas, t_cg=15.0863, pos=1, conse=6,
                b_output_cm=False, starting_date=0, n_cm=0, cm_output_interval=0, b_c2=False, gap_days=365.25):
    """
    pixel-based COLD algorithm.
    Zhu, Z., Zhang, J., Yang, Z., Aljaddani, A. H., Cohen, W. B., Qiu, S., &
    Zhou, C. (2020). Continuous monitoring of land disturbance based on Landsat time series.
    Remote Sensing of Environment, 38, 111116.
    Parameters
    ----------
    dates: 1d array of shape(observation numbers), list of ordinal dates
    ts_b: 1d array of shape(observation numbers), time series of blue band.
    ts_g: 1d array of shape(observation numbers), time series of green band
    ts_r: 1d array of shape(observation numbers), time series of red band
    ts_n: 1d array of shape(observation numbers), time series of nir band
    ts_s1: 1d array of shape(observation numbers), time series of swir1 band
    ts_s2: 1d array of shape(observation numbers), time series of swir2 band
    ts_t: 1d array of shape(observation numbers), time series of thermal band
    qas: 1d array, the QA cfmask bands. '0' - clear; '1' - water; '2' - shadow; '3' - snow; '4' - cloud
    t_cg: threshold of change magnitude, default is chi2.ppf(0.99,5)
    pos: position id of the pixel
    conse: consecutive observation number
    b_output_cm: bool, 'True' means outputting change magnitude and change magnitude dates, only for object-based COLD
    starting_date: the starting date of the whole dataset to enable reconstruct CM_date, all pixels for a tile
                    should have the same date, only for b_output_cm is True. Only b_output_cm == 'True'
    n_cm: the length of outputted change magnitude. Only b_output_cm == 'True'
    cm_output_interval: the temporal interval of outputting change magnitudes. Only b_output_cm == 'True'
    b_c2: bool, a temporal parameter to indicate if collection 2. C2 needs ignoring thermal band for valid pixel
          test due to the current low quality
    gap_days: define the day number of the gap year for determining i_dense. Setting a large value (e.g., 1500)
                if the gap year in the middle of the time range

    Returns
    ----------
    change records: the COLD outputs that characterizes each temporal segment if b_output_cm==False
    Or
    [change records, cm_outputs, cm_outputs_date] if b_output_cm==True
    """
    _validate_params(func_name='cold_detect', t_cg=t_cg, pos=pos, conse=conse, b_output_cm=b_output_cm,
                     starting_date=starting_date, n_cm=n_cm, cm_output_interval=cm_output_interval, b_c2=b_c2,
                     gap_days=gap_days)
    # make sure it is c contiguous array and 64 bit
    dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas = _validate_data(dates, ts_b, ts_g, ts_r, ts_n, ts_s1,
                                                                            ts_s2, ts_t, qas)

    return _cold_detect(dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas, t_cg, pos, conse, b_output_cm,
                        starting_date, n_cm, cm_output_interval, b_c2, gap_days)



def obcold_reconstruct(dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas, break_dates, pos=1, conse=6, b_c2=False):
    """
    re-contructructing change records using break dates.
    Parameters
    ----------
    dates: 1d array of shape(observation numbers), list of ordinal dates
    ts_b: 1d array of shape(observation numbers), time series of blue band.
    ts_g: 1d array of shape(observation numbers), time series of green band
    ts_r: 1d array of shape(observation numbers), time series of red band
    ts_n: 1d array of shape(observation numbers), time series of nir band
    ts_s1: 1d array of shape(observation numbers), time series of swir1 band
    ts_s2: 1d array of shape(observation numbers), time series of swir2 band
    ts_t: 1d array of shape(observation numbers), time series of thermal band
    qas: 1d array, the QA cfmask bands. '0' - clear; '1' - water; '2' - shadow; '3' - snow; '4' - cloud
    break_dates: 1d array, the break dates obtained from other procedures such as obia
    conse: consecutive observation number (for calculating change magnitudes)
    b_c2: bool, a temporal parameter to indicate if collection 2. C2 needs ignoring thermal band for valid pixel test
          due to its current low quality

    Returns
    ----------
    change records: the COLD outputs that characterizes each temporal segment
    """
    _validate_params(func_name='sccd_detect', pos=pos, conse=conse, b_c2=b_c2)
    dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas, break_dates = _validate_data(dates, ts_b, ts_g, ts_r, ts_n,
                                                                                         ts_s1, ts_s2, ts_t, qas,
                                                                                         break_dates)

    return _obcold_reconstruct(dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas, break_dates, pos, conse, b_c2)


def sccd_detect(dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas, t_cg=15.0863, pos=1, conse=6, b_c2=False,
                b_pinpoint=False, gate_tcg=9.236, b_monitor_init=False):
    """
    pixel-based offline SCCD algorithm.
    Ye, S., Rogan, J., Zhu, Z., & Eastman, J. R. (2021). A near-real-time approach for monitoring forest
    disturbance using Landsat time series: Stochastic continuous change detection.
    Remote Sensing of Environment, 252, 112167.

    Parameters
    ----------
    dates: 1d array of shape(observation numbers), list of ordinal dates
    ts_b: 1d array of shape(observation numbers), time series of blue band.
    ts_g: 1d array of shape(observation numbers), time series of green band
    ts_r: 1d array of shape(observation numbers), time series of red band
    ts_n: 1d array of shape(observation numbers), time series of nir band
    ts_s1: 1d array of shape(observation numbers), time series of swir1 band
    ts_s2: 1d array of shape(observation numbers), time series of swir2 band
    ts_t: 1d array of shape(observation numbers), time series of thermal band
    qas: 1d array, the QA cfmask bands. '0' - clear; '1' - water; '2' - shadow; '3' - snow; '4' - cloud
    t_cg: threshold of change magnitude, default is chi2.ppf(0.99,5)
    pos: position id of the pixel, i.e.,  (row -1) * ncols + col, row and col starts from 1
    conse: consecutive observation number
    b_c2: bool, a temporal parameter to indicate if collection 2. C2 needs ignoring thermal band for valid pixel
                test due to its current low quality
    b_pinpoint: bool, output pinpoint break where pinpoint is an overdetection of break using conse =3
                        and threshold = gate_tcg, which are used to simulate the situation of NRT scenario and
                        for training a machine-learning model
    gate_tcg: the gate change magnitude threshold for defining anomaly
    b_monitor_init: bool, output change metrics during the initialization stage
    Note that passing 2-d array to c as 2-d pointer does not work, so have to pass separate bands
    Returns
    ----------
    SccdOutput: the SCCD outputs that characterizes each temporal segment if b_pinpoint==False
    Or
    [SccdOutput, SccdReccgPinpoint]: b_pinpoint==True
    """
    _validate_params(func_name='sccd_detect', t_cg=t_cg, pos=pos, conse=conse, b_c2=b_c2, b_pinpoint=b_pinpoint,
                     gate_tcg=gate_tcg, b_monitor_init=b_monitor_init)
    # make sure it is c contiguous array and 64 bit
    dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas = _validate_data(dates, ts_b, ts_g, ts_r, ts_n, ts_s1,
                                                                            ts_s2, ts_t, qas)

    sccd_wrapper = SccdWrapper()
    tmp = copy.deepcopy(sccd_wrapper.sccd_detect(dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas, t_cg,
                                    pos, conse, b_c2, b_pinpoint, gate_tcg, b_monitor_init))
    return tmp


def sccd_update(sccd_pack, dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas, t_cg=15.0863, pos=1, conse=6,
                b_c2=False, b_pinpoint=False, gate_tcg=9.236, b_monitor_init=False):
    """
    SCCD online update for new observations
    Ye, S., Rogan, J., Zhu, Z., & Eastman, J. R. (2021). A near-real-time approach for monitoring forest
    disturbance using Landsat time series: Stochastic continuous change detection.
    Remote Sensing of Environment, 252, 112167.

    Parameters
    ----------
    sccd_pack:
    dates: 1d array of shape(observation numbers), list of ordinal dates
    ts_b: 1d array of shape(observation numbers), time series of blue band.
    ts_g: 1d array of shape(observation numbers), time series of green band
    ts_r: 1d array of shape(observation numbers), time series of red band
    ts_n: 1d array of shape(observation numbers), time series of nir band
    ts_s1: 1d array of shape(observation numbers), time series of swir1 band
    ts_s2: 1d array of shape(observation numbers), time series of swir2 band
    ts_t: 1d array of shape(observation numbers), time series of thermal band
    qas: 1d array, the QA cfmask bands. '0' - clear; '1' - water; '2' - shadow; '3' - snow; '4' - cloud
    t_cg: threshold of change magnitude, default is chi2.ppf(0.99,5)
    pos: position id of the pixel, i.e.,  (row -1) * ncols + col, row and col starts from 1
    conse: consecutive observation number
    b_c2: bool, a temporal parameter to indicate if collection 2. C2 needs ignoring thermal band for valid pixel
                test due to its current low quality
    b_pinpoint: bool, output pinpoint break where pinpoint is an overdetection of break using conse =3
                        and threshold = gate_tcg, which are used to simulate the situation of NRT scenario and
                        for training a machine-learning model
    gate_tcg: the gate change magnitude threshold for defining anomaly
    b_monitor_init: bool, output change metrics during the initialization stage
    Note that passing 2-d array to c as 2-d pointer does not work, so have to pass separate bands
    Returns
    ----------
    change records: the SCCD outputs that characterizes each temporal segment
    """
    if not isinstance(sccd_pack, SccdOutput):
        raise ValueError("The type of sccd_pack has to be namedtuple 'SccdOutput'!")

    _validate_params(func_name='sccd_update', t_cg=t_cg, pos=pos, conse=conse, b_c2=b_c2, b_pinpoint=b_pinpoint,
                     gate_tcg=gate_tcg, b_monitor_init=b_monitor_init)

    dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas = _validate_data(dates, ts_b, ts_g, ts_r, ts_n, ts_s1,
                                                                            ts_s2, ts_t, qas)

    return _sccd_update(sccd_pack, dates, ts_b, ts_g, ts_r, ts_n, ts_s1, ts_s2, ts_t, qas, t_cg, pos, conse, b_c2,
                        gate_tcg, b_monitor_init)

