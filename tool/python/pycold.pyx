from libc.stdlib cimport malloc
from libc.string cimport strcpy, strlen
import numpy as np
import math
from libcpp cimport bool
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free


cdef int NUM_FC = 40  # define the maximum number of outputted curves
cdef int SPATIAL_OUTPUT_INTERVAL = 30  # for object-based COLD only

reccg_dt = np.dtype([('t_start', np.int32),  # time when series model gets started
                     ('t_end', np.int32),  # time when series model gets ended
                     ('t_break', np.int32),  # time when the first break (change) is observed
                     ('pos', np.int32),  # the location of each time series model
                     ('num_obs', np.int32),  # the number of "good" observations used for model estimation
                     ('category', np.short),  # the quality of the model estimation (what model is used, what process is used)
                     ('change_prob', np.short),  # the probability of a pixel that have undergone change (between 0 and 100)
                     ('coefs', np.double, (7, 8)),  # coefficients for each time series model for each spectral band
                     ('rmse', np.double, 7),  # RMSE for each time series model for each spectral band
                     ('magnitude', np.double, 7)])  # the magnitude of change difference between model prediction
                                          # and observation for each spectral band)


cdef extern from "../../src/output.h":
    ctypedef struct Output_t:
        int t_start
        int t_end
        int t_break
        int pos
        int num_obs
        short int category
        short int change_prob
        double coefs[7][8]
        double rmse[7]
        double magnitude[7]

cdef extern from "../../src/cold.h":
    cdef int cold(short *buf_b, short *buf_g, short *buf_r, short *buf_n, short *buf_s1, short *buf_s2,
                  short *buf_t, short *fmask_buf, int *valid_date_array, int valid_num_scenes, int num_samples,
                  int col_pos, int row_pos, double tcg, int conse, bool b_outputCM, int starting_date, Output_t *rec_cg,
                  int *num_fc, short int *CM_outputs, char *CM_outputs_date);



def pycold(int[:] dates, short[:] ts_b, short[:] ts_g, short[:] ts_r, short[:] ts_n,
           short[:] ts_s1, short[:] ts_s2, short[:] ts_t, short[:] qas, double t_cg = 15.0863, int conse=6,
           int num_samples=5000, int col_pos=1, int row_pos=1, bint b_outputCM=False, int starting_date=0):
    """
    dates: 1d array of shape(observation numbers), list of ordinal dates
    ts_b: 1d array of shape(observation numbers), the time series of blue band.
    ts_g: 1d array of shape(observation numbers), the time series of green band
    ts_r: 1d array of shape(observation numbers), the time series of red band
    ts_n: 1d array of shape(observation numbers), the time series of nir band
    ts_s1: 1d array of shape(observation numbers), the time series of swir1 band
    ts_s2: 1d array of shape(observation numbers), the time series of swir2 band
    ts_t: 1d array of shape(observation numbers), the time series of thermal band
    qas: 1d array, the QA bands
    t_cg: threshold of change magnitude, default is chi2.ppf(0.99,5)
    conse: consecutive observation number
    num_samples: column number per scanline, used to save pixel position
    col_pos: column position of current processing pixel, used to save pixel position
    row_pos: raw position of current processing pixel, used to save pixel position
    b_outputCM: bool, 'True' means outputting change magnitude and change magnitude dates, only for object-based COLD
    starting_date: the starting date of the whole dataset to enable reconstruct CM_date,
                   all pixels for a tile should have the same date, only for b_outputCM is True
    Note: passing 2-d array to c as 2-d pointer does not work, so have to pass separate bands
    """

    cdef int valid_num_scenes = qas.shape[0]
    # allocate memory for rec_cg
    cdef int num_fc = 0
    cdef Output_t t
    cdef Output_t* rec_cg = <Output_t*> PyMem_Malloc(NUM_FC * sizeof(t))
    # rec_cg = np.ndarray((NUM_FC, dtype=reccg_dt)
    # cdef Output_t [:] rec_cg_view = rec_cg  # memory view, ::1 means C contiguous

    # CM_outputs and CM_outputs_date are not used so far, but left for object-based cold (under development)
    if b_outputCM == True:
        CM_outputs = np.full(math.ceil((dates[valid_num_scenes-1] - starting_date) / SPATIAL_OUTPUT_INTERVAL) + 1, -9999, dtype=np.short)
        CM_outputs_date = np.full(math.ceil((dates[valid_num_scenes-1] - starting_date) / SPATIAL_OUTPUT_INTERVAL) + 1, 255, dtype=np.byte)
    else:  # set the length to 1 to save memory, as they won't be assigned values
        CM_outputs = np.full(1, -9999, dtype=np.short)
        CM_outputs_date = np.full(1, 255, dtype=np.byte)
    cdef short [:] CM_outputs_view = CM_outputs  # memory view
    cdef char [:] CM_outputs_date_view = CM_outputs_date  # memory view

    assert ts_b.shape[0] == dates.shape[0]
    assert ts_g.shape[0] == dates.shape[0]
    assert ts_r.shape[0] == dates.shape[0]
    assert ts_n.shape[0] == dates.shape[0]
    assert ts_s1.shape[0] == dates.shape[0]
    assert ts_s2.shape[0] == dates.shape[0]
    assert ts_t.shape[0] == dates.shape[0]
    assert qas.shape[0] == dates.shape[0]

    # examine if the qa input has filled value
    for x in range(valid_num_scenes):
        if qas[x] == 255:
            raise ValueError("qa array has filled values (255); please remove the rows of filled values")

    ret = cold(&ts_b[0], &ts_g[0], &ts_r[0], &ts_n[0], &ts_s1[0], &ts_s2[0], &ts_t[0], &qas[0], &dates[0],
               valid_num_scenes, num_samples, col_pos, row_pos, t_cg, conse, b_outputCM, starting_date, rec_cg,
               &num_fc, &CM_outputs_view[0], &CM_outputs_date_view[0])

    if b_outputCM == False:
        return np.asarray(<Output_t[:num_fc]>rec_cg) # np.asarray uses also the buffer-protocol and is able to construct
                                                     # a dtype-object from cython's array
    else:  # for object-based COLD
        return [np.asarray(<Output_t[:num_fc]>rec_cg), CM_outputs_view, CM_outputs_date_view]
