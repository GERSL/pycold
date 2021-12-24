from libc.stdlib cimport malloc
from libc.string cimport strcpy, strlen
import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
import math
from libcpp cimport bool
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free


cdef int NUM_FC = 40  # define the maximum number of outputted curves

DTYPE_date = np.float64
ctypedef np.float64_t DTYPE_t

reccg_dt = np.dtype([('t_start', np.int32),  # time when series model gets started
                     ('t_end', np.int32),  # time when series model gets ended
                     ('t_break', np.int32),  # time when the first break (change) is observed
                     ('pos', np.int32),  # the location of each time series model
                     ('num_obs', np.int32),  # the number of "good" observations used for model estimation
                     ('category', np.short),  # the quality of the model estimation (what model is used, what process is used)
                     ('change_prob', np.short),  # the probability of a pixel that have undergone change (between 0 and 100)
                     ('coefs', np.float, (7, 8)),  # coefficients for each time series model for each spectral band
                     ('rmse', np.float, 7),  # RMSE for each time series model for each spectral band
                     ('magnitude', np.float, 7)])  # the magnitude of change difference between model prediction
                                          # and observation for each spectral band)


cdef extern from "../../cxx/output.h":
    ctypedef struct Output_t:
        int t_start
        int t_end
        int t_break
        int pos
        int num_obs
        short int category
        short int change_prob
        float coefs[7][8]
        float rmse[7]
        float magnitude[7]

cdef extern from "../../cxx/cold.h":
    cdef int cold(long *buf_b, long *buf_g, long *buf_r, long *buf_n, long *buf_s1, long *buf_s2,
                  long *buf_t, long *fmask_buf, long *valid_date_array, int valid_num_scenes, int num_samples,
                  int col_pos, int row_pos, double tcg, int conse, bool b_outputCM, int starting_date, Output_t *rec_cg,
                  int *num_fc, int CM_OUTPUT_INTERVAL, short int *CM_outputs, unsigned char *CMdirection_outputs, unsigned char *CM_outputs_date);



def cold_detect(np.ndarray[np.int64_t, ndim=1] dates, np.ndarray[np.int64_t, ndim=1] ts_b, np.ndarray[np.int64_t, ndim=1] ts_g,
                np.ndarray[np.int64_t, ndim=1] ts_r, np.ndarray[np.int64_t, ndim=1] ts_n, np.ndarray[np.int64_t, ndim=1] ts_s1,
                np.ndarray[np.int64_t, ndim=1] ts_s2, np.ndarray[np.int64_t, ndim=1] ts_t, np.ndarray[np.int64_t, ndim=1] qas,
                double t_cg = 15.0863, int conse=6, int num_samples=5000, int col_pos=1, int row_pos=1, bint b_outputCM=False,
                int starting_date=0, int CM_OUTPUT_INTERVAL=32):
    """
    Helper function to do COLD algorithm.

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
    	conse: consecutive observation number
    	num_samples: column number per scanline, used to save pixel position
    	col_pos: column position of current processing pixel, used to save pixel position
    	row_pos: raw position of current processing pixel, used to save pixel position
    	b_outputCM: bool, 'True' means outputting change magnitude and change magnitude dates, only for object-based COLD
    	starting_date: the starting date of the whole dataset to enable reconstruct CM_date,
                   	all pixels for a tile should have the same date, only for b_outputCM is True
    	CM_OUTPUT_INTERVAL: the temporal interval of outputting change magnitudes
    	Note that passing 2-d array to c as 2-d pointer does not work, so have to pass separate bands
    	Returns
    	----------
    	change records: the COLD outputs that characterizes each temporal segment
    """

    cdef int valid_num_scenes = qas.shape[0]
    # allocate memory for rec_cg
    cdef int num_fc = 0
    cdef Output_t t
    cdef Output_t* rec_cg = <Output_t*> PyMem_Malloc(NUM_FC * sizeof(t))
    # rec_cg = np.ndarray((NUM_FC, dtype=reccg_dt)
    # cdef Output_t [:] rec_cg_view = rec_cg  # memory view, ::1 means C contiguous

    cdef long [:] dates_view = dates
    cdef long [:] ts_b_view = ts_b
    cdef long [:] ts_g_view = ts_g
    cdef long [:] ts_r_view = ts_r
    cdef long [:] ts_n_view = ts_n
    cdef long [:] ts_s1_view = ts_s1
    cdef long [:] ts_s2_view = ts_s2
    cdef long [:] ts_t_view = ts_t
    cdef long [:] qas_view = qas

    assert ts_b_view.shape[0] == dates_view.shape[0]
    assert ts_g_view.shape[0] == dates_view.shape[0]
    assert ts_r_view.shape[0] == dates_view.shape[0]
    assert ts_n_view.shape[0] == dates_view.shape[0]
    assert ts_s1_view.shape[0] == dates_view.shape[0]
    assert ts_s2_view.shape[0] == dates_view.shape[0]
    assert ts_t_view.shape[0] == dates_view.shape[0]
    assert qas_view.shape[0] == dates_view.shape[0]

    # CM_outputs and CM_outputs_date are not used so far, but left for object-based cold (under development)
    if b_outputCM == True:
        CM_outputs = np.full(math.ceil((dates[valid_num_scenes-1] - starting_date) / CM_OUTPUT_INTERVAL) + 1, -9999, dtype=np.short)
        CM_direction_outputs = np.full(math.ceil((dates[valid_num_scenes-1] - starting_date) / CM_OUTPUT_INTERVAL) + 1, 255, dtype=np.uint8)
        CM_outputs_date = np.full(math.ceil((dates[valid_num_scenes-1] - starting_date) / CM_OUTPUT_INTERVAL) + 1, 255, dtype=np.uint8)
    else:  # set the length to 1 to save memory, as they won't be assigned values
        CM_outputs = np.full(1, -9999, dtype=np.short)
        CM_direction_outputs = np.full(1, 255, dtype=np.uint8)
        CM_outputs_date = np.full(1, 255, dtype=np.uint8)
    cdef short [:] CM_outputs_view = CM_outputs  # memory view
    cdef unsigned char [:] CM_direction_outputs_view = CM_direction_outputs  # memory view
    cdef unsigned char [:] CM_outputs_date_view = CM_outputs_date  # memory view

    result = cold(&ts_b_view[0], &ts_g_view[0], &ts_r_view[0], &ts_n_view[0], &ts_s1_view[0], &ts_s2_view[0], &ts_t_view[0],
                 &qas_view[0], &dates_view[0], valid_num_scenes, num_samples, col_pos, row_pos, t_cg, conse, b_outputCM,
                 starting_date, rec_cg, &num_fc, CM_OUTPUT_INTERVAL, &CM_outputs_view[0], &CM_direction_outputs_view[0], &CM_outputs_date_view[0])
    if result != 0:
        raise RuntimeError("cold function fails for col_pos = {} and row_pos = {}".format(col_pos, row_pos))
    else:
        if b_outputCM == False:
            return np.asarray(<Output_t[:num_fc]>rec_cg) # np.asarray uses also the buffer-protocol and is able to construct
                                                         # a dtype-object from cython's array
        else:  # for object-based COLD
            return [np.asarray(<Output_t[:num_fc]>rec_cg), CM_outputs, CM_direction_outputs, CM_outputs_date]
