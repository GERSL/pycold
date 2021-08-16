from libc.stdlib cimport malloc
from libc.string cimport strcpy, strlen
from cython.view cimport array as cvarray
import numpy as np
import math


cdef int NUM_FC 40  # define the maximum number of outputted curves
cdef int SPATIAL_OUTPUT_INTERVAL 30  # for object-based COLD only

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

cdef extern from "../../src/cold.h":
    cdef int cold
    (
        short int **buf,            # I: Landsat spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed
        short int *fmask_buf,       # I: the time series of cfmask values. 0 - clear; 1 - water; 2 - shadow; 3 - snow; 4 - cloud
        int *valid_date_array,      # I: valid date as matlab serial date form (counting from Jan 0, 0000). Note ordinal date in python is from (Jan 1th, 0001)
        int valid_num_scenes,       # I: number of valid scenes
        int num_samples,            # I: column number per scanline, used to save pixel position
        int col_pos,                # I: column position of current processing pixel, used to save pixel position
        int row_pos,                # I: raw position of current processing pixel, used to save pixel position
        double tcg,                 # I: threshold of change threshold
        int conse,                  # I: consecutive observation number
        bool b_outputCM,            # I: indicate if outputting change magnitudes for object-based cold, for cold only, it is the false
        int starting_date,          # I: (optional) the starting date of the whole dataset to enable reconstruct CM_date, all pixels for a tile should have the same date, only for b_outputCM is True
        Output_t *rec_cg,           # O: outputted structure for CCDC results
        int *num_fc,                # O: number of fitting curves
        short int* CM_outputs,      # I/O: (optional) maximum change magnitudes at every SPATIAL_OUTPUT_INTERVAL days, only for b_outputCM is True
        char* CM_outputs_date       # I/O: (optional) dates for maximum change magnitudes at every SPATIAL_OUTPUT_INTERVAL days, only for b_outputCM is True
    );


def pycold(int[:] dates, short[:,::1] spectral_inputs, short[:] qas, double t_cg = 15.086, int conse=6,
           int num_samples=5000, int col_pos=0, int row_pos=0, bint b_outputCM=False, int starting_date=0) -> Output_t:
    """
    dates: 1d array, list of ordinal dates
    spectral_inputs: 2d array of shape(observation numbers, bands), note that the band order is blue, green,
                    red, NIR, SWIR1, SWIR2, TIR. short[:,::1] give Cython the information that the data is C-contiguous
    qas: 1d array, the QA bands
    t_cg: threshold of change magnitude, default is chi2.ppf(0.99,5)
    conse: consecutive observation number
    num_samples: column number per scanline, used to save pixel position
    col_pos: column position of current processing pixel, used to save pixel position
    row_pos: raw position of current processing pixel, used to save pixel position
    b_outputCM: bool, 'True' means outputting change magnitude and change magnitude dates, only for object-based COLD
    starting_date: the starting date of the whole dataset to enable reconstruct CM_date,
                   all pixels for a tile should have the same date, only for b_outputCM is True
    """

    cdef int valid_num_scenes = qas.shape[0]
    # allocate memory for rec_cg
    cdef int* num_fc = <int*>malloc(1 * sizeof(int))

    rec_cg = np.empty(NUM_FC, reccg_dt)
    cdef reccg_dt rec_cg_view[:] = rec_cg
    # CM_outputs and CM_outputs_date are not used so far, but left for object-based cold (under development)
    CM_outputs = np.full(math.ceil((dates[valid_num_scenes-1] - starting_date) / SPATIAL_OUTPUT_INTERVAL) + 1, -9999, dtype=np.short)
    cdef short [:] CM_outputs_view = CM_outputs  # memory view
    CM_outputs_date = np.full(math.ceil((dates[valid_num_scenes-1] - starting_date) / SPATIAL_OUTPUT_INTERVAL) + 1, -9999, dtype=np.byte)
    cdef char [:] CM_outputs_date_view = CM_outputs_date  # memory view


    if spectral_inputs.shape[0] == 0:
        raise ValueError("No time series inputs!")
    if spectral_inputs.shape[1] != 7:
        raise ValueError("The band number of spectral_inputs is required to be 7 (the last column is thermal band)!")
    assert spectral_inputs.shape[0] == dates.shape[0]:
    asset spectral_inputs.shape[0] == qas.shape[0]:

    # examine if the qa input has filled value
    for x in range(valid_num_scenes):
        if qas[x] == 255:
            raise ValueError("qa array has filled values (255); please remove the rows of filled values and continue")


    ret = cold(spectral_inputs, qas, dates, valid_num_scenes, num_samples, col_pos, row_pos, tcg, conse,
               b_outputCM, starting_date, rec_cg_view, num_fc, CM_outputs_view, CM_outputs_date_view);

    free(num_fc)
    if b_outputCM == False:
        return rec_cg_view[0:*num_fc]
    else:  # for object-based COLD
        return [rec_cg_view[0:*num_fc], CM_outputs_view, CM_outputs_date_view]



