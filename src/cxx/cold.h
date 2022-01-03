#ifndef CCD_H
#define CCD_H
#include <stdbool.h>
#include "output.h"
// #include <xgboost/c_api.h>

int preprocessing
(
    long *buf_b,            /* I:  Landsat blue spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_g,            /* I:  Landsat green spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_r,            /* I:  Landsat red spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_n,            /* I:  Landsat NIR spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_s1,           /* I:  Landsat swir1 spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_s2,           /* I:  Landsat swir2 spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_t,            /* I:  Landsat thermal spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *fmask_buf,        /* I:   mask time series  */
    int *valid_num_scenes, /* I/O: * number of scenes after cfmask counts and  */
    int *id_range,
    int *clear_sum,      /* I/O: Total number of clear cfmask pixels          */
    int *water_sum,      /* I/O: counter for cfmask water pixels.             */
    int *shadow_sum,     /* I/O: counter for cfmask shadow pixels.            */
    int *sn_sum,         /* I/O: Total number of snow cfmask pixels           */
    int *cloud_sum      /* I/O: counter for cfmask cloud pixels.             */
);


int cold
(
    long *buf_b,            /* I:  Landsat blue spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_g,            /* I:  Landsat green spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_r,            /* I:  Landsat red spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_n,            /* I:  Landsat NIR spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_s1,           /* I:  Landsat swir1 spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_s2,           /* I:  Landsat swir2 spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_t,            /* I:  Landsat thermal spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *fmask_buf,       /* I:  the time series of cfmask values. 0 - clear; 1 - water; 2 - shadow; 3 - snow; 4 - cloud  */
    long *valid_date_array,      /* I:  valid date as matlab serial date form (counting from Jan 0, 0000). Note ordinal date in python is from (Jan 1th, 0001) */
    int valid_num_scenes,       /* I: number of valid scenes  */
    int pos,                     /* I: the position id of pixel */
    double tcg,                 /* I: threshold of change threshold  */
    int conse,                  /* I: consecutive observation number   */
    bool b_outputCM,            /* I: indicate if outputting change magnitudes for object-based cold, for cold only, it is the false */
    int starting_date,          /* I: (optional) the starting date of the whole dataset to enable reconstruct CM_date, all pixels for a tile should have the same date, only for b_outputCM is True */
    Output_t *rec_cg,           /* O: outputted structure for CCDC results    */
    int *num_fc,                /* O: number of fitting curves                   */
    int CM_OUTPUT_INTERVAL,     /* I: the interval of days of outputting change maganitude                   */
    short int* CM_outputs,      /* I/O: (optional) maximum change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
    unsigned char* CMdirection_outputs,      /* I/O: direction of change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
    unsigned char* CM_outputs_date       /* I/O: (optional) dates for maximum change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
);


int stand_procedure
(
    int valid_num_scenes,             /* I:  number of valid scenes  */
    long *valid_date_array,            /* I: valid date time series  */
    long *buf_b,            /* I:  Landsat blue spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_g,            /* I:  Landsat green spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_r,            /* I:  Landsat red spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_n,            /* I:  Landsat NIR spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_s1,           /* I:  Landsat swir1 spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_s2,           /* I:  Landsat swir2 spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_t,            /* I:  Landsat thermal spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *fmask_buf,       /* I:  mask-based time series  */
    int *id_range,
    double tcg,                 /* I: threshold of change threshold  */
    int conse,                  /* I: consecutive observation number   */
    bool b_outputCM,              /* I: indicate if outputting change magnitudes object-based cold*/
    int starting_date,           /* I: the starting date of the whole dataset to enable reconstruct CM_date, all pixels for a tile should have the same date, only for b_outputCM is True */
    Output_t *rec_cg,           /* O: outputted structure for CCDC results     */
    int *num_fc,                /* O: number of fitting curves                       */
    int CM_OUTPUT_INTERVAL,
    short int* CM_outputs,      /* I/O: maximum change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
    unsigned char* CMdirection_outputs,      /* I/O: direction of change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
    unsigned char* CM_outputs_date      /* I/O: dates for maximum change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
);

int inefficientobs_procedure
(
    int valid_num_scenes,             /* I:  number of scenes  */
    long *valid_date_array,    /* I: valid date time series  */
    long *buf_b,            /* I:  Landsat blue spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_g,            /* I:  Landsat green spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_r,            /* I:  Landsat red spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_n,            /* I:  Landsat NIR spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_s1,           /* I:  Landsat swir1 spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_s2,           /* I:  Landsat swir2 spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_t,            /* I:  Landsat thermal spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *fmask_buf,      /* I:  mask-based time series  */
    int *id_range,
    float sn_pct,
    Output_t *rec_cg,
    int *num_fc
);

const char *check_parameter(
        int mode,
        char* in_path,
        char* out_path,
        int n_cores,
        int row,
        int col,
        int task,
        char* mask_path,
        double probability_threshold,
        int conse,
        int min_days_conse,
        int output_mode,
        int verbose
);


int tsalgorithm_executor(
    int mode,
    char* in_path,
    char* out_path,
    int n_cores,
    int row,
    int col,
    int task,
    char* mask_path,
    float probability_threshold,
    int conse,
    int min_days_conse,
    int output_mode,
    int verbose
);

int obcold_reconstruction_procedure
(
    short int **buf,            /* I:  pixel-based time series  */
    short int *fmask_buf,      /* I:  mask-based time series  */
    int *valid_date_array,    /* I: valid date time series  */
    int valid_num_scenes,             /* I:  number of scenes  */
    Output_t *rec_cg,    /* O: Initialize NUM of Functional Curves    */
    int *num_fc,
    int num_year,       /*I: the number of focused years */
    int *break_dates, /*an array of break dates with a fixed length of num_year, '0' means no breaks */
    int num_samples,            /* I: column number per scanline                    */
    int col_pos,                /* I: column position of current processing pixel   */
    int row_pos,
    int conse
);

double angle_decaying(
    double input,
    double lowbound,
    double highbound
);


int ccd_scanline
(
    int row,                 /* I:   input row. The beginning number is 1, not 0!   */
    char *in_path,           /* I:   Landsat ARD directory  */
    char **scene_list,       /* I:   current scene name in list of sceneIDs       */
    //char *mask_path,       /* I:   mask_path */
    char *user_mask_path,     /* I:   the path of inputted mask       */
    int  num_samples,        /* I:   number of image samples (X width)      */
    int  num_scenes,         /* I:   current num. in list of scenes to read */
    int *sdate,              /* I:   Original array of julian date values         */
    char *out_path,
    int method,
    double probability_threshold,
    int min_days_conse,
    int n_focus_variable,
    int n_total_variable,
    int* focus_blist,
    bool NDVI_INCLUDED,
    bool NBR_INCLUDED,
    bool RGI_INCLUDED,
    bool TCTWETNESS_INCLUDED,
    bool TCTGREENNESS_INCLUDED,
    bool EVI_INCLUDED,
    bool DI_INCLUDED,
    bool NDMI_INCLUDED,
//    BoosterHandle booster,
    bool b_landspecific_mode,
    char *auxiliary_var_path,
    int conse,
    bool b_outputCM,
    bool b_outputCM_reconstruction,
    int CM_OUTPUT_INTERVAL
);

int stand_procedure_fixeddays
(
    int valid_num_scenes,             /* I:  number of valid scenes  */
    long *valid_date_array,            /* I: valid date time series  */
    long *buf_b,            /* I:  Landsat blue spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_g,            /* I:  Landsat green spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_r,            /* I:  Landsat red spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_n,            /* I:  Landsat NIR spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_s1,           /* I:  Landsat swir1 spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_s2,           /* I:  Landsat swir2 spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *buf_t,            /* I:  Landsat thermal spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    long *fmask_buf,       /* I:  mask-based time series  */
    int *id_range,
    double tcg,                 /* I: threshold of change threshold  */
    int conse,                  /* I: consecutive observation number   */
    bool b_outputCM,              /* I: indicate if cold is running as the first step of object-based cold*/
    int starting_date,           /* I: the starting date of the whole dataset to enable reconstruct CM_date, all pixels for a tile should have the same date, only for b_outputCM is True */
    Output_t *rec_cg,           /* O: outputted structure for CCDC results     */
    int *num_fc,                /* O: number of fitting curves                       */
    int CM_OUTPUT_INTERVAL,
    short int* CM_outputs,      /* I/O: maximum change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
    unsigned char* CMdirection_outputs,      /* I/O: direction of change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
    unsigned char* CM_outputs_date,      /* I/O: dates for maximum change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
    int min_days_conse
);

#endif // CCD_H
