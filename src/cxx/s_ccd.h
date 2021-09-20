#ifndef CCD_STOCHASTIC_H
#define CCD_STOCHASTIC_H

#endif // CCD_STOCHASTIC_H
#include "KFAS.h"

int sccd
(
    short int **buf,            /* I/O:  pixel-based time series           */
    short int *fmask_buf,       /* I:  mask-based time series              */
    int *valid_date_array,      /* I: valid date time series               */
    int valid_num_scenes,       /* I: number of valid scenes under cfmask fill counts  */
    Output_t_sccd *rec_cg,           /* O: outputted structure for CCDC results    */
    int *num_fc,                /* O: number of fitting curves                       */
    int num_samples,            /* I: column number per scanline                    */
    int col_pos,                /* I: column position of current processing pixel   */
    int row_pos,                 /* I: raw position of current processing pixel  */
    bool b_fastmode,
    char* states_output_dir,
    double probability_threshold,
    int min_days_conse,
    int training_type, /* for training process*/
    int monitorwindow_lowerlin,
    int monitorwindow_upperlim,
    short int *sensor_buf,
    int n_focus_variable,
    int total_variable,
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
    bool b_landspecific,
    short int auxval,
    int conse
);


int sccd_stand_procedure
(
    int valid_num_scenes,             /* I:  number of scenes  */
    int *valid_date_array,           /* I: valid date time series  */
    short int *buf_b,            /* I:  Landsat blue spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    short int *buf_g,            /* I:  Landsat green spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    short int *buf_r,            /* I:  Landsat red spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    short int *buf_n,            /* I:  Landsat NIR spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    short int *buf_s1,           /* I:  Landsat swir1 spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    short int *buf_s2,           /* I:  Landsat swir2 spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    short int *buf_t,            /* I:  Landsat thermal spectral time series.The dimension is (n_obs, 7). Invalid (qa is filled value (255)) must be removed */
    short int *fmask_buf,           /* I:  mask-based time series  */
    int *id_range,
    Output_t_sccd *rec_cg,
    int *num_curve,                 /* Intialize NUM of Functional Curves    */
    char *states_output_dir,
    bool b_fastmode,
    double user_probability_threshold,
    int min_days_conse,
    int training_type, /* for training process*/
    int distdays_lowerlim,  /* for training process*/
    int distdays_upperlim,  /* for training process*/
    short int *sensor_buf,
    int user_n_focus_variable,
    int user_n_total_variable,
    int* user_focus_blist,
    bool NDVI_INCLUDED,
    bool NBR_INCLUDED,
    bool RGI_INCLUDED,
    bool TCTWETNESS_INCLUDED,
    bool TCTGREENNESS_INCLUDED,
    bool EVI_INCLUDED,
    bool DI_INCLUDED,
    bool NDMI_INCLUDED,
//    BoosterHandle booster,
    bool b_landspecific,
    short int auxval,
    int conse
);

int sccd_inefficientobs_procedure
(
    int valid_num_scenes,             /* I:  number of scenes  */
    int *valid_date_array,    /* I: valid date time series  */
    short int **buf,            /* I:  pixel-based time series  */
    short int *fmask_buf,      /* I:  mask-based time series  */
    int *id_range,
    double sn_pct,
    Output_t_sccd *rec_cg,
    int *num_curve
);

int step1_ccd_initialize
(
    int conse,              /* I: adjusted consecutive observation number               */
    float* adj_rmse,           /* I: the adjusted RMS                        */
    int n_clr,                 /* I: number of clear observations                         */
    double reg_TCG,               /* I: the record of fitting curve                        */
    int* i_dense,               /* I: used to count i for dense time point check          */
    int* num_curve,            /* I/O: the number of fitting curve                        */
    int *clrx,                  /* I/O: clear pixel curve in X direction (date)             */
    float **clry,               /* I/O: clear pixel curve in Y direction (spectralbands)    */
    int* cur_i,                /* I/O: the current number of monitoring observation          */
    int* i_start,              /* I/O: the start number of current curve                   */
    int* i_start_copy,
    int* end,                  /* I/O: the end of remaining observations.               */
    double **fit_cft,           /* I/O: Fitted coefficients 2-D array.                   */
    Output_t_sccd *rec_cg,            /* I/O: records of change points info                    */
    int i_span_min,                /* I: the minimum value for i_span                    */
    int* prev_i_break,              /*I : the i_break of the last curve                    */
    double *rmse,                     /* O: Root Mean Squared Error array.        */
    int n_focus_variable,
    int n_total_variable,
    int* focus_blist,
    int min_days_conse
);

int step1_update_cft
(
    int adj_conse,              /* I: adjusted consecutive observation number               */
    float* adj_rmse,           /* I: the adjusted RMS                        */
    int n_clr,                 /* I: number of clear observations                         */
    double reg_TCG,               /* I: the record of fitting curve                        */
    int *clrx,                  /* I/O: clear pixel curve in X direction (date)             */
    float **clry,               /* I/O: clear pixel curve in Y direction (spectralbands)    */
    int cur_i,                /* I/O: the current number of monitoring observation          */
    int i_start,              /* I/O: the start number of current curve                   */
    double **fit_cft,           /* I/O: Fitted coefficients 2-D array.                   */
    Output_t_sccd *rec_cg,            /* I/O: records of change points info                    */
    double *rmse,                     /* O: Root Mean Squared Error array.        */
    int *num_curve,
    int end,
    int *prev_i_break,
    int n_focus_variable,
    int n_total_variable,
    int* focus_blist,
    double t_cg_outelier
);

int step1_ssm_initialize(
    ssmodel* instance,          /* I/O: the outputted initial SS model                */
    int *clrx,                  /* I: clear pixel curve in X direction (date)             */
    float *clry,                /* I: clear pixel curve in Y direction (spectralbands)    */
    int stable_start,           /* I:  the start of the stable stage  */
    int stable_end,             /* I:  the start of the stable stage  */
    gsl_vector *next_a,               /* I:  initial a1  */
    gsl_matrix *next_P,          /* I:  initial P1  */
    int i_b,
    double **level_state_records,
    double **annual_state_records,
    double **semi_state_records,
    double **third_state_records,
    double **rmse_records,
    int starting_date,
    int prev_break_date,             /*I: the i_break of the previous curve*/
    int num_curve,
    Output_t_sccd *rec_cg,
    double rmse,
    double adjust_rmse,
    bool b_fastmode,
    double *fit_cft,
    int* valid_count,
    double *sum_square_smooth_Q,
    double *sum_square_smooth_H,
    double *sum_smooth_Q,
    double *sum_smooth_H,
    double *sum_square_vt,
    double *sum_vt,
    double *sum_kalman_coef,
    double **temporal_rmse_square_sum,
    double **temporal_rmse_sum,
    int *temporal_rmse_count
);

int step2_KF_ChangeDetection
(
    ssmodel* instance,
    gsl_matrix** vec_next_P,
    gsl_vector** vec_next_a,
    int *clrx,
    float **clry,
    int cur_i,
    double **level_state_records,
    double **annual_state_records,
    double **semi_state_records,
    double **third_state_records,
    double **rmse_records,
    Output_t_sccd* rec_cg,
    int *num_curve,
    int *end,
    int conse,
    int starting_date,
    int i_start,
    int *prev_i_break,             /*I: the i_break of the last curve*/
    bool b_fastmode,
    double** fit_cft1,
    double** fit_cft2,
    double** fit_cft3,
    double** fit_cft4,
    int *clrx_record1,
    int *clrx_record2,
    int *clrx_record3,
    int *clrx_record4,
    int *valid_count,
    double *sum_square_smooth_Q,
    double *sum_square_smooth_H,
    double *sum_smooth_Q,
    double *sum_smooth_H,
    double *sum_square_vt,
    double *sum_vt,
    double *sum_kalman_coef,
    int* i_count,
    float* adj_rmse,
    double t_cg_adjust,
    double t_cg_min,
    double t_cg_gradual,
    double **temporal_rmse_square_sum,
    double **temporal_rmse_sum,
    int *temporal_rmse_count,
    double **yearly_temporal_rmse_square_sum,
    double **yearly_temporal_rmse_sum,
    int *yearly_temporal_rmse_count,
    double probability_threshold,
    double *change_magnitude, /* for training process*/
    int training_type, /* for training process*/
    int n_focus_variable,
    int n_total_variable,
    int* focus_blist
);

int step2_chow_test
(
    int *clrx,
    float **clry,
    int first_seg_start,
    int first_seg_end,
    int second_seg_start,
    int second_seg_end,
    double f_prob,
    int n_focus_variable,
    int n_total_variable,
    int* focus_blist
);

int KF_ts_predict_conse
(
    ssmodel *instance,         /* i: the inputted ssm instance   */
    int *clrx,               /* i: the inputted dates   */
    gsl_matrix* P_ini,         /* i: a m x m matrix of the covariance matrix for pred_start */
    gsl_vector* at_ini,        /* i: a m vector of a for pred_start */
    int pred_start, /* i: close, included for prediction */
    int pred_end,   /* i: close, included for prediction */
    double *pred_y,   /* O: the predicted obs values */
    double *pred_y_f,  /*O: the predicted f (RMSE) values */
    bool b_fastmode,
    double *fit_cft,
    bool b_foutput
);



int KF_ts_predict_single_conse
(
    ssmodel *instance,         /* i: the inputted ssm instance   */
    int clrx,
    gsl_matrix* P_ini,         /* i: a m x m matrix of the covariance matrix for pred_start */
    gsl_vector* at_ini,        /* i: a m vector of a for pred_start */
    double *pred_y,   /* O: the predicted obs values */
    double *pred_y_f,  /*O: the predicted f (RMSE) values */
    bool b_fastmode,
    double *fit_cft
);

int KF_ts_filter_falsechange
(
    ssmodel *instance,         /* i: the inputted ssm instance   */
    int *clrx,               /* i: the inputted dates   */
    gsl_matrix* P_ini,         /* i/O: a m x m matrix of the covariance matrix for pred_start */
    gsl_vector* at_ini,        /* i/O: a m vector of a for pred_start */
    int cur_i,
    int i_band,                /* i: the band number */
    double **level_state_records,
    double **annual_state_records,
    double **semi_state_records,
    double **third_state_records,
    int starting_date,
    bool b_fastmode
);

int KF_ts_filter_regular
(
    ssmodel *instance,         /* i: the inputted ssm instance   */
    int *clrx,               /* i: the inputted dates   */
    float *clry,               /* i: the inputted observations   */
    gsl_matrix* P_ini,         /* i/O: a m x m matrix of the covariance matrix for pred_start */
    gsl_vector* at_ini,        /* i/O: a m vector of a for pred_start */
    int cur_i,
    int i_b,                /* i: the band number */
    double **rmse_records,   /* O: the predicted obs values */
    double **level_state_records,
    double **annual_state_records,
    double **semi_state_records,
    double **third_state_records,
    int starting_date,   /* starting_date is used to compute record index based on clrx*/
    bool b_fastmode,
    double *fit_cft,
    int *valid_count,
    double *sum_square_smooth_Q,
    double *sum_square_smooth_H,
    double *sum_smooth_Q,
    double *sum_smooth_H,
    double *sum_square_vt,
    double *sum_vt,
    double *sum_kalman_coef,
    bool b_start,
    double **temporal_rmse_square_sum,
    double **temporal_rmse_sum,
    int *temporal_rmse_count
);

int KF_ts_filter_norecord
(
    ssmodel *instance,         /* i: the inputted ssm instance   */
    int *clrx,               /* i: the inputted dates   */
    double *clry,               /* i: the inputted observations   */
    gsl_matrix* P_ini,         /* i/O: a m x m matrix of the covariance matrix for pred_start */
    gsl_vector* at_ini,        /* i/O: a m vector of a for pred_start */
    int cur_i,
    int i_band,                /* i: the band number */
    int starting_date   /* starting_date is used to compute record index based on clrx*/
);

/************************************************************************
FUNCTION: step3_processingend

PURPOSE:
Step 3 of S-CCD: processing the end of time series.
RETURN VALUE:
Type = int (SUCCESS OR FAILURE)

Programmer: Su Ye
**************************************************************************/
int step3_processingend
(
    ssmodel* instance,
    gsl_matrix** vec_next_P,
    gsl_vector** vec_next_a,
    int *clrx,
    float **clry,
    int cur_i,
    double **level_state_records,
    double **annual_state_records,
    double **semi_state_records,
    double **third_state_records,
    double **rmse_records,
    Output_t_sccd* rec_cg,
    int *num_curve,
    int *end,
    int conse,
    int starting_date,
    int bl_train_complete,
    int i_start,
    int prev_i_break,             /*I: the i_break of the last curve*/
    float* adj_rmse,           /* I: the adjusted RMS                        */
    bool b_fastmode,
    double **fit_cft1,
    double **fit_cft2,
    double **fit_cft3,
    double **fit_cft4,
    int *clrx_record1,
    int *clrx_record2,
    int *clrx_record3,
    int *clrx_record4,
    double probability_threshold,
    int* valid_count,
    double *sum_square_smooth_Q,
    double *sum_square_smooth_H,
    double *sum_smooth_Q,
    double *sum_smooth_H,
    double *sum_square_vt,
    double *sum_vt,
    double *sum_kalman_coef,
    double **temporal_rmse_square_sum,
    double **temporal_rmse_sum,
    int *temporal_rmse_count,
    int n_focus_variable,
    int n_total_variable,
    int* focus_blist,
    int min_days_conse
);

int grid_searching
(
     ssmodel *instance,
     double rmse,
     double adj_rmse,
     double mean_y,
     double *best_h,
     double *best_q,
     double interval,
     int actual_m,
     double ini_a,
     double *best_p
);

int Jazwinski_searching
(
     ssmodel *instance,
     double rmse,
     double mean_y,
     double *best_h,
     double *best_q,
     double interval,
     int actual_m
);

int getcategorycurve_old
(
    Output_t_sccd* t,
    int num_curve
);
