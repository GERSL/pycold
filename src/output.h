#ifndef OUTPUT_H
#define OUTPUT_H
#include <stdbool.h>
#include "defines.h"
#include "input.h"

typedef struct {
  int row;
  int col;
} Position_t;

/* Structure for the 'output' data type */
typedef struct
{
    int t_start;           /* time when series model gets started */
    int t_end;             /* time when series model gets ended */
    int t_break;           /* time when the first break (change) is observed */
    int pos;        /* the location of each time series model */
    int num_obs;           /* the number of "good" observations used for model
                              estimation */
    short int category;          /* the quality of the model estimation (what model
                              is used, what process is used) */
    /*
    The current category in output structure:
    first digit:
    0: normal model (no change)
    1: change at the beginning of time series model
    2: change at the end of time series model
    3: disturbance change in the middle
    4: fmask fail scenario
    5: permanent snow scenario
    6: outside user mask
    second digit:
    1: model has only constant term
    4: model has 3 coefs + 1 const
    6: model has 5 coefs + 1 const
    8: model has 7 coefs + 1 const*/
    short int change_prob;     /* the probability of a pixel that have undergone
                                  change (between 0 and 100) */
    double coefs[TOTAL_IMAGE_BANDS][NUM_COEFFS];
                           /*  coefficients for each time series model for each 
                               spectral band*/    
    double rmse[TOTAL_IMAGE_BANDS];
                           /*  RMSE for each time series model for each 
                               spectral band*/
    double magnitude[TOTAL_IMAGE_BANDS];/* the magnitude of change (difference between model
                                  prediction and observation for each spectral band)*/
} Output_t;


typedef struct
{
    int t_start;           /* time when series model gets started */
    int t_end;             /* time when series model gets ended */
    int t_break;           /* time when the first break (change) is observed */
    int pos;        /* the location of each time series model */
    // float change_prob;     /* the probability of a pixel that have undergone */
                              /* change (between 0 and 100) */
    int num_obs;           /* the number of "good" observations used for model
                              estimation */
    short int category;          /* the quality of the model estimation (what model
                              is used, what process is used) */
    /*
    The current category in output structure:
    first digit:
    0: normal model (abrupt)
    1: change at the beginning of time series model
    2: change at the end of time series model
    3: disturbance change in the extended period of initialization
    4: fmask fail scenario
    5: permanent snow scenario
    6: outside user mask
    7: normal model (short term)
    second digit:
    0: undefined
    1: land disturbance
    2: reforestation/afforestation
    3: regrowth */

    short int land_type;
    int t_confirmed;          /* lag for near real-time */
    int change_prob;     /* the probability of a pixel that have undergone
                              change (between 0 and 100) */
    double coefs[TOTAL_IMAGE_BANDS][SCCD_MAX_NUM_C];
                           /*  coefficients for each time series model for each
                               spectral band*/
//    double obs_disturb[TOTAL_IMAGE_BANDS];
//    double state_disturb[TOTAL_IMAGE_BANDS][SCCD_MAX_NUM_C - 1];
    double rmse[TOTAL_IMAGE_BANDS];

    double magnitude[TOTAL_IMAGE_BANDS];/* the magnitude of change (difference between model
                                  prediction and observation for each spectral band)*/


//    double coefs[TOTAL_IMAGE_BANDS+TOTAL_INDICES][SCCD_MAX_NUM_C];
//                           /*  coefficients for each time series model for each
//                               spectral band*/
//    double obs_disturb[TOTAL_IMAGE_BANDS+TOTAL_INDICES];
//    double state_disturb[TOTAL_IMAGE_BANDS+TOTAL_INDICES][SCCD_MAX_NUM_C - 1];
//    double rmse[TOTAL_IMAGE_BANDS+TOTAL_INDICES];

//    double magnitude[TOTAL_IMAGE_BANDS+TOTAL_INDICES];/* the magnitude of change (difference between model
//                                  prediction and observation for each spectral band)*/


} Output_t_sccd;

/* SY 09242018 */
int firstDegradationYear
(
    Output_t* t,  /*I: outputted structure */
    int num_fc,   /*I: the number of inputted models */
    int n_lassoband,  /*I:the number of lasso band to focus on*/
    bool updirection, /*I: updirection or downdirection for the inputted band */
    int *dyear        /*O: outputed degradation year */
);

int write_raw_binary
(
    FILE *rb_fptr,      /* I: pointer to the raw binary file */
    int nlines,         /* I: number of lines to write to the file */
    int nsamps,         /* I: number of samples to write to the file */
    int size,           /* I: number of bytes per pixel (ex. sizeof(uint8)) */
    int *img_array     /* I: array of nlines * nsamps * size to be written
                              to the raw binary file */
);

void write_envi_header
(
    FILE *fptr,            /* I: pointer to the raw binary file */
    Input_meta_t *meta     /* I: saved header file info */
);

void write_output_csv
(
    FILE *csv_fptr,      /* I: pointer to the csv file */
    Output_t t           /* I: outputted structure     */
);

int write_output_binary
(
    FILE *fptr,      /* I: pointer to the binary file */
    Output_t t          /* I: outputted structure     */
);

int write_output_binary_sccd
(
    FILE *fptr,      /* I: pointer to the binary file */
    Output_t_sccd t          /* I: outputted structure     */
);

#endif
