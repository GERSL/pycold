#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <stdbool.h>
#include <unistd.h>
#include <math.h>
#include "defines.h"
#include "cold.h"
#include "const.h"
#include "utilities.h"
#include "2d_array.h"
#include "input.h"
#include "output.h"
#include "misc.h"
#include "s_ccd.h"
#include "distribution_math.h"
#include <sys/stat.h>


int lasso_blist[NUM_LASSO_BANDS] = {1, 2, 3, 4, 5};

#define safe_xgboost(call) {                                            \
int err = (call);                                                       \
if (err != 0) {                                                         \
  fprintf(stderr, "%s:%d: error in %s: %s\n", __FILE__, __LINE__, #call, XGBGetLastError()); \
  exit(1);                                                              \
}                                                                       \
}

// convert NLCD category index for grand county study
int getlabelfromNLCD(int maskval){
    if(maskval == 21 || maskval == 22 || maskval == 23 || maskval == 24)
        return 0;
    else if(maskval == 81 || maskval == 82 )
        return 1;
    else if(maskval == 41 || maskval == 42 || maskval == 43)
        return 2;
    else if(maskval == 71 || maskval == 52)
        return 3;
    else if(maskval == 90 || maskval == 95)
        return 4;
    else if(maskval == 11)
        return 5;
    else if(maskval==31)
        return 6;
    else
        return NA_VALUE;
}

long getMicrotime(){
    struct timeval currentTime;
    gettimeofday(&currentTime, NULL);
    return currentTime.tv_sec * (int)1e6 + currentTime.tv_usec;
}

void substr(char dest[], char src[], int offset, int len)
{
    int i;
    for(i = 0; i < len && src[offset + i] != '\0'; i++)
        dest[i] = src[i + offset];
    dest[i] = '\0';
}




/******************************************************************************
MODULE:  cold

PURPOSE:  main function for cold

RETURN VALUE:
Type = int (SUCCESS OR FAILURE)

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
11/14/2018   Su Ye         Original Development
******************************************************************************/
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
    bool b_c2,                  /* I: a temporal parameter to indicate if collection 2. C2 needs ignoring thermal band due to the current low quality  */
    Output_t *rec_cg,           /* O: outputted structure for CCDC results    */
    int *num_fc,                /* O: number of fitting curves                   */
    int CM_OUTPUT_INTERVAL,
    short int* CM_outputs,      /* I/O: (optional) maximum change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
//    unsigned char* CMdirection_outputs,      /* I/O: direction of change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
    short int* CM_outputs_date,        /* I/O: (optional) dates for maximum change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
    double gap_days                /* I: the day number of gap to define i_dense; it is useful for the cases that gap is in the middle of time series      */
)
{
    int clear_sum = 0;      /* Total number of clear cfmask pixels          */
    int water_sum = 0;      /* counter for cfmask water pixels.             */
    int shadow_sum = 0;     /* counter for cfmask shadow pixels.            */
    int sn_sum = 0;         /* Total number of snow cfmask pixels           */
    int cloud_sum = 0;      /* counter for cfmask cloud pixels.             */
    float sn_pct;           /* Percent of snow pixels.                      */
    int status;
    int *id_range;
    int i;
    char FUNC_NAME[] = "cold";
    int result;

    if(valid_num_scenes == 0){
        return (SUCCESS);
    }

    id_range = (int*)calloc(valid_num_scenes, sizeof(int));

    // check if valid_date_array is sorted
    for (i = 0; i < valid_num_scenes - 1; i++)
        if (valid_date_array[i] > valid_date_array[i+1]){
            RETURN_ERROR("The inputted data does not follow an ascending order!", FUNC_NAME, ERROR);
        }

    status = preprocessing(buf_b, buf_g, buf_r, buf_n, buf_s1, buf_s2, buf_t,
                           fmask_buf, &valid_num_scenes, id_range, &clear_sum,
                           &water_sum, &shadow_sum, &sn_sum, &cloud_sum, b_c2);
    // printf("preprocessing finished \n");

    if (status != SUCCESS)
    {
        RETURN_ERROR("Error for preprocessing.", FUNC_NAME, ERROR);
    }

    // clear_pct is not used anymore in V13.01
    //clr_pct = (double) clear_sum / (double) (valid_num_scenes);

    sn_pct = (float) sn_sum/ (float) (sn_sum + clear_sum + 0.01);

    /******************************************************************/
    /*************** rec_cg initialization ****************************/
    /******************************************************************/
//    for(k = 0; k < NUM_FC; k++)
//    {
//       rec_cg[k].pos = -1;
//       rec_cg[k].category = -1;
//       rec_cg[k].t_start = -1;
//       rec_cg[k].t_end = -1;
//       rec_cg[k].t_break = -1;
//       rec_cg[k].num_obs = -1;
//       rec_cg[k].change_prob = -1;

//        for (i = 0; i < TOTAL_IMAGE_BANDS; i++)
//            rec_cg[k].magnitude[i] = -1;

//        for (i = 0; i < TOTAL_IMAGE_BANDS; i++)
//            for(j = 0; j < NUM_COEFFS; j++)
//                rec_cg[k].coefs[i][j] = -1;

//        for (i = 0; i < TOTAL_IMAGE_BANDS; i++)
//            rec_cg[k].rmse[i] = -1;
//    }

    // if ((clr_pct < T_CLR)||(clear_sum < N_TIMES * MAX_NUM_C)){
    if (clear_sum < N_TIMES * MAX_NUM_C){
        result = inefficientobs_procedure(valid_num_scenes,valid_date_array, buf_b, buf_g, buf_r, buf_n, buf_s1, buf_s2, buf_t,
                                 fmask_buf, id_range,sn_pct,rec_cg, num_fc, b_c2);
    }
    else{

        /**************************************************************/
        /*                                                            */
        /* standard_procedure for CCD                                 */
        /*                                                            */
        /**************************************************************/
         result = stand_procedure(valid_num_scenes, valid_date_array, buf_b, buf_g, buf_r, buf_n, buf_s1, buf_s2, buf_t, fmask_buf, id_range,
                                 tcg, conse, b_outputCM, starting_date, rec_cg, num_fc, CM_OUTPUT_INTERVAL, CM_outputs,
                                 CM_outputs_date, b_c2, gap_days);
//        result = stand_procedure_fixeddays(valid_num_scenes, valid_date_array, buf_b, buf_g, buf_r, buf_n, buf_s1, buf_s2, buf_t, fmask_buf, id_range,
//                                 tcg, conse, b_outputCM, starting_date, rec_cg, num_fc, CM_OUTPUT_INTERVAL, CM_outputs,
//                                 CM_outputs_date, (conse - 1) * 16);


       // printf("stand procedure finished \n");
    }

    for(i = 0; i < *num_fc; i++)
    {
        rec_cg[i].pos = pos;
    }

    free(id_range);

    // for debug
    // printf("free stage 6 \n");

    if (result == SUCCESS)
    {
        return (SUCCESS);
    }
    else
    {
        return (FAILURE);
    }
}

/****************************************************************************
MODULE:  stand_procedure for fixed days peek window
PURPOSE:  standard procedure when having enough clear pixels
RETURN VALUE:
Type = int (SUCCESS OR FAILURE)
HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
11/14/2018   Su Ye          Modification from main function in original CCDC.c
******************************************************************************/

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
    short int* CM_outputs_date,      /* I/O: dates for maximum change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
    int min_days_conse
)
{
    int status;
    int i, j, k, k_new, b;
    int m;                       /*the number of ID to be removed */
    char FUNC_NAME[] = "stand_procedure";
    //char msg_str[MAX_STR_LEN];       /* Input data scene name                 */
    int *rm_ids;
    int *cpx;                        /* nunber of clear pixels X ?            */
    float **cpy;                     /* nunber of clear pixels Y ?            */
    int i_rec;                       /* start of model before noise removal   */

    int end;
    float **fit_cft;                 /* Fitted coefficients 2-D array.        */
    float **rec_v_dif;
    float **rec_v_dif_copy;
    float v_dif_norm = 0.0;
    float mini_rmse;                 /* Mimimum RMSE                          */

    float v_dif_mean;
    float *rmse;                     /* Root Mean Squared Error array.        */
    int i_count;                     /* Count difference of i each iteration  */
    int i_break;                     /* for recording break points, i is index*/
    float **v_diff;
    float *vec_mag; /* what is the differece */ /* they are used in 2 different branches */
    float *vec_magg;/* these two?            */ /* this one is never freed */
    float vec_magg_min;
    float **v_dif_mag;               /* vector for magnitude of differences.  */
    float **temp_v_dif;              /* for the thermal band.......           */
    int rec_fc;                      /* Record num. of functional curves      */
    int i_conse;
    int i_b;
    int ini_conse;
    int i_ini;                       /* for recording begin of time, i is index*/
    int i_dense;
    float ts_pred_temp;
    int i_span, rm_ids_len;
    float time_span;                 /* Span of time in no. of years.         */
    int *bl_ids;
    int *ids;
    int bl_train;                    /* Flag for which way to train the model.*/
    int i_start;                     /* The first observation for TSFit       */
    int num_c = 8;                   /* Max number of coefficients for model  */
    int update_num_c = 8;            /* Number of coefficients to update      */
    int n_rmse;                      /* number of RMSE values                 */
    float *d_yr;
    int d_rt;
    float* v_start;  /* Vector for start of observation(s)    */
    float* v_end;    /* Vector for end of observastion(s)     */
    float* v_slope;  /* Vector for anormalized slope values   */
    float* v_dif;    /* Vector for difference values          */
    float* tmpcg_rmse; /* to temporarily change RMSE          */

    int ids_old_len;
    int *ids_old;
    int id_last;                     /* The last stable id.                   */

    int ids_len;                     /* number of ids, incremented continuously*/
    float break_mag;
    int adj_conse;
    float adj_TCG;
    float min_rmse[TOTAL_IMAGE_BANDS]; /* Adjusted RMSE for all bands          */
    float date_vario;           /* I: median date                                          */
    float max_date_difference;   /* I: maximum difference between two neighbor dates        */
    float** v_diff_tmp;

    int n_clr;                  /* I: the number of clear pixels                          */
    int *clrx;                  /* I: clear pixel curve in X direction (date)             */
    float **clry;               /* I: clear pixel curve in Y direction (spectralbands)    */
    float mean_angle;           /* I: mean angle of vec_diff                              */
    int pre_end;
    // short int tmp_max_prob = 0;
    short int tmp_CM = 0;
    int tmp;
    int current_CM_n;
    float prob_angle; // change probability for angle
    int i_span_skip = 0;
    unsigned char tmp_direction;  // tmp variable for storing the directionality of change magnitude
    int conse_min;
    int adj_conse_end;

    fit_cft = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, LASSO_COEFFS,
                                         sizeof (float));
    if (fit_cft == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft memory", FUNC_NAME, FAILURE);
    }

    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
    {
            min_rmse[k] = 0.0;
    }



    rmse = (float *)calloc(TOTAL_IMAGE_BANDS, sizeof(float));
    if (rmse == NULL)
    {
        RETURN_ERROR ("Allocating rmse memory", FUNC_NAME, FAILURE);
    }

    n_clr = 0;
    clrx = (int* )malloc(valid_num_scenes * sizeof(int));
    clry = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, valid_num_scenes,
                                         sizeof (float));

//        for (i = 0; i < valid_num_scenes; i++)
//        {
//            printf("clrx: %d\n", (int)clrx[i]);
//            printf("id: %d\n", (int)id_range[i]);
//            for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
//            {
//               printf("buf: %d\n", (int)buf[k][i]);
//               printf("clry: %f\n", (double)clry[k][i]);
//            }

//        }

    for (i = 0; i < valid_num_scenes; i++)
    {
        if ((fmask_buf[i] < 2) && (id_range[i] == 1))
        {
            // remain the first element for replicated date
            if((n_clr > 0) && (valid_date_array[i] == clrx[n_clr - 1]))
                continue;
            else
            {
                clrx[n_clr] = valid_date_array[i];
                //printf("%d is %d\n", n_clr + 1, clrx[n_clr]);
                for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                {
                    if (k == 0)
                        clry[k][n_clr] = (float)buf_b[i];
                    else if (k == 1)
                        clry[k][n_clr] = (float)buf_g[i];
                    else if (k == 2)
                        clry[k][n_clr] = (float)buf_r[i];
                    else if (k == 3)
                        clry[k][n_clr] = (float)buf_n[i];
                    else if (k == 4)
                        clry[k][n_clr] = (float)buf_s1[i];
                    else if (k == 5)
                        clry[k][n_clr] = (float)buf_s2[i];
                    else if (k == 6)
                        clry[k][n_clr] = (float)(buf_t[i] * 10 - 27320);
                    //printf("%3.2f\n", clry[k][n_clr]);
                }
                n_clr++;
            }
        }
    }


//    for (k = 0; k < n_clr; k++)
//    {
//        printf("clrx %d: %d\n", k+1, (int)clrx[k]);
//    }


    temp_v_dif = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, n_clr,
                                     sizeof (float));
    if (temp_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating temp_v_dif memory",FUNC_NAME, FAILURE);
    }

    rec_v_dif = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, n_clr,
                                     sizeof (float));
    if (rec_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating rec_v_dif memory",FUNC_NAME, FAILURE);
    }

    rec_v_dif_copy = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, n_clr,
                                     sizeof (float));
    if (rec_v_dif_copy == NULL)
    {
        RETURN_ERROR ("Allocating rec_v_dif_copy memory",FUNC_NAME, FAILURE);
    }


    ids = (int *)calloc(n_clr, sizeof(int));
    if (ids == NULL)
    {
        RETURN_ERROR("ERROR allocating ids memory", FUNC_NAME, FAILURE);
    }
    ids_old = (int *)calloc(n_clr, sizeof(int));
    if (ids_old == NULL)
    {
        RETURN_ERROR("ERROR allocating ids_old memory", FUNC_NAME, FAILURE);
    }

    bl_ids = (int *)calloc(n_clr, sizeof(int));
    if (bl_ids == NULL)
    {
        RETURN_ERROR("ERROR allocating bl_ids memory", FUNC_NAME, FAILURE);
    }

    rm_ids = (int *)calloc(n_clr, sizeof(int));
    if (rm_ids == NULL)
    {
        RETURN_ERROR("ERROR allocating rm_ids memory", FUNC_NAME, FAILURE);
    }

    v_start = (float *)calloc(NUM_LASSO_BANDS, sizeof(float));
    if (v_start == NULL)
    {
        RETURN_ERROR("ERROR allocating v_start memory", FUNC_NAME, FAILURE);
    }

    v_end = (float *)calloc(NUM_LASSO_BANDS, sizeof(float));
    if (v_end == NULL)
    {
        RETURN_ERROR("ERROR allocating v_end memory", FUNC_NAME, FAILURE);
    }

    v_slope = (float *)calloc(NUM_LASSO_BANDS, sizeof(float));
    if (v_slope == NULL)
    {
        RETURN_ERROR("ERROR allocating v_slope memory", FUNC_NAME, FAILURE);
    }

    v_dif = (float *)calloc(NUM_LASSO_BANDS, sizeof(float));
    if (v_dif == NULL)
    {
        RETURN_ERROR("ERROR allocating v_dif memory", FUNC_NAME, FAILURE);
    }

    tmpcg_rmse = (float *)calloc(NUM_LASSO_BANDS, sizeof(float));
    if (tmpcg_rmse == NULL)
    {
        RETURN_ERROR("ERROR allocating tmpcg_rmse memory", FUNC_NAME, FAILURE);
    }



    end = n_clr;

//    /**************************************************************/
//    /*                                                            */
//    /* Remove repeated ids.                                       */
//    /*                                                            */
//    /**************************************************************/

//    matlab_unique(clrx, clry, n_clr, &end);
//    n_clr = end;
    /**************************************************************/
    /*                                                            */
    /* calculate variogram for each band and dates.                           */
    /*                                                            */
    /**************************************************************/
    status = adjust_median_variogram(clrx, clry, TOTAL_IMAGE_BANDS, 0, end-1, &date_vario,
                                     &max_date_difference, min_rmse, 1);
    if (status != SUCCESS)
    {
            RETURN_ERROR("ERROR calling median_variogram routine", FUNC_NAME,
                         FAILURE);
    }

    /* adjust T_cg based delta days*/
//    adj_conse = round (conse * 16 / (double)date_vario);
//    if (adj_conse < conse)
//        adj_conse = conse;

//     /* adjust conse based delta days*/
//    if(adj_conse > conse)
//    {
//        // adj_TCG = chi2inv(1 - pow(1 - PROB_T_CG, (double)conse / (double)adj_conse), NUM_LASSO_BANDS);
//        adj_TCG =X2(NUM_LASSO_BANDS, 1 - pow(1 - probability_threshold, (double)conse / (double)adj_conse));
//    }
//    else
//    {
//        adj_TCG = tcg;
//    }
    conse_min = conse;
    adj_TCG = tcg;

    v_dif_mag = (float **) allocate_2d_array(TOTAL_IMAGE_BANDS, conse_min * 2,
                sizeof (float));
    if (v_dif_mag == NULL)
    {
        RETURN_ERROR ("Allocating v_dif_mag memory",
                                 FUNC_NAME, FAILURE);
    }

    vec_mag = (float *)calloc(conse_min * 2, sizeof(float));
    if (vec_mag == NULL)
    {
        RETURN_ERROR ("Allocating vec_mag memory", FUNC_NAME, FAILURE);
    }

    v_diff = (float **) allocate_2d_array(NUM_LASSO_BANDS,
                      conse_min * 2, sizeof (float));
    if (v_diff == NULL)
    {
        RETURN_ERROR ("Allocating v_diff memory",
                      FUNC_NAME, FAILURE);
    }

    /**************************************************************/
    /*                                                            */
    /* Start with mininum requirement of clear obs.               */
    /*                                                            */
    /**************************************************************/

    i = N_TIMES * MIN_NUM_C;

    /**************************************************************/
    /*                                                            */
    /* The first observation for TSFit.                           */
    /*                                                            */
    /**************************************************************/

    i_start = 1;

    i_dense = 1;

    /**************************************************************/
    /*                                                            */
    /* Record the start of the model initialization               */
    /*     (0=>initial;1=>done)                                   */
    /*                                                            */
    /**************************************************************/

    bl_train = 0;

//    /**************************************************************/
//    /*                                                            */
//    /* initialize number of the functional curves                 */
//    /*                                                            */
//    /**************************************************************/


//    *num_fc = *num_fc +1;

    /**************************************************************/
    /*                                                            */
    /* Record the *num_fc at the beginning of each pixel.          */
    /*                                                            */
    /**************************************************************/
    rec_fc = *num_fc;


    /**************************************************************/
    /*                                                            */
    /* While loop - process til the last clear observation - adj_conse*/
    /*                                                            */
    /**************************************************************/
    // printf("dd_step1\n");
    while ((i <= end - conse_min)  && (clrx[end - 1] - clrx[i - 1 + 1] >= min_days_conse))  // i started from 1
    {
        // printf("i = %d; end = %d; clrx[end - 1] = %d; clrx[i] = %d. \n", i, end, clrx[end-1], clrx[i]);
        /**********************************************************/
        /*                                                        */
        /* span of "i"                                            */
        /*                                                        */
        /**********************************************************/

        i_span = i - i_start + 1;

        /**********************************************************/
        /*                                                        */
        /* span of time (num of years)                            */
        /*                                                        */
        /**********************************************************/
        time_span = (float)(clrx[i-1] - clrx[i_start-1]) / NUM_YEARS;

        /**********************************************************/
        /*                                                        */
        /* basic requrirements: 1) enough observations;           */
        /*                      2) enough time                    */
        /*                                                        */
        /**********************************************************/


//        tmp_direction = 255;
//        posi_count = 0;
//        nega_count = 0;
        if ((i_span >= N_TIMES * MIN_NUM_C) && (time_span >= (float)MIN_YEARS))
        {
            /******************************************************/
            /*                                                    */
            /* Initializing model.                                */
            /*                                                    */
            /******************************************************/

            if (bl_train == 0)
            {
                /*************************************************  */
                /*                                                  */
                /* Step 1: noise removal.                           */
                /*                                                  */
                /****************************************************/

                /*************************************************  */
                /*                                                  */
                /*  check maximum time gap as first                 */
                /*                                                  */
                /****************************************************/

                int max_date_diff = 0;
                for (k = i_start - 1; k < i - 1; k++)
                {
                    if (clrx[k + 1] - clrx[k] > max_date_diff)
                    {
                        max_date_diff = clrx[k + 1] - clrx[k];
                    }
                    //printf("%d \n", clrx[k]);
                }

                if (max_date_diff > NUM_YEARS)      //SY 09192018
                {
                    i++;
                    i_start++;
                    i_dense = i_start;
                    continue; //SY 02122019
                }

                status = auto_mask(clrx, clry, i_start-1, i+conse_min-1,
                               (float)(clrx[i+conse_min-1]-clrx[i_start-1]) / NUM_YEARS,
                               min_rmse[1], min_rmse[4], (float)T_CONST, bl_ids);
                // printf("ddstep2 auto_mask finished \n");
                if (status != SUCCESS)
                {
                    RETURN_ERROR("ERROR calling auto_mask during model initilization",
                                  FUNC_NAME, FAILURE);
                }

                /**************************************************/
                /*                                                */
                /* Clear the IDs buffers.                         */
                /*                                                */
                /**************************************************/

                for (k = 0; k < n_clr; k++)
                    ids[k] = 0;

                /**************************************************/
                /*                                                */
                /* IDs to be removed.                             */
                /*                                                */
                /**************************************************/

                for (k = i_start-1; k < i; k++)
                {
                    ids[k-i_start+1] = k;
                }
                m= 0;
                i_span = 0;
                for (k = 0; k < i-i_start + 1; k++) /** 02282019 SY **/
                {
                    if (bl_ids[k] == 1)
                    {
                        rm_ids[m] = ids[k];
                        //printf("%d \n", ids[k]);
                        m++;
                    }
                    else
                        i_span++;  /* update i_span after noise removal */
                }


                rm_ids_len = m;

                /**************************************************/
                /*                                                */
                /* Check if there are enough observation.         */
                /*                                                */
                /**************************************************/

                if (i_span < (N_TIMES * MIN_NUM_C))
                {
                    /**********************************************/
                    /*                                            */
                    /* Move forward to the i+1th clear observation*/
                    /*                                            */
                    /**********************************************/

                    i++;

                    /**********************************************/
                    /*                                            */
                    /* Not enough clear observations.             */
                    /*                                            */
                    /**********************************************/

                    continue;
                }

                if (end == 0)
                    RETURN_ERROR("No available data point", FUNC_NAME, FAILURE);

                /**************************************************/
                /*                                                */
                /* Allocate memory for cpx, cpy.                  */
                /*                                                */
                /**************************************************/

                cpx = (int* )malloc(end * sizeof(int));
                if (cpx == NULL)
                    RETURN_ERROR("ERROR allocating cpx memory", FUNC_NAME, FAILURE);

                cpy = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, end,
                                 sizeof (float));
                if (cpy == NULL)
                {
                    RETURN_ERROR ("Allocating cpy memory", FUNC_NAME, FAILURE);
                }

                /**************************************************/
                /*                                                */
                /* Remove noise pixels between i_start & i.       */
                /*                                                */
                /**************************************************/

                m = 0;
                for (k = 0, k_new=0; k < end; k++)
                {
                    if (m < rm_ids_len && k == rm_ids[m])
                    {
                        m++;
                        continue;
                    }
                    cpx[k_new] = clrx[k];
                    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                    {
                        cpy[b][k_new] = clry[b][k];
                    }
                    k_new++;
                }


                /**************************************************/
                /*                                                */
                /* Record i before noise removal.                 */
                /* This is very important, ie model is not yet    */
                /* initialized.   The multitemporal masking shall */
                /* be done again instead of removing outliers  In */
                /* every masking.                                 */
                /*                                                */
                /**************************************************/

                i_rec = i;

                /**************************************************/
                /*                                                */
                /* Update i afer noise removal.                   */
                /*     (i_start stays the same).                  */
                /*                                                */
                /**************************************************/

                i = i_start + i_span - 1;

                /**************************************************/
                /*                                                */
                /* Update span of time (num of years).            */
                /*                                                */
                /**************************************************/

                time_span=(cpx[i-1] - cpx[i_start-1]) / NUM_YEARS;

                /**************************************************/
                /*                                                */
                /* Check if there is enough time.                 */
                /*                                                */
                /**************************************************/

                if (time_span < MIN_YEARS)
                {
                    i = i_rec;   /* keep the original i */

                    /**********************************************/
                    /*                                            */
                    /* Move forward to the i+1th clear observation*/
                    /*                                            */
                    /**********************************************/

                    i++;
                    free(cpx);
                    status = free_2d_array ((void **) cpy);
                    if (status != SUCCESS)
                    {
                          RETURN_ERROR ("Freeing memory: cpy\n",
                                FUNC_NAME, FAILURE);
                    }
                    continue;    /* not enough time span */
                }

                //SY 09272018
                /**************************************************/
                /*                                                */
                /* updated end after checking if enought time_span*/
                /*                                                */
                /**************************************************/
                end = k_new;

                /**************************************************/
                /*                                                */
                /* Remove noise in original arrays.               */
                /*                                                */
                /**************************************************/

                for (k = 0; k < end; k++)
                {
                    clrx[k] = cpx[k];
                    for (m = 0; m < TOTAL_IMAGE_BANDS; m++)
                    {
                        clry[m][k] = cpy[m][k];
                    }
                }

//                    for (k = 0; k < n_clr; k++)
//                    {
//                        printf("clrx %d: %d\n", k+1, (int)clrx[k]);
//                    }

                free(cpx);
                status = free_2d_array ((void **) cpy);
                if (status != SUCCESS)
                {
                    RETURN_ERROR ("Freeing memory: cpy\n",
                         FUNC_NAME, FAILURE);
                }


                /**************************************************/
                /*                                                */
                /* Step 2) model fitting: initialize model testing*/
                /*         variables defining computed variables. */
                /*                                                */
                /**************************************************/

                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    /**********************************************/
                    /*                                            */
                    /* Initial model fit.                         */
                    /*                                            */
                    /**********************************************/
                    // printf("dd_step2 auto_ts_fit starts \n");
                    status = auto_ts_fit_float(clrx, clry,  i_b, i_b, i_start-1, i-1,
                             MIN_NUM_C, fit_cft, &rmse[i_b], rec_v_dif);
                    // printf("dd_step2 first using auto_ts_fit successed \n");
//                        for (k = 0; k < MAX_NUM_C; k++)
//                        {

//                               printf("%f\n", (float)fit_cft[i_b][k]);

//                        }
//                    printf("rmse for band %d is %f\n", i_b, rmse[i_b]);
                    // printf("auto_ts_fit finished \n");
                    if (status != SUCCESS)
                    {
                        RETURN_ERROR ("Calling auto_ts_fit_float during model initilization\n",
                             FUNC_NAME, FAILURE);
                    }



                }

                v_dif_norm = 0.0;

                for (i_b = 0; i_b < NUM_LASSO_BANDS; i_b++)
                {

                    /**********************************************/
                    /*                                            */
                    /* Calculate min. rmse.                       */
                    /*                                            */
                    /**********************************************/

                    mini_rmse = max((float)min_rmse[lasso_blist[i_b]], rmse[lasso_blist[i_b]]);

                    /**********************************************/
                    /*                                            */
                    /* Compare the first observation.             */
                    /*                                            */
                    /**********************************************/

                    v_start[i_b] = rec_v_dif[lasso_blist[i_b]][0]
                            / mini_rmse;

                    /**********************************************/
                    /*                                            */
                    /* Compare the last clear observation.        */
                    /*                                            */
                    /**********************************************/

                    v_end[i_b] = rec_v_dif[lasso_blist[i_b]][i-i_start]
                                            / mini_rmse;

                    /**********************************************/
                    /*                                            */
                    /* Anormalized slope values.                  */
                    /*                                            */
                    /**********************************************/
                    v_slope[i_b] = fit_cft[lasso_blist[i_b]][1] *
                                    (clrx[i-1]-clrx[i_start-1])/mini_rmse/SLOPE_SCALE;

                    /**********************************************/
                    /*                                            */
                    /* Difference in model intialization.         */
                    /*                                            */
                    /**********************************************/

                    v_dif[i_b] = fabs(v_slope[i_b]) + max(fabs(v_start[i_b]), fabs(v_end[i_b]));
                    v_dif_norm += v_dif[i_b] * v_dif[i_b];
                    //printf("%f \n", v_dif[i_b]);

                }

//                for(b = 0; b < TOTAL_IMAGE_BANDS; b++)
//                    printf("%.10f, %.10f, %.10f, %.10f,%.10f, %.10f, %.10f, %.10f\n", fit_cft[b][0],
//                            fit_cft[b][1], fit_cft[b][2], fit_cft[b][3], fit_cft[b][4], fit_cft[b][5], fit_cft[b][6], fit_cft[b][7]);


                /**************************************************/
                /*                                                */
                /* Find stable start for each curve.              */
                /*                                                */
                /**************************************************/
                if (v_dif_norm > adj_TCG)
                {
                    /**********************************************/
                    /*                                            */
                    /* Start from next clear observation.         */
                    /*                                            */
                    /**********************************************/

                    i_start++;

                    /**********************************************/
                    /*                                            */
                    /* Move forward to the i+1th clear observation*/
                    /*                                            */
                    /**********************************************/

                    i++;

                    /**********************************************/
                    /*                                            */
                    /* Keep all data and move to the next obs.    */
                    /*                                            */
                    /**********************************************/

                    continue;
                }

                /**************************************************/
                /*                                                */
                /* Model is ready.                                */
                /*                                                */
                /**************************************************/

                bl_train = 1;

                /**************************************************/
                /*                                                */
                /* Count difference of i for each iteration.      */
                /*                                                */
                /**************************************************/

                i_count = 0;

                /**************************************************/
                /*                                                */
                /* Find the previous break point.                 */
                /*                                                */
                /**************************************************/

                if (*num_fc == rec_fc)
                {
                    i_break = 1; /* first curve */
                }
                else
                {
                    /**********************************************/
                    /*                                            */
                    /* After the first curve, compare rmse to     */
                    /* determine which curve to determine t_break.*/
                    /*                                            */
                    /**********************************************/

                    for (k = 0; k < end; k++)
                    {
                        if (clrx[k] >= rec_cg[*num_fc-1].t_break)
                        {
                            i_break = k + 1;
                            break;
                        }
                    }
                }

                if (i_start > i_break)
                {
                    /**********************************************/
                    /*                                            */
                    /* Model fit at the beginning of the time     */
                    /* series.                                    */
                    /*                                            */
                    /**********************************************/

                    for(i_ini = i_start - 2; i_ini >= i_break - 1; i_ini--) // SY 09192018
                    {
                        if ((i_ini - (i_break - 1) + 1) < conse_min)
                        {
                            ini_conse = i_ini - (i_break - 1) + 1;
                        }
                        else
                        {
                            ini_conse = conse_min;
                        }

                        if (ini_conse == 0)
                        {
                            RETURN_ERROR ("No data point for model fit at "
                                  "the begining", FUNC_NAME, FAILURE);
                        }

                        /******************************************/
                        /*                                        */
                        /* Allocate memory for model_v_dif,       */
                        /* v_diff, vec_magg for the non-stdin     */
                        /* branch here.                           */
                        /*                                        */
                        /******************************************/

                        v_diff_tmp = (float **) allocate_2d_array(NUM_LASSO_BANDS,
                                    ini_conse, sizeof (float));
                        if (v_diff_tmp == NULL)
                        {
                            RETURN_ERROR ("Allocating v_diff_tmp memory",
                                           FUNC_NAME, FAILURE);
                        }

                        vec_magg = (float *) malloc(ini_conse * sizeof (float));
                        if (vec_magg == NULL)
                        {
                            RETURN_ERROR ("Allocating vec_magg memory",
                                          FUNC_NAME, FAILURE);
                        }

                        /******************************************/
                        /*                                        */
                        /* Detect change.                         */
                        /* value of difference for adj_conse      */
                        /* observations                          */
                        /* Record the magnitude of change.        */
                        /*                                        */
                        /******************************************/

                        vec_magg_min = 9999.0;
                        for (i_conse = 1; i_conse < ini_conse + 1; i_conse++) // SY 09192018
                        {
                            v_dif_norm = 0.0;
                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {

                                /**********************************/
                                /*                                */
                                /* Absolute differences.          */
                                /*                                */
                                /**********************************/

                                // SY 09192018 moving fitting into (i_b == lasso_blist[b])to save time //
                                // SY 02/13/2019 delete these speed-up modification as non-lasso bands
                                // are important for change agent classification
                                auto_ts_predict_float(clrx, fit_cft, MIN_NUM_C, i_b, i_ini-i_conse+1,
                                                i_ini-i_conse+1, &ts_pred_temp);
                                v_dif_mag[i_b][i_conse-1] = (float)clry[i_b][i_ini-i_conse+1] -
                                                   ts_pred_temp;// SY 09192018
                                // printf("auto_ts_predict finished \n");
                                /**********************************/
                                /*                                */
                                /* Normalize to z-score.          */
                                /*                                */
                                /**********************************/

                                for (b = 0; b < NUM_LASSO_BANDS; b++)
                                {
                                   if (i_b == lasso_blist[b])
                                   {
                                        /**************************/
                                        /*                        */
                                        /* Minimum rmse.          */
                                        /*                        */
                                        /**************************/

                                        mini_rmse = max((float)min_rmse[i_b], rmse[i_b]);

                                        /**************************/
                                        /*                        */
                                        /* z-scores.              */
                                        /*                        */
                                        /**************************/

                                        v_diff_tmp[b][i_conse-1] = v_dif_mag[i_b][i_conse-1] // SY 09192018
                                                                      / mini_rmse;
                                        v_dif_norm += v_diff_tmp[b][i_conse-1] * v_diff_tmp[b][i_conse-1]; // SY 09192018
            //                                        printf("i_b: %d\n",i_b);
//                                        printf("clry: %f\n", clry[i_b][i_ini-i_conse]);

//                                        printf("ts_pred_temp: %f\n",ts_pred_temp);
//                                        printf("v_dif_norm: %f\n",v_dif_norm);
//                                        printf("mini_rmse: %f\n",mini_rmse);
//                                        printf("v_dif_mag: %f\n",v_dif_mag[i_b][i_conse]);
                                    }
                                  }
                                }
                                vec_magg[i_conse-1] = v_dif_norm; // SY 09192018

                                if (vec_magg_min > vec_magg[i_conse-1])
                                {
                                    vec_magg_min =  vec_magg[i_conse-1]; // SY 09192018
                                }
                        }

                        /******************************************/
                        /*                                        */
                        /* Change angle.                      */
                        /*                                        */
                        /******************************************/

                        mean_angle = MeanAngl_float(v_diff_tmp, NUM_LASSO_BANDS, ini_conse);
                        /******************************************/
                        /*                                        */
                        /* Change detection.                      */
                        /*                                        */
                        /******************************************/
//                        for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
//                        {
//                            printf("%f, %f, %f, %f\n", fit_cft[b][0], fit_cft[b][1], fit_cft[b][2], fit_cft[b][3]);

//                        }
                        if (b_outputCM == TRUE){
                            if((ini_conse >= conse_min) & (clrx[i_ini] - clrx[i_break-1] >= min_days_conse)){  // only save CM based on >= conse obs
                            // if ((ini_conse >= adj_conse) && (*num_fc == rec_fc)){ // only save out the first curve
                                prob_angle = angle_decaying(mean_angle, (double)NSIGN, 135.0);
                                current_CM_n = (clrx[i_ini+1] - starting_date) / CM_OUTPUT_INTERVAL;  // looking back, t_break = current_i + 1
                                tmp = round(prob_angle * vec_magg_min * 100);
                                if (tmp > MAX_SHORT) // MAX_SHORT is upper limit of short 16
                                    tmp = MAX_SHORT;
                                tmp_CM = (short int) (tmp);
                                if(tmp_CM > CM_outputs[current_CM_n])
                                {
                                    CM_outputs[current_CM_n] = tmp_CM;
                                    CM_outputs_date[current_CM_n] = (short int)(clrx[i_ini+1] - ORDINAL_LANDSAT4_LAUNCH);
                                    // printf("date = %d\n", clrx[i_ini+1]);
                                }
                            }
                        }

                        if ((vec_magg_min > adj_TCG) && (mean_angle < NSIGN)) /* change detected */
                        {
                            free(vec_magg);
                            status = free_2d_array ((void **) v_diff_tmp);
                            if (status != SUCCESS)
                            {
                                RETURN_ERROR ("Freeing memory: v_diff_tmp\n",
                                              FUNC_NAME, FAILURE);
                            }
                            break;
                        }
                        else if (vec_magg[0] > T_MAX_CG) /* false change */
                        {
                            for (k = i_ini; k < end - 1; k++)
                            {
                                clrx[k] = clrx[k+1];
                                for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                                {
                                    clry[b][k] = clry[b][k+1];
                                }
                            }
                            i--;
                            end--;

                        }

                        /******************************************/
                        /*                                        */
                        /* Free the temporary memory.             */
                        /*                                        */
                        /******************************************/

                        free(vec_magg);
                        status = free_2d_array ((void **) v_diff_tmp);
                        if (status != SUCCESS)
                        {
                            RETURN_ERROR ("Freeing memory: v_diff_tmp\n",
                                          FUNC_NAME, FAILURE);
                        }

                        /**************************************/
                        /*                                    */
                        /* Update i_start if i_ini is not a   */
                        /* confirmed break.                   */
                        /*                                    */
                        /**************************************/

                        i_start = i_ini + 1;

                    } // end for (i_ini = i_start-1; i_ini >= i_break; i_ini--)

                }    // end for if (i_start > i_break)

                /**************************************************/
                /*                                                */
                /* Enough to fit simple model and confirm a break.*/
                /*                                                */
                /**************************************************/

                if ((*num_fc == rec_fc) && ((i_start - i_dense) >= LASSO_MIN))
                {
                    /**********************************************/
                    /*                                            */
                    /* Defining computed variables.               */
                    /*                                            */
                    /**********************************************/

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {

                        status = auto_ts_fit_float(clrx, clry, i_b, i_b,  i_dense-1, i_start-2,
                                 MIN_NUM_C, fit_cft, &rmse[i_b], temp_v_dif);//SY 02132019
                        if (status != SUCCESS)
                        {
                              RETURN_ERROR ("Calling auto_ts_fit_float with enough observations\n",
                                         FUNC_NAME, FAILURE);
                        }

                    }

                    /**********************************************/
                    /*                                            */
                    /* Record time of curve end,                  */
                    /* postion of the pixels.                     */
                    /*                                            */
                    /**********************************************/


                    rec_cg[*num_fc].t_end = clrx[i_start-2];
                    //rec_cg[*num_fc].pos.row = row;
                    //rec_cg[*num_fc].pos.col = col;

                    /**********************************************/
                    /*                                            */
                    /* Record break time, fit category, change    */
                    /* probability, time of curve start, number   */
                    /* of observations, change magnitude.         */
                    /*                                            */
                    /**********************************************/

                    rec_cg[*num_fc].t_break = clrx[i_start -1];
                    rec_cg[*num_fc].category = 10 + MIN_NUM_C;
                    rec_cg[*num_fc].change_prob = 100;
                    rec_cg[*num_fc].t_start = clrx[0];
                    rec_cg[*num_fc].num_obs = i_start - i_dense + 1;  //SY 09182018

//                    if ((i_start - 1 + adj_conse) < end)
//                        rec_cg[*num_fc].t_confirmed = clrx[i_start + adj_conse - 1];
//                    else
//                        rec_cg[*num_fc].t_confirmed = clrx[end - 1];

                    for (i_b = 0; i_b <  TOTAL_IMAGE_BANDS; i_b++)
                    {
                        quick_sort_float(v_dif_mag[i_b], 0, ini_conse-1);
                        matlab_2d_float_median(v_dif_mag, i_b, ini_conse,
                                              &v_dif_mean);
                        rec_cg[*num_fc].magnitude[i_b] = -v_dif_mean;
                    }

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        for (k = 0; k < LASSO_COEFFS; k++)
                        {
                            /**************************************/
                            /*                                    */
                            /* Record fitted coefficients.        */
                            /*                                    */
                            /**************************************/
                            if(k < MIN_NUM_C)
                                rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                            else
                                rec_cg[*num_fc].coefs[i_b][k] = 0;
                        }

                        /******************************************/
                        /*                                        */
                        /* Record rmse of the pixel.              */
                        /*                                        */
                        /******************************************/

                        rec_cg[*num_fc].rmse[i_b] = rmse[i_b];


                     }

                    /**********************************************/
                    /*                                            */
                    /* Identified and move on for the next        */
                    /* functional curve.                          */
                    /*                                            */
                    /**********************************************/

                    *num_fc = *num_fc + 1;
                    if (*num_fc >= NUM_FC)
                    {
                        /******************************************/
                        /*                                        */
                        /* Reallocate memory for rec_cg.          */
                        /*                                        */
                        /******************************************/

                        rec_cg = realloc(rec_cg, (*num_fc + 1) * sizeof(Output_t));
                        if (rec_cg == NULL)
                        {
                            RETURN_ERROR("ERROR allocating rec_cg memory",
                                         FUNC_NAME, FAILURE);
                        }
                    }
                }
            } /* end of initializing model */

            /******************************************************/
            /*                                                    */
            /* Allocate memory for v_diff for the non-stdin branch*/
            /*                                                    */
            /******************************************************/

//            for (k = 0; k < n_clr; k++)
//            {
//                printf("clrx %d: %d\n", k+1, (int)clrx[k]);
//            }


            /******************************************************/
            /*                                                    */
            /* Continuous monitoring started!!!                   */
            /*                                                    */
            /******************************************************/
            // printf("dd_step3 initialization finished\n");
            if (bl_train == 1)
            {
                // printf("processing %d obs finished \n", i);
                /**************************************************/
                /*                                                */
                /* Clears the IDs buffers.                        */
                /*                                                */
                /**************************************************/

                if(clrx[i + conse_min - 1] - clrx[i] < min_days_conse) //peek window starts from i (current pixel + 1) and ended with i + conse -1
                {
                    adj_conse = conse_min;
                    while(adj_conse < 2 * conse_min)
                    {
                        if (i + adj_conse - 1 >= end)
                            break;
                        if(clrx[i + adj_conse - 1] - clrx[i] >= min_days_conse)
                            break;

                        adj_conse = adj_conse + 1;
                    }
                    //adj_conse_record = adj_conse;
                    //adj_conse = round((double)adj_conse / CONSE) * CONSE;
                    //t_cg_adjust = X2(n_focus_variable, 1 - pow(1 - probability_threshold, (double)CONSE / (double)adj_conse));
                }
                else
                {
                    adj_conse = conse_min;
                    //t_cg_adjust = s_tcg;
                    //t_cg_outelier = S_TCG_MIN;
                }
                // printf("adj_conse = %d \n", adj_conse)

                for (k = 0; k < n_clr; k++)
                {
                    ids[k] = 0;
                }

                /**************************************************/
                /*                                                */
                /* All IDs.                                       */
                /*                                                */
                /**************************************************/

                ids_len = 0;
                for (k = i_start-1; k < i; k++)
                {
                    ids[k-i_start+1] = k;
                    ids_len++;
                }
                i_span = i - i_start +1;

                if(i_span_skip == 0)
                    i_span_skip = i_span;
                /**************************************************/
                /*                                                */
                /* Determine the time series model.               */
                /*                                                */
                /**************************************************/

                update_cft(i_span, N_TIMES, MIN_NUM_C, MID_NUM_C, MAX_NUM_C,
                          num_c, &update_num_c);

                /************************************************************/
                /*                                                          */
                /* initial model fit when there are not many observations.  */
                /* if (i_count == 0 || ids_old_len < (N_TIMES * MAX_NUM_C)) */
                /*                                                          */
                /************************************************************/

                if (i_count == 0 || i_span <= (N_TIMES * MAX_NUM_C))
                {
                    /**********************************************/
                    /*                                            */
                    /* update i_count at each iteration.          */
                    /*                                            */
                    /**********************************************/

                    i_count = clrx[i-1] - clrx[i_start-1];

                    pre_end = i;

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {

                        status = auto_ts_fit_float(clrx, clry, i_b, i_b, i_start-1, i-1, update_num_c,
                                             fit_cft, &rmse[i_b], rec_v_dif);
                        if (status != SUCCESS)
                        {
                            RETURN_ERROR ("Calling auto_ts_fit_float during continuous monitoring\n",
                                          FUNC_NAME, FAILURE);
                        }

                    }
                    // printf("auto_ts_fit_float finished \n", i);

                    /**********************************************/
                    /*                                            */
                    /* Updating information for the first         */
                    /* iteration.  Record time of curve start and */
                    /* time of curve end.                         */
                    /*                                            */
                    /**********************************************/

                    rec_cg[*num_fc].t_start = clrx[i_start-1];
                    rec_cg[*num_fc].t_end = clrx[i-1];

                    /**********************************************/
                    /*                                            */
                    /* No break at the moment.                    */
                    /*                                            */
                    /**********************************************/

                    rec_cg[*num_fc].t_break = 0;


                    /**********************************************/
                    /*                                            */
                    /* Record change probability, number of       */
                    /* observations, fit category.                */
                    /*                                            */
                    /**********************************************/

                    rec_cg[*num_fc].change_prob = 0;
                    rec_cg[*num_fc].num_obs = i-i_start+1;
                    rec_cg[*num_fc].category = 0 + update_num_c;

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {

                        /******************************************/
                        /*                                        */
                        /* Record rmse of the pixel.              */
                        /*                                        */
                        /******************************************/

                        rec_cg[*num_fc].rmse[i_b] = rmse[i_b];

                        /******************************************/
                        /*                                        */
                        /* Record change magnitude.               */
                        /*                                        */
                        /******************************************/
                        rec_cg[*num_fc].magnitude[i_b] = 0.0;

                        for (k = 0; k < LASSO_COEFFS; k++)
                        {
                            /**************************************/
                            /*                                    */
                            /* Record fitted coefficients.        */
                            /*                                    */
                            /**************************************/

                            rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                         }


                     }
                    // printf("auto_ts_fit2 finished \n", i);


                    /**********************************************/
                    /*                                            */
                    /* Detect change, value of difference for     */
                    /* adj_conse observations.                        */
                    /*                                            */
                    /**********************************************/
                    for (m = 0; m < conse_min * 2; m++)
                    {
                        vec_mag[m] = 0;
                        for (b = 0; b < NUM_LASSO_BANDS; b++)
                            v_diff[b][m] = 0;
                        for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                            v_dif_mag[b][m] = 0;
                    }


                    for (i_conse = 1; i_conse < conse_min+1; i_conse++)//SY 09192018
                    {
                        v_dif_norm = 0.0;
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        {
                            /**************************************/
                            /*                                    */
                            /* Absolute differences.              */
                            /*                                    */
                            /**************************************/
                            // printf("auto_ts_predict started finished \n", i);
                            auto_ts_predict_float(clrx, fit_cft, update_num_c, i_b, i+i_conse-1, i+i_conse-1,
                                            &ts_pred_temp); //SY 09192018
                            // printf("auto_ts_predict finished \n", i);
                            v_dif_mag[i_b][i_conse-1] = (float)clry[i_b][i+i_conse-1] - ts_pred_temp;//SY 09192018

                            /**************************************/
                            /*                                    */
                            /* Normalize to z-score.              */
                            /*                                    */
                            /**************************************/

                            for (b = 0; b < NUM_LASSO_BANDS; b++)
                            {
                                if (i_b == lasso_blist[b])
                                {
                                    /******************************/
                                    /*                            */
                                    /* Minimum rmse,              */
                                    /* z-scores.                  */
                                    /*                            */
                                    /******************************/

                                    mini_rmse = max(min_rmse[i_b], rmse[i_b]);
                                    v_diff[b][i_conse-1] = v_dif_mag[i_b][i_conse-1] / mini_rmse;
                                    v_dif_norm += v_diff[b][i_conse-1] * v_diff[b][i_conse-1];
                                }
                            }
                        }
                        vec_mag[i_conse-1] = v_dif_norm;//SY 09192018
                    }

                    /**********************************************/
                    /*                                            */
                    /* Clears the IDs_old buffers.                */
                    /*                                            */
                    /**********************************************/

                    for (k = 0; k < ids_len; k++)
                    {
                        ids_old[k] = 0;
                    }

                    /**********************************************/
                    /*                                            */
                    /* IDs that have not been updated.            */
                    /*                                            */
                    /**********************************************/

                    for (k = 0; k < ids_len; k++)
                    {
                        ids_old[k] = ids[k];
                    }
                    ids_old_len = ids_len;

                } // end for if (i_count == 0 || i_span <= (N_TIMES * MAX_NUM_C))
                else // update frequency
                {
                    if((i - pre_end >= UPDATE_FREQ) && (i - pre_end) > (int)(i_span_skip * SKIP_PERCENTAGE))
                    {
                        /******************************************/
                        /*                                        */
                        /* Update coefficent at each iteration year. */
                        /*                                        */
                        /******************************************/

                        i_count = clrx[i-1] - clrx[i_start-1];
                        pre_end = i;
                        // update i_span_skip
                        i_span_skip = i_span;
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        {
                            status = auto_ts_fit_float(clrx, clry, i_b, i_b, i_start-1, i-1, update_num_c,
                                             fit_cft, &rmse[i_b], rec_v_dif);
                            // printf("auto_ts_fit2 finished \n", i);
                            if (status != SUCCESS)
                            {
                                RETURN_ERROR ("Calling auto_ts_fit_float for change detection with "
                                     "enough observations\n", FUNC_NAME, FAILURE);
                            }


                        }

                        /******************************************/
                        /*                                        */
                        /* Record fitted coefficients.            */
                        /*                                        */
                        /******************************************/
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        {

                            for (k = 0; k < LASSO_COEFFS; k++)
                            {
                                /**********************************/
                                /*                                */
                                /* Record fitted coefficients.    */
                                /*                                */
                                /**********************************/

                                rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                            }
                            /**************************************/
                            /*                                    */
                            /* Record rmse of the pixel.          */
                            /*                                    */
                            /**************************************/

                            rec_cg[*num_fc].rmse[i_b] = rmse[i_b];

                            /******************************************/
                            /*                                        */
                            /* Record number of observations, fit     */
                            /* category.                              */
                            /*                                        */
                            /******************************************/

                            rec_cg[*num_fc].num_obs = i - i_start + 1;
                            rec_cg[*num_fc].category = 0 + update_num_c;
                        }

                        /******************************************/
                        /*                                        */
                        /* Clears the IDs_Old buffers.            */
                        /*                                        */
                        /******************************************/

                        for (k = 0; k < ids_len; k++)
                        {
                            ids_old[k] = 0;
                        }

                        /******************************************/
                        /*                                        */
                        /* IDs that have not been updated.        */
                        /*                                        */
                        /******************************************/

                        for (k = 0; k < ids_len; k++)
                        {
                            ids_old[k] = ids[k];
                        }
                        ids_old_len = ids_len;

                    } //  if(i_update == UPDATE_FREQ || i_update > UPDATE_FREQ)


                    /**********************************************/
                    /*                                            */
                    /* Record time of curve end.                  */
                    /*                                            */
                    /**********************************************/

                    rec_cg[*num_fc].t_end = clrx[i-1];

                    /**********************************************/
                    /*                                            */
                    /* Use fixed number for RMSE computing.       */
                    /*                                            */
                    /**********************************************/

                    n_rmse = N_TIMES * rec_cg[*num_fc].category;

                    /**********************************************/
                    /*                                            */
                    /* Better days counting for RMSE calculating  */
                    /* relative days distance.                    */
                    /*                                            */
                    /**********************************************/

                    if (ids_old_len == 0)
                    {
                        RETURN_ERROR ("No data points for RMSE calculating",
                                     FUNC_NAME, FAILURE);
                    }

                    d_yr = malloc(ids_old_len * sizeof(float));
                    if (d_yr == NULL)
                    {
                        RETURN_ERROR ("Allocating d_yr memory",
                                     FUNC_NAME, FAILURE);
                    }

                    for(m = 0; m < ids_old_len; m++)
                    {
                        d_rt = clrx[ids_old[m]] - clrx[i+conse_min-1];
                        d_yr[m] = fabs(round((float)d_rt / NUM_YEARS) * NUM_YEARS - (float)d_rt);
                    }

                    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                    {
                        for (m = 0; m < ids_old_len; m++)
                            rec_v_dif_copy[b][m] = rec_v_dif[b][m];
                    }

                    /**********************************************/
                    /*                                            */
                    /* Sort the rec_v_dif based on d_yr.          */
                    /*                                            */
                    /**********************************************/

                    quick_sort_2d_float(d_yr, rec_v_dif_copy, 0, ids_old_len-1, TOTAL_IMAGE_BANDS);
                    for(b = 0; b < NUM_LASSO_BANDS; b++)
                        tmpcg_rmse[b] = 0.0;

                    /**********************************************/
                    /*                                            */
                    /* Temporarily changing RMSE.                 */
                    /*                                            */
                    /**********************************************/

                    for (b = 0; b < NUM_LASSO_BANDS; b++)
                    {
                        matlab_2d_array_norm_float(rec_v_dif_copy, lasso_blist[b], n_rmse,
                                         &tmpcg_rmse[b]);
                        tmpcg_rmse[b] /= sqrtf((float)(n_rmse - rec_cg[*num_fc].category));
                    }

                    /**********************************************/
                    /*                                            */
                    /* Free allocated memories.                   */
                    /*                                            */
                    /**********************************************/
                    free(d_yr);

                    /**********************************************/
                    /*                                            */
                    /* Move the ith col to i-1th col.             */
                    /*                                            */
                    /**********************************************/

                    for (m = 0; m < conse_min-1; m++)
                    {
                        vec_mag[m] = vec_mag[m+1];
                        for (b = 0; b < NUM_LASSO_BANDS; b++)
                            v_diff[b][m] = v_diff[b][m+1];
                        for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                            v_dif_mag[b][m] = v_dif_mag[b][m+1];
                    }

                    for (b = 0; b < NUM_LASSO_BANDS; b++)
                        v_diff[b][conse_min-1] = 0.0;
                    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                        v_dif_mag[b][conse_min-1] = 0.0;
                    vec_mag[conse_min-1] = 0.0;

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {

                        /******************************************/
                        /*                                        */
                        /* .     */
                        /*                                        */
                        /******************************************/
    //                        if (i == 45){
    //                            i = 45;
    //                       }
                        auto_ts_predict_float(clrx, fit_cft, update_num_c, i_b, i+conse_min-1,
                             i+conse_min-1, &ts_pred_temp);
                        v_dif_mag[i_b][conse_min-1] = (float)clry[i_b][i+conse_min-1] - ts_pred_temp;

                        /******************************************/
                        /*                                        */
                        /* Normalized to z-scores.                */
                        /*                                        */
                        /******************************************/
                        for (b = 0; b < NUM_LASSO_BANDS; b++)
                        {
                            if (i_b == lasso_blist[b])
                            {
                                /**********************************/
                                /*                                */
                                /* Minimum rmse.                  */
                                /*                                */
                                /**********************************/

                                mini_rmse = max((double)min_rmse[i_b], tmpcg_rmse[b]);

                                /**********************************/
                                /*                                */
                                /* Z-score.                       */
                                /*                                */
                                /**********************************/

                                v_diff[b][conse_min-1] = v_dif_mag[i_b][conse_min-1] / mini_rmse;
                                vec_mag[conse_min-1] += v_diff[b][conse_min-1] * v_diff[b][conse_min-1]; // SY 02132014
                            }
                        }
                     }

                    /* padding procedure */
                    if(adj_conse > conse_min)
                    {
                        d_yr = malloc(ids_old_len * sizeof(float));
                        if (d_yr == NULL)
                        {
                            RETURN_ERROR ("Allocating d_yr memory",
                                         FUNC_NAME, FAILURE);
                        }

                        for (i_conse = conse_min; i_conse < adj_conse; i_conse++)
                        {
                            for(m = 0; m < ids_old_len; m++)
                            {
                                d_rt = clrx[ids_old[m]] - clrx[i+ i_conse];
                                d_yr[m] = fabs(round((float)d_rt / NUM_YEARS) * NUM_YEARS - (float)d_rt);
                            }

                            for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                            {
                                for (m = 0; m < ids_old_len; m++)
                                    rec_v_dif_copy[b][m] = rec_v_dif[b][m];
                            }

                            /**********************************************/
                            /*                                            */
                            /* Sort the rec_v_dif based on d_yr.          */
                            /*                                            */
                            /**********************************************/

                            quick_sort_2d_float(d_yr, rec_v_dif_copy, 0, ids_old_len-1, TOTAL_IMAGE_BANDS);

                            /**********************************************/
                            /*                                            */
                            /* Temporarily changing RMSE.                 */
                            /*                                            */
                            /**********************************************/

                            for (b = 0; b < NUM_LASSO_BANDS; b++)
                            {
                                matlab_2d_array_norm_float(rec_v_dif_copy, lasso_blist[b], n_rmse,
                                                 &tmpcg_rmse[b]);
                                tmpcg_rmse[b] /= sqrt((float)(n_rmse - rec_cg[*num_fc].category));
                            }

                            vec_mag[i_conse] = 0.0;
                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {
                                auto_ts_predict_float(clrx, fit_cft, update_num_c, i_b, i+i_conse,
                                     i+i_conse, &ts_pred_temp);
                                v_dif_mag[i_b][i_conse] = (float)clry[i_b][i+i_conse] - ts_pred_temp;

                                for (b = 0; b < NUM_LASSO_BANDS; b++)
                                {
                                    if (i_b == lasso_blist[b])
                                    {

                                        mini_rmse = max((float)min_rmse[i_b], tmpcg_rmse[b]);

                                        v_diff[b][i_conse] = v_dif_mag[i_b][i_conse] / mini_rmse;
                                        vec_mag[i_conse] += v_diff[b][i_conse] * v_diff[b][i_conse]; // SY 02132014
                                     }
                                 }
                             }
                        } // for (i_conse = conse_min; i_conse < adj_conse; i_conse++)
                        /**********************************************/
                        /*                                            */
                        /* Free allocated memories.                   */
                        /*                                            */
                        /**********************************************/
                        free(d_yr);
                    }
                } // else update frequency

                mean_angle = MeanAngl_float(v_diff, NUM_LASSO_BANDS, adj_conse);


                break_mag = 9999.0;
                for (m = 0; m < adj_conse; m++)
                {
                    if (break_mag > vec_mag[m])
                    {
                        break_mag = vec_mag[m];
                    }
                }

//                if(i == 532 - 4){
//                  for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
//                    for(j = 0; j < MAX_NUM_C; j++)
//                        printf("%f\n", fit_cft[k][j]);

//                  for(m = 0; m < n_rmse; m++)
//                  {

//                      printf("%f\n", d_yr[m]);
//                      printf("%f\n", rec_v_dif_copy[0][m]);

//                  }
//                }

                if (b_outputCM == TRUE)
                {
                    prob_angle = angle_decaying(mean_angle, (double)NSIGN, 135.0);
                    // prob_MCM = Chi_Square_Distribution(break_mag, NUM_LASSO_BANDS);
                    current_CM_n = (clrx[i] - starting_date) / CM_OUTPUT_INTERVAL;
                    tmp = round(prob_angle * break_mag * 100);
                    if (tmp > MAX_SHORT) // MAX_SHORT is upper limit of short 16
                        tmp = MAX_SHORT;
                    tmp_CM = (short int) (tmp);

                    /*********************************************/
                    /*      change direction by majority vote    */
                    /*********************************************/
                    if(tmp_CM > CM_outputs[current_CM_n])
                    {
                        CM_outputs[current_CM_n] = tmp_CM;
                        CM_outputs_date[current_CM_n] = (short int)(clrx[i] - ORDINAL_LANDSAT4_LAUNCH);
                    }
                }

//                if (clrx[i] > 731426 - 1){
//                    printf("clry[4][i + adj_conse] = %f\n", clry[4][i + adj_conse - 1]);
//                }
                if ((break_mag > adj_TCG) && (mean_angle < NSIGN))
                {

//                    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
//                        for(j = 0; j < MAX_NUM_C; j++)
//                            printf("%f\n", fit_cft[k][j]);
//                    if (verbose)
//                    {
//                        printf("Change Magnitude = %.2f\n", break_mag - adj_TCG);
//                    }

                    /**********************************************/
                    /*                                            */
                    /* Record break time.                        */
                    /*                                            */
                    /**********************************************/

                    rec_cg[*num_fc].t_break = clrx[i];
                    rec_cg[*num_fc].change_prob = 100;
//                    if ((i_start - 1 + adj_conse) < end)
//                        rec_cg[*num_fc].t_confirmed = clrx[i + adj_conse - 1];
//                    else
//                        rec_cg[*num_fc].t_confirmed = clrx[end - 1];

                    //rec_cg[*num_fc].t_confirmed = clrx[i + adj_conse - 1];

                    /* check if it is a lasso band               */
                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        quick_sort_float(v_dif_mag[i_b], 0, adj_conse-1);
                        matlab_2d_float_median(v_dif_mag, i_b, adj_conse,
                                               &rec_cg[*num_fc].magnitude[i_b]);
                    }
                    /**********************************************/
                    /*                                            */
                    /* Identified and move on for the next        */
                    /* functional curve.                          */
                    /*                                            */
                    /**********************************************/

                    *num_fc = *num_fc + 1;

                    if (*num_fc >= NUM_FC)
                    {
                        /******************************************/
                        /*                                        */
                        /* Reallocate memory for rec_cg.          */
                        /*                                        */
                        /******************************************/

                        rec_cg = realloc(rec_cg, (*num_fc + 1) * sizeof(Output_t));
                        if (rec_cg == NULL)
                        {
                            RETURN_ERROR("ERROR allocating rec_cg memory",
                                                 FUNC_NAME, FAILURE);
                        }
                    }

                    /**********************************************/
                    /*                                            */
                    /* Start from i+1 for the next functional     */
                    /* curve.                                     */
                    /*                                            */
                    /**********************************************/

                    i_start = i + 1;

                    /**********************************************/
                    /*                                            */
                    /* Start training again.                      */
                    /*                                            */
                    /**********************************************/

                    bl_train = 0;

                    i_span_skip = 0;

                }

                else if (vec_mag[0] > T_MAX_CG)  /*false change*/
                {
                    /**********************************************/
                    /*                                            */
                    /* Remove noise.                              */
                    /*                                            */
                    /**********************************************/

                    for (m = i; m < end -1; m++)
                    {
                        clrx[m] = clrx[m+1];
                        for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                            clry[b][m] = clry[b][m+1];
                    }

                    i--;   /* stay & check again after noise removal */
                    end--;
                }

           } /* end of continuous monitoring */
        }  /* end of checking basic requrirements */

        /**********************************************************/
        /*                                                        */
        /* Move forward to the i+1th clear observation.           */
        /*                                                        */
        /**********************************************************/

        i++;

    } /* end of "while (i <= end - adj_conse) */

    /**************************************************************/
    /*                                                            */
    /* Two ways for processing the end of the time series.        */
    /*                                                            */
    /**************************************************************/
    // printf("main part finished \n");
    if (bl_train == 1)
    {

        /**********************************************************/
        /*                                                        */
        /* If no break, find at the end of the time series,       */
        /* define probability of change based on adj_conse.           */
        /*                                                        */
        /**********************************************************/
        adj_conse_end = end - i + 1;

        /* padding peek window observation */
        if(adj_conse_end > conse_min)
        {
            d_yr = malloc(ids_old_len * sizeof(float));
            if (d_yr == NULL)
            {
                RETURN_ERROR ("Allocating d_yr memory",
                             FUNC_NAME, FAILURE);
            }

            for (i_conse = conse_min; i_conse < adj_conse_end; i_conse++)
            {
                for(m = 0; m < ids_old_len; m++)
                {
                    d_rt = clrx[ids_old[m]] - clrx[i+ i_conse];
                    d_yr[m] = fabs(round((float)d_rt / NUM_YEARS) * NUM_YEARS - (float)d_rt);
                }

                for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                {
                    for (m = 0; m < ids_old_len; m++)
                        rec_v_dif_copy[b][m] = rec_v_dif[b][m];
                }

                /**********************************************/
                /*                                            */
                /* Sort the rec_v_dif based on d_yr.          */
                /*                                            */
                /**********************************************/

                quick_sort_2d_float(d_yr, rec_v_dif_copy, 0, ids_old_len-1, TOTAL_IMAGE_BANDS);
                for(b = 0; b < NUM_LASSO_BANDS; b++)
                    tmpcg_rmse[b] = 0.0;

                /**********************************************/
                /*                                            */
                /* Temporarily changing RMSE.                 */
                /*                                            */
                /**********************************************/

                for (b = 0; b < NUM_LASSO_BANDS; b++)
                {
                    matlab_2d_array_norm_float(rec_v_dif_copy, lasso_blist[b], n_rmse,
                                     &tmpcg_rmse[b]);
                    tmpcg_rmse[b] /= sqrt((float)(n_rmse - rec_cg[*num_fc].category));
                }

                vec_mag[i_conse] = 0;
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    auto_ts_predict_float(clrx, fit_cft, update_num_c, i_b, i+i_conse,
                         i+i_conse, &ts_pred_temp);
                    v_dif_mag[i_b][i_conse] = (float)clry[i_b][i_conse] - ts_pred_temp;

                    for (b = 0; b < NUM_LASSO_BANDS; b++)
                    {
                        if (i_b == lasso_blist[b])
                        {

                            mini_rmse = max((double)min_rmse[i_b], tmpcg_rmse[b]);

                            v_diff[b][i_conse] = v_dif_mag[i_b][i_conse] / mini_rmse;
                            vec_mag[i_conse] += v_diff[b][i_conse] * v_diff[b][i_conse]; // SY 02132014
                         }
                     }
                 }
            } // for (i_conse = conse_min; i_conse < adj_conse; i_conse++)
            /**********************************************/
            /*                                            */
            /* Free allocated memories.                   */
            /*                                            */
            /**********************************************/
            free(d_yr);
        }



        id_last = adj_conse_end;
        for (i_conse = adj_conse_end - 1; i_conse >= 0; i_conse--)
        {
            v_diff_tmp =(float **) allocate_2d_array(NUM_LASSO_BANDS, adj_conse_end - i_conse, sizeof (float));
            for (k = 0; k < NUM_LASSO_BANDS; k++)
                for(j = 0; j < adj_conse - i_conse; j++)
                  v_diff_tmp[k][j] = v_diff[k][i_conse+j];

            mean_angle = MeanAngl_float(v_diff_tmp, NUM_LASSO_BANDS, adj_conse_end - i_conse);
            //float tt = vec_mag[i_conse];
            if ((vec_mag[i_conse] <= adj_TCG)||(mean_angle >= NSIGN))
            {
                /**************************************************/
                /*                                                */
                /* The last stable ID.                            */
                /*                                                */
                /**************************************************/

                id_last = i_conse + 1;
                status = free_2d_array ((void **) v_diff_tmp);
                if (status != SUCCESS)
                {
                    RETURN_ERROR ("Freeing memory: v_diff_tmp\n",
                                  FUNC_NAME, FAILURE);
                }
                break;
            }
            status = free_2d_array ((void **) v_diff_tmp);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: v_diff_tmp\n",
                              FUNC_NAME, FAILURE);
            }
        }

        /**********************************************************/
        /*                                                        */
        /* Update change probability, end time of the curve.      */
        /*                                                        */
        /**********************************************************/

        rec_cg[*num_fc].change_prob = (short int)((double)(adj_conse_end - id_last) * 100.0 / (double)adj_conse_end) ;
//        rec_cg[*num_fc].t_confirmed = 0;
        rec_cg[*num_fc].t_end = clrx[end -1 - adj_conse_end + id_last];    // 11/18/2018 SY

        /**********************************************************/
        /*                                                        */
        /* Mean value fit for the rest of the pixels < adj_conse & > 1*/
        /*                                                        */
        /**********************************************************/

        if (adj_conse_end > id_last)
        {
            /******************************************************/
            /*                                                    */
            /* Update time of the probable change.                */
            /*                                                    */
            /******************************************************/

            rec_cg[*num_fc].t_break = clrx[end - adj_conse_end+id_last];

            /******************************************************/
            /*                                                    */
            /* Update magnitude of change.                        */
            /*                                                    */
            /******************************************************/

            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
                quick_sort_float(v_dif_mag[i_b], id_last, adj_conse_end-1);
                //printf("%f\n", v_dif_mag[i_b][adj_conse-1]);
                matlab_float_2d_partial_median(v_dif_mag, i_b, id_last, adj_conse_end-1,
                                     &rec_cg[*num_fc].magnitude[i_b]);
            }
         }
    }

    else if (bl_train == 0)
    {
        /**********************************************************/
        /*                                                        */
        /* If break found close to the end of the time series,    */
        /* use [adj_conse,MIN_NUM_C*N_TIMES+adj_conse) to fit curve.      */
        /*                                                        */
        /* Update i_start.                                        */
        /*                                                        */
        /**********************************************************/
        if (*num_fc == rec_fc)
        {
            /******************************************************/
            /*                                                    */
            /* First curve.                                       */
            /*                                                    */
            /******************************************************/

            i_start = 1;
        }
        else
        {
            for (k = 0; k < n_clr; k++)
            {
                 if (clrx[k] >= rec_cg[*num_fc-1].t_break)
                 {
                       i_start = k + 1;
                       break;
                 }
             }
        }

        for (m = 0; m < n_clr; m++)
        {
            bl_ids[m] = 0;
        }

        if ((end - i_start + 1) >= LASSO_MIN) //04/02/2019 change adj_conse to CONSE_END
        {
            /******************************************************/
            /*                                                    */
            /* Multitemporal cloud mask.                          */
            /*                                                    */
            /******************************************************/

            status = auto_mask(clrx, clry, i_start-1, end-1,
                           (float)(clrx[end-1]-clrx[i_start-1]) / NUM_YEARS,
                           min_rmse[1], min_rmse[4], (float)T_CONST, bl_ids);
            if (status != SUCCESS)
                RETURN_ERROR("ERROR calling auto_mask at the end of time series",
                              FUNC_NAME, FAILURE);

            /******************************************************/
            /*                                                    */
            /* Clears the IDs buffers.                            */
            /*                                                    */
            /******************************************************/

            for (m = 0; m < n_clr-1; m++)
            {
                ids[m] = 0;
            }

            /******************************************************/
            /*                                                    */
            /* IDs to be removed.                                 */
            /*                                                    */
            /******************************************************/

            for (k = i_start-1; k < end; k++)
            {
                ids[k-i_start+1] = k;
            }
            m= 0;
            i_span = 0;
            for (k = 0; k < end-i_start+1; k++)
            {
                if (bl_ids[k] == 1)
                {
                    rm_ids[m] = ids[k];
                    m++;
                }
                else
                    i_span++;  /* update i_span after noise removal */
            }
            rm_ids_len = m;

            /******************************************************/
            /*                                                    */
            /* Remove noise pixels between i_start & i .          */
            /*                                                    */
            /******************************************************/

            m = 0;
            for (k = 0,k_new=0; k < end; k++) /* 03192019 SY*/
            {
                if (m < rm_ids_len && k == rm_ids[m])
                {
                    m++;
                    continue;
                }
                clrx[k_new] = clrx[k];
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    clry[i_b][k_new] = clry[i_b][k];
                k_new++;
            }
            end = k_new;
         }

         if ((end - i_start + 1) >= LASSO_MIN)     // 09/28/2018 SY delete equal sign //11/15/2018 put back equal sign
         {
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
//                    for(k=i_start-1;k<end;k++)
//                    {
//                        printf("clrx %d: %d\n", k,clrx[k]);
//                        printf("clry %d: %f\n", k,clry[i_b][k]);
//                    }
                status = auto_ts_fit_float(clrx, clry, i_b, i_b, i_start-1, end-1, MIN_NUM_C,
                                     fit_cft, &rmse[i_b], temp_v_dif);
                if (status != SUCCESS)
                {
                     RETURN_ERROR ("Calling auto_ts_fit_float at the end of time series\n",
                            FUNC_NAME, FAILURE);
                }

            }

            /******************************************************/
            /*                                                    */
            /* Record time of curve start, time of curve end,     */
            /* break time, postion of the pixel.                  */
            /*                                                    */
            /******************************************************/

            if (*num_fc == rec_fc)
            {
                    rec_cg[*num_fc].t_start = clrx[0];
            }
            else
            {
                    rec_cg[*num_fc].t_start = rec_cg[*num_fc-1].t_break;
            }
            rec_cg[*num_fc].t_end = clrx[end-1];
            rec_cg[*num_fc].t_break = 0;
            //rec_cg[*num_fc].pos.row = row;
            //rec_cg[*num_fc].pos.col = col;




            /******************************************************/
            /*                                                    */
            /* Record change probability, number of observations, */
            /* fit category.                                      */
            /*                                                    */
            /******************************************************/

            rec_cg[*num_fc].change_prob = 0;
//            rec_cg[*num_fc].t_confirmed = 0;
            rec_cg[*num_fc].num_obs = end - i_start + 1;
            rec_cg[*num_fc].category = 20 + MIN_NUM_C; /* simple model fit at the end */

            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {

                /******************************************************/
                /*                                                    */
                /* Record fitted coefficients.                        */
                /*                                                    */
                /******************************************************/

                for (k = 0; k < LASSO_COEFFS; k++)
                {
                    rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                }
                rec_cg[*num_fc].rmse[i_b] = rmse[i_b];

                /******************************************************/
                /*                                                    */
                /* Record change magnitude.                           */
                /*                                                    */
                /******************************************************/

                rec_cg[*num_fc].magnitude[i_b] = 0.0;


             }

            //*num_fc = *num_fc + 1;
            if (*num_fc >= NUM_FC)
            {
                /**************************************************/
                /*                                                */
                /* Reallocate memory for rec_cg.                  */
                /*                                                */
                /**************************************************/

                rec_cg = realloc(rec_cg, (*num_fc + 1) * sizeof(Output_t));
                if (rec_cg == NULL)
                {
                    RETURN_ERROR("ERROR allocating rec_cg memory",
                                 FUNC_NAME, FAILURE);
                }
            }
        }else{
           *num_fc = *num_fc - 1;
        }
    }

    *num_fc = *num_fc + 1;


    /******************************************************************/
    /*                                                                */
    /* Free memory allocations for this section.                      */
    /*                                                                */
    /******************************************************************/

    // for debug
    // printf("free stage 1 \n");
    status = free_2d_array ((void **) fit_cft);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fit_cft\n", FUNC_NAME,
                      FAILURE);
    }

    free(rmse);

    free(clrx);

    // for debug
    // printf("free stage 1 \n");

    status = free_2d_array ((void **) clry);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: clry\n", FUNC_NAME,
                      FAILURE);
    }

    status = free_2d_array ((void **) temp_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: temp_v_dif\n",
                      FUNC_NAME, FAILURE);
    }

    // for debug
    // printf("free stage 2 \n");

    status = free_2d_array ((void **) rec_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: rec_v_dif\n",
                      FUNC_NAME, FAILURE);
    }
    status = free_2d_array ((void **) rec_v_dif_copy);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: \n",
                      FUNC_NAME, FAILURE);
    }


    // for debug
    // printf("free stage 3 \n");

    free(ids);
    free(ids_old);
    free(bl_ids);
    free(rm_ids);


    free(v_start);
    free(v_end);
    free(v_slope);
    free(v_dif);

    free(tmpcg_rmse);

    // for debug
    // printf("free stage 4 \n");

    status = free_2d_array ((void **) v_dif_mag);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: v_dif_mag\n",
                   FUNC_NAME, FAILURE);
    }
    free(vec_mag);

    status = free_2d_array ((void **) v_diff);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: v_diff\n",
                      FUNC_NAME, FAILURE);
    }

    // for debug
    // printf("free stage 5 \n");

    /******************************************************************/
    /*                                                                */
    /* Output rec_cg structure to the output file.                    */
    /* Note: can use fread to read out the structure from the output  */
    /* file.                                                          */
    /* If output was stdout, skip this step.                          */
    /*                                                                */
    /******************************************************************/


    return (SUCCESS);
}



/******************************************************************************
MODULE:  stand_procedure
PURPOSE:  standard procedure when having enough clear pixels
RETURN VALUE:
Type = int (SUCCESS OR FAILURE)
HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
11/14/2018   Su Ye          Modification from main function in original CCDC.c
******************************************************************************/

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
    bool b_outputCM,              /* I: indicate if cold is running as the first step of object-based cold*/
    int starting_date,           /* I: the starting date of the whole dataset to enable reconstruct CM_date, all pixels for a tile should have the same date, only for b_outputCM is True */
    Output_t *rec_cg,           /* O: outputted structure for CCDC results     */
    int *num_fc,                /* O: number of fitting curves                       */
    int CM_OUTPUT_INTERVAL,
    short int* CM_outputs,      /* I/O: maximum change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
    short int* CM_outputs_date,      /* I/O: dates for maximum change magnitudes at every CM_OUTPUT_INTERVAL days, only for b_outputCM is True*/
    bool b_c2,
    double gap_days
)
{
    int status;
    int i, j, k, k_new, b;
    int m;                       /*the number of ID to be removed */
    char FUNC_NAME[] = "stand_procedure";
    //char msg_str[MAX_STR_LEN];       /* Input data scene name                 */
    int *rm_ids;
    int *cpx;                        /* nunber of clear pixels X ?            */
    float **cpy;                     /* nunber of clear pixels Y ?            */
    int i_rec;                       /* start of model before noise removal   */

    int end;
    float **fit_cft;                 /* Fitted coefficients 2-D array.        */
    float **rec_v_dif;
    float **rec_v_dif_copy;
    float v_dif_norm = 0.0;
    float mini_rmse;                 /* Mimimum RMSE                          */

    float v_dif_mean;
    float *rmse;                     /* Root Mean Squared Error array.        */
    int i_count;                     /* Count difference of i each iteration  */
    int i_break;                     /* for recording break points, i is index*/
    float **v_diff;
    float *vec_mag; /* what is the differece */ /* they are used in 2 different branches */
    float *vec_magg;/* these two?            */ /* this one is never freed */
    float vec_magg_min;
    float **v_dif_mag;               /* vector for magnitude of differences.  */
    float **temp_v_dif;              /* for the thermal band.......           */
    int rec_fc;                      /* Record num. of functional curves      */
    int i_conse;
    int i_b;
    int ini_conse;
    int i_ini;                       /* for recording begin of time, i is index*/
    int i_dense;
    float ts_pred_temp;
    int i_span, rm_ids_len;
    float time_span;                 /* Span of time in no. of years.         */
    int *bl_ids;
    int *ids;
    int bl_train;                    /* Flag for which way to train the model.*/
    int i_start;                     /* The first observation for TSFit       */
    int num_c = 8;                   /* Max number of coefficients for model  */
    int update_num_c = 8;            /* Number of coefficients to update      */
    int n_rmse;                      /* number of RMSE values                 */
    float *d_yr;
    int d_rt;
    float* v_start;  /* Vector for start of observation(s)    */
    float* v_end;    /* Vector for end of observastion(s)     */
    float* v_slope;  /* Vector for anormalized slope values   */
    float* v_dif;    /* Vector for difference values          */
    float* tmpcg_rmse; /* to temporarily change RMSE          */

    int ids_old_len;
    int *ids_old;
    int id_last;                     /* The last stable id.                   */

    int ids_len;                     /* number of ids, incremented continuously*/
    float break_mag;
    int adj_conse;
    float adj_TCG;
    float min_rmse[TOTAL_IMAGE_BANDS]; /* Adjusted RMSE for all bands          */
    float date_vario;           /* I: median date                                          */
    float max_date_difference;   /* I: maximum difference between two neighbor dates        */
    float** v_diff_tmp;

    int n_clr;                  /* I: the number of clear pixels                          */
    int *clrx;                  /* I: clear pixel curve in X direction (date)             */
    float **clry;               /* I: clear pixel curve in Y direction (spectralbands)    */
    float mean_angle;           /* I: mean angle of vec_diff                              */
    int pre_end;
    // short int tmp_max_prob = 0;
    short int tmp_CM = 0;
    int tmp;
    int current_CM_n;
    float prob_angle; // change probability for angle
    int i_span_skip = 0;

    fit_cft = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, LASSO_COEFFS,
                                         sizeof (float));
    if (fit_cft == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft memory", FUNC_NAME, FAILURE);
    }

    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
    {
            min_rmse[k] = 0.0;
    }



    rmse = (float *)calloc(TOTAL_IMAGE_BANDS, sizeof(float));
    if (rmse == NULL)
    {
        RETURN_ERROR ("Allocating rmse memory", FUNC_NAME, FAILURE);
    }

    n_clr = 0;
    clrx = (int* )malloc(valid_num_scenes * sizeof(int));
    clry = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, valid_num_scenes,
                                         sizeof (float));

//        for (i = 0; i < valid_num_scenes; i++)
//        {
//            printf("clrx: %d\n", (int)clrx[i]);
//            printf("id: %d\n", (int)id_range[i]);
//            for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
//            {
//               printf("buf: %d\n", (int)buf[k][i]);
//               printf("clry: %f\n", (double)clry[k][i]);
//            }

//        }

    for (i = 0; i < valid_num_scenes; i++)
    {
        if ((fmask_buf[i] < 2) && (id_range[i] == 1))
        {
            // remain the first element for replicated date
            if((n_clr > 0) && (valid_date_array[i] == clrx[n_clr - 1]))
                continue;
            else
            {
                clrx[n_clr] = valid_date_array[i];
                //printf("%d is %d\n", n_clr + 1, clrx[n_clr]);
                for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                {
                    if (k == 0)
                        clry[k][n_clr] = (float)buf_b[i];
                    else if (k == 1)
                        clry[k][n_clr] = (float)buf_g[i];
                    else if (k == 2)
                        clry[k][n_clr] = (float)buf_r[i];
                    else if (k == 3)
                        clry[k][n_clr] = (float)buf_n[i];
                    else if (k == 4)
                        clry[k][n_clr] = (float)buf_s1[i];
                    else if (k == 5)
                        clry[k][n_clr] = (float)buf_s2[i];
                    else if (k == 6){
                        if ((b_c2 == TRUE) && (buf_t[i] == 0))
                            clry[k][n_clr] = 0;
                        else
                            clry[k][n_clr] = (float)(buf_t[i] * 10 - 27320);
                    }

                    //printf("%3.2f\n", clry[k][n_clr]);
                }
                n_clr++;
            }
        }
    }


//    for (k = 0; k < n_clr; k++)
//    {
//        printf("clrx %d: %d\n", k+1, (int)clrx[k]);
//    }


    temp_v_dif = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, n_clr,
                                     sizeof (float));
    if (temp_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating temp_v_dif memory",FUNC_NAME, FAILURE);
    }

    rec_v_dif = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, n_clr,
                                     sizeof (float));
    if (rec_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating rec_v_dif memory",FUNC_NAME, FAILURE);
    }

    rec_v_dif_copy = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, n_clr,
                                     sizeof (float));
    if (rec_v_dif_copy == NULL)
    {
        RETURN_ERROR ("Allocating rec_v_dif_copy memory",FUNC_NAME, FAILURE);
    }


    ids = (int *)calloc(n_clr, sizeof(int));
    if (ids == NULL)
    {
        RETURN_ERROR("ERROR allocating ids memory", FUNC_NAME, FAILURE);
    }
    ids_old = (int *)calloc(n_clr, sizeof(int));
    if (ids_old == NULL)
    {
        RETURN_ERROR("ERROR allocating ids_old memory", FUNC_NAME, FAILURE);
    }

    bl_ids = (int *)calloc(n_clr, sizeof(int));
    if (bl_ids == NULL)
    {
        RETURN_ERROR("ERROR allocating bl_ids memory", FUNC_NAME, FAILURE);
    }

    rm_ids = (int *)calloc(n_clr, sizeof(int));
    if (rm_ids == NULL)
    {
        RETURN_ERROR("ERROR allocating rm_ids memory", FUNC_NAME, FAILURE);
    }

    v_start = (float *)calloc(NUM_LASSO_BANDS, sizeof(float));
    if (v_start == NULL)
    {
        RETURN_ERROR("ERROR allocating v_start memory", FUNC_NAME, FAILURE);
    }

    v_end = (float *)calloc(NUM_LASSO_BANDS, sizeof(float));
    if (v_end == NULL)
    {
        RETURN_ERROR("ERROR allocating v_end memory", FUNC_NAME, FAILURE);
    }

    v_slope = (float *)calloc(NUM_LASSO_BANDS, sizeof(float));
    if (v_slope == NULL)
    {
        RETURN_ERROR("ERROR allocating v_slope memory", FUNC_NAME, FAILURE);
    }

    v_dif = (float *)calloc(NUM_LASSO_BANDS, sizeof(float));
    if (v_dif == NULL)
    {
        RETURN_ERROR("ERROR allocating v_dif memory", FUNC_NAME, FAILURE);
    }

    tmpcg_rmse = (float *)calloc(NUM_LASSO_BANDS, sizeof(float));
    if (tmpcg_rmse == NULL)
    {
        RETURN_ERROR("ERROR allocating tmpcg_rmse memory", FUNC_NAME, FAILURE);
    }



    end = n_clr;

//    /**************************************************************/
//    /*                                                            */
//    /* Remove repeated ids.                                       */
//    /*                                                            */
//    /**************************************************************/

//    matlab_unique(clrx, clry, n_clr, &end);
//    n_clr = end;
    /**************************************************************/
    /*                                                            */
    /* calculate variogram for each band and dates.                           */
    /*                                                            */
    /**************************************************************/
    status = adjust_median_variogram(clrx, clry, TOTAL_IMAGE_BANDS, 0, end-1, &date_vario,
                                     &max_date_difference, min_rmse, 1);
    if (status != SUCCESS)
    {
            RETURN_ERROR("ERROR calling median_variogram routine", FUNC_NAME,
                         FAILURE);
    }

    /* adjust T_cg based delta days*/
//    adj_conse = round (conse * 16 / (double)date_vario);
//    if (adj_conse < conse)
//        adj_conse = conse;

//     /* adjust conse based delta days*/
//    if(adj_conse > conse)
//    {
//        // adj_TCG = chi2inv(1 - pow(1 - PROB_T_CG, (double)conse / (double)adj_conse), NUM_LASSO_BANDS);
//        adj_TCG =X2(NUM_LASSO_BANDS, 1 - pow(1 - probability_threshold, (double)conse / (double)adj_conse));
//    }
//    else
//    {
//        adj_TCG = tcg;
//    }
    adj_conse = conse;
    adj_TCG = tcg;

    v_dif_mag = (float **) allocate_2d_array(TOTAL_IMAGE_BANDS, adj_conse,
                sizeof (float));
    if (v_dif_mag == NULL)
    {
        RETURN_ERROR ("Allocating v_dif_mag memory",
                                 FUNC_NAME, FAILURE);
    }

    vec_mag = (float *)calloc(adj_conse, sizeof(float));
    if (vec_mag == NULL)
    {
        RETURN_ERROR ("Allocating vec_mag memory", FUNC_NAME, FAILURE);
    }

    v_diff = (float **) allocate_2d_array(NUM_LASSO_BANDS,
                      adj_conse, sizeof (float));
    if (v_diff == NULL)
    {
        RETURN_ERROR ("Allocating v_diff memory",
                      FUNC_NAME, FAILURE);
    }

    /**************************************************************/
    /*                                                            */
    /* Start with mininum requirement of clear obs.               */
    /*                                                            */
    /**************************************************************/

    i = N_TIMES * MIN_NUM_C;

    /**************************************************************/
    /*                                                            */
    /* The first observation for TSFit.                           */
    /*                                                            */
    /**************************************************************/

    i_start = 1;

    i_dense = 1;

    /**************************************************************/
    /*                                                            */
    /* Record the start of the model initialization               */
    /*     (0=>initial;1=>done)                                   */
    /*                                                            */
    /**************************************************************/

    bl_train = 0;

//    /**************************************************************/
//    /*                                                            */
//    /* initialize number of the functional curves                 */
//    /*                                                            */
//    /**************************************************************/


//    *num_fc = *num_fc +1;

    /**************************************************************/
    /*                                                            */
    /* Record the *num_fc at the beginning of each pixel.          */
    /*                                                            */
    /**************************************************************/
    rec_fc = *num_fc;


    /**************************************************************/
    /*                                                            */
    /* While loop - process til the last clear observation - adj_conse*/
    /*                                                            */
    /**************************************************************/
    // printf("dd_step1\n");
    while (i <= end - adj_conse)
    {
        /**********************************************************/
        /*                                                        */
        /* span of "i"                                            */
        /*                                                        */
        /**********************************************************/

        i_span = i - i_start + 1;

        /**********************************************************/
        /*                                                        */
        /* span of time (num of years)                            */
        /*                                                        */
        /**********************************************************/
        time_span = (float)(clrx[i-1] - clrx[i_start-1]) / NUM_YEARS;

        /**********************************************************/
        /*                                                        */
        /* basic requrirements: 1) enough observations;           */
        /*                      2) enough time                    */
        /*                                                        */
        /**********************************************************/


//        tmp_direction = 255;
//        posi_count = 0;
//        nega_count = 0;
        if ((i_span >= N_TIMES * MIN_NUM_C) && (time_span >= (float)MIN_YEARS))
        {
            /******************************************************/
            /*                                                    */
            /* Initializing model.                                */
            /*                                                    */
            /******************************************************/

            if (bl_train == 0)
            {
                /*************************************************  */
                /*                                                  */
                /* Step 1: noise removal.                           */
                /*                                                  */
                /****************************************************/

                /*************************************************  */
                /*                                                  */
                /*  check maximum time gap as first                 */
                /*                                                  */
                /****************************************************/

                int max_date_diff = 0;
                for (k = i_start - 1; k < i - 1; k++)
                {
                    if (clrx[k + 1] - clrx[k] > max_date_diff)
                    {
                        max_date_diff = clrx[k + 1] - clrx[k];
                    }
                    //printf("%d \n", clrx[k]);
                }

                if (max_date_diff > gap_days)      //SY 09192018
                {
                    i++;
                    i_start++;
                    i_dense = i_start;
                    continue; //SY 02122019
                }

                status = auto_mask(clrx, clry, i_start-1, i+adj_conse-1,
                               (float)(clrx[i+adj_conse-1]-clrx[i_start-1]) / NUM_YEARS,
                               min_rmse[1], min_rmse[4], (float)T_CONST, bl_ids);
                // printf("ddstep2 auto_mask finished \n");
                if (status != SUCCESS)
                {
                    RETURN_ERROR("ERROR calling auto_mask during model initilization",
                                  FUNC_NAME, FAILURE);
                }

                /**************************************************/
                /*                                                */
                /* Clear the IDs buffers.                         */
                /*                                                */
                /**************************************************/

                for (k = 0; k < n_clr; k++)
                    ids[k] = 0;

                /**************************************************/
                /*                                                */
                /* IDs to be removed.                             */
                /*                                                */
                /**************************************************/

                for (k = i_start-1; k < i; k++)
                {
                    ids[k-i_start+1] = k;
                }
                m= 0;
                i_span = 0;
                for (k = 0; k < i-i_start + 1; k++) /** 02282019 SY **/
                {
                    if (bl_ids[k] == 1)
                    {
                        rm_ids[m] = ids[k];
                        //printf("%d \n", ids[k]);
                        m++;
                    }
                    else
                        i_span++;  /* update i_span after noise removal */
                }


                rm_ids_len = m;

                /**************************************************/
                /*                                                */
                /* Check if there are enough observation.         */
                /*                                                */
                /**************************************************/

                if (i_span < (N_TIMES * MIN_NUM_C))
                {
                    /**********************************************/
                    /*                                            */
                    /* Move forward to the i+1th clear observation*/
                    /*                                            */
                    /**********************************************/

                    i++;

                    /**********************************************/
                    /*                                            */
                    /* Not enough clear observations.             */
                    /*                                            */
                    /**********************************************/

                    continue;
                }

                if (end == 0)
                    RETURN_ERROR("No available data point", FUNC_NAME, FAILURE);

                /**************************************************/
                /*                                                */
                /* Allocate memory for cpx, cpy.                  */
                /*                                                */
                /**************************************************/

                cpx = (int* )malloc(end * sizeof(int));
                if (cpx == NULL)
                    RETURN_ERROR("ERROR allocating cpx memory", FUNC_NAME, FAILURE);

                cpy = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, end,
                                 sizeof (float));
                if (cpy == NULL)
                {
                    RETURN_ERROR ("Allocating cpy memory", FUNC_NAME, FAILURE);
                }

                /**************************************************/
                /*                                                */
                /* Remove noise pixels between i_start & i.       */
                /*                                                */
                /**************************************************/

                m = 0;
                for (k = 0, k_new=0; k < end; k++)
                {
                    if (m < rm_ids_len && k == rm_ids[m])
                    {
                        m++;
                        continue;
                    }
                    cpx[k_new] = clrx[k];
                    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                    {
                        cpy[b][k_new] = clry[b][k];
                    }
                    k_new++;
                }


                /**************************************************/
                /*                                                */
                /* Record i before noise removal.                 */
                /* This is very important, ie model is not yet    */
                /* initialized.   The multitemporal masking shall */
                /* be done again instead of removing outliers  In */
                /* every masking.                                 */
                /*                                                */
                /**************************************************/

                i_rec = i;

                /**************************************************/
                /*                                                */
                /* Update i afer noise removal.                   */
                /*     (i_start stays the same).                  */
                /*                                                */
                /**************************************************/

                i = i_start + i_span - 1;

                /**************************************************/
                /*                                                */
                /* Update span of time (num of years).            */
                /*                                                */
                /**************************************************/

                time_span=(cpx[i-1] - cpx[i_start-1]) / NUM_YEARS;

                /**************************************************/
                /*                                                */
                /* Check if there is enough time.                 */
                /*                                                */
                /**************************************************/

                if (time_span < MIN_YEARS)
                {
                    i = i_rec;   /* keep the original i */

                    /**********************************************/
                    /*                                            */
                    /* Move forward to the i+1th clear observation*/
                    /*                                            */
                    /**********************************************/

                    i++;
                    free(cpx);
                    status = free_2d_array ((void **) cpy);
                    if (status != SUCCESS)
                    {
                          RETURN_ERROR ("Freeing memory: cpy\n",
                                FUNC_NAME, FAILURE);
                    }
                    continue;    /* not enough time span */
                }

                //SY 09272018
                /**************************************************/
                /*                                                */
                /* updated end after checking if enought time_span*/
                /*                                                */
                /**************************************************/
                end = k_new;

                /**************************************************/
                /*                                                */
                /* Remove noise in original arrays.               */
                /*                                                */
                /**************************************************/

                for (k = 0; k < end; k++)
                {
                    clrx[k] = cpx[k];
                    for (m = 0; m < TOTAL_IMAGE_BANDS; m++)
                    {
                        clry[m][k] = cpy[m][k];
                    }
                }

//                    for (k = 0; k < n_clr; k++)
//                    {
//                        printf("clrx %d: %d\n", k+1, (int)clrx[k]);
//                    }

                free(cpx);
                status = free_2d_array ((void **) cpy);
                if (status != SUCCESS)
                {
                    RETURN_ERROR ("Freeing memory: cpy\n",
                         FUNC_NAME, FAILURE);
                }


                /**************************************************/
                /*                                                */
                /* Step 2) model fitting: initialize model testing*/
                /*         variables defining computed variables. */
                /*                                                */
                /**************************************************/

                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    /**********************************************/
                    /*                                            */
                    /* Initial model fit.                         */
                    /*                                            */
                    /**********************************************/
                    // printf("dd_step2 auto_ts_fit starts \n");
                    status = auto_ts_fit_float(clrx, clry,  i_b, i_b, i_start-1, i-1,
                             MIN_NUM_C, fit_cft, &rmse[i_b], rec_v_dif);
                    // printf("dd_step2 first using auto_ts_fit successed \n");
//                        for (k = 0; k < MAX_NUM_C; k++)
//                        {

//                               printf("%f\n", (float)fit_cft[i_b][k]);

//                        }
//                    printf("rmse for band %d is %f\n", i_b, rmse[i_b]);
                    // printf("auto_ts_fit finished \n");
                    if (status != SUCCESS)
                    {
                        RETURN_ERROR ("Calling auto_ts_fit_float during model initilization\n",
                             FUNC_NAME, FAILURE);
                    }



                }

                v_dif_norm = 0.0;

                for (i_b = 0; i_b < NUM_LASSO_BANDS; i_b++)
                {

                    /**********************************************/
                    /*                                            */
                    /* Calculate min. rmse.                       */
                    /*                                            */
                    /**********************************************/

                    mini_rmse = max((float)min_rmse[lasso_blist[i_b]], rmse[lasso_blist[i_b]]);

                    /**********************************************/
                    /*                                            */
                    /* Compare the first observation.             */
                    /*                                            */
                    /**********************************************/

                    v_start[i_b] = rec_v_dif[lasso_blist[i_b]][0]
                            / mini_rmse;

                    /**********************************************/
                    /*                                            */
                    /* Compare the last clear observation.        */
                    /*                                            */
                    /**********************************************/

                    v_end[i_b] = rec_v_dif[lasso_blist[i_b]][i-i_start]
                                            / mini_rmse;

                    /**********************************************/
                    /*                                            */
                    /* Anormalized slope values.                  */
                    /*                                            */
                    /**********************************************/
                    v_slope[i_b] = fit_cft[lasso_blist[i_b]][1] *
                                    (clrx[i-1]-clrx[i_start-1])/mini_rmse/SLOPE_SCALE;

                    /**********************************************/
                    /*                                            */
                    /* Difference in model intialization.         */
                    /*                                            */
                    /**********************************************/

                    v_dif[i_b] = fabs(v_slope[i_b]) + max(fabs(v_start[i_b]), fabs(v_end[i_b]));
                    v_dif_norm += v_dif[i_b] * v_dif[i_b];
                    //printf("%f \n", v_dif[i_b]);

                }

//                for(b = 0; b < TOTAL_IMAGE_BANDS; b++)
//                    printf("%.10f, %.10f, %.10f, %.10f,%.10f, %.10f, %.10f, %.10f\n", fit_cft[b][0],
//                            fit_cft[b][1], fit_cft[b][2], fit_cft[b][3], fit_cft[b][4], fit_cft[b][5], fit_cft[b][6], fit_cft[b][7]);


                /**************************************************/
                /*                                                */
                /* Find stable start for each curve.              */
                /*                                                */
                /**************************************************/
                if (v_dif_norm > adj_TCG)
                {
                    /**********************************************/
                    /*                                            */
                    /* Start from next clear observation.         */
                    /*                                            */
                    /**********************************************/

                    i_start++;

                    /**********************************************/
                    /*                                            */
                    /* Move forward to the i+1th clear observation*/
                    /*                                            */
                    /**********************************************/

                    i++;

                    /**********************************************/
                    /*                                            */
                    /* Keep all data and move to the next obs.    */
                    /*                                            */
                    /**********************************************/

                    continue;
                }

                /**************************************************/
                /*                                                */
                /* Model is ready.                                */
                /*                                                */
                /**************************************************/

                bl_train = 1;

                /**************************************************/
                /*                                                */
                /* Count difference of i for each iteration.      */
                /*                                                */
                /**************************************************/

                i_count = 0;

                /**************************************************/
                /*                                                */
                /* Find the previous break point.                 */
                /*                                                */
                /**************************************************/

                if (*num_fc == rec_fc)
                {
                    i_break = i_dense; /* first curve */
                }
                else
                {
                    /**********************************************/
                    /*                                            */
                    /* After the first curve, compare rmse to     */
                    /* determine which curve to determine t_break.*/
                    /*                                            */
                    /**********************************************/

                    for (k = 0; k < end; k++)
                    {
                        if (clrx[k] >= rec_cg[*num_fc-1].t_break)
                        {
                            i_break = k + 1;
                            break;
                        }
                    }
                }

                if (i_start > i_break)
                {
                    /**********************************************/
                    /*                                            */
                    /* Model fit at the beginning of the time     */
                    /* series.                                    */
                    /*                                            */
                    /**********************************************/

                    for(i_ini = i_start - 2; i_ini >= i_break - 1; i_ini--) // SY 09192018
                    {
                        if ((i_ini - (i_break - 1) + 1) < adj_conse)
                        {
                            ini_conse = i_ini - (i_break - 1) + 1;
                        }
                        else
                        {
                            ini_conse = adj_conse;
                        }

                        if (ini_conse == 0)
                        {
                            RETURN_ERROR ("No data point for model fit at "
                                  "the begining", FUNC_NAME, FAILURE);
                        }

                        /******************************************/
                        /*                                        */
                        /* Allocate memory for model_v_dif,       */
                        /* v_diff, vec_magg for the non-stdin     */
                        /* branch here.                           */
                        /*                                        */
                        /******************************************/

                        v_diff_tmp = (float **) allocate_2d_array(NUM_LASSO_BANDS,
                                    ini_conse, sizeof (float));
                        if (v_diff_tmp == NULL)
                        {
                            RETURN_ERROR ("Allocating v_diff_tmp memory",
                                           FUNC_NAME, FAILURE);
                        }

                        vec_magg = (float *) malloc(ini_conse * sizeof (float));
                        if (vec_magg == NULL)
                        {
                            RETURN_ERROR ("Allocating vec_magg memory",
                                          FUNC_NAME, FAILURE);
                        }

                        /******************************************/
                        /*                                        */
                        /* Detect change.                         */
                        /* value of difference for adj_conse      */
                        /* observations                          */
                        /* Record the magnitude of change.        */
                        /*                                        */
                        /******************************************/

                        vec_magg_min = 9999.0;
                        for (i_conse = 1; i_conse < ini_conse + 1; i_conse++) // SY 09192018
                        {
                            v_dif_norm = 0.0;
                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {

                                /**********************************/
                                /*                                */
                                /* Absolute differences.          */
                                /*                                */
                                /**********************************/

                                // SY 09192018 moving fitting into (i_b == lasso_blist[b])to save time //
                                // SY 02/13/2019 delete these speed-up modification as non-lasso bands
                                // are important for change agent classification
                                auto_ts_predict_float(clrx, fit_cft, MIN_NUM_C, i_b, i_ini-i_conse+1,
                                                i_ini-i_conse+1, &ts_pred_temp);
                                v_dif_mag[i_b][i_conse-1] = (float)clry[i_b][i_ini-i_conse+1] -
                                                   ts_pred_temp;// SY 09192018
                                // printf("auto_ts_predict finished \n");
                                /**********************************/
                                /*                                */
                                /* Normalize to z-score.          */
                                /*                                */
                                /**********************************/

                                for (b = 0; b < NUM_LASSO_BANDS; b++)
                                {
                                   if (i_b == lasso_blist[b])
                                   {
                                        /**************************/
                                        /*                        */
                                        /* Minimum rmse.          */
                                        /*                        */
                                        /**************************/

                                        mini_rmse = max((float)min_rmse[i_b], rmse[i_b]);

                                        /**************************/
                                        /*                        */
                                        /* z-scores.              */
                                        /*                        */
                                        /**************************/

                                        v_diff_tmp[b][i_conse-1] = v_dif_mag[i_b][i_conse-1] // SY 09192018
                                                                      / mini_rmse;
                                        v_dif_norm += v_diff_tmp[b][i_conse-1] * v_diff_tmp[b][i_conse-1]; // SY 09192018
            //                                        printf("i_b: %d\n",i_b);
//                                        printf("clry: %f\n", clry[i_b][i_ini-i_conse]);

//                                        printf("ts_pred_temp: %f\n",ts_pred_temp);
//                                        printf("v_dif_norm: %f\n",v_dif_norm);
//                                        printf("mini_rmse: %f\n",mini_rmse);
//                                        printf("v_dif_mag: %f\n",v_dif_mag[i_b][i_conse]);
                                    }
                                  }
                                }
                                vec_magg[i_conse-1] = v_dif_norm; // SY 09192018

                                if (vec_magg_min > vec_magg[i_conse-1])
                                {
                                    vec_magg_min =  vec_magg[i_conse-1]; // SY 09192018
                                }
                        }

                        /******************************************/
                        /*                                        */
                        /* Change angle.                      */
                        /*                                        */
                        /******************************************/

                        mean_angle = MeanAngl_float(v_diff_tmp, NUM_LASSO_BANDS, ini_conse);
                        /******************************************/
                        /*                                        */
                        /* Change detection.                      */
                        /*                                        */
                        /******************************************/
//                        for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
//                        {
//                            printf("%f, %f, %f, %f\n", fit_cft[b][0], fit_cft[b][1], fit_cft[b][2], fit_cft[b][3]);

//                        }
                        if (b_outputCM == TRUE){
                            if(ini_conse >= adj_conse){  // only save CM based on >= conse obs
                            // if ((ini_conse >= adj_conse) && (*num_fc == rec_fc)){ // only save out the first curve
                                prob_angle = angle_decaying(mean_angle, (double)NSIGN, 90.0);
                                current_CM_n = (clrx[i_ini+1] - starting_date) / CM_OUTPUT_INTERVAL;  // looking back, t_break = current_i + 1
                                tmp = round(prob_angle * vec_magg_min * 100);
                                if (tmp > MAX_SHORT) // MAX_SHORT is upper limit of short 16
                                    tmp = MAX_SHORT;
                                tmp_CM = (short int) (tmp);
                                if(tmp_CM > CM_outputs[current_CM_n])
                                {
                                    /*********************************************/
                                    /*      change direction by majority vote    */
                                    /*********************************************/
//                                    tmp_direction = 0;
//                                    for (b = 0; b < NUM_LASSO_BANDS; b++)
//                                    {
//                                        posi_count = 0;
//                                        nega_count = 0;
//                                        for(j = 0; j < adj_conse; j++){
//                                            if (v_diff_tmp[b][j] > 0){
//                                                posi_count++;
//                                            }else{
//                                                nega_count++;
//                                            }
//                                        }

//                                        if (posi_count > nega_count){
//                                            tmp_direction = tmp_direction + pow(2, b);
//                                        }
//                                    }
                                    CM_outputs[current_CM_n] = tmp_CM;
//                                    CMdirection_outputs[current_CM_n] = tmp_direction;
                                    CM_outputs_date[current_CM_n] = (short int)(clrx[i_ini+1] - ORDINAL_LANDSAT4_LAUNCH);
                                    // printf("date = %d\n", clrx[i_ini+1]);
                                }
                            }
                        }

                        if ((vec_magg_min > adj_TCG) && (mean_angle < NSIGN)) /* change detected */
                        {
                            free(vec_magg);
                            status = free_2d_array ((void **) v_diff_tmp);
                            if (status != SUCCESS)
                            {
                                RETURN_ERROR ("Freeing memory: v_diff_tmp\n",
                                              FUNC_NAME, FAILURE);
                            }
                            break;
                        }
                        else if (vec_magg[0] > T_MAX_CG) /* false change */
                        {
                            for (k = i_ini; k < end - 1; k++)
                            {
                                clrx[k] = clrx[k+1];
                                for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                                {
                                    clry[b][k] = clry[b][k+1];
                                }
                            }
                            i--;
                            end--;

                        }

                        /******************************************/
                        /*                                        */
                        /* Free the temporary memory.             */
                        /*                                        */
                        /******************************************/

                        free(vec_magg);
                        status = free_2d_array ((void **) v_diff_tmp);
                        if (status != SUCCESS)
                        {
                            RETURN_ERROR ("Freeing memory: v_diff_tmp\n",
                                          FUNC_NAME, FAILURE);
                        }

                        /**************************************/
                        /*                                    */
                        /* Update i_start if i_ini is not a   */
                        /* confirmed break.                   */
                        /*                                    */
                        /**************************************/

                        i_start = i_ini + 1;

                    } // end for (i_ini = i_start-1; i_ini >= i_break; i_ini--)

                }    // end for if (i_start > i_break)

                /**************************************************/
                /*                                                */
                /* Enough to fit simple model and confirm a break.*/
                /*                                                */
                /**************************************************/

                if ((*num_fc == rec_fc) && ((i_start - i_dense) >= LASSO_MIN))
                {
                    /**********************************************/
                    /*                                            */
                    /* Defining computed variables.               */
                    /*                                            */
                    /**********************************************/

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {

                        status = auto_ts_fit_float(clrx, clry, i_b, i_b,  i_dense-1, i_start-2,
                                 MIN_NUM_C, fit_cft, &rmse[i_b], temp_v_dif);//SY 02132019
                        if (status != SUCCESS)
                        {
                              RETURN_ERROR ("Calling auto_ts_fit_float with enough observations\n",
                                         FUNC_NAME, FAILURE);
                        }

                    }

                    /**********************************************/
                    /*                                            */
                    /* Record time of curve end,                  */
                    /* postion of the pixels.                     */
                    /*                                            */
                    /**********************************************/


                    rec_cg[*num_fc].t_end = clrx[i_start-2];
                    //rec_cg[*num_fc].pos.row = row;
                    //rec_cg[*num_fc].pos.col = col;

                    /**********************************************/
                    /*                                            */
                    /* Record break time, fit category, change    */
                    /* probability, time of curve start, number   */
                    /* of observations, change magnitude.         */
                    /*                                            */
                    /**********************************************/

                    rec_cg[*num_fc].t_break = clrx[i_start -1];
                    rec_cg[*num_fc].category = 10 + MIN_NUM_C;
                    rec_cg[*num_fc].change_prob = 100;
                    rec_cg[*num_fc].t_start = clrx[i_dense-1];
                    rec_cg[*num_fc].num_obs = i_start - i_dense + 1;  //SY 09182018

//                    if ((i_start - 1 + adj_conse) < end)
//                        rec_cg[*num_fc].t_confirmed = clrx[i_start + adj_conse - 1];
//                    else
//                        rec_cg[*num_fc].t_confirmed = clrx[end - 1];

                    for (i_b = 0; i_b <  TOTAL_IMAGE_BANDS; i_b++)
                    {
                        quick_sort_float(v_dif_mag[i_b], 0, ini_conse-1);
                        matlab_2d_float_median(v_dif_mag, i_b, ini_conse,
                                              &v_dif_mean);
                        rec_cg[*num_fc].magnitude[i_b] = -v_dif_mean;
                    }

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        for (k = 0; k < LASSO_COEFFS; k++)
                        {
                            /**************************************/
                            /*                                    */
                            /* Record fitted coefficients.        */
                            /*                                    */
                            /**************************************/
                            if(k < MIN_NUM_C)
                                rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                            else
                                rec_cg[*num_fc].coefs[i_b][k] = 0;

                        }

                        /******************************************/
                        /*                                        */
                        /* Record rmse of the pixel.              */
                        /*                                        */
                        /******************************************/

                        rec_cg[*num_fc].rmse[i_b] = rmse[i_b];


                     }

                    /**********************************************/
                    /*                                            */
                    /* Identified and move on for the next        */
                    /* functional curve.                          */
                    /*                                            */
                    /**********************************************/

                    *num_fc = *num_fc + 1;
                    if (*num_fc >= NUM_FC)
                    {
                        /******************************************/
                        /*                                        */
                        /* Reallocate memory for rec_cg.          */
                        /*                                        */
                        /******************************************/

                        rec_cg = realloc(rec_cg, (*num_fc + 1) * sizeof(Output_t));
                        if (rec_cg == NULL)
                        {
                            RETURN_ERROR("ERROR allocating rec_cg memory",
                                         FUNC_NAME, FAILURE);
                        }
                    }
                }
            } /* end of initializing model */

            /******************************************************/
            /*                                                    */
            /* Allocate memory for v_diff for the non-stdin branch*/
            /*                                                    */
            /******************************************************/

//            for (k = 0; k < n_clr; k++)
//            {
//                printf("clrx %d: %d\n", k+1, (int)clrx[k]);
//            }


            /******************************************************/
            /*                                                    */
            /* Continuous monitoring started!!!                   */
            /*                                                    */
            /******************************************************/
            // printf("dd_step3 initialization finished\n");
            if (bl_train == 1)
            {
                // printf("processing %d obs finished \n", i);
                /**************************************************/
                /*                                                */
                /* Clears the IDs buffers.                        */
                /*                                                */
                /**************************************************/

                for (k = 0; k < n_clr; k++)
                {
                    ids[k] = 0;
                }

                /**************************************************/
                /*                                                */
                /* All IDs.                                       */
                /*                                                */
                /**************************************************/

                ids_len = 0;
                for (k = i_start-1; k < i; k++)
                {
                    ids[k-i_start+1] = k;
                    ids_len++;
                }
                i_span = i - i_start +1;

                if(i_span_skip == 0)
                    i_span_skip = i_span;
                /**************************************************/
                /*                                                */
                /* Determine the time series model.               */
                /*                                                */
                /**************************************************/

                update_cft(i_span, N_TIMES, MIN_NUM_C, MID_NUM_C, MAX_NUM_C,
                          num_c, &update_num_c);

                /************************************************************/
                /*                                                          */
                /* initial model fit when there are not many observations.  */
                /* if (i_count == 0 || ids_old_len < (N_TIMES * MAX_NUM_C)) */
                /*                                                          */
                /************************************************************/

                if (i_count == 0 || i_span <= (N_TIMES * MAX_NUM_C))
                {
                    /**********************************************/
                    /*                                            */
                    /* update i_count at each iteration.          */
                    /*                                            */
                    /**********************************************/

                    i_count = clrx[i-1] - clrx[i_start-1];

                    pre_end = i;

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {

                        status = auto_ts_fit_float(clrx, clry, i_b, i_b, i_start-1, i-1, update_num_c,
                                             fit_cft, &rmse[i_b], rec_v_dif);
                        if (status != SUCCESS)
                        {
                            RETURN_ERROR ("Calling auto_ts_fit_float during continuous monitoring\n",
                                          FUNC_NAME, FAILURE);
                        }

                    }
                    // printf("auto_ts_fit_float finished \n", i);

                    /**********************************************/
                    /*                                            */
                    /* Updating information for the first         */
                    /* iteration.  Record time of curve start and */
                    /* time of curve end.                         */
                    /*                                            */
                    /**********************************************/

                    rec_cg[*num_fc].t_start = clrx[i_start-1];
                    rec_cg[*num_fc].t_end = clrx[i-1];

                    /**********************************************/
                    /*                                            */
                    /* No break at the moment.                    */
                    /*                                            */
                    /**********************************************/

                    rec_cg[*num_fc].t_break = 0;


                    /**********************************************/
                    /*                                            */
                    /* Record change probability, number of       */
                    /* observations, fit category.                */
                    /*                                            */
                    /**********************************************/

                    rec_cg[*num_fc].change_prob = 0;
                    rec_cg[*num_fc].num_obs = i-i_start+1;
                    rec_cg[*num_fc].category = 0 + update_num_c;

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {

                        /******************************************/
                        /*                                        */
                        /* Record rmse of the pixel.              */
                        /*                                        */
                        /******************************************/

                        rec_cg[*num_fc].rmse[i_b] = rmse[i_b];

                        /******************************************/
                        /*                                        */
                        /* Record change magnitude.               */
                        /*                                        */
                        /******************************************/
                        rec_cg[*num_fc].magnitude[i_b] = 0.0;

                        for (k = 0; k < LASSO_COEFFS; k++)
                        {
                            /**************************************/
                            /*                                    */
                            /* Record fitted coefficients.        */
                            /*                                    */
                            /**************************************/

                            rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                         }


                     }
                    // printf("auto_ts_fit2 finished \n", i);


                    /**********************************************/
                    /*                                            */
                    /* Detect change, value of difference for     */
                    /* adj_conse observations.                        */
                    /*                                            */
                    /**********************************************/
                    for (m = 0; m < adj_conse; m++)
                    {
                        vec_mag[m] = 0;
                        for (b = 0; b < NUM_LASSO_BANDS; b++)
                            v_diff[b][m] = 0;
                        for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                            v_dif_mag[b][m] = 0;
                    }


                    for (i_conse = 1; i_conse < adj_conse+1; i_conse++)//SY 09192018
                    {
                        v_dif_norm = 0.0;
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        {
                            /**************************************/
                            /*                                    */
                            /* Absolute differences.              */
                            /*                                    */
                            /**************************************/
                            // printf("auto_ts_predict started finished \n", i);
                            auto_ts_predict_float(clrx, fit_cft, update_num_c, i_b, i+i_conse-1, i+i_conse-1,
                                            &ts_pred_temp); //SY 09192018
                            // printf("auto_ts_predict finished \n", i);
                            v_dif_mag[i_b][i_conse-1] = (float)clry[i_b][i+i_conse-1] - ts_pred_temp;//SY 09192018

                            /**************************************/
                            /*                                    */
                            /* Normalize to z-score.              */
                            /*                                    */
                            /**************************************/

                            for (b = 0; b < NUM_LASSO_BANDS; b++)
                            {
                                if (i_b == lasso_blist[b])
                                {
                                    /******************************/
                                    /*                            */
                                    /* Minimum rmse,              */
                                    /* z-scores.                  */
                                    /*                            */
                                    /******************************/

                                    mini_rmse = max(min_rmse[i_b], rmse[i_b]);
                                    v_diff[b][i_conse-1] = v_dif_mag[i_b][i_conse-1] / mini_rmse;
                                    v_dif_norm += v_diff[b][i_conse-1] * v_diff[b][i_conse-1];
                                }
                            }
                        }
                        vec_mag[i_conse-1] = v_dif_norm;//SY 09192018
                    }

                    /**********************************************/
                    /*                                            */
                    /* Clears the IDs_old buffers.                */
                    /*                                            */
                    /**********************************************/

                    for (k = 0; k < ids_len; k++)
                    {
                        ids_old[k] = 0;
                    }

                    /**********************************************/
                    /*                                            */
                    /* IDs that have not been updated.            */
                    /*                                            */
                    /**********************************************/

                    for (k = 0; k < ids_len; k++)
                    {
                        ids_old[k] = ids[k];
                    }
                    ids_old_len = ids_len;

                } // end for if (i_count == 0 || i_span <= (N_TIMES * MAX_NUM_C))
                else // update frequency
                {
                    if((i - pre_end >= UPDATE_FREQ) && (i - pre_end) > (int)(i_span_skip * SKIP_PERCENTAGE))
                    // if ((float)(clrx[i-1] - clrx[i_start-1]) >= (1.33*(float)i_count))
                    {
                        /******************************************/
                        /*                                        */
                        /* Update coefficent at each iteration year. */
                        /*                                        */
                        /******************************************/

                        i_count = clrx[i-1] - clrx[i_start-1];
                        pre_end = i;
                        // update i_span_skip
                        i_span_skip = i_span;
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        {
                            status = auto_ts_fit_float(clrx, clry, i_b, i_b, i_start-1, i-1, update_num_c,
                                             fit_cft, &rmse[i_b], rec_v_dif);
                            // printf("auto_ts_fit2 finished \n", i);
                            if (status != SUCCESS)
                            {
                                RETURN_ERROR ("Calling auto_ts_fit_float for change detection with "
                                     "enough observations\n", FUNC_NAME, FAILURE);
                            }


                        }

                        /******************************************/
                        /*                                        */
                        /* Record fitted coefficients.            */
                        /*                                        */
                        /******************************************/
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        {

                            for (k = 0; k < LASSO_COEFFS; k++)
                            {
                                /**********************************/
                                /*                                */
                                /* Record fitted coefficients.    */
                                /*                                */
                                /**********************************/

                                rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                            }
                            /**************************************/
                            /*                                    */
                            /* Record rmse of the pixel.          */
                            /*                                    */
                            /**************************************/

                            rec_cg[*num_fc].rmse[i_b] = rmse[i_b];

                            /******************************************/
                            /*                                        */
                            /* Record number of observations, fit     */
                            /* category.                              */
                            /*                                        */
                            /******************************************/

                            rec_cg[*num_fc].num_obs = i - i_start + 1;
                            rec_cg[*num_fc].category = 0 + update_num_c;
                        }

                        /******************************************/
                        /*                                        */
                        /* Clears the IDs_Old buffers.            */
                        /*                                        */
                        /******************************************/

                        for (k = 0; k < ids_len; k++)
                        {
                            ids_old[k] = 0;
                        }

                        /******************************************/
                        /*                                        */
                        /* IDs that have not been updated.        */
                        /*                                        */
                        /******************************************/

                        for (k = 0; k < ids_len; k++)
                        {
                            ids_old[k] = ids[k];
                        }
                        ids_old_len = ids_len;

                    } //  if(i_update == UPDATE_FREQ || i_update > UPDATE_FREQ)


                    /**********************************************/
                    /*                                            */
                    /* Record time of curve end.                  */
                    /*                                            */
                    /**********************************************/

                    rec_cg[*num_fc].t_end = clrx[i-1];

                    /**********************************************/
                    /*                                            */
                    /* Use fixed number for RMSE computing.       */
                    /*                                            */
                    /**********************************************/

                    n_rmse = N_TIMES * rec_cg[*num_fc].category;

                    /**********************************************/
                    /*                                            */
                    /* Better days counting for RMSE calculating  */
                    /* relative days distance.                    */
                    /*                                            */
                    /**********************************************/

                    if (ids_old_len == 0)
                    {
                        RETURN_ERROR ("No data points for RMSE calculating",
                                     FUNC_NAME, FAILURE);
                    }

                    d_yr = malloc(ids_old_len * sizeof(float));
                    if (d_yr == NULL)
                    {
                        RETURN_ERROR ("Allocating d_yr memory",
                                     FUNC_NAME, FAILURE);
                    }

                    for(m = 0; m < ids_old_len; m++)
                    {
                        d_rt = clrx[ids_old[m]] - clrx[i+adj_conse-1];
                        d_yr[m] = fabs(round((float)d_rt / NUM_YEARS) * NUM_YEARS - (float)d_rt);
                    }

                    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                    {
                        for (m = 0; m < ids_old_len; m++)
                            rec_v_dif_copy[b][m] = rec_v_dif[b][m];
                    }

                    /**********************************************/
                    /*                                            */
                    /* Sort the rec_v_dif based on d_yr.          */
                    /*                                            */
                    /**********************************************/

                    quick_sort_2d_float(d_yr, rec_v_dif_copy, 0, ids_old_len-1, TOTAL_IMAGE_BANDS);
                    for(b = 0; b < NUM_LASSO_BANDS; b++)
                        tmpcg_rmse[b] = 0.0;

                    /**********************************************/
                    /*                                            */
                    /* Temporarily changing RMSE.                 */
                    /*                                            */
                    /**********************************************/

                    for (b = 0; b < NUM_LASSO_BANDS; b++)
                    {
                        matlab_2d_array_norm_float(rec_v_dif_copy, lasso_blist[b], n_rmse,
                                         &tmpcg_rmse[b]);
                        tmpcg_rmse[b] /= sqrtf((float)(n_rmse - rec_cg[*num_fc].category));
                    }

                    /**********************************************/
                    /*                                            */
                    /* Free allocated memories.                   */
                    /*                                            */
                    /**********************************************/
                    free(d_yr);

                    /**********************************************/
                    /*                                            */
                    /* Move the ith col to i-1th col.             */
                    /*                                            */
                    /**********************************************/

                    for (m = 0; m < adj_conse-1; m++)
                    {
                        vec_mag[m] = vec_mag[m+1];
                        for (b = 0; b < NUM_LASSO_BANDS; b++)
                            v_diff[b][m] = v_diff[b][m+1];
                        for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                            v_dif_mag[b][m] = v_dif_mag[b][m+1];
                    }

                    for (b = 0; b < NUM_LASSO_BANDS; b++)
                        v_diff[b][adj_conse-1] = 0.0;
                    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                        v_dif_mag[b][adj_conse-1] = 0.0;
                    vec_mag[adj_conse-1] = 0.0;

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {

                        /******************************************/
                        /*                                        */
                        /* .     */
                        /*                                        */
                        /******************************************/
    //                        if (i == 45){
    //                            i = 45;
    //                       }
                        auto_ts_predict_float(clrx, fit_cft, update_num_c, i_b, i+adj_conse-1,
                             i+adj_conse-1, &ts_pred_temp);
                        v_dif_mag[i_b][adj_conse-1] = (float)clry[i_b][i+adj_conse-1] - ts_pred_temp;

                        /******************************************/
                        /*                                        */
                        /* Normalized to z-scores.                */
                        /*                                        */
                        /******************************************/
                        for (b = 0; b < NUM_LASSO_BANDS; b++)
                        {
                            if (i_b == lasso_blist[b])
                            {
                                /**********************************/
                                /*                                */
                                /* Minimum rmse.                  */
                                /*                                */
                                /**********************************/

                                mini_rmse = max((double)min_rmse[i_b], tmpcg_rmse[b]);

                                /**********************************/
                                /*                                */
                                /* Z-score.                       */
                                /*                                */
                                /**********************************/

                                v_diff[b][adj_conse-1] = v_dif_mag[i_b][adj_conse-1] / mini_rmse;
                                vec_mag[adj_conse-1] += v_diff[b][adj_conse-1] * v_diff[b][adj_conse-1]; // SY 02132014
                            }
                        }
                     }
                } // else update frequency

                mean_angle = MeanAngl_float(v_diff, NUM_LASSO_BANDS, adj_conse);


                break_mag = 9999.0;
                for (m = 0; m < adj_conse; m++)
                {
                    if (break_mag > vec_mag[m])
                    {
                        break_mag = vec_mag[m];
                    }
                }

//                if(i == 532 - 4){
//                  for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
//                    for(j = 0; j < MAX_NUM_C; j++)
//                        printf("%f\n", fit_cft[k][j]);

//                  for(m = 0; m < n_rmse; m++)
//                  {

//                      printf("%f\n", d_yr[m]);
//                      printf("%f\n", rec_v_dif_copy[0][m]);

//                  }
//                }

                if (b_outputCM == TRUE)
                {
                    prob_angle = angle_decaying(mean_angle, (double)NSIGN, 90.0);
                    // prob_MCM = Chi_Square_Distribution(break_mag, NUM_LASSO_BANDS);
                    current_CM_n = (clrx[i] - starting_date) / CM_OUTPUT_INTERVAL;
                    tmp = round(prob_angle * break_mag * 100);
                    if (tmp > MAX_SHORT) // MAX_SHORT is upper limit of short 16
                        tmp = MAX_SHORT;
                    tmp_CM = (short int) (tmp);

                    /*********************************************/
                    /*      change direction by majority vote    */
                    /*********************************************/
                    if(tmp_CM > CM_outputs[current_CM_n])
                    {
//                        tmp_direction = 0;
//                        for (b = 0; b < NUM_LASSO_BANDS; b++)
//                        {
//                            posi_count = 0;
//                            nega_count = 0;
//                            for(j = 0; j < adj_conse; j++){
//                                if (v_diff[b][j] > 0){
//                                    posi_count++;
//                                }else{
//                                    nega_count++;
//                                }
//                            }
//                            if (posi_count > nega_count){
//                                tmp_direction = tmp_direction + pow(2, b);
//                            }
//                        }
                        CM_outputs[current_CM_n] = tmp_CM;
//                        CMdirection_outputs[current_CM_n] = tmp_direction;
                        CM_outputs_date[current_CM_n] = (short int)(clrx[i] - ORDINAL_LANDSAT4_LAUNCH);
                    }
                }

//                if (clrx[i] > 731426 - 1){
//                    printf("clry[4][i + adj_conse] = %f\n", clry[4][i + adj_conse - 1]);
//                }
                if ((break_mag > adj_TCG) && (mean_angle < NSIGN))
                {

//                    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
//                        for(j = 0; j < MAX_NUM_C; j++)
//                            printf("%f\n", fit_cft[k][j]);
//                    if (verbose)
//                    {
//                        printf("Change Magnitude = %.2f\n", break_mag - adj_TCG);
//                    }

                    /**********************************************/
                    /*                                            */
                    /* Record break time.                        */
                    /*                                            */
                    /**********************************************/

                    rec_cg[*num_fc].t_break = clrx[i];
                    rec_cg[*num_fc].change_prob = 100;
//                    if ((i_start - 1 + adj_conse) < end)
//                        rec_cg[*num_fc].t_confirmed = clrx[i + adj_conse - 1];
//                    else
//                        rec_cg[*num_fc].t_confirmed = clrx[end - 1];

                    //rec_cg[*num_fc].t_confirmed = clrx[i + adj_conse - 1];

                    /* check if it is a lasso band               */
                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        quick_sort_float(v_dif_mag[i_b], 0, adj_conse-1);
                        matlab_2d_float_median(v_dif_mag, i_b, adj_conse,
                                               &rec_cg[*num_fc].magnitude[i_b]);
                    }
                    /**********************************************/
                    /*                                            */
                    /* Identified and move on for the next        */
                    /* functional curve.                          */
                    /*                                            */
                    /**********************************************/

                    *num_fc = *num_fc + 1;

                    if (*num_fc >= NUM_FC)
                    {
                        /******************************************/
                        /*                                        */
                        /* Reallocate memory for rec_cg.          */
                        /*                                        */
                        /******************************************/

                        rec_cg = realloc(rec_cg, (*num_fc + 1) * sizeof(Output_t));
                        if (rec_cg == NULL)
                        {
                            RETURN_ERROR("ERROR allocating rec_cg memory",
                                                 FUNC_NAME, FAILURE);
                        }
                    }

                    /**********************************************/
                    /*                                            */
                    /* Start from i+1 for the next functional     */
                    /* curve.                                     */
                    /*                                            */
                    /**********************************************/

                    i_start = i + 1;

                    /**********************************************/
                    /*                                            */
                    /* Start training again.                      */
                    /*                                            */
                    /**********************************************/

                    bl_train = 0;

                    i_span_skip = 0;

                }

                else if (vec_mag[0] > T_MAX_CG)  /*false change*/
                {
                    /**********************************************/
                    /*                                            */
                    /* Remove noise.                              */
                    /*                                            */
                    /**********************************************/

                    for (m = i; m < end -1; m++)
                    {
                        clrx[m] = clrx[m+1];
                        for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                            clry[b][m] = clry[b][m+1];
                    }

                    i--;   /* stay & check again after noise removal */
                    end--;
                }

           } /* end of continuous monitoring */
        }  /* end of checking basic requrirements */

        /**********************************************************/
        /*                                                        */
        /* Move forward to the i+1th clear observation.           */
        /*                                                        */
        /**********************************************************/

        i++;

    } /* end of "while (i <= end - adj_conse) */

    /**************************************************************/
    /*                                                            */
    /* Two ways for processing the end of the time series.        */
    /*                                                            */
    /**************************************************************/
    // printf("main part finished \n");
    if (bl_train == 1)
    {

        /**********************************************************/
        /*                                                        */
        /* If no break, find at the end of the time series,       */
        /* define probability of change based on adj_conse.           */
        /*                                                        */
        /**********************************************************/

        id_last = adj_conse;
        for (i_conse = adj_conse - 1; i_conse >= 0; i_conse--)
        {
            v_diff_tmp =(float **) allocate_2d_array(NUM_LASSO_BANDS, adj_conse - i_conse, sizeof (float));
            for (k = 0; k < NUM_LASSO_BANDS; k++)
                for(j = 0; j < adj_conse - i_conse; j++)
                  v_diff_tmp[k][j] = v_diff[k][i_conse+j];

            mean_angle = MeanAngl_float(v_diff_tmp, NUM_LASSO_BANDS, adj_conse - i_conse);
            //float tt = vec_mag[i_conse];
            if ((vec_mag[i_conse] <= adj_TCG)||(mean_angle >= NSIGN))
            {
                /**************************************************/
                /*                                                */
                /* The last stable ID.                            */
                /*                                                */
                /**************************************************/

                id_last = i_conse + 1;
                status = free_2d_array ((void **) v_diff_tmp);
                if (status != SUCCESS)
                {
                    RETURN_ERROR ("Freeing memory: v_diff_tmp\n",
                                  FUNC_NAME, FAILURE);
                }
                break;
            }
            status = free_2d_array ((void **) v_diff_tmp);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: v_diff_tmp\n",
                              FUNC_NAME, FAILURE);
            }
        }

        /**********************************************************/
        /*                                                        */
        /* Update change probability, end time of the curve.      */
        /*                                                        */
        /**********************************************************/

        rec_cg[*num_fc].change_prob = (short int)((double)(adj_conse - id_last) * 100.0 / (double)adj_conse) ;
//        rec_cg[*num_fc].t_confirmed = 0;
        rec_cg[*num_fc].t_end = clrx[end -1 - adj_conse + id_last];    // 11/18/2018 SY

        /**********************************************************/
        /*                                                        */
        /* Mean value fit for the rest of the pixels < adj_conse & > 1*/
        /*                                                        */
        /**********************************************************/

        if (adj_conse > id_last)
        {
            /******************************************************/
            /*                                                    */
            /* Update time of the probable change.                */
            /*                                                    */
            /******************************************************/

            rec_cg[*num_fc].t_break = clrx[end - adj_conse+id_last];

            /******************************************************/
            /*                                                    */
            /* Update magnitude of change.                        */
            /*                                                    */
            /******************************************************/

            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
                quick_sort_float(v_dif_mag[i_b], id_last, adj_conse-1);
                //printf("%f\n", v_dif_mag[i_b][adj_conse-1]);
                matlab_float_2d_partial_median(v_dif_mag, i_b, id_last, adj_conse-1,
                                     &rec_cg[*num_fc].magnitude[i_b]);
            }
         }
    }
    else if (bl_train == 0)
    {
        /**********************************************************/
        /*                                                        */
        /* If break found close to the end of the time series,    */
        /* use [adj_conse,MIN_NUM_C*N_TIMES+adj_conse) to fit curve.      */
        /*                                                        */
        /* Update i_start.                                        */
        /*                                                        */
        /**********************************************************/
        if (*num_fc == rec_fc)
        {
            /******************************************************/
            /*                                                    */
            /* First curve.                                       */
            /*                                                    */
            /******************************************************/

            i_start = 1;
        }
        else
        {
            for (k = 0; k < n_clr; k++)
            {
                 if (clrx[k] >= rec_cg[*num_fc-1].t_break)
                 {
                       i_start = k + 1;
                       break;
                 }
             }
        }

        for (m = 0; m < n_clr; m++)
        {
            bl_ids[m] = 0;
        }

        if ((end - i_start + 1) >= LASSO_MIN) //04/02/2019 change adj_conse to CONSE_END
        {
            /******************************************************/
            /*                                                    */
            /* Multitemporal cloud mask.                          */
            /*                                                    */
            /******************************************************/

            status = auto_mask(clrx, clry, i_start-1, end-1,
                           (float)(clrx[end-1]-clrx[i_start-1]) / NUM_YEARS,
                           min_rmse[1], min_rmse[4], (float)T_CONST, bl_ids);
            if (status != SUCCESS)
                RETURN_ERROR("ERROR calling auto_mask at the end of time series",
                              FUNC_NAME, FAILURE);

            /******************************************************/
            /*                                                    */
            /* Clears the IDs buffers.                            */
            /*                                                    */
            /******************************************************/

            for (m = 0; m < n_clr-1; m++)
            {
                ids[m] = 0;
            }

            /******************************************************/
            /*                                                    */
            /* IDs to be removed.                                 */
            /*                                                    */
            /******************************************************/

            for (k = i_start-1; k < end; k++)
            {
                ids[k-i_start+1] = k;
            }
            m= 0;
            i_span = 0;
            for (k = 0; k < end-i_start+1; k++)
            {
                if (bl_ids[k] == 1)
                {
                    rm_ids[m] = ids[k];
                    m++;
                }
                else
                    i_span++;  /* update i_span after noise removal */
            }
            rm_ids_len = m;

            /******************************************************/
            /*                                                    */
            /* Remove noise pixels between i_start & i .          */
            /*                                                    */
            /******************************************************/

            m = 0;
            for (k = 0,k_new=0; k < end; k++) /* 03192019 SY*/
            {
                if (m < rm_ids_len && k == rm_ids[m])
                {
                    m++;
                    continue;
                }
                clrx[k_new] = clrx[k];
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    clry[i_b][k_new] = clry[i_b][k];
                k_new++;
            }
            end = k_new;
         }

         if ((end - i_start + 1) >= LASSO_MIN)     // 09/28/2018 SY delete equal sign //11/15/2018 put back equal sign
         {
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
//                    for(k=i_start-1;k<end;k++)
//                    {
//                        printf("clrx %d: %d\n", k,clrx[k]);
//                        printf("clry %d: %f\n", k,clry[i_b][k]);
//                    }
                status = auto_ts_fit_float(clrx, clry, i_b, i_b, i_start-1, end-1, MIN_NUM_C,
                                     fit_cft, &rmse[i_b], temp_v_dif);
                if (status != SUCCESS)
                {
                     RETURN_ERROR ("Calling auto_ts_fit_float at the end of time series\n",
                            FUNC_NAME, FAILURE);
                }

            }

            /******************************************************/
            /*                                                    */
            /* Record time of curve start, time of curve end,     */
            /* break time, postion of the pixel.                  */
            /*                                                    */
            /******************************************************/

            if (*num_fc == rec_fc)
            {
                    rec_cg[*num_fc].t_start = clrx[0];
            }
            else
            {
                    rec_cg[*num_fc].t_start = rec_cg[*num_fc-1].t_break;
            }
            rec_cg[*num_fc].t_end = clrx[end-1];
            rec_cg[*num_fc].t_break = 0;
            //rec_cg[*num_fc].pos.row = row;
            //rec_cg[*num_fc].pos.col = col;




            /******************************************************/
            /*                                                    */
            /* Record change probability, number of observations, */
            /* fit category.                                      */
            /*                                                    */
            /******************************************************/

            rec_cg[*num_fc].change_prob = 0;
//            rec_cg[*num_fc].t_confirmed = 0;
            rec_cg[*num_fc].num_obs = end - i_start + 1;
            rec_cg[*num_fc].category = 20 + MIN_NUM_C; /* simple model fit at the end */

            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {

                /******************************************************/
                /*                                                    */
                /* Record fitted coefficients.                        */
                /*                                                    */
                /******************************************************/

                for (k = 0; k < LASSO_COEFFS; k++)
                {
                    rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                }
                rec_cg[*num_fc].rmse[i_b] = rmse[i_b];

                /******************************************************/
                /*                                                    */
                /* Record change magnitude.                           */
                /*                                                    */
                /******************************************************/

                rec_cg[*num_fc].magnitude[i_b] = 0.0;


             }

            //*num_fc = *num_fc + 1;
            if (*num_fc >= NUM_FC)
            {
                /**************************************************/
                /*                                                */
                /* Reallocate memory for rec_cg.                  */
                /*                                                */
                /**************************************************/

                rec_cg = realloc(rec_cg, (*num_fc + 1) * sizeof(Output_t));
                if (rec_cg == NULL)
                {
                    RETURN_ERROR("ERROR allocating rec_cg memory",
                                 FUNC_NAME, FAILURE);
                }
            }
        }else{
           *num_fc = *num_fc - 1;
        }
    }

    *num_fc = *num_fc + 1;


    /******************************************************************/
    /*                                                                */
    /* Free memory allocations for this section.                      */
    /*                                                                */
    /******************************************************************/

    // for debug
    // printf("free stage 1 \n");
    status = free_2d_array ((void **) fit_cft);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fit_cft\n", FUNC_NAME,
                      FAILURE);
    }

    free(rmse);

    free(clrx);

    // for debug
    // printf("free stage 1 \n");

    status = free_2d_array ((void **) clry);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: clry\n", FUNC_NAME,
                      FAILURE);
    }

    status = free_2d_array ((void **) temp_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: temp_v_dif\n",
                      FUNC_NAME, FAILURE);
    }

    // for debug
    // printf("free stage 2 \n");

    status = free_2d_array ((void **) rec_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: rec_v_dif\n",
                      FUNC_NAME, FAILURE);
    }
    status = free_2d_array ((void **) rec_v_dif_copy);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: \n",
                      FUNC_NAME, FAILURE);
    }


    // for debug
    // printf("free stage 3 \n");

    free(ids);
    free(ids_old);
    free(bl_ids);
    free(rm_ids);


    free(v_start);
    free(v_end);
    free(v_slope);
    free(v_dif);

    free(tmpcg_rmse);

    // for debug
    // printf("free stage 4 \n");

    status = free_2d_array ((void **) v_dif_mag);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: v_dif_mag\n",
                   FUNC_NAME, FAILURE);
    }
    free(vec_mag);

    status = free_2d_array ((void **) v_diff);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: v_diff\n",
                      FUNC_NAME, FAILURE);
    }

    // for debug
    // printf("free stage 5 \n");

    /******************************************************************/
    /*                                                                */
    /* Output rec_cg structure to the output file.                    */
    /* Note: can use fread to read out the structure from the output  */
    /* file.                                                          */
    /* If output was stdout, skip this step.                          */
    /*                                                                */
    /******************************************************************/


    return (SUCCESS);
}

/******************************************************************************
MODULE:  inefficientobs_procedure

PURPOSE:  the procedure for inefficient clear pixels

RETURN VALUE:
Type = int (SUCCESS OR FAILURE)

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
11/14/2018   Su Ye          Modification from main function in original CCDC.c
******************************************************************************/
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
    int *num_fc,
    bool b_c2
)
{
    int n_sn = 0;
    int i, k;
    int end;
    int i_start;
    int i_span;
    int *clrx;                       /* clear pixel curve in X direction (date)    */
    float **clry;                    /* clear pixel curve in Y direction (spectralbands)    */
    float **fit_cft;                 /* Fitted coefficients 2-D array.        */
    float *rmse;                     /* Root Mean Squared Error array.        */
    float **temp_v_dif;              /* for the thermal band.......           */
    char FUNC_NAME[] = "inefficientobs_procedure";
    int n_clr = 0;
    int i_b;
    int status;

    clrx = malloc(valid_num_scenes * sizeof(int));
    if (clrx == NULL)
    {
        RETURN_ERROR("ERROR allocating clrx memory", FUNC_NAME, FAILURE);
    }

    clry = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, valid_num_scenes,
                                         sizeof (float));
    if (clry == NULL)
    {
        RETURN_ERROR ("Allocating clry memory", FUNC_NAME, FAILURE);
    }

    temp_v_dif = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, valid_num_scenes,
                                     sizeof (float));
    if (temp_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating temp_v_dif memory",FUNC_NAME, FAILURE);
    }

    fit_cft = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, LASSO_COEFFS,
                                         sizeof (float));
    if (fit_cft == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft memory", FUNC_NAME, FAILURE);
    }

    rmse = (float *)calloc(TOTAL_IMAGE_BANDS, sizeof(float));
    if (rmse == NULL)
    {
        RETURN_ERROR ("Allocating rmse memory", FUNC_NAME, FAILURE);
    }

    if (sn_pct > T_SN)
    {
        /**********************************************************/
        /*                                                        */
        /* Snow observations are "good" now.                      */
        /*                                                        */
        /**********************************************************/

        for (i = 0; i < valid_num_scenes; i++)
        {
        if ((fmask_buf[i] == CFMASK_SNOW) || (fmask_buf[i] < 2))
            {
                clrx[n_sn] = valid_date_array[i];
                for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                {
                    if (k == 0)
                        clry[k][n_sn] = (float)buf_b[i];
                    else if(k == 1)
                        clry[k][n_sn] = (float)buf_g[i];
                    else if(k == 2)
                        clry[k][n_sn] = (float)buf_r[i];
                    else if(k == 3)
                        clry[k][n_sn] = (float)buf_n[i];
                    else if(k == 4)
                        clry[k][n_sn] = (float)buf_s1[i];
                    else if(k == 5)
                        clry[k][n_sn] = (float)buf_s2[i];
                    else if(k == 6){
                        if (b_c2 == TRUE){
                            clry[k][n_sn] = 0;
                        }
                        else{
                            clry[k][n_sn] = (float)(buf_t[i] * 10 - 27320);
                        }
                   }
                }
                n_sn++;
            }
        }
        end = n_sn;

//        for (i = 0; i < n_sn; i++)
//        {
//           printf("thermal: %f\n", clry[TOTAL_IMAGE_BANDS-1][i]);
//           //printf("clry: %f\n", (float)clry[k][i]);
//        }

//        /**************************************************************/
//        /*                                                            */
//        /* Remove repeated ids.                                       */
//        /*                                                            */
//        /**************************************************************/

//        matlab_unique(clrx, clry, n_clr, &end);

        if (n_sn < N_TIMES * MIN_NUM_C) // not enough snow pixels
        {
            free(clrx);
            status = free_2d_array((void **) clry);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: clry\n",
                              FUNC_NAME, FAILURE);
            }
            free(rmse);
            status = free_2d_array ((void **) temp_v_dif);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: temp_v_dif\n",
                              FUNC_NAME, FAILURE);
            }
            status = free_2d_array ((void **) fit_cft);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: fit_cft\n", FUNC_NAME,
                              FAILURE);
            }
            return (SUCCESS);
        }

        /**********************************************************/
        /*                                                        */
        /* Start model fit for snow persistent pixels.            */
        /*                                                        */
        /**********************************************************/

//        if (verbose)
//            printf ("Fit permanent snow observations, now pixel = %f\n",
//                   100.0 * sn_pct);

        i_start = 1; /* the first observation for TSFit */

        /**********************************************************/
        /*                                                        */
        /* Treat saturated and unsaturated pixels differently.    */
        /*                                                        */
        /**********************************************************/


        for (k = 0; k < TOTAL_IMAGE_BANDS; k++)  //
        {
            i_span = 0;
            if (k != TOTAL_IMAGE_BANDS - 1) // for optical bands sy 10/01/2018
            {
                for (i = 0; i < end; i++)
                {
                    if (clry[k][i] > 0.0 && clry[k][i] < 10000.0)
                    {
                        clrx[i_span] = clrx[i];
                        clry[k][i_span] = clry[k][i];
                        i_span++;
                    }
                 }

                if (i_span < MIN_NUM_C * N_TIMES){
                    fit_cft[k][0] = 10000; // fixed value for saturated pixels
                    for(i = 1; i < MAX_NUM_C; i++)
                        fit_cft[k][i] = 0;
                }
                else
                {
                    status = auto_ts_fit_float(clrx, clry, k, k, 0, i_span-1, MIN_NUM_C,
                             fit_cft, &rmse[k], temp_v_dif);

                    if (status != SUCCESS)
                        RETURN_ERROR ("Calling auto_ts_fit_float\n",
                               FUNC_NAME, EXIT_FAILURE);
                }
             }
             else
             {
                for (i = 0; i < end; i++)
                {
                    //printf("thermal: %f\n", clry[k][i]);
                    //printf("clry: %f\n", (float)clry[k][i]);
                    if (clry[k][i] > -9320 && clry[k][i] < 7070)
                    {
                        clrx[i_span] = clrx[i];
                        clry[k][i_span] = clry[k][i];
                        i_span++;
                    }
                 }

                status = auto_ts_fit_float(clrx, clry, k, k, 0, i_span-1, MIN_NUM_C,
                         fit_cft, &rmse[k], temp_v_dif);

                if (status != SUCCESS)
                    RETURN_ERROR ("Calling auto_ts_fit_float\n",
                           FUNC_NAME, EXIT_FAILURE);

              }
        }


        /**********************************************************/
        /*                                                        */
        /*                                                        */
        /**********************************************************/

        rec_cg[*num_fc].t_start = clrx[i_start-1];
        rec_cg[*num_fc].t_end = clrx[end-1];

        /**********************************************************/
        /*                                                        */
        /* No break at the moment.                                */
        /*                                                        */
        /**********************************************************/

        rec_cg[*num_fc].t_break = 0;

        /**********************************************************/
        /*                                                        */
        /* Record postion of the pixel.                           */
        /*                                                        */
        /**********************************************************/

        //rec_cg[*num_fc].pos.row = row;
        //rec_cg[*num_fc].pos.col = col;

        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
        {
            for (k = 0; k < LASSO_COEFFS; k++)
            {
                /**************************************************/
                /*                                                */
                /* Record fitted coefficients.                    */
                /*                                                */
                /**************************************************/

                rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
            }

            /******************************************************/
            /*                                                    */
            /* Record rmse of the pixel.                          */
            /*                                                    */
            /******************************************************/

            rec_cg[*num_fc].rmse[i_b] = rmse[i_b];
        }

        /**********************************************************/
        /*                                                        */
        /* Record change probability, number of observations.     */
        /*                                                        */
        /**********************************************************/

        rec_cg[*num_fc].change_prob = 0;
//        rec_cg[*num_fc].t_confirmed = 0;
        rec_cg[*num_fc].num_obs = n_sn;
        rec_cg[*num_fc].category = 50 + MIN_NUM_C; /* snow pixel */

        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
        {
            /******************************************************/
            /*                                                    */
            /* Record change magnitude.                           */
            /*                                                    */
            /******************************************************/

            rec_cg[*num_fc].magnitude[i_b] = 0.0;
        }


        if (*num_fc >= NUM_FC)
        {

            /******************************************************/
            /*                                                    */
            /* Reallocate memory for rec_cg.                      */
            /*                                                    */
            /******************************************************/

            rec_cg = realloc(rec_cg, (*num_fc + 1) * sizeof(Output_t));
            if (rec_cg == NULL)
            {
                RETURN_ERROR("ERROR allocating rec_cg memory",
                             FUNC_NAME, FAILURE);
            }
        }
    }  // if sn_pct > T_SN

    else
    {

        /**********************************************************/
        /*                                                        */
        /* normal inefficient observation procedure.              */
        /*                                                        */
        /**********************************************************/

        n_clr = 0;

        for (i = 0; i < valid_num_scenes; i++)
        {
            if (id_range[i] == 1)
            {
                clrx[n_clr] = valid_date_array[i];
                for (k = 0; k < TOTAL_IMAGE_BANDS; k++){
                    if (k == 0)
                        clry[k][n_clr] = (float)buf_b[i];
                    else if(k == 1)
                        clry[k][n_clr] = (float)buf_g[i];
                    else if(k == 2)
                        clry[k][n_clr] = (float)buf_r[i];
                    else if(k == 3)
                        clry[k][n_clr] = (float)buf_n[i];
                    else if(k == 4)
                        clry[k][n_clr] = (float)buf_s1[i];
                    else if(k == 5)
                        clry[k][n_clr] = (float)buf_s2[i];
                    else if(k == 6){
                        if (b_c2 == TRUE){
                            clry[k][n_sn] = 0;
                        }
                        else{
                            clry[k][n_sn] = (float)(buf_t[i] * 10 - 27320);
                        }
                   }
                }
                n_clr++;
            }
        }
        end = n_clr;

        /**************************************************************/
        /*                                                            */
        /* Remove repeated ids.                                       */
        /*                                                            */
        /**************************************************************/

        matlab_unique(clrx, clry, n_clr, &end);


        n_clr = 0;
        float band2_median; // probably not good practice to declare here....
        quick_sort_float(clry[1], 0, end - 1);
        matlab_2d_float_median(clry, 1, end, &band2_median);

        n_clr = 0;
        for (i = 0; i < end; i++)
        {
            if (clry[1][i] < (band2_median + 400.0))
            {
                clrx[n_clr] = clrx[i];
                for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                {
                     clry[k][n_clr] = clry[k][i];
                }
                n_clr++;
            }
        }
        end = n_clr;


        /**********************************************************/
        /*                                                        */
        /* The first observation for TSFit.                       */
        /*                                                        */
        /**********************************************************/

        i_start = 1; /* the first observation for TSFit */

        if (n_clr < N_TIMES * MIN_NUM_C)
        {
        // num_fc = 0, so won't output any curve
            *num_fc = *num_fc - 1;
            WARNING_MESSAGE("Not enough good clear observations (<12obs)\n", FUNC_NAME);
        }
        else
        {
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
                status = auto_ts_fit_float(clrx, clry, i_b, i_b, 0, end-1, MIN_NUM_C,
                                     fit_cft, &rmse[i_b], temp_v_dif);
                if (status != SUCCESS)
                {
                    RETURN_ERROR ("Calling auto_ts_fit_float for clear persistent pixels\n",
                                  FUNC_NAME, FAILURE);
                }
            }
            /**********************************************************/
            /*                                                        */
            /* Update information at each iteration.                  */
            /* Record time of curve start, time of curve end.         */
            /*                                                        */
            /**********************************************************/

            rec_cg[*num_fc].t_start = clrx[i_start-1];
            rec_cg[*num_fc].t_end = clrx[end-1];

            /**********************************************************/
            /*                                                        */
            /* No break at the moment.                                */
            /*                                                        */
            /**********************************************************/

            rec_cg[*num_fc].t_break = 0;

            /**********************************************************/
            /*                                                        */
            /* Record postion of the pixel.                           */
            /*                                                        */
            /**********************************************************/

            //rec_cg[*num_fc].pos.row = row;
            //rec_cg[*num_fc].pos.col = col;

            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
                for (k = 0; k < LASSO_COEFFS; k++)
                {
                    /**************************************************/
                    /*                                                */
                    /* Record fitted coefficients.                    */
                    /*                                                */
                    /**************************************************/

                    rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                }

                /******************************************************/
                /*                                                    */
                /* Record rmse of the pixel.                          */
                /*                                                    */
                /******************************************************/

                rec_cg[*num_fc].rmse[i_b] = rmse[i_b];
            }

            /**********************************************************/
            /*                                                        */
            /* Record change probability, number of observations,     */
            /* fit category.                                          */
            /*                                                        */
            /**********************************************************/
            rec_cg[*num_fc].change_prob = 0;
//            rec_cg[*num_fc].t_confirmed = 0;
            rec_cg[*num_fc].num_obs = n_clr;
            rec_cg[*num_fc].category = 40 + MIN_NUM_C;

            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
                /******************************************************/
                /*                                                    */
                /* Record change magnitude.                           */
                /*                                                    */
                /******************************************************/
                rec_cg[*num_fc].magnitude[i_b] = 0.0;
            }
        }

    }
    *num_fc = *num_fc + 1;

    free(clrx);
    status = free_2d_array((void **) clry);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: clry\n",
                      FUNC_NAME, FAILURE);
    }
    free(rmse);
    status = free_2d_array ((void **) temp_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: temp_v_dif\n",
                      FUNC_NAME, FAILURE);
    }
    status = free_2d_array ((void **) fit_cft);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fit_cft\n", FUNC_NAME,
                      FAILURE);
    }

    return (SUCCESS);

}




/******************************************************************************
MODULE:  obcold_reconstruction_procedure
PURPOSE:  standard procedure when having enough clear pixels
RETURN VALUE:
Type = int (SUCCESS OR FAILURE)
HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
02/14/2021   Su Ye          rescontructing rec_cg using breakdate
******************************************************************************/
int obcold_reconstruction_procedure
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
    long *break_dates, /*an array of break dates with a fixed length of num_year, '0' means no breaks */
    int break_date_len,       /*I: the length of break_dates */
    int pos,              /*I: the position of the pixel */
    bool b_c2,                  /* I: a temporal parameter to indicate if collection 2. C2 needs ignoring thermal band due to the current low quality  */
    int conse,
    Output_t *rec_cg,    /* O: Initialize NUM of Functional Curves    */
    int *num_fc
)
{
    int status;
    int i, j, k, k_new, b;
    int m;                       /*the number of ID to be removed */
    char FUNC_NAME[] = "obcold_reconstruction_procedure";
    //char msg_str[MAX_STR_LEN];       /* Input data scene name                 */
    int *rm_ids;

    int end;
    float **fit_cft;                 /* Fitted coefficients 2-D array.        */
    float **rec_v_dif;
    float v_dif_norm = 0.0;
    float mini_rmse;                 /* Mimimum RMSE                          */

    int i_b;

    float ts_pred_temp;
    int *bl_ids;
    int num_c = 8;                   /* Max number of coefficients for model  */
    int update_num_c = 8;            /* Number of coefficients to update      */

    float min_rmse[TOTAL_IMAGE_BANDS]; /* Adjusted RMSE for all bands          */


    int n_clr;                  /* I: the number of clear pixels                          */
    int *clrx;                  /* I: clear pixel curve in X direction (date)             */
    float **clry;               /* I: clear pixel curve in Y direction (spectralbands)    */
    int *clrx_tmp;                        /* nunber of clear pixels X ?            */
    float **clry_tmp;                     /* nunber of clear pixels Y ?            */
    int rm_ids_len;
    float* v_dif;
    float* v_dif_allbands;
    int *id_range;
    float date_vario;           /* I: median date                                          */
    float max_date_difference;   /* I: maximum difference between two neighbor dates        */
    int clear_sum = 0;      /* Total number of clear cfmask pixels          */
    int water_sum = 0;      /* counter for cfmask water pixels.             */
    int shadow_sum = 0;     /* counter for cfmask shadow pixels.            */
    int sn_sum = 0;         /* Total number of snow cfmask pixels           */
    int cloud_sum = 0;      /* counter for cfmask cloud pixels.             */
    float sn_pct;
    float* rmse;
    int end_copy;
    int result;
    int *break_list;
    int i_break;
    int i_break_tmp;

    int break_num = 0;
    int i_last_break = 0;

    int adj_conse = conse;
    float **v_dif_mag;
    int last_break_record = 0;  // used to record the year of the last break

    if(valid_num_scenes == 0){
        return (SUCCESS);
    }
    id_range = (int*)calloc(valid_num_scenes, sizeof(int));

    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
    {
            min_rmse[k] = 0.0;
    }

    fit_cft = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, LASSO_COEFFS,
                                         sizeof (float));
    if (fit_cft == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft memory", FUNC_NAME, FAILURE);
    }

    rmse = (float *)calloc(TOTAL_IMAGE_BANDS, sizeof(float));
    if (rmse == NULL)
    {
        RETURN_ERROR ("Allocating rmse memory", FUNC_NAME, FAILURE);
    }

    v_dif = (float *)calloc(NUM_LASSO_BANDS, sizeof(float));
    if (v_dif == NULL)
    {
        RETURN_ERROR("ERROR allocating v_dif memory", FUNC_NAME, FAILURE);
    }

    v_dif_allbands = (float *)calloc(TOTAL_BANDS, sizeof(float));
    if (v_dif_allbands == NULL)
    {
        RETURN_ERROR("ERROR allocating v_dif memory", FUNC_NAME, FAILURE);
    }

    n_clr = 0;

    // check if validate_date_array follow a ascending order
    for (i = 0; i < valid_num_scenes - 1; i++)
        if (valid_date_array[i] > valid_date_array[i+1]){
            RETURN_ERROR("The inputted data does not follow an ascending order!", FUNC_NAME, ERROR);
        }


    clrx = (int* )malloc(valid_num_scenes * sizeof(int));
    clry = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, valid_num_scenes,
                                         sizeof (float));


    clrx_tmp = (int* )malloc(valid_num_scenes * sizeof(int));
    clry_tmp = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, valid_num_scenes,
                                         sizeof (float));



    status = preprocessing(buf_b, buf_g, buf_r, buf_n, buf_s1, buf_s2, buf_t,
                           fmask_buf, &valid_num_scenes, id_range, &clear_sum,
                           &water_sum, &shadow_sum, &sn_sum, &cloud_sum, b_c2);
    if (status != SUCCESS)
    {
        RETURN_ERROR("Error for preprocessing.", FUNC_NAME, ERROR);
    }

    for (i = 0; i < valid_num_scenes; i++)
    {
        if ((fmask_buf[i] < 2) && (id_range[i] == 1))
        {
            // remain the first element for replicated date
            if((n_clr > 0) && (valid_date_array[i] == clrx[n_clr - 1]))
                continue;
            else
            {
                clrx[n_clr] = valid_date_array[i];
                //printf("%d is %d\n", n_clr + 1, clrx[n_clr]);
                for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                {
                    if (k == 0)
                        clry[k][n_clr] = (float)buf_b[i];
                    else if (k == 1)
                        clry[k][n_clr] = (float)buf_g[i];
                    else if (k == 2)
                        clry[k][n_clr] = (float)buf_r[i];
                    else if (k == 3)
                        clry[k][n_clr] = (float)buf_n[i];
                    else if (k == 4)
                        clry[k][n_clr] = (float)buf_s1[i];
                    else if (k == 5)
                        clry[k][n_clr] = (float)buf_s2[i];
                    else if (k == 6)
                        clry[k][n_clr] = (float)(buf_t[i] * 10 - 27320);
                    //printf("%3.2f\n", clry[k][n_clr]);
                }
                n_clr++;
            }
        }
    }

    rec_v_dif = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, n_clr,
                                     sizeof (float));
    if (rec_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating rec_v_dif memory",FUNC_NAME, FAILURE);
    }
    bl_ids = (int *)calloc(n_clr, sizeof(int));
    if (bl_ids == NULL)
    {
        RETURN_ERROR("ERROR allocating bl_ids memory", FUNC_NAME, FAILURE);
    }

    rm_ids = (int *)calloc(n_clr, sizeof(int));
    if (rm_ids == NULL)
    {
        RETURN_ERROR("ERROR allocating rm_ids memory", FUNC_NAME, FAILURE);
    }


    end = n_clr;

    break_list = (int* )malloc(break_date_len * sizeof(int));

    if (clear_sum < N_TIMES * MAX_NUM_C){
        result = inefficientobs_procedure(valid_num_scenes,valid_date_array, buf_b, buf_g, buf_r, buf_n, buf_s1, buf_s2, buf_t,
                                 fmask_buf, id_range,sn_pct,rec_cg, num_fc, b_c2);
    }else
    {
        /**************************************************************/
        /*                                                            */
        /* calculate variogram for each band and dates.                           */
        /*                                                            */
        /**************************************************************/
        status = adjust_median_variogram(clrx, clry, TOTAL_IMAGE_BANDS, 0, end-1, &date_vario,
                                         &max_date_difference, min_rmse, 1);
        if (status != SUCCESS)
        {
                RETURN_ERROR("ERROR calling median_variogram routine", FUNC_NAME,
                             FAILURE);
        }


        /* sort the date array first in case that they are not sorted */
        quick_sort_long(break_dates, 0, break_date_len-1);

        // guarantee that the interval of two breaks has at least 365 days
        for(i = 0; i< break_date_len; i++){
            if((break_dates[i] > 0) & (break_dates[i] - last_break_record >= AVE_DAYS_IN_A_YEAR)){
               break_list[break_num] = (int)break_dates[i];
               last_break_record = break_dates[i];
               break_num ++;
            }
        }


        for(j = 0; j < break_num + 1; j++) // break_num + 1 is because the segment = break_num + 1
        {
            if (j == break_num) // last break, using last_obs + 1 as break
                i_break_tmp = end;
            else
                i_break_tmp = find_index_clrx(clrx, end, break_list[j]);


            if(i_break_tmp - 1 - i_last_break + 1 < LASSO_MIN)
                continue;

            // if 6 <=number of observations < 12, no needs to filter out outliers, direct fitting
            if(i_break_tmp - 1 - i_last_break + 1 >= MIN_NUM_C * N_TIMES && j < break_num)
            {
                end_copy = end;

                /*****************************************************/
                /*         test 1:    t-mask                        */
                /****************************************************/
                status = auto_mask(clrx, clry, i_last_break, i_break_tmp - 1,
                                    (float)(clrx[i_break_tmp - 1]-clrx[i_last_break] + 1) / NUM_YEARS,
                                    min_rmse[1], min_rmse[4], (float)T_CONST, bl_ids);
                // printf("auto_mask finished \n");
                if (status != SUCCESS)
                {
                    RETURN_ERROR("ERROR calling auto_mask during model initilization",
                                  FUNC_NAME, FAILURE);
                }

                /**************************************************/
                /*                                                */
                /* IDs to be removed.                             */
                /*                                                */
                /**************************************************/
                rm_ids_len = 0;

                for (k = i_last_break; k < i_break_tmp; k++) /** 02282019 SY **/
                {
                    if (bl_ids[k - i_last_break] == 1)
                    {
                        rm_ids[rm_ids_len] = k;
                        rm_ids_len++;
                    }
                }

                if(i_break_tmp - 1 - i_last_break + 1 - rm_ids_len < LASSO_MIN) // first test not passsed
                {
                    end = end_copy;
                    continue;
                }

                // update clrx and clry and save them to tmp array
                m = 0;
                for (k = 0, k_new=0; k < end; k++)
                {
                    if (m < rm_ids_len && k == rm_ids[m])
                    {
                        m++;
                        continue;
                    }
                    clrx_tmp[k_new] = clrx[k];
                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        clry_tmp[i_b][k_new] = clry[i_b][k];
                    }
                    k_new++;
                }

                // update other parameters
                end = k_new;
                i_break_tmp = i_break_tmp - rm_ids_len;


                /**************************************************/
                /* Test 2:        COLD test analysis           */
                /**************************************************/
                update_cft(i_break_tmp - 1 - i_last_break + 1, N_TIMES, MIN_NUM_C, MID_NUM_C, MAX_NUM_C,
                          num_c, &update_num_c);

                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    // no need to include blue and BT for here
                    for (b = 0; b < NUM_LASSO_BANDS; b++)
                    {
                        /**********************************************/
                        /*                                            */
                        /* Initial model fit.                         */
                        /*                                            */
                        /**********************************************/
                        if (i_b == lasso_blist[b])
                        {
                            status = auto_ts_fit_float(clrx_tmp, clry_tmp,  lasso_blist[b], lasso_blist[b], i_last_break,
                                                       i_break_tmp - 1,
                                     update_num_c, fit_cft, &rmse[i_b], rec_v_dif);
                            if (status != SUCCESS)
                            {
                                RETURN_ERROR ("Calling auto_ts_fit_float during model initilization\n",
                                     FUNC_NAME, FAILURE);
                            }
                        }

                    }
                }

                rm_ids_len = 0;
                for (k = i_last_break; k < i_break_tmp; k++)//SY 09192018
                {
                    v_dif_norm = 0.0;
                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {

                        for (b = 0; b < NUM_LASSO_BANDS; b++)
                        {
                            if (i_b == lasso_blist[b])
                            {
                                auto_ts_predict_float(clrx_tmp, fit_cft, update_num_c, i_b, k, k,
                                                &ts_pred_temp); //SY 09192018
                                mini_rmse = max(min_rmse[i_b], rmse[i_b]);
                                v_dif[b] = (clry_tmp[i_b][k] - ts_pred_temp) / mini_rmse;
                                v_dif_norm += v_dif[b]* v_dif[b];
                            }
                        }
                     }
                     if (v_dif_norm > T_MAX_CG_reconstruct){
                        rm_ids[rm_ids_len] = k;
                        rm_ids_len++;
                     }
                 }

                if(i_break_tmp - 1 - i_last_break + 1 - rm_ids_len < LASSO_MIN)
                {
                    end = end_copy;
                    continue;
                }

                m = 0;
                for (k = 0, k_new=0; k < end; k++)
                {
                    if (m < rm_ids_len && k == rm_ids[m])
                    {
                        m++;
                        continue;
                    }
                    clrx[k_new] = clrx_tmp[k];
                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        clry[i_b][k_new] = clry_tmp[i_b][k];
                    }
                    k_new++;
                }

                // update end and i_break_tmp
                end = k_new;
           }

            // the above is the ultimate test, if passed, we think it is the 'final' break point
            if (j == break_num) // last break, using last_obs + 1 as break
                i_break = end;
            else
                i_break = find_index_clrx(clrx, end, break_list[j]);


           update_cft(i_break - 1 - i_last_break + 1, N_TIMES, MIN_NUM_C, MID_NUM_C, MAX_NUM_C,
                      num_c, &update_num_c);


            /*final fitting*/
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
                status = auto_ts_fit_float(clrx, clry, i_b, i_b, i_last_break, i_break - 1, update_num_c,
                                 fit_cft, &rmse[i_b], rec_v_dif);
                // printf("auto_ts_fit2 finished \n", i);
                if (status != SUCCESS)
                {
                    RETURN_ERROR ("Calling auto_ts_fit_float for change detection with "
                         "enough observations\n", FUNC_NAME, FAILURE);
                }
            }

            rec_cg[*num_fc].t_start = clrx[i_last_break];
            rec_cg[*num_fc].t_end = clrx[i_break - 1];

            // rec_cg[*num_fc].t_confirmed = 0;
            rec_cg[*num_fc].num_obs = i_break - i_last_break + 1;
            rec_cg[*num_fc].category = 0 + update_num_c;

            /* note that each time we calculate change magnitude for the last curve, not current curve*/
//                        if (*num_fc > 0){
//                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
//                            {
//                                auto_ts_predict(clrx, fit_cft, update_num_c, i_b, i_last_break, i_last_break, &ts_pred_temp);
//                                auto_ts_predict(clrx, fit_cft_last, update_num_c, i_b, i_last_break, i_last_break, &ts_pred_temp_last);
//                                rec_cg[*num_fc - 1].magnitude[i_b] = ts_pred_temp - ts_pred_temp_last;
//                                rec_cg[*num_fc - 1].change_prob = 100;
//                            }
//                        }

            // if segment is not the last segment
            if (j != (break_num + 1) - 1)
            {

                if (end - i_break  < conse){
                    adj_conse = end - i_break;
                }else{
                    adj_conse = conse;
                }
                v_dif_mag = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, adj_conse,
                                                         sizeof (float));
                if (v_dif_mag == NULL)
                {
                    RETURN_ERROR ("Allocating v_dif_mag memory",
                                             FUNC_NAME, FAILURE);
                }

                for(i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    for (k = 0; k < adj_conse; k++){
                        auto_ts_predict_float(clrx, fit_cft, update_num_c, i_b, i_break+k,
                             i_break+k, &ts_pred_temp);
                        v_dif_mag[i_b][k] = (float)clry[i_b][i_break+k] - ts_pred_temp;
                    }
                    quick_sort_float(v_dif_mag[i_b], 0, adj_conse-1);
                    matlab_2d_float_median(v_dif_mag, i_b, adj_conse,
                                           &rec_cg[*num_fc].magnitude[i_b]);
                }


                rec_cg[*num_fc].t_break = clrx[i_break];
                rec_cg[*num_fc].change_prob = 100;
                status = free_2d_array((void **)v_dif_mag);
                if (status != SUCCESS)
                {
                    RETURN_ERROR ("Freeing memory: v_dif_mag\n",
                                  FUNC_NAME, FAILURE);
                }
            }
            // note for the last segment, we didn't assign magnitude, change probability and t_break so far

            /******************************************/
            /*                                        */
            /* Record fitted coefficients.            */
            /*                                        */
            /******************************************/
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {

                for (k = 0; k < LASSO_COEFFS; k++)
                {
                    /**********************************/
                    /*                                */
                    /* Record fitted coefficients.    */
                    /*                                */
                    /**********************************/

                    rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                }
                /**************************************/
                /*                                    */
                /* Record rmse of the pixel.          */
                /*                                    */
                /**************************************/

                rec_cg[*num_fc].rmse[i_b] = rmse[i_b];

            }

            *num_fc = *num_fc + 1;

            // finally, change i_last_break; if the last segment, no need to update
            if(j < break_num){
                i_last_break = i_break;
            }


        } // end  for(j = 0; j < break_num; j++)

        if(*num_fc == 0){ // apply inefficient procedure to force to have a model
            result = inefficientobs_procedure(valid_num_scenes,valid_date_array, buf_b, buf_g, buf_r, buf_n, buf_s1, buf_s2, buf_t,
                                     fmask_buf, id_range,sn_pct,rec_cg, num_fc, b_c2);
        }
    }

    // the last step, force to update the last segment
    if(*num_fc - 1 >= 0){
        rec_cg[*num_fc - 1].change_prob = 0;
        rec_cg[*num_fc - 1].t_break = 0;
        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
        {
            rec_cg[*num_fc - 1].magnitude[i_b] = 0;
        }
    }


    for(i = 0; i < *num_fc; i++)
    {
        rec_cg[i].pos = pos;
    }



    free(id_range);

    status = free_2d_array((void **) fit_cft);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fit_cft\n", FUNC_NAME,
                      FAILURE);
    }

    free(rmse);
    free(v_dif);
    free(v_dif_allbands);
    free(clrx);
    status = free_2d_array ((void **) clry);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: clry\n", FUNC_NAME,
                      FAILURE);
    }

    free(clrx_tmp);
    status = free_2d_array ((void **) clry_tmp);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: cpy\n", FUNC_NAME,
                      FAILURE);
    }

    status = free_2d_array ((void **) rec_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: rec_v_dif\n",
                      FUNC_NAME, FAILURE);
    }

    free(bl_ids);
    free(rm_ids);
    free(break_list);
    return (SUCCESS);
}




