#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "math.h"
#include "defines.h"
#include "misc.h"
#include "const.h"
#include "utilities.h"
#include "2d_array.h"
#include "output.h"
#include "s_ccd.h"
#include "distribution_math.h"
#define N_Columns 8
#define Max_DOF 30

int DEFAULT_FOCUS_BLIST[DEFAULT_N_FOCUS_VARIABLE] = {1, 2, 3, 4, 5};

//#define safe_xgboost(call) {                                            \
//int err = (call);                                                       \
//if (err != 0) {                                                         \
//  fprintf(stderr, "%s:%d: error in %s: %s\n", __FILE__, __LINE__, #call, XGBGetLastError()); \
//  exit(1);                                                              \
//}                                                                       \
//}

/******************************************************************************
MODULE:  sccd

PURPOSE:  main function for stochastic ccd

RETURN VALUE:
Type = int (SUCCESS OR FAILURE)

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
11/14/2018   Su Ye         Original Development
******************************************************************************/

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
    int monitorwindow_lowerlin,  /* for training process*/
    int monitorwindow_upperlim,   /* for training process*/
    short int *sensor_buf,
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
    bool b_landspecific,
    short int auxval,
    int conse
)
{
    int clear_sum = 0;      /* Total number of clear cfmask pixels          */
    int water_sum = 0;      /* counter for cfmask water pixels.             */
    int shadow_sum = 0;     /* counter for cfmask shadow pixels.            */
    int sn_sum = 0;         /* Total number of snow cfmask pixels           */
    int cloud_sum = 0;      /* counter for cfmask cloud pixels.             */
    double sn_pct;           /* Percent of snow pixels.                      */
    int status;
    int *id_range;
    int i, j, k;
    char FUNC_NAME[] = "sccd";
    int result;

    id_range = (int*)calloc(valid_num_scenes, sizeof(int));

    status = preprocessing(buf, fmask_buf, &valid_num_scenes, id_range, &clear_sum,
                           &water_sum, &shadow_sum, &sn_sum, &cloud_sum);


    if (status != SUCCESS)
    {
        RETURN_ERROR("Error for preprocessing.", FUNC_NAME, ERROR);
    }

    // clear_pct is not used anymore in V13.01
    //clr_pct = (double) clear_sum / (double) (valid_num_scenes);

    sn_pct = (double) sn_sum/ (double) (sn_sum + clear_sum + 0.01);


    for(i = 0; i < NUM_FC; i++)
    {
        rec_cg[i].pos = num_samples * (row_pos - 1) + col_pos;

        //the below is to initialize
        rec_cg[i].category = -9999;
        //rec_cg[i].land_type = 0;

        for (j = 0; j < TOTAL_IMAGE_BANDS+TOTAL_INDICES; j++){
//            rec_cg[i].obs_disturb[j] = NA_VALUE;
            rec_cg[i].rmse[j] = NA_VALUE;
            rec_cg[i].magnitude[j] = NA_VALUE;
            for(k = 0; k < SCCD_MAX_NUM_C - 1; k++){
//                rec_cg[i].state_disturb[j][k] = NA_VALUE;
                rec_cg[i].coefs[j][k] = NA_VALUE;
            }
            rec_cg[i].coefs[j][SCCD_MAX_NUM_C - 1] = NA_VALUE;
        }
    }

    if (clear_sum < N_TIMES * MAX_NUM_C){
        result = sccd_inefficientobs_procedure(valid_num_scenes,valid_date_array, buf,
                                 fmask_buf, id_range,sn_pct,rec_cg, num_fc);
    }
    else{

        /**************************************************************/
        /*                                                            */
        /* standard_procedure for CCD                                 */
        /*                                                            */
        /**************************************************************/

//       result = sccd_stand_procedure(valid_num_scenes, valid_date_array, buf,
//                                     fmask_buf, id_range, rec_cg, num_fc,
//                                     states_output_dir, b_fastmode, probability_threshold,
//                                     min_days_conse, training_type, monitorwindow_lowerlin,
//                                     monitorwindow_upperlim, sensor_buf, n_focus_variable,
//                                     n_total_variable, focus_blist, NDVI_INCLUDED,
//                                     NBR_INCLUDED, RGI_INCLUDED,TCTWETNESS_INCLUDED,TCTGREENNESS_INCLUDED,
//                                     EVI_INCLUDED, DI_INCLUDED, NDMI_INCLUDED, booster, b_landspecific, auxval);
       result = sccd_stand_procedure(valid_num_scenes, valid_date_array, buf,
                                     fmask_buf, id_range, rec_cg, num_fc,
                                     states_output_dir, b_fastmode, probability_threshold,
                                     min_days_conse, training_type, monitorwindow_lowerlin,
                                     monitorwindow_upperlim, sensor_buf, n_focus_variable,
                                     n_total_variable, focus_blist, NDVI_INCLUDED,
                                     NBR_INCLUDED, RGI_INCLUDED,TCTWETNESS_INCLUDED,TCTGREENNESS_INCLUDED,
                                     EVI_INCLUDED, DI_INCLUDED, NDMI_INCLUDED, b_landspecific, auxval, conse);

    }

    free(id_range);
    id_range = NULL;
    if (result == SUCCESS)
    {
        return (SUCCESS);
    }
    else
    {
        return (FAILURE);
    }
}




// Solve for f so that Pr[F < f] = pr where F has an F distribution with
// v1, v2 degrees of freedom.

double F(int v1, int v2, double pr) {
   double f0 = 0.0;
   double f = 2.0;
   double diff;
   double delta = 0.1;
   double p;
   int n = 0;

   f = 1.0;
   while ( (p = F_Distribution(f,v1,v2)) < pr) {
      f += delta;
      if (p > 0.999) {
         delta /= 10.0;
         f = (f - delta) / 10;
      }
   }

   f0 = f - delta;
   while ( (p = F_Distribution(f0,v1,v2)) > pr) {
      f0 -= delta;
      if (p < 0.001) {
         delta /= 10.0;
         f0 = (f0 + delta) / 10;
      }
   }

   while (fabs(f - f0) > 1.e-6) {
      diff = F_Distribution(f, v1, v2) - pr;
      diff /= F_Density(f,v1,v2);
      diff /= 2.0;
      f0 = f;
      f -= diff;
      n++;
      if (n > 40) exit(0);
   }
   return f;
}

// standard observation day as doy between i_start and i_end
int convert_standardordinal(
    int i_start,
    int i_end
)
{
    int start_year = (int)((i_start - ORDINALDAY_19710101) / AVE_DAYS_IN_A_YEAR) + 1971;
    int end_year = (int)((i_end - ORDINALDAY_19710101) / AVE_DAYS_IN_A_YEAR) + 1971;
    int ordinaldays_standard;
    convert_year_doy_to_ordinal(start_year, JULY1ST_DOY, &ordinaldays_standard);
    if(ordinaldays_standard < i_start) // meaning i_start is at behind july 1st of that year
        convert_year_doy_to_ordinal(start_year + 1, JULY1ST_DOY, &ordinaldays_standard);
    return ordinaldays_standard;
}

int step1_ssm_initialize(
    ssmodel *instance,          /* I/O: the outputted initial SS model                */
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
    double rmse,                     /*I: the rmse from lasso process */
    double adj_rmse,                 /*I: the rmse from temporal variogram */
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
)
{
    int i, j, k;
    int status;
    double mean_y = 0;
    double interval;
    double best_q;
    double best_H;
    double tmp_q;
    double tmp_lik;
    double tmp_f_rmse = 0;
    double tmp_v_rmse = 0;
    gsl_vector* kt;
    double variogram = 0;
    char FUNC_NAME[] = "step1_ssm_initialize";
    double var_smooth_Q;
    double var_smooth_H;
    bool b_start;
    double* state_sum;
    double max_rmse;
    double ini_p;
    double ini_a;
    int z_sum;
    double unadjusted_rmse;

//    printf("fit_cft[0][0]: %f\n", fit_cft[0]);
//    printf("fit_cft[0][1]: %f\n", fit_cft[1]);
//    printf("fit_cft[0][2]: %f\n", fit_cft[2]);
//    printf("fit_cft[0][3]: %f\n", fit_cft[3]);
//    printf("fit_cft[0][4]: %f\n", fit_cft[4]);
//    printf("fit_cft[0][5]: %f\n", fit_cft[5]);
//    fit_cft2vec_a(fit_cft, next_a, clrx[stable_start]);
//    vec_a2fit_cft(next_a, fit_cft, clrx[stable_start]);
//    printf("fit_cft[0][0]: %f\n", fit_cft[0]);
//    printf("fit_cft[0][1]: %f\n", fit_cft[1]);
//    printf("fit_cft[0][2]: %f\n", fit_cft[2]);
//    printf("fit_cft[0][3]: %f\n", fit_cft[3]);
//    printf("fit_cft[0][4]: %f\n", fit_cft[4]);
//    printf("fit_cft[0][5]: %f\n", fit_cft[5]);

    /* assign memory and initialize them to zero */
    state_sum = (double *) calloc(instance->m, sizeof(double));
    if (state_sum == NULL)
    {
        RETURN_ERROR ("Allocating state_sum memory",FUNC_NAME, FAILURE);
    }

    kt = gsl_vector_alloc(instance->m);

    /********************************************************************************************/
    /* first step: update a and p between prev_i_break - 1(included) and stable_start(not included) */
    /********************************************************************************************/
//    if(b_fastmode == FALSE)
//    {
//        for(i = prev_break_date; i < clrx[stable_start]; i++)
//        {
//                level_state_records[i_b][i - starting_date] = NA_VALUE;
//                annual_state_records[i_b][i - starting_date] = NA_VALUE;
//                semi_state_records[i_b][i - starting_date] = NA_VALUE;
//                third_state_records[i_b][i - starting_date] = NA_VALUE;
//        }
//    }



    /*
           level       cycle     cycle*       cycle     cycle*      cycle     cycle*
    level      1    0.00000000 0.00000000  0.00000000 0.00000000  0.00000000  0.00000000
    cycle      0    0.99985204 0.01720158  0.00000000 0.00000000  0.00000000  0.00000000
    cycle*     0   -0.01720158 0.99985204  0.00000000 0.00000000  0.00000000  0.00000000
    cycle      0    0.00000000 0.00000000  0.99940821 0.03439806  0.00000000  0.00000000
    cycle*     0    0.00000000 0.00000000 -0.03439806 0.99940821  0.00000000  0.00000000
    cycle      0    0.00000000 0.00000000  0.00000000 0.00000000  0.99866864  0.05158437
    cycle*     0    0.00000000 0.00000000  0.00000000 0.00000000 -0.05158437  0.99866864

    */
    for(i = 0; i < instance->m; i++)
    {
        for (j = 0; j < instance->m; j++)
        {
            if((i == 0)&&(j == 0))
            {
              gsl_matrix_set(instance->T, i, j, 1.0);
              continue;
            }

            if (instance ->structure % 10 == 1)
            {
                if((i == 1)&&(j == 1))
                {
                  gsl_matrix_set(instance->T, i, j, cos((double)TWO_PI / (double)NUM_YEARS));
                  continue;
                }

                if((i == 1)&&(j == 2))
                {
                  gsl_matrix_set(instance->T, i, j, sin((double)TWO_PI / (double)NUM_YEARS));
                  continue;
                }

                if((i == 2)&&(j == 2))
                {
                  gsl_matrix_set(instance->T, i, j, cos((double)TWO_PI / (double)NUM_YEARS));
                  continue;
                }

                if((i == 2)&&(j == 1))
                {
                  gsl_matrix_set(instance->T, i, j, -sin((double)TWO_PI / (double)NUM_YEARS));
                  continue;
                }
            }

            if (instance ->structure / 10 == 1)
            {
                if((i == 3)&&(j == 3))
                {
                  gsl_matrix_set(instance->T, i, j, cos((double)TWO_PI / (double)NUM_YEARS * 2.0));
                  continue;
                }

                if((i == 4)&&(j == 4))
                {
                  gsl_matrix_set(instance->T, i, j, cos((double)TWO_PI / (double)NUM_YEARS * 2.0));
                  continue;
                }

                if((i == 3)&&(j == 4))
                {
                  gsl_matrix_set(instance->T, i, j, sin((double)TWO_PI / (double)NUM_YEARS * 2.0));
                  continue;
                }

                if((i == 4)&&(j == 3))
                {
                  gsl_matrix_set(instance->T, i, j, -sin((double)TWO_PI / (double)NUM_YEARS * 2.0));
                  continue;
                }

            }


            if(instance ->structure / 100 == 1)
            {
                if((i == 5)&&(j == 5))
                {
                  gsl_matrix_set(instance->T, i, j, cos((double)TWO_PI / (double)NUM_YEARS * 3.0));
                  continue;
                }

                if((i == 6)&&(j == 6))
                {
                  gsl_matrix_set(instance->T, i, j, cos((double)TWO_PI / (double)NUM_YEARS * 3.0));
                  continue;
                }

                if((i == 5)&&(j == 6))
                {
                  gsl_matrix_set(instance->T, i, j, sin((double)TWO_PI / (double)NUM_YEARS * 3.0));
                  continue;
                }

                if((i == 6)&&(j == 5))
                {
                  gsl_matrix_set(instance->T, i, j, -sin((double)TWO_PI / (double)NUM_YEARS * 3.0));
                  continue;
                }

            }

        }
    }


//    printf("%f, %f, %f\n", fit_cft[0], fit_cft[2], fit_cft[3]);
    /* initialize a1 and p1*/
    // printf("success 1 and m = %d \n", instance->m);
    fit_cft2vec_a(fit_cft, next_a, clrx[stable_start], instance->m, instance->structure);
//    vec_a2fit_cft(next_a, fit_cft, clrx[stable_start]);

//    printf("%f, %f, %f\n", fit_cft[0], fit_cft[2], fit_cft[3]);
    if(instance->m == 1)
    {
        /*   initialize Z     */
        gsl_vector_set(instance->Z, 0, 1.0);

    }
    else if(instance->m == 3)
    {
        /*   initialize Z     */
        gsl_vector_set(instance->Z, 0, 1.0);
        gsl_vector_set(instance->Z, 1, 1.0);
        gsl_vector_set(instance->Z, 2, 0.0);

    }
    else if(instance->m == 5)
    {
        /*   initialize Z     */
        gsl_vector_set(instance->Z, 0, 1.0);
        gsl_vector_set(instance->Z, 1, 1.0);
        gsl_vector_set(instance->Z, 2, 0.0);
        gsl_vector_set(instance->Z, 3, 1.0);
        gsl_vector_set(instance->Z, 4, 0.0);

    }
    else
    {
        /*   initialize Z     */
        gsl_vector_set(instance->Z, 0, 1.0);
        gsl_vector_set(instance->Z, 1, 1.0);
        gsl_vector_set(instance->Z, 2, 0.0);
        gsl_vector_set(instance->Z, 3, 1.0);
        gsl_vector_set(instance->Z, 4, 0.0);
        gsl_vector_set(instance->Z, 5, 1.0);
        gsl_vector_set(instance->Z, 6, 0.0);
    }

//    if((num_curve == 0)||(num_curve == 1 && (rec_cg[0].category / 10) == 1))
//    {
    //double t1, t2;

    gsl_vector_memcpy(instance->a1, next_a);
    gsl_matrix_memcpy(instance->P1, next_P);

   for(k = stable_start; k < stable_end + 1; k++)
   {
       instance->clear_date_array[k - stable_start] = clrx[k];
       instance->yt[k - stable_start] = clry[k];
   }
   instance->n = clrx[stable_end] - clrx[stable_start] + 1;            /* number of time points */

   // array_1d_mean(clry, stable_start, stable_end, &mean_y);

   //singleband_minvariogram(clrx, clry, stable_start, stable_end, &variogram);

   interval = (double)(clrx[stable_end] - clrx[stable_start]) / (stable_end - stable_start);
   unadjusted_rmse = rmse * rmse * (double)(stable_end - stable_start + 1 - SCCD_MAX_NUM_C) / (double)(stable_end - stable_start + 1);
   /* grid searching method */
   //grid_searching(instance, (double)rmse, (double)variogram, (double)mean_y, &best_H, &best_q, interval, instance->m);
   //Jazwinski_searching(instance, (double)rmse, (double)mean_y, &best_H, &best_q, interval,instance->m);

   /* initialize p based on a intensity*/
   if(instance->m == 1)
       ini_a = gsl_vector_get(next_a, 0);
   else if(instance->m == 3)
       ini_a = gsl_vector_get(next_a, 0) + gsl_vector_get(next_a, 1);
   else if(instance->m == 5)
       ini_a = gsl_vector_get(next_a, 0) + gsl_vector_get(next_a, 1) + gsl_vector_get(next_a, 3);
   else
       ini_a = gsl_vector_get(next_a, 0) + gsl_vector_get(next_a, 1) + gsl_vector_get(next_a, 3) + gsl_vector_get(next_a, 5);
   //ini_a = gsl_vector_get(next_a, 0);

   // calculate z_sum
   z_sum = 0;
   for (k = 0; k < instance -> m; k++)
       z_sum = z_sum + gsl_vector_get(instance->Z, k);

   ini_p = pow((double)ini_a * INITIAL_P_RATIO, 2) / z_sum;

   for(k = 0; k < instance->m; k++)
   {
       //t1 = gsl_vector_get(next_a, k) * gsl_vector_get(next_a, k)  * INITIAL_P_RATIO;
       //t2 = rmse * rmse / (instance ->m + 1) * 2;
       // gsl_matrix_set(next_P, k, k, gsl_vector_get(next_a, k) * gsl_vector_get(next_a, k)  * 0.005);
       //gsl_matrix_set(next_P, k, k, gsl_vector_get(next_a, k) * gsl_vector_get(next_a, k)  * INITIAL_P_RATIO);
       //gsl_matrix_set(next_P, k, k, rmse * rmse / (instance ->m + 1) / 2 * 0.1);
       //gsl_matrix_set(next_P, k, k, pow(clry[stable_start]/ (instance ->m + 1) * 2, 2) * INITIAL_P_RATIO);
       gsl_matrix_set(next_P, k, k, ini_p);
    }
//   best_H = rmse * rmse * 0.9;
   //rmse = rmse * sqrt((double)(stable_end - stable_start + 1 - SCCD_MAX_NUM_C) / (double)(stable_end - stable_start + 1));

//   if(ini_p < rmse * rmse)
//       best_H = rmse * rmse - ini_p;
//   else
//       best_H = INI_P00;
   //best_H = rmse * rmse;
   best_H = unadjusted_rmse;

   best_q = INI_Q00;
   // grid_searching(instance, (double)rmse, (double)variogram, (double)mean_y, &best_H, &best_q, interval, instance->m, ini_a, &ini_p);


// Jazwinski method
//   tmp_lik = ssmloglik_gridsearching(instance, best_H, best_q, &tmp_f_rmse, &tmp_v_rmse);
//   if( tmp_v_rmse > tmp_f_rmse)
//       best_q = (tmp_v_rmse - tmp_f_rmse) / ((instance->m + 1)/2  * interval);

   //printf("best_q is %f\n", best_q);

   for (i = 0; i < instance->m; i++)
       gsl_matrix_set(instance->Q, i, i, best_q);

   instance->H = best_H;


//        printf("instance->Q[1] = %f\n", gsl_matrix_get(instance->Q, 1, 1));
//        printf("instance->Q[2] = %f\n", gsl_matrix_get(instance->Q, 2, 2));
//        printf("instance->Q[3] = %f\n", gsl_matrix_get(instance->Q, 3, 3));
//        printf("instance->Q[4] = %f\n", gsl_matrix_get(instance->Q, 4, 4));
//        printf("instance->Q[5] = %f\n", gsl_matrix_get(instance->Q, 5, 5));
//        printf("instance->H = %f\n", instance->H);

    if((status = ! LBFGS_STOP) && (status = ! LBFGS_SUCCESS) )
    {
        RETURN_ERROR ("Fitting state-space model fails",FUNC_NAME, FAILURE);
    }

    /******************************************************/
    /*     fourth step: update states from stable_start   */
    /******************************************************/

    for(i = stable_start; i < stable_end  + 1; i++)
    {
        if (i < stable_start + SCCD_MAX_NUM_C)
            b_start = TRUE;
        else
        {
            for (j = 0; j < instance->m; j++)
                state_sum[j] = state_sum[j] + fabsf(gsl_vector_get(next_a, j));
                //state_sum[j] = state_sum[j] + gsl_vector_get(next_a, j) * gsl_vector_get(next_a, j);
            b_start = FALSE;
        }


        KF_ts_filter_regular(instance, clrx, clry, next_P, next_a,  i,
                             i_b, rmse_records, level_state_records, annual_state_records,
                              semi_state_records, third_state_records,starting_date, b_fastmode,
                             fit_cft, valid_count, sum_square_smooth_Q,
                             sum_square_smooth_H, sum_smooth_Q, sum_smooth_H,
                             sum_square_vt, sum_vt, sum_kalman_coef, b_start, temporal_rmse_square_sum, temporal_rmse_sum,
                             temporal_rmse_count);

    }
    //*sum_square_vt = rmse * rmse * (*valid_count);

    var_smooth_Q = sum_square_smooth_Q[i_b] / (valid_count[i_b] - 1);
    for (i = 0; i < instance->m; i++)
    {
        //printf("%f\n", state_sum[i]);
        gsl_matrix_set(instance->Q, i, i, var_smooth_Q * state_sum[i] / state_sum[0]);
    }

    //var_smooth_H = sum_square_smooth_H[i_b] / (valid_count[i_b] - 1);
    //instance->H = var_smooth_H;



//    for(i = 1; i < 6; i++)
//        printf("instance->Q: %f\n", gsl_matrix_get(instance->Q, i, i));

//    printf("instance->H: %f\n", instance->H);



    /* free memory*/
    gsl_vector_free(kt);
    free(state_sum);
    return SUCCESS;

}


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
    double *rmse,                      /* O: Root Mean Squared Error array.        */
    int n_focus_variable,
    int n_total_variable,
    int* focus_blist,
    int min_days_conse
)
{
    int status;
    int k, m, b;
    int i_b;
    int *bl_ids, *ids;
    int *rm_ids;
    int i_span;
    int rm_ids_len;
    int *cpx;
    float **cpy;
    int k_new;
    int i_rec;
    double v_dif_norm;

    double time_span;
    double mini_rmse;                 /* Mimimum RMSE                          */
    double* v_start;                  /* Vector for start of observation(s)    */
    double* v_end;                    /* Vector for end of observastion(s)     */
    double* v_slope;                 /* Vector for anormalized slope values   */
    double* v_dif;                    /* Vector for difference values          */
    double **v_diff;
    double *vec_magg;           /* this one is never freed */
    double vec_magg_min;
    double **v_dif_mag;               /* vector for magnitude of differences.  */
    double **v_dif_magg;
    int i_conse;
    int ini_conse;
    int i_ini;
    double ts_pred_temp;
    double **tmp_v_dif;              /* temporal regression difference       */
    double **rec_v_dif;
    double v_dif_mean;
    double **fit_cft_tmp;
    double mean_angle;           /* I: mean angle of vec_diff                              */
    char FUNC_NAME[] = "step1_ccd_initialize";
    int update_num_c;
    int real_end;

    // i_dense = *i_start;

    /****************************************************/
    /*                                                  */
    /*     allocac memory for variables                 */
    /*                                                  */
    /****************************************************/
    rec_v_dif = (double **)allocate_2d_array(n_total_variable, n_clr,
                                     sizeof (double));
    if (rec_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating rec_v_dif memory",FUNC_NAME, FAILURE);
    }

    ids = (int *)malloc(n_clr * sizeof(int));
    if (ids == NULL)
    {
        RETURN_ERROR("ERROR allocating ids memory", FUNC_NAME, FAILURE);
    }

    bl_ids = (int *)malloc(n_clr * sizeof(int));
    if (bl_ids == NULL)
    {
        RETURN_ERROR("ERROR allocating bl_ids memory", FUNC_NAME, FAILURE);
    }

    rm_ids = (int *)malloc(n_clr * sizeof(int));
    if (rm_ids == NULL)
    {
        RETURN_ERROR("ERROR allocating rm_ids memory", FUNC_NAME, FAILURE);
    }

    v_start = (double *)malloc(n_focus_variable * sizeof(double));
    if (v_start == NULL)
    {
        RETURN_ERROR("ERROR allocating v_start memory", FUNC_NAME, FAILURE);
    }

    v_end = (double *)malloc(n_focus_variable * sizeof(double));
    if (v_end == NULL)
    {
        RETURN_ERROR("ERROR allocating v_end memory", FUNC_NAME, FAILURE);
    }

    v_slope = (double *)malloc(n_focus_variable * sizeof(double));
    if (v_slope == NULL)
    {
        RETURN_ERROR("ERROR allocating v_slope memory", FUNC_NAME, FAILURE);
    }

    v_dif = (double *)malloc(n_focus_variable * sizeof(double));
    if (v_dif == NULL)
    {
        RETURN_ERROR("ERROR allocating v_dif memory", FUNC_NAME, FAILURE);
    }

    tmp_v_dif = (double **)allocate_2d_array(n_total_variable, n_clr,
                                     sizeof (double));
    if (tmp_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating tmp_v_dif memory",FUNC_NAME, FAILURE);
    }

    v_dif_mag = (double **) allocate_2d_array(n_total_variable, conse,
                sizeof (double));
    if (v_dif_mag == NULL)
    {
        RETURN_ERROR ("Allocating v_dif_mag memory",
                                 FUNC_NAME, FAILURE);
    }

//    vec_mag = (double *)calloc(conse, sizeof(double));
//    if (vec_mag == NULL)
//    {
//        RETURN_ERROR ("Allocating vec_mag memory", FUNC_NAME, FAILURE);
//    }
    /*************************************************  */
    /*                                                  */
    /*  check maximum time gap as first                 */
    /*                                                  */
    /****************************************************/

    int max_date_diff = 0;
    for (k = *i_start; k < *cur_i; k++)
    {
        if (clrx[k + 1] - clrx[k] > max_date_diff)
        {
            max_date_diff = clrx[k + 1] - clrx[k];
        }
    }

    if (max_date_diff > NUM_YEARS)      //SY 09192018
    {
        // *cur_i = *cur_i + 1;
        *i_start = *i_start + 1;
        *i_start_copy = *i_start_copy + 1;
        *i_dense = *i_dense + 1;


        status = free_2d_array ((void **) rec_v_dif);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: rec_v_dif\n",
                          FUNC_NAME, FAILURE);
        }

        free(ids);
        free(bl_ids);
        free(rm_ids);
    //    free(rmse);
    //    rmse = NULL;
        free(v_start);
        free(v_end);
        free(v_slope);
        free(v_dif);

        status = free_2d_array ((void **) tmp_v_dif);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: tmp_v_dif\n",
                          FUNC_NAME, FAILURE);
        }

        status = free_2d_array ((void **) v_dif_mag);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: v_dif_mag\n",
                       FUNC_NAME, FAILURE);
        }

        return INCOMPLETE;
    }
    /*************************************************  */
    /*                                                  */
    /* Step 2: noise removal.                           */
    /*                                                  */
    /****************************************************/
    real_end = min(*cur_i+conse, *end - 1);
    status = auto_mask(clrx, clry, *i_start, real_end,
                   (double)(clrx[real_end]-clrx[*i_start]) / NUM_YEARS,
                   adj_rmse[1], adj_rmse[4], SCCD_T_CONST, bl_ids);
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

    for (k = *i_start; k < *cur_i + 1; k++)
    {
        ids[k- *i_start] = k;
    }

    m= 0;

    i_span = 0;

    for (k = 0; k < *cur_i - *i_start + 1; k++)
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

    /**************************************************/
    /*                                                */
    /* Check if there are enough observation.         */
    /*                                                */
    /**************************************************/

    if (i_span < i_span_min)
    {
        /**********************************************/
        /*                                            */
        /* Move forward to the i+1th clear observation*/
        /*                                            */
        /**********************************************/

        // *cur_i= *cur_i + 1;

        /**********************************************/
        /*                                            */
        /* Not enough clear observations.             */
        /*                                            */
        /**********************************************/
        status = free_2d_array ((void **) rec_v_dif);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: rec_v_dif\n",
                          FUNC_NAME, FAILURE);
        }

        free(ids);
        free(bl_ids);
        free(rm_ids);
        free(v_start);
        free(v_end);
        free(v_slope);
        free(v_dif);

        status = free_2d_array ((void **) tmp_v_dif);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: tmp_v_dif\n",
                          FUNC_NAME, FAILURE);
        }

        status = free_2d_array ((void **) v_dif_mag);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: v_dif_mag\n",
                       FUNC_NAME, FAILURE);
        }

        return INCOMPLETE;
    }

    if (*end == 0)
        RETURN_ERROR("No available data point", FUNC_NAME, FAILURE);

    /**************************************************/
    /*                                                */
    /* Allocate memory for cpx, cpy.                  */
    /*                                                */
    /**************************************************/

    cpx = malloc((*end) * sizeof(int));
    if (cpx == NULL)
        RETURN_ERROR("ERROR allocating cpx memory", FUNC_NAME, FAILURE);

    cpy = (float **) allocate_2d_array (n_total_variable, *end,
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
    for (k = 0, k_new=0; k < *end; k++)
    {
        if (m < rm_ids_len && k == rm_ids[m])
        {
            m++;
            continue;
        }
        cpx[k_new] = clrx[k];
        for (b = 0; b < n_total_variable; b++)
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

    i_rec = *cur_i;

    /**************************************************/
    /*                                                */
    /* Update i afer noise removal.                   */
    /*     (i_start stays the same).                  */
    /*                                                */
    /**************************************************/

    *cur_i = *i_start + i_span - 1;

    /**************************************************/
    /*                                                */
    /* Update span of time (num of years).            */
    /*                                                */
    /**************************************************/

    time_span=(cpx[*cur_i] - cpx[*i_start]) / NUM_YEARS;

    /**************************************************/
    /*                                                */
    /* Check if there is enough time.                 */
    /*                                                */
    /**************************************************/

    if (time_span < MIN_YEARS)
    {
        *cur_i = i_rec;   /* keep the original i */

        /**********************************************/
        /*                                            */
        /* Move forward to the i+1th clear observation*/
        /*                                            */
        /**********************************************/

        // *cur_i= *cur_i + 1;

        status = free_2d_array ((void **) rec_v_dif);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: rec_v_dif\n",
                          FUNC_NAME, FAILURE);
        }

        free(ids);
        free(bl_ids);
        free(rm_ids);
    //    free(rmse);
    //    rmse = NULL;
        free(v_start);
        free(v_end);
        free(v_slope);
        free(v_dif);

        status = free_2d_array ((void **) tmp_v_dif);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: tmp_v_dif\n",
                          FUNC_NAME, FAILURE);
        }

        status = free_2d_array ((void **) v_dif_mag);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: v_dif_mag\n",
                       FUNC_NAME, FAILURE);
        }

        free(cpx);
        status = free_2d_array ((void **) cpy);
        if (status != SUCCESS)
        {
              RETURN_ERROR ("Freeing memory: cpy\n",
                    FUNC_NAME, FAILURE);
        }

        return INCOMPLETE;    /* not enough time */
    }

    //SY 09272018
    /**************************************************/
    /*                                                */
    /* updated end after checking if enought time_span*/
    /*                                                */
    /**************************************************/
    *end = k_new;

    /**************************************************/
    /*                                                */
    /* Remove noise in original arrays.               */
    /*                                                */
    /**************************************************/

    for (k = 0; k < *end; k++)
    {
        clrx[k] = cpx[k];
        for (m = 0; m < n_total_variable; m++)
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
    /* Step 3) model fitting: initialize model testing*/
    /*         variables defining computed variables. */
    /*                                                */
    /**************************************************/
    if (*cur_i - *i_start + 1 > N_TIMES * SCCD_MAX_NUM_C - 1)
       update_num_c = SCCD_MAX_NUM_C;
    else if (*cur_i - *i_start + 1 > N_TIMES * MID_NUM_C - 1)
       update_num_c = MID_NUM_C;
    else
       update_num_c = MIN_NUM_C;

    // update_num_c = MIN_NUM_C;

    for (b = 0; b < n_total_variable; b++)
    {
        /**********************************************/
        /*                                            */
        /* Initial model fit.                         */
        /*                                            */
        /**********************************************/

        status = auto_ts_fit(clrx, clry,  b, b, *i_start, *cur_i,
                 update_num_c, fit_cft, &rmse[b], rec_v_dif);
//                        for (k = 0; k < SCCD_MAX_NUM_C; k++)
//                        {

//                               printf("%f\n", (double)fit_cft[b][k]);

//                        }


        if (status != SUCCESS)
        {
            RETURN_ERROR ("Calling auto_ts_fit during model initilization\n",
                 FUNC_NAME, FAILURE);
        }
    }

    /*** testing codes ***/
//    FILE* fptr_save;
//    char states_output_dir_full[MAX_STR_LEN];
//    char states_output_dir[MAX_STR_LEN];
//    sprintf(states_output_dir, "/home/su/Documents/Jupyter/source/LandsatARD/Stable_recf_B");
//    for(i_b = 0; i_b < NUM_BANDS; i_b++)
//    {
//        sprintf(states_output_dir_full, "%s%d.csv", states_output_dir, i_b + 1);
//        fptr_save= fopen(states_output_dir_full, "w");

//        fprintf(fptr_save, "%s, %s, %s\n", "clrx", "clry", "rec_v_dif");

//        for (int i = *i_start; i < *cur_i + 1; i++)
//        {

//            fprintf(fptr_save, "%i, %f, %f\n",
//                    clrx[i], clry[i_b][i],rec_v_dif[i_b][i]);

//        }
//        fclose(fptr_save);
//    }

    v_dif_norm = 0.0;

    for (i_b = 0; i_b < n_focus_variable; i_b++)
    {

        /**********************************************/
        /*                                            */
        /* Calculate min. rmse.                       */
        /*                                            */
        /**********************************************/

        mini_rmse = max(adj_rmse[focus_blist[i_b]], rmse[focus_blist[i_b]]);
        //mini_rmse = rmse[focus_blist[i_b]];
        /**********************************************/
        /*                                            */
        /* Compare the first observation.             */
        /*                                            */
        /**********************************************/

        v_start[i_b] = (rec_v_dif[focus_blist[i_b]][0]) / mini_rmse ;

        /**********************************************/
        /*                                            */
        /* Compare the last clear observation.        */
        /*                                            */
        /**********************************************/

        v_end[i_b] = rec_v_dif[focus_blist[i_b]][*cur_i-*i_start]/ mini_rmse;

        /**********************************************/
        /*                                            */
        /* Anormalized slope values.                  */
        /*                                            */
        /**********************************************/
        v_slope[i_b] = fit_cft[focus_blist[i_b]][1] *
                        (clrx[*cur_i]-clrx[*i_start])/mini_rmse;

        /**********************************************/
        /*                                            */
        /* Difference in model intialization.         */
        /*                                            */
        /**********************************************/

        v_dif[i_b] = fabs(v_slope[i_b]) + max(fabs(v_start[i_b]), fabs(v_end[i_b]));
        //v_dif[i_b] = fabs(v_slope[i_b]) + 1;
        v_dif_norm += v_dif[i_b] * v_dif[i_b];

    }



    /**************************************************/
    /*                                                */
    /* Find stable start for each curve.              */
    /*                                                */
    /**************************************************/
    if (v_dif_norm > reg_TCG)
    {
        /**********************************************/
        /*                                            */
        /* Start from next clear observation.         */
        /*                                            */
        /**********************************************/

        *i_start = *i_start + 1;

        *i_start_copy = *i_start_copy + 1;
        /**********************************************/
        /*                                            */
        /* Move forward to the i+1th clear observation*/
        /*                                            */
        /**********************************************/

        //*cur_i = *cur_i + 1;

        /**********************************************/
        /*                                            */
        /* Keep all data and move to the next obs.    */
        /*                                            */
        /**********************************************/

        status = free_2d_array ((void **) rec_v_dif);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: rec_v_dif\n",
                          FUNC_NAME, FAILURE);
        }

        free(ids);
        free(bl_ids);
        free(rm_ids);
    //    free(rmse);
    //    rmse = NULL;
        free(v_start);
        free(v_end);
        free(v_slope);
        free(v_dif);

        status = free_2d_array ((void **) tmp_v_dif);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: tmp_v_dif\n",
                          FUNC_NAME, FAILURE);
        }

        status = free_2d_array ((void **) v_dif_mag);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: v_dif_mag\n",
                       FUNC_NAME, FAILURE);
        }

        return INCOMPLETE;
    }


    /**************************************************/
    /*                                                */
    /* Find the previous break point.                 */
    /*                                                */
    /**************************************************/

//    if (*num_curve == 0)
//    {
//        i_break = 0; /* first curve */
//    }
//    else
//    {
//        /**********************************************/
//        /*                                            */
//        /* After the first curve, compare rmse to     */
//        /* determine which curve to determine t_break.*/
//        /*                                            */
//        /**********************************************/

//        for (k = 0; k < *end; k++)
//        {
//            if (clrx[k] >= rec_cg[*num_curve-1].t_break)
//            {
//                i_break = k;
//                break;
//            }
//        }
//    }

    /**************************************************/
    /*                                                */
    /* Step 4) look back: to include fitting points   */
    /*         and find change points.                */
    /*                                                */
    /**************************************************/
    if(*num_curve == 0)
    {
        *prev_i_break = *i_dense;
    }

    //if ( *num_curve == 0 && *i_start > *prev_i_break)
    //if (*i_start < 0) //a simple way to omit looking back procedure
    if (*i_start > *prev_i_break)
    {
        /**********************************************/
        /*                                            */
        /* Model fit at the beginning of the time     */
        /* series.                                    */
        /*                                            */
        /**********************************************/

        for(i_ini = *i_start - 1; i_ini >= *prev_i_break; i_ini--) // SY 09192018
        {
            // printf("ini is %d\n", i_ini);
            if (i_ini - *prev_i_break + 1 < conse)
            {
                ini_conse = i_ini - *prev_i_break + 1;
            }
            else
            {
                ini_conse = conse;
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

            v_diff = (double **) allocate_2d_array(n_focus_variable,
                        ini_conse, sizeof (double));
            if (v_diff == NULL)
            {
                RETURN_ERROR ("Allocating v_diff memory",
                               FUNC_NAME, FAILURE);
            }

            vec_magg = (double *) malloc(ini_conse * sizeof (double));
            if (vec_magg == NULL)
            {
                RETURN_ERROR ("Allocating vec_magg memory",
                              FUNC_NAME, FAILURE);
            }

            /******************************************/
            /*                                        */
            /* Detect change.                         */
            /* value of difference for conse          */
            /* Record the magnitude of change.        */
            /*                                        */
            /******************************************/

            vec_magg_min = 9999.0;
            for (i_conse = 1; i_conse < ini_conse + 1; i_conse++) // SY 09192018
            {
                v_dif_norm = 0.0;
                for (i_b = 0; i_b < n_total_variable; i_b++)
                {

                    /**********************************/
                    /*                                */
                    /* Absolute differences.          */
                    /*                                */
                    /**********************************/

                    // SY 09192018 moving fitting into (i_b == focus_blist[b])to save time //
                    // SY 02/13/2019 delete these speed-up modification as non-lasso bands
                    // are important for change agent classification
                    // printf("i_b = %d\n", i_b);
                    auto_ts_predict(clrx, fit_cft, update_num_c, i_b, i_ini-i_conse+1,
                                    i_ini-i_conse+1, &ts_pred_temp);
                    v_dif_mag[i_b][i_conse-1] = (double)clry[i_b][i_ini-i_conse+1] -
                                       ts_pred_temp;// SY 09192018

                    /**********************************/
                   /*                                */
                   /* Normalize to z-score.          */
                   /*                                */
                   /**********************************/

                   for (b = 0; b < n_focus_variable; b++)
                   {
                       if (i_b == focus_blist[b])
                       {
                            /**************************/
                            /*                        */
                            /* Minimum rmse.          */
                            /*                        */
                            /**************************/

                            mini_rmse = max(adj_rmse[i_b], rmse[i_b]);

                            /**************************/
                            /*                        */
                            /* z-scores.              */
                            /*                        */
                            /**************************/

                            v_diff[b][i_conse-1] = v_dif_mag[i_b][i_conse-1] // SY 09192018
                                                          / mini_rmse;
                            v_dif_norm += v_diff[b][i_conse-1] * v_diff[b][i_conse-1]; // SY 09192018
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

            mean_angle = MeanAngl(v_diff, n_focus_variable, ini_conse);

            if ((vec_magg_min > reg_TCG) && (mean_angle < NSIGN)) /* i_start found*/
            {
                free(vec_magg);
                status = free_2d_array ((void **) v_diff);
                if (status != SUCCESS)
                {
                    RETURN_ERROR ("Freeing memory: v_diff\n",
                                  FUNC_NAME, FAILURE);
                }
                break;
            }
            else if ((vec_magg[0] > T_MAX_CG) && (i_ini > *prev_i_break)) /* false change */
            {
                for (k = i_ini; k < *end - 1; k++)
                {
                    clrx[k] = clrx[k+1];
                    for (b = 0; b < n_total_variable; b++)
                    {
                        clry[b][k] = clry[b][k+1];
                    }
                }
                *cur_i = *cur_i - 1;
                *end = *end - 1;
            }

            /**************************************/
            /*                                    */
            /* Update i_start if i_ini is not a   */
            /* confirmed break.                   */
            /*                                    */
            /**************************************/
            *i_start_copy = *i_start_copy - (*i_start - i_ini -1);

            *i_start = i_ini + 1;

            /******************************************/
            /*                                        */
            /* Free the temporary memory.             */
            /*                                        */
            /******************************************/

            free(vec_magg);
            status = free_2d_array ((void **) v_diff);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: v_diff\n",
                              FUNC_NAME, FAILURE);
            }
        } // end for (i_ini = i_start-1; i_ini >= prev_i_break; i_ini--)
    }    // end for if (i_start > prev_i_break)

    /**************************************************/
    /*                                                */
    /* Enough to fit simple model and confirm a break.*/
    /*                                                */
    /**************************************************/
    /* fit all curve 09102019 SY */
    //if ((*num_curve == 0 && *i_start - *i_dense >= conse)||(*i_start - *prev_i_break >= MIN_NUM_C * N_TIMES && clrx[*i_start]- clrx[*prev_i_break] > NUM_YEARS)) // change conse to MIN_NUM_C*times
    //if ((*num_curve == 0 && *i_start - *i_dense >= SCCD_MAX_NUM_C * N_TIMES) ||(*num_curve > 0 && *i_start  - *prev_i_break >= SCCD_MAX_NUM_C * N_TIMES && clrx[*i_start]- clrx[*prev_i_break] > 1.5 * NUM_YEARS)) // change conse to MIN_NUM_C*time
    if(*num_curve == 0 && *i_start - *i_dense >= conse && clrx[*i_start] - clrx[*i_dense] >= min_days_conse)
    {
        /**********************************************/
        /*                                            */
        /* Defining computed variables.               */
        /*                                            */
        /**********************************************/
        // printf("%d\n", conse);
        fit_cft_tmp = (double **) allocate_2d_array (n_total_variable, LASSO_COEFFS, sizeof (double));
        if (fit_cft_tmp == NULL)
        {
            RETURN_ERROR ("Allocating fit_cft_tmp memory", FUNC_NAME, FAILURE);
        }
        for (i_b = 0; i_b < n_total_variable; i_b++)
        {
            if(*num_curve == 0)
                status = auto_ts_fit(clrx, clry, i_b, i_b,  *i_dense, *i_start,
                         MIN_NUM_C, fit_cft_tmp, &rmse[i_b], tmp_v_dif);
            else
                status = auto_ts_fit(clrx, clry, i_b, i_b,  *prev_i_break, *i_start,
                         MIN_NUM_C, fit_cft_tmp, &rmse[i_b], tmp_v_dif);//SY 09182018
            if (status != SUCCESS)
            {
                  RETURN_ERROR ("Calling auto_ts_fit with enough observations\n",
                             FUNC_NAME, FAILURE);
            }

        }

        ini_conse = conse;

        v_dif_magg = (double **) allocate_2d_array(n_total_variable,
                    ini_conse, sizeof (double));
        if (v_diff == NULL)
        {
            RETURN_ERROR ("Allocating v_diff memory",
                           FUNC_NAME, FAILURE);
        }


        for (i_conse = 1; i_conse < ini_conse + 1; i_conse++) // SY 09192018
        {
            v_dif_norm = 0.0;
            for (i_b = 0; i_b < n_total_variable; i_b++)
            {

                /**********************************/
                /*                                */
                /* Absolute differences.          */
                /*                                */
                /**********************************/

                // SY 09192018 moving fitting into (i_b == focus_blist[b])to save time //
                // SY 02/13/2019 delete these speed-up modification as non-lasso bands
                // are important for change agent classification
                auto_ts_predict(clrx, fit_cft, update_num_c, i_b, *i_start-i_conse,
                                *i_start-i_conse, &ts_pred_temp);
                v_dif_magg[i_b][i_conse - 1] = (double)clry[i_b][*i_start-i_conse] -
                                   ts_pred_temp;// SY 09192018

            }
       }



        rec_cg[*num_curve].t_end = clrx[*i_start-1];
        //rec_cg[*num_curve].pos.row = row;
        //rec_cg[*num_curve].pos.col = col;

        /**********************************************/
        /*                                            */
        /* Record break time, fit category, change    */
        /* probability, time of curve start, number   */
        /* of observations, change magnitude.         */
        /*                                            */
        /**********************************************/

        rec_cg[*num_curve].t_break = clrx[*i_start];

        rec_cg[*num_curve].category = 10;
        // rec_cg[*num_curve].change_prob = 100;
        rec_cg[*num_curve].t_start = clrx[0];
        rec_cg[*num_curve].num_obs = *i_start - *prev_i_break;  //SY 09182018
        rec_cg[*num_curve].t_confirmed = clrx[*i_start + i_conse - 1];
        rec_cg[*num_curve].change_prob = 100;
        *prev_i_break = *i_start;

        for (i_b = 0; i_b < n_total_variable; i_b++)
        {
            quick_sort_double(v_dif_magg[i_b], 0, ini_conse-1);
            matlab_2d_double_median(v_dif_magg, i_b, ini_conse,
                                  &v_dif_mean);
            mini_rmse = max((double)adj_rmse[i_b], rmse[i_b]);
            rec_cg[*num_curve].magnitude[i_b] = -v_dif_mean / mini_rmse;
            //rec_cg[*num_curve].magnitude[i_b] = -v_dif_mean;
        }

        for (i_b = 0; i_b < n_total_variable; i_b++)
        {

            /**********************************************/
            /*                                            */
            /* Record time of curve end,                  */
            /* postion of the pixels.                     */
            /*                                            */
            /**********************************************/

            for (k = 0; k < SCCD_MAX_NUM_C; k++)
            {
                /**************************************/
                /*                                    */
                /* Record fitted coefficients.        */
                /*                                    */
                /**************************************/
                if(k < MIN_NUM_C)
                    rec_cg[*num_curve].coefs[i_b][k] = (double)fit_cft_tmp[i_b][k];
                else
                    rec_cg[*num_curve].coefs[i_b][k] = 0;
            }

            /******************************************/
            /*                                        */
            /* Record rmse of the pixel.              */
            /*                                        */
            /******************************************/

            rec_cg[*num_curve].rmse[i_b] = rmse[i_b];


         }

        /**********************************************/
        /*                                            */
        /* Identified and move on for the next        */
        /* functional curve.                          */
        /*                                            */
        /**********************************************/
        status = free_2d_array ((void **) fit_cft_tmp);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: fit_cft_tmp\n", FUNC_NAME,FAILURE);
        }


        *num_curve = *num_curve + 1;
        if (*num_curve >= NUM_FC)
        {
            /******************************************/
            /*                                        */
            /* Reallocate memory for rec_cg.          */
            /*                                        */
            /******************************************/

            rec_cg = realloc(rec_cg, (*num_curve + 1) * sizeof(Output_t_sccd));
            if (rec_cg == NULL)
            {
                RETURN_ERROR("ERROR allocating rec_cg memory",
                             FUNC_NAME, FAILURE);
            }
        }

        status = free_2d_array ((void **) v_dif_magg);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: v_diff\n",
                          FUNC_NAME, FAILURE);
        }
    }

    status = free_2d_array ((void **) rec_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: rec_v_dif\n",
                      FUNC_NAME, FAILURE);
    }

    free(ids);
    free(bl_ids);
    free(rm_ids);
//    free(rmse);
//    rmse = NULL;
    free(v_start);
    free(v_end);
    free(v_slope);
    free(v_dif);


    status = free_2d_array ((void **) tmp_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: tmp_v_dif\n",
                      FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **) v_dif_mag);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: v_dif_mag\n",
                   FUNC_NAME, FAILURE);
    }

    return (SUCCESS);
}

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
)
{
    int num_c = SCCD_MAX_NUM_C;
    int i_span;
    i_span = cur_i - i_start +1 - 1;
    double **rec_v_dif;
    char FUNC_NAME[] = "step1_update_cft";
    int status;
    double v_dif_norm;
    int i_b, b;
    int m;
    double *vec_mag;
    double **v_diff;
    double **v_dif_mag;
    double ts_pred_temp;
    double mean_angle;
    int RETURN_VALUE;
    double break_mag;
    double mini_rmse;
    double tmp;
    int i, k;
    int i_count;
    int update_num_c;


    rec_v_dif = (double **)allocate_2d_array(n_total_variable, n_clr,
                                     sizeof (double));
    if (rec_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating rec_v_dif memory",FUNC_NAME, FAILURE);
    }

    vec_mag = (double *)calloc(adj_conse, sizeof(double));
    if (vec_mag == NULL)
    {
        RETURN_ERROR ("Allocating vec_mag memory", FUNC_NAME, FAILURE);
    }

    v_diff = (double **) allocate_2d_array(n_focus_variable,
                      adj_conse, sizeof (double));
    if (v_diff == NULL)
    {
        RETURN_ERROR ("Allocating v_diff memory",
                      FUNC_NAME, FAILURE);
    }

    v_dif_mag = (double **) allocate_2d_array(n_total_variable, adj_conse,
                sizeof (double));
    if (v_dif_mag == NULL)
    {
        RETURN_ERROR ("Allocating v_dif_mag memory",
                                 FUNC_NAME, FAILURE);
    }

    update_cft(i_span, N_TIMES, MIN_NUM_C, MID_NUM_C, SCCD_MAX_NUM_C,
              num_c, &update_num_c);


    for (i_b = 0; i_b < n_total_variable; i_b++)
    {

        status = auto_ts_fit(clrx, clry, i_b, i_b, i_start, cur_i - 1, update_num_c,
                             fit_cft, &rmse[i_b], rec_v_dif);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Calling auto_ts_fit during continuous monitoring\n",
                          FUNC_NAME, FAILURE);
        }

    }



    /**********************************************/
    /*                                            */
    /* Detect change, value of difference for     */
    /* adj_conse observations.                    */
    /*                                            */
    /**********************************************/
    for (m = 0; m < adj_conse; m++)
    {
        vec_mag[m] = 0;
        for (b = 0; b < n_focus_variable; b++)
            v_diff[b][m] = 0;
        for (b = 0; b < n_total_variable; b++)
            v_dif_mag[b][m] = 0;
    }

    for (i = 0; i < adj_conse; i++)//SY 09192018
    {
        v_dif_norm = 0.0;
        for (i_b = 0; i_b < n_total_variable; i_b++)
        {
            /**************************************/
            /*                                    */
            /* Absolute differences.              */
            /*                                    */
            /**************************************/

            auto_ts_predict(clrx, fit_cft, update_num_c, i_b, cur_i+i, cur_i+i,
                            &ts_pred_temp); //SY 09192018
            v_dif_mag[i_b][i] = (double)clry[i_b][cur_i+i] - ts_pred_temp;//SY 09192018

            /**************************************/
            /*                                    */
            /* Normalize to z-score.              */
            /*                                    */
            /**************************************/

            for (b = 0; b < n_focus_variable; b++)
            {
                if (i_b == focus_blist[b])
                {
                    /******************************/
                    /*                            */
                    /* Minimum rmse,              */
                    /* z-scores.                  */
                    /*                            */
                    /******************************/

                    mini_rmse = max(adj_rmse[i_b], rmse[i_b]);
                    v_diff[b][i] = v_dif_mag[i_b][i] / mini_rmse;
                    v_dif_norm += v_diff[b][i] * v_diff[b][i];
                }
            }
        }
        vec_mag[i] = v_dif_norm;//SY 09192018
    }

    mean_angle = MeanAngl(v_diff, n_focus_variable, adj_conse);
    break_mag = 9999.0;
    for (m = 0; m < adj_conse; m++)
    {
        if (break_mag > vec_mag[m])
        {
            break_mag = vec_mag[m];
        }
    }

    if ((break_mag > reg_TCG) && (mean_angle < NSIGN))
    {

        /**********************************************/
        /*                                            */
        /* Record break time.                        */
        /*                                            */
        /**********************************************/

        rec_cg[*num_curve].t_break = clrx[cur_i];
        for (i_b = 0; i_b < n_total_variable; i_b++)
        {
            quick_sort_double(v_dif_mag[i_b], 0, adj_conse-1);
            matlab_2d_double_median(v_dif_mag, i_b, adj_conse,
                                   &tmp);
            mini_rmse = max(adj_rmse[i_b], rmse[i_b]);
            rec_cg[*num_curve].magnitude[i_b] = tmp / mini_rmse;
            //rec_cg[*num_curve].magnitude[i_b] = tmp;
        }

        /**********************************************/
        /*                                            */
        /* Updating information for the first         */
        /* iteration.  Record time of curve start and */
        /* time of curve end.                         */
        /*                                            */
        /**********************************************/

        rec_cg[*num_curve].t_start = clrx[i_start];
        rec_cg[*num_curve].t_end = clrx[cur_i-1];
        rec_cg[*num_curve].num_obs = cur_i -i_start;
        rec_cg[*num_curve].t_confirmed = clrx[cur_i + adj_conse - 1];
        rec_cg[*num_curve].change_prob = 100;
        /**********************************************/
        /*                                            */
        /* No break at the moment.                    */
        /*                                            */
        /**********************************************/

        rec_cg[*num_curve].t_break = clrx[cur_i];

        *prev_i_break = cur_i;
        /**********************************************/
        /*                                            */
        /* Record change probability, number of       */
        /* observations, fit category.                */
        /*                                            */
        /**********************************************/

        rec_cg[*num_curve].category = 30;

//        for(i_b = 0; i_b < n_total_variable; i_b++)
//        {
//            rec_cg[*num_curve].obs_disturb[i_b] = 0;
//            for (k = 0; k < SCCD_MAX_NUM_C - 1; k++)
//            {
//                rec_cg[*num_curve].state_disturb[i_b][k] = 0;
//            }
//        }

        for (i_b = 0; i_b < n_total_variable; i_b++)
        {

            /******************************************/
            /*                                        */
            /* Record rmse of the pixel.              */
            /*                                        */
            /******************************************/

            rec_cg[*num_curve].rmse[i_b] = rmse[i_b];


            for (k = 0; k < SCCD_MAX_NUM_C; k++)
            {
                /**************************************/
                /*                                    */
                /* Record fitted coefficients.        */
                /*                                    */
                /**************************************/

                rec_cg[*num_curve].coefs[i_b][k] = fit_cft[i_b][k];
             }


         }
        /**********************************************/
        /*                                            */
        /* Identified and move on for the next        */
        /* functional curve.                          */
        /*                                            */
        /**********************************************/

        *num_curve = *num_curve + 1;

        if (*num_curve >= NUM_FC)
        {
            /******************************************/
            /*                                        */
            /* Reallocate memory for rec_cg.          */
            /*                                        */
            /******************************************/

            rec_cg = realloc(rec_cg, (*num_curve + 1) * sizeof(Output_t_sccd));
            if (rec_cg == NULL)
            {
                RETURN_ERROR("ERROR allocating rec_cg memory",
                                     FUNC_NAME, FAILURE);
            }
        }

        RETURN_VALUE = CHANGEDETECTED;

    }
    else if (vec_mag[0] > t_cg_outelier)  /*false change*/
    {
        /**********************************************/
        /*                                            */
        /* Remove noise.                              */
        /*                                            */
        /**********************************************/

        for (m = i; m < end -1; m++)
        {
            clrx[m] = clrx[m+1];
            for (b = 0; b < n_total_variable; b++)
                clry[b][m] = clry[b][m+1];
        }
        RETURN_VALUE = FALSECHANGE;
    }

    else
        RETURN_VALUE = REGULAREND;



    status = free_2d_array ((void **) rec_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: rec_v_dif\n",
                      FUNC_NAME, FAILURE);
    }
    free(vec_mag);

    status = free_2d_array ((void **) v_diff);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: v_diff\n",
                      FUNC_NAME, FAILURE);
    }
    status = free_2d_array ((void **) v_dif_mag);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: v_dif_mag\n",
                   FUNC_NAME, FAILURE);
    }

    return (RETURN_VALUE);

}


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
)
{
    int i, j;
    gsl_vector* kt;
    gsl_vector* kt_tmp;
    double yt_tmp;
    double ft_tmp;
//    gsl_matrix* pt;         /*  A m x m matrix containing the covariance matrix for last_pred_loc */
//    gsl_vector* at;
    gsl_matrix* mm;
    mm = gsl_matrix_alloc(instance->m, instance->m);

    kt = gsl_vector_alloc(instance->m);
    kt_tmp = gsl_vector_alloc(instance->m);


    /* make copy so that recursion won't really change values*/
//    gsl_matrix_memcpy(pt, P_ini);
//    gsl_vector_memcpy(at, at_ini);

    /* the predict from clrx[pred_start] to clrx[pred_end] as if they are all missing obs */
    for (i = pred_start; i < pred_end + 1; i++)
    {
        /* this loop predicts every values between two observation dates */
        for (j = 0; j < clrx[i + 1] - clrx[i]; j++) /* predict ith observation */
        {
            /**for observation date, we need make predication
             * but without updating KF parameters */
            if (j == 0)
            {
                /* first get at from fit_cft*/
                //if (b_fastmode == TRUE)
                fit_cft2vec_a(fit_cft, at_ini, clrx[i], instance->m, instance->structure);

                /* predict y without updating */
                gsl_blas_ddot(instance->Z, at_ini, &yt_tmp);

                pred_y[i - pred_start] = yt_tmp;
                if(b_foutput == TRUE)
                {
                    /* kt = pt*zt */
                    gsl_blas_dsymv(CblasUpper, 1.0, P_ini,
                                   instance->Z, 0.0, kt_tmp);

                    /* ft = kt *ztt + ht */
                    gsl_blas_ddot(instance->Z, kt_tmp, &ft_tmp);
                    ft_tmp = ft_tmp + instance->H;

                    pred_y_f[i - pred_start] = ft_tmp;
                }
                else
                    pred_y_f[i - pred_start] = 0;

                //printf("ft for %d is %f: \n", clrx[i] + j, ft_tmp);

                /* update fit_cft using new at*/
//                if (b_fastmode == TRUE)
//                    vec_a2fit_cft(at_ini, fit_cft, clrx[i]);
            }

            if(b_foutput == TRUE)
            {
                filter1step_missingobs(instance->Z, instance->H, instance->T,
                                   instance->Q, at_ini, &ft_tmp, P_ini, kt, instance->m, b_fastmode);

//                gsl_blas_dsymm(CblasRight, CblasUpper, 1.0, P_ini, instance->T, 0.0, mm);

//                /* pt = mm * tt^T */
//                gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, P_ini, instance->T, 0.0, P_ini);
            }


             //printf("ft for %d is %f: \n", clrx[i] + j, ft_tmp);


        }
    }


    gsl_vector_free(kt);
//    gsl_matrix_free(pt);
//    gsl_vector_free(at);
    gsl_vector_free(kt_tmp);
    gsl_matrix_free(mm);
    return SUCCESS;

}


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
)
{
    int i;
    gsl_vector* kt_tmp;
    double yt_tmp;
    double ft_tmp;
    kt_tmp = gsl_vector_alloc(instance->m);

    /* the predict from clrx[pred_start] to clrx[pred_end] as if they are all missing obs */


    /**for observation date, we need make predication
     * but without updating KF parameters */

        /* predict y without updating */
    if (b_fastmode == TRUE)
        fit_cft2vec_a(fit_cft, at_ini, clrx, instance->m, instance->structure);

    gsl_blas_ddot(instance->Z, at_ini, &yt_tmp);

    *pred_y = (double)yt_tmp;
    /*kt = pt*zt*/
    gsl_blas_dsymv(CblasUpper, 1.0, P_ini,
                   instance->Z, 0.0, kt_tmp);

    /*ft = kt *ztt + ht*/
    gsl_blas_ddot(instance->Z, kt_tmp, &ft_tmp);
    *pred_y_f = (double)ft_tmp + (double)instance->H;

        //printf("ft for %d is %f: \n", clrx[i] + j, ft_tmp);

             //printf("ft for %d is %f: \n", clrx[i] + j, ft_tmp);

    gsl_vector_free(kt_tmp);
    return SUCCESS;

}

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
)
{
    int j;
    gsl_vector* kt;
    gsl_vector* att; // filtered states
    double vt;
    double ft;
    double tmp_smooth_Q;
    double tmp_smooth_H;
    double old_state;
    int d_yr;

    kt = gsl_vector_alloc(instance->m);
    att = gsl_vector_alloc(instance->m);

    /* this loop predicts every missing values between two observation dates */
    for (j = 0; j < clrx[cur_i + 1] - clrx[cur_i]; j++) /* predict ith observation */
    {
//        level_state_records[i_b][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(at_ini, 0);
//        if (b_fastmode == FALSE)
//        {
//            /* first, save state info for each date for current i*/
//            annual_state_records[i_b][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(at_ini, 1);
//            semi_state_records[i_b][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(at_ini, 3);
//            if(instance->m == 7)
//                third_state_records[i_b][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(at_ini, 5);
//            else
//                third_state_records[i_b][clrx[cur_i] + j - starting_date] = 0;
//        }

//        if (cur_i == 93)
//        {
//            printf("%d\n", clrx[cur_i] + j);
//            printf("%f\n", level_state_records[0][clrx[cur_i] + j - starting_date]);
//            printf("%f\n", annual_state_records[0][clrx[cur_i] + j - starting_date]);
//            printf("%f\n", semi_state_records[0][clrx[cur_i] + j - starting_date]);
//        }
        if (0 == j)
        {
            /* first get at from fit_cft*/

            fit_cft2vec_a(fit_cft, at_ini, clrx[cur_i], instance->m, instance->structure);

            // printf("i = %d: %f, %f, %f\n", clrx[cur_i], gsl_vector_get(at_ini, 0), gsl_vector_get(at_ini, 1), gsl_vector_get(at_ini, 3));

            /* record old state */
            old_state = gsl_vector_get(at_ini, 0);

            /* predicts valid obs values between two observation dates */
            filter1step_validobs((double)clry[cur_i], instance->Z, &instance->H, instance->T, instance->Q,
                            at_ini, P_ini, &vt, &ft, kt, instance->m, b_fastmode, FALSE, att);

            level_state_records[i_b][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(att, 0);
            if (b_fastmode == FALSE)
            {
                if(instance->structure % 10 == 1)
                    annual_state_records[i_b][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(att, 1);
                else
                    annual_state_records[i_b][clrx[cur_i] + j - starting_date] = 0;

                /* first, save state info for each date for current i*/
                if(instance->structure / 10 == 1)
                    semi_state_records[i_b][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(att, 3);
                else
                    semi_state_records[i_b][clrx[cur_i] + j - starting_date] = 0;

                if(instance->structure / 100 == 1)
                    third_state_records[i_b][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(att, 5);
                else
                    third_state_records[i_b][clrx[cur_i] + j - starting_date] = 0;
            }

            if(b_start == FALSE)
            {
                tmp_smooth_Q = (gsl_vector_get(att, 0) - old_state) / (sqrtf((double)abs(clrx[cur_i] - clrx[cur_i - 1])));
                // tmp_smooth_Q = (gsl_vector_get(att, 0) - old_state) / (abs(clrx[cur_i] - clrx[cur_i - 1]));
                // tmp_smooth_Q = vt * (1 - (double)instance ->H / ft) / (sqrtf((double)abs(clrx[cur_i] - clrx[cur_i - 1])));
//                if(cur_i == 92)
//                {
//                    printf("%d ,%d\n", clrx[cur_i], clrx[cur_i - 1]);
//                }
                // sum_square_vt[i_b] = sum_square_vt[i_b] + vt * vt;
                if((instance ->H > ft)||(ft == 0))
                {
                    tmp_smooth_H = vt;
                    sum_square_vt[i_b] = sum_square_vt[i_b] + vt * vt;
                    sum_vt[i_b] = sum_vt[i_b] + vt;
                    sum_kalman_coef[i_b] = sum_kalman_coef[i_b] + 1;
                }
                else
                {
                    tmp_smooth_H = vt * (double)instance ->H / ft;
                    //sum_square_vt[i_b] = sum_square_vt[i_b] + (double) vt * vt * (1 - (double)instance ->H / ft) * (1 - (double)instance ->H / ft);
                    sum_square_vt[i_b] = sum_square_vt[i_b] + vt * vt;
                    sum_vt[i_b] = sum_vt[i_b] + vt;
                    sum_kalman_coef[i_b] = sum_kalman_coef[i_b] + (1 - (double)instance ->H / ft);
                }
                sum_smooth_Q[i_b] = sum_smooth_Q[i_b] + tmp_smooth_Q;
                sum_square_smooth_Q[i_b] = sum_square_smooth_Q[i_b] + tmp_smooth_Q * tmp_smooth_Q;
                sum_smooth_H[i_b] = sum_smooth_H[i_b] + tmp_smooth_H;
                sum_square_smooth_H[i_b] = sum_square_smooth_H[i_b] + tmp_smooth_H * tmp_smooth_H;
                valid_count[i_b] = valid_count[i_b] + 1;
                d_yr = clrx[cur_i] - (int)((int)(clrx[cur_i] / NUM_YEARS) * NUM_YEARS);

                temporal_rmse_square_sum[i_b][d_yr/TEMPORAL_RMSE_BIN] = temporal_rmse_square_sum[i_b][d_yr/TEMPORAL_RMSE_BIN] + vt * vt;
                temporal_rmse_sum[i_b][d_yr/TEMPORAL_RMSE_BIN] = temporal_rmse_sum[i_b][d_yr/TEMPORAL_RMSE_BIN] + vt;
                if(i_b == 0) // using i_b == 0 to avoid repeated count
                   temporal_rmse_count[d_yr/TEMPORAL_RMSE_BIN] = temporal_rmse_count[d_yr/TEMPORAL_RMSE_BIN] + 1;

            }

            /* update fit_cft using new at*/

            vec_a2fit_cft(at_ini, fit_cft, clrx[cur_i]+1, instance->m, instance->structure);

//            if((instance ->H > ft)||(ft == 0))
//                rmse_records[i_b][cur_i] = (double) (vt);
//            else
//                rmse_records[i_b][cur_i] = (double) (vt) * sqrt((1 - instance ->H / ft));

            //rmse_records[i_b][cur_i] = (double) (vt) / sqrtf((double)ft);
             rmse_records[i_b][cur_i] = (double) (vt);

            // printf("rmse for %d time point for band %d: %f\n", cur_i, i_band, rmse[i_b][cur_i]);
        }
        else
        {
            level_state_records[i_b][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(at_ini, 0);
            if (b_fastmode == FALSE)
            {
                /* first, save state info for each date for current i*/
                annual_state_records[i_b][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(at_ini, 1);
                semi_state_records[i_b][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(at_ini, 3);
                if(instance->m == 7)
                    third_state_records[i_b][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(at_ini, 5);
                else
                    third_state_records[i_b][clrx[cur_i] + j - starting_date] = 0;
            }
            filter1step_missingobs(instance->Z, instance->H, instance->T, instance->Q,
                            at_ini, &ft, P_ini,  kt, instance->m, b_fastmode);

        }
        // printf("ft for %d is %f: \n", clrx[cur_i] + j, ft);
    }

    gsl_vector_free(kt);
    gsl_vector_free(att);

    return SUCCESS;
}

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
)
{
    int j;
    gsl_vector* kt;
    gsl_vector* att;
    double vt;
    double ft;

    kt = gsl_vector_alloc(instance->m);
    att = gsl_vector_alloc(instance->m);
    /* this loop predicts every missing values between two observation dates */
    for (j = 0; j < clrx[cur_i + 1] - clrx[cur_i]; j++) /* predict ith observation */
    {
        if (0 == j)
        {
            /* predicts valid obs values between two observation dates */
            filter1step_validobs(clry[cur_i], instance->Z, &instance->H, instance->T, instance->Q,
                            at_ini, P_ini, &vt, &ft, kt, instance->m, FALSE, FALSE, att);

            // printf("rmse for %d time point for band %d: %f\n", cur_i, i_band, rmse[i_b][cur_i]);
        }
        else
        {
            filter1step_missingobs(instance->Z, instance->H, instance->T, instance->Q,
                            at_ini, &ft, P_ini,  kt, instance->m, FALSE);
        }

        // printf("ft for %d is %f: \n", clrx[cur_i] + j, ft);
    }

    gsl_vector_free(kt);
    gsl_vector_free(att);

    return SUCCESS;
}

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
    int starting_date,        /* starting_date is used to compute record index based on clrx*/
    bool b_fastmode
)
{
    int j;
    gsl_vector* kt;
    double ft;

    kt = gsl_vector_alloc(instance->m);

    /* this loop predicts every missing values between two observation dates */
    for (j = 0; j < clrx[cur_i+1] - clrx[cur_i]; j++) /* predict ith observation */
    {
        if(b_fastmode == FALSE)
        {
            /* first, save state info for each date for current i*/
            level_state_records[i_band][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(at_ini, 0);
            annual_state_records[i_band][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(at_ini, 1);
            semi_state_records[i_band][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(at_ini, 3);
            if(instance->m == 7)
                third_state_records[i_band][clrx[cur_i] + j - starting_date] = (double)gsl_vector_get(at_ini, 5);
            else
                third_state_records[i_band][clrx[cur_i] + j - starting_date] = 0;
        }


        filter1step_missingobs(instance->Z, instance->H, instance->T, instance->Q,
                        at_ini, &ft, P_ini,  kt, instance->m, b_fastmode);

        /* next, save f as f updates after filtering*/
        // f_records[i_band][clrx[cur_i]  + j - starting_date] = (double) ft;
        // printf("ft for %d is %f: \n", clrx[cur_i] + j, ft);
    }

    gsl_vector_free(kt);

    return SUCCESS;

}

/************************************************************************
FUNCTION: step2_KF_ChangeDetection

PURPOSE:
Step 2 of S-CCD: change detection using kalman filter.
RETURN VALUE:
Type = int (SUCCESS OR FAILURE)

Programmer: Su Ye
**************************************************************************/

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
    double **fit_cft1,
    double **fit_cft2,
    double **fit_cft3,
    double **fit_cft4,
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
    int *i_count,
    float *adj_rmse,
    double t_cg_adjust,
    double t_cg_outelier,
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
)
{
    int i_b, b, m, i, k;
    int center_date;
    int status;
    int i_conse;
    // int j_conse;
    double *pred_y;
    double *pred_y_f;
    char FUNC_NAME[] = "step2_KF_ChangeDetection";
    double pred_y_single;   /* O: the predicted obs values */
    double pred_f_single;  /*O: the predicted f (RMSE) values */
    double v_dif_mag_standard;           /* normalized magnitude, =v_diff in ccd.c */
    double **v_dif_mag;
    //double **v_dif_gradiant;
    double **v_dif;
    double *v_dif_mag_norm;
    double max_rmse_tmp;
    double *medium_v_dif;
    // double medium_v_dif_gradient[n_focus_variable];
    // bool grad_change_flag;
    bool abrupt_change_flag;
    double mean_angle;
    int RETURN_VALUE;
    // double c0, c1;
    double **fit_cft_average;
    gsl_matrix** copy_vec_next_P;
    gsl_vector** copy_vec_next_a;


    double* tmpcg_rmse_square; /* to temporarily change RMSE          */
    double* tmpcg_rmse;
    double* max_rmse;

    double Q_tmp[SCCD_MAX_NUM_C - 1];

    int n_rmse;
    double *d_yr;

    int d_rt;
    double **rec_v_dif_copy;
    double sum_median_tmp;
    double prob_conse = 0.0;
    double prob_conse_avg = 0.0;
    double mediam_mag;
    int center_dyr;
    double weight1, weight2, weight3, weight4;
    int max_num;
    double rmse_tmp;

    v_dif = (double **) allocate_2d_array (n_focus_variable, conse, sizeof (double));
    if (v_dif == NULL)
    {
        RETURN_ERROR ("Allocating v_dif memory",
                      FUNC_NAME, FAILURE);
    }

    v_dif_mag = (double **) allocate_2d_array(n_total_variable, conse,
                sizeof (double));
    if (v_dif_mag == NULL)
    {
        RETURN_ERROR ("Allocating v_dif_mag memory",
                                 FUNC_NAME, FAILURE);
    }

//    v_dif_gradiant = (double **) allocate_2d_array(n_total_variable, (conse - 1) * (conse - 2) /2,
//                sizeof (double));
//    if (v_dif_gradiant == NULL)
//    {
//        RETURN_ERROR ("Allocating v_dif_gradiant memory",
//                                 FUNC_NAME, FAILURE);
//    }


    copy_vec_next_a = (gsl_vector **) allocate_2d_array(n_total_variable, 1, sizeof(gsl_vector));
    if (copy_vec_next_a == NULL)
    {
        RETURN_ERROR ("Allocating copy_vec_next_a memory", FUNC_NAME, FAILURE);
    }

    copy_vec_next_P = (gsl_matrix **)allocate_2d_array(n_total_variable, 1, sizeof(gsl_matrix));
    if (copy_vec_next_P == NULL)
    {
        RETURN_ERROR ("Allocating copy_vec_next_P memory", FUNC_NAME, FAILURE);
    }

    fit_cft_average = (double **) allocate_2d_array (n_total_variable, LASSO_COEFFS, sizeof (double));
    if (fit_cft_average == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft_average memory", FUNC_NAME, FAILURE);
    }

    pred_y_f = (double *)malloc(conse * sizeof(double));
    if (pred_y_f == NULL)
    {
        RETURN_ERROR ("Allocating pred_y_f memory", FUNC_NAME, FAILURE);
    }

    pred_y = (double *)malloc(conse * sizeof(double));
    if (pred_y == NULL)
    {
        RETURN_ERROR ("Allocating pred_y memory", FUNC_NAME, FAILURE);
    }

    v_dif_mag_norm = (double *)malloc(conse * sizeof(double));
    if (v_dif_mag_norm == NULL)
    {
        RETURN_ERROR ("Allocating v_dif_mag_norm memory", FUNC_NAME, FAILURE);
    }

    medium_v_dif = (double *)malloc(n_focus_variable * sizeof(double));
    if (medium_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating v_dif_mag_norm memory", FUNC_NAME, FAILURE);
    }

    rec_v_dif_copy = (double **)allocate_2d_array(n_total_variable, cur_i - i_start,
                                     sizeof (double));
    if (rec_v_dif_copy == NULL)
    {
        RETURN_ERROR ("Allocating rec_v_dif_copy memory",FUNC_NAME, FAILURE);
    }

    d_yr = (double*)malloc((valid_count[0] + 1) * sizeof(double));
    if (d_yr == NULL)
    {
        RETURN_ERROR ("Allocating d_yr memory",
                     FUNC_NAME, FAILURE);
    }

    tmpcg_rmse_square= (double*)malloc(n_total_variable * sizeof(double)); /* to temporarily change RMSE          */
    if (tmpcg_rmse_square == NULL)
    {
        RETURN_ERROR ("Allocating tmpcg_rmse_square memory",
                     FUNC_NAME, FAILURE);
    }

    tmpcg_rmse = (double*)malloc(n_total_variable * sizeof(double));
    if (tmpcg_rmse == NULL)
    {
        RETURN_ERROR ("Allocating tmpcg_rmse memory",
                     FUNC_NAME, FAILURE);
    }

    max_rmse = (double*)malloc(n_total_variable * sizeof(double));
    if (max_rmse == NULL)
    {
        RETURN_ERROR ("Allocating max_rmse memory",
                     FUNC_NAME, FAILURE);
    }
    /* make a copy for fit_cft, initial P and at */
    for(i_b = 0; i_b < n_total_variable; i_b++)
    {
        copy_vec_next_a[i_b] =  gsl_vector_alloc (instance[i_b].m);
        copy_vec_next_P[i_b] = gsl_matrix_alloc (instance[i_b].m, instance[i_b].m);
        gsl_matrix_memcpy(copy_vec_next_P[i_b], vec_next_P[i_b]);
        gsl_vector_memcpy(copy_vec_next_a[i_b], vec_next_a[i_b]);
        // var_smooth_Q[i_b] = sum_square_smooth_Q[i_b] / (valid_count[i_b] - 1);
        // var_smooth_H[i_b] = sum_square_smooth_H[i_b] / (valid_count[i_b] - 1);

    }


//    weight2 = 0;
//    weight3 = 0;
//    weight4 = 0;

    center_date = (int)(clrx[cur_i + conse - 1] + clrx[cur_i]) / 2;
    //center_date = clrx[cur_i + conse - 1];
//    weight1 = (double)(center_date - *clrx_record1);
//    weight2 = (double)(center_date - *clrx_record2);
//    weight3 = (double)(center_date - *clrx_record3);
//    weight4 = (double)(center_date - *clrx_record4);

//    weight1 = 1.0 / (double)(center_date - *clrx_record1);
//    weight2 = 1.0 / (double)(center_date - *clrx_record2);
//    weight3 = 1.0 / (double)(center_date - *clrx_record3);
//    weight4 = 0;
    weight1 = WEIGHT1;
    weight2 = WEIGHT2;
    weight3 = WEIGHT3;
    weight4 = WEIGHT4;
//    if(cur_i == 185){
//        cur_i = 185;
//    }
    for(i_b = 0; i_b < n_total_variable; i_b++)
    {
        for(i = 0; i < SCCD_MAX_NUM_C; i++)
            fit_cft_average[i_b][i] = fit_cft1[i_b][i] * weight1 / (weight1 + weight2 + weight3 + weight4)
                    + fit_cft2[i_b][i] * weight2 / (weight1 + weight2 + weight3 + weight4)
                    + fit_cft3[i_b][i] * weight3 / (weight1 + weight2 + weight3 + weight4)
                    + fit_cft4[i_b][i] * weight4 / (weight1 + weight2 + weight3 + weight4);
    }

    n_rmse = valid_count[0];

    if(n_rmse > LASSO_COEFFS * N_TIMES)
    {

//            for(i = 0; i < TEMPORAL_RMSE_BIN_NUM; i++)
//                printf("%d, %f, %f, %f, %f, %f, %f, %f\n", temporal_rmse_count[i],
//                       temporal_rmse_square_sum[0][i], temporal_rmse_square_sum[1][i],
//                        temporal_rmse_square_sum[2][i], temporal_rmse_square_sum[3][i],
//                        temporal_rmse_square_sum[4][i], temporal_rmse_square_sum[5][i], temporal_rmse_square_sum[6][i]);
//            for(i = 0; i < TEMPORAL_RMSE_BIN_NUM; i++)
//                printf("%d, %f, %f, %f, %f, %f, %f, %f\n", temporal_rmse_count[i],
//                       temporal_rmse_sum[0][i], temporal_rmse_sum[1][i],
//                        temporal_rmse_sum[2][i], temporal_rmse_sum[3][i],
//                        temporal_rmse_sum[4][i], temporal_rmse_sum[5][i], temporal_rmse_sum[6][i]);

        i = 1; /* step length */
        //center_date = (int) clrx[cur_i];
        center_dyr = center_date - (int)((int)(center_date / NUM_YEARS) * NUM_YEARS);
        center_dyr = center_dyr / TEMPORAL_RMSE_BIN;
        for(i_b = 0; i_b < n_total_variable; i_b++)
        {
            tmpcg_rmse_square[i_b] = temporal_rmse_square_sum[i_b][center_dyr];
            tmpcg_rmse[i_b] = temporal_rmse_sum[i_b][center_dyr];
        }

        n_rmse = temporal_rmse_count[center_dyr];

        //max_num = max(LASSO_COEFFS * N_TIMES, valid_count[0] / 4);
        while(n_rmse <= LASSO_COEFFS * N_TIMES)
        {
            if(center_dyr - i < 0)
            {
                if(temporal_rmse_count[TEMPORAL_RMSE_BIN_NUM + center_dyr - i] > 0)
                {
                    n_rmse = n_rmse + temporal_rmse_count[TEMPORAL_RMSE_BIN_NUM + center_dyr - i];
                    for(i_b = 0; i_b < n_total_variable; i_b++)
                    {
                        tmpcg_rmse_square[i_b] = tmpcg_rmse_square[i_b] + temporal_rmse_square_sum[i_b][TEMPORAL_RMSE_BIN_NUM + center_dyr - i];
                        tmpcg_rmse[i_b] = tmpcg_rmse[i_b] + temporal_rmse_sum[i_b][TEMPORAL_RMSE_BIN_NUM + center_dyr - i];
                    }
                }
            }
            else
            {
                if(temporal_rmse_count[center_dyr - i] > 0)
                {
                    n_rmse = n_rmse + temporal_rmse_count[center_dyr - i];
                    for(i_b = 0; i_b < n_total_variable; i_b++)
                    {
                        tmpcg_rmse_square[i_b] = tmpcg_rmse_square[i_b] + temporal_rmse_square_sum[i_b][center_dyr - i];
                        tmpcg_rmse[i_b] = tmpcg_rmse[i_b] + temporal_rmse_sum[i_b][center_dyr - i];
                    }
                }
            }

            if(center_dyr + i > TEMPORAL_RMSE_BIN_NUM - 1)
            {
                if(temporal_rmse_count[center_dyr + i - TEMPORAL_RMSE_BIN_NUM] > 0)
                {
                    n_rmse = n_rmse + temporal_rmse_count[center_dyr + i - TEMPORAL_RMSE_BIN_NUM];
                    for(i_b = 0; i_b < n_total_variable; i_b++)
                    {
                        tmpcg_rmse_square[i_b] = tmpcg_rmse_square[i_b] + temporal_rmse_square_sum[i_b][center_dyr + i - TEMPORAL_RMSE_BIN_NUM];
                        tmpcg_rmse[i_b] = tmpcg_rmse[i_b] + temporal_rmse_sum[i_b][center_dyr + i - TEMPORAL_RMSE_BIN_NUM];
                    }
                }
            }
            else
            {
                if(temporal_rmse_count[center_dyr + i] > 0)
                {
                    n_rmse = n_rmse + temporal_rmse_count[center_dyr + i];
                    for(i_b = 0; i_b < n_total_variable; i_b++)
                    {
                        tmpcg_rmse_square[i_b] = tmpcg_rmse_square[i_b] + temporal_rmse_square_sum[i_b][center_dyr + i];
                        tmpcg_rmse[i_b] = tmpcg_rmse[i_b] + temporal_rmse_sum[i_b][center_dyr + i];
                    }
                }
            }
            i++;
        }

        for (i_b = 0; i_b < n_total_variable; i_b++)
        {
            tmpcg_rmse_square[i_b] = sqrtf(tmpcg_rmse_square[i_b]) / sqrtf((double)n_rmse);
            tmpcg_rmse[i_b] = tmpcg_rmse[i_b] / (double)n_rmse;
//            tmpcg_rmse_square[i_b] = sqrtf(tmpcg_rmse_square[i_b]) / sqrtf((double)n_rmse - SCCD_NUM_C);
//            tmpcg_rmse[i_b] = tmpcg_rmse[i_b] / (double)(n_rmse - SCCD_NUM_C);
        }

    }
    else
    {
        for (i_b = 0; i_b < n_total_variable; i_b++)
        {
            tmpcg_rmse_square[i_b] = sqrtf(sum_square_vt[i_b]) / sqrtf((double)n_rmse);
            tmpcg_rmse[i_b] = sum_vt[i_b] / (double)n_rmse;
//            tmpcg_rmse_square[i_b] = sqrtf(sum_square_vt[i_b]) / sqrtf((double)n_rmse - SCCD_NUM_C);
//            tmpcg_rmse[i_b] = sum_vt[i_b] / (double)(n_rmse - SCCD_NUM_C);
        }
    }


    /* previous rmse number is valid_count+1, cuz we didn't count the first obs*/
//    for(m = 0; m < valid_count[0] + 1; m++)
//    {
//        d_rt =  clrx[cur_i+ conse - 1] - clrx[cur_i - 1 - m];
//        d_yr[m] = fabs((double)d_rt - round((double)d_rt / NUM_YEARS) * NUM_YEARS);
//    }

//    for (i_b = 0; i_b < n_total_variable; i_b++)
//    {
//        for (m = 0; m < valid_count[0] + 1; m++)
//            rec_v_dif_copy[i_b][m] = rmse_records[i_b][cur_i - 1 - m];
//    }

//    /**********************************************/
//    /*                                            */
//    /* Sort the rec_v_dif based on d_yr.          */
//    /*                                            */
//    /**********************************************/
//    n_rmse = valid_count[0] + 1; // plus 1 is because valid_count does not include the first obs
//    if(n_rmse > SCCD_MAX_NUM_C * N_TIMES * 2)
//       n_rmse = SCCD_MAX_NUM_C * N_TIMES * 2;

//    quick_sort_2d_double(d_yr, rec_v_dif_copy, 0, valid_count[0], n_total_variable);

//    for(i_b = 0; i_b < n_total_variable; i_b++)
//        tmpcg_rmse_square[i_b] = 0.0;

//    /**********************************************/
//    /*                                            */
//    /* Temporarily changing RMSE.                 */
//    /*                                            */
//    /**********************************************/

//    for (i_b = 0; i_b < n_total_variable; i_b++)
//    {
//        matlab_2d_array_norm(rec_v_dif_copy, i_b, n_rmse,
//                         &tmpcg_rmse_square[i_b]);
//        tmpcg_rmse_square[i_b] /= sqrt(n_rmse);
//    }


    for (k = 0; k < conse; k++)
    {
        v_dif_mag_norm[k] = 0;
    }

    // v_dif_mag_norm_ini = 0;
    /* eliminate gradual change */
    abrupt_change_flag = TRUE;

    for(i_conse = 0; i_conse < conse; i_conse++)
    {
        for(i_b = 0; i_b < n_total_variable; i_b++)
        {
//            if(i_conse == 0)
//            {
//                KF_ts_predict_conse(&instance[i_b], clrx, copy_vec_next_P[i_b] , copy_vec_next_a[i_b], cur_i + i_conse,
//                            cur_i + i_conse, &pred_y_single, &pred_f_single, b_fastmode, fit_cft_average[i_b], TRUE);
//                max_rmse[i_b] = sqrtf(pred_f_single);
//            }
//            else
                /* kalman predict "conse" consecutive points without updates */
//                KF_ts_predict_conse(&instance[i_b], clrx, copy_vec_next_P[i_b] , copy_vec_next_a[i_b], cur_i + i_conse,
//                            cur_i + i_conse, &pred_y_single, &pred_f_single, b_fastmode, fit_cft_average[i_b], FALSE);
                KF_ts_predict_conse(&instance[i_b], clrx, copy_vec_next_P[i_b] , copy_vec_next_a[i_b], cur_i + i_conse,
                            cur_i + i_conse, &pred_y_single, &pred_f_single, b_fastmode, fit_cft_average[i_b], FALSE);

             //max_rmse[i_b] = max(adj_rmse[i_b], sqrtf(sum_square_vt[i_b] / (valid_count[i_b] - 1)));
//            if(pred_f_single < instance[i_b].H)
//                kalman_ratio = 1;
//            else
            //kalman_ratio[i_b] = sum_kalman_coef[i_b] / valid_count[i_b];
            //kalman_ratio[i_b] = 1 - (double)instance[i_b].H / pred_f_single;
            //kalman_ratio[i_b] = sum_kalman_coef[i_b] / valid_count[i_b];
            //max_rmse[i_b] = max(adj_rmse[i_b] * kalman_ratio[i_b], sqrtf(sum_square_vt[i_b] / (double)(valid_count[i_b]) - SCCD_MAX_NUM_C + 1));
            //max_rmse[i_b] = max(adj_rmse[i_b], sqrtf(pred_f_single));
            //max_rmse[i_b] = max(adj_rmse[i_b], sqrtf(sum_square_vt[i_b] / (double) (valid_count[i_b])));
            //max_rmse[i_b] = max(sqrtf(pred_f_single), tmpcg_rmse_square[i_b]);
            max_rmse[i_b] = max(adj_rmse[i_b], tmpcg_rmse_square[i_b]);
            //max_rmse[i_b] = tmpcg_rmse_square[i_b];
            //max_rmse[i_b] = sqrtf(pred_f_single);
            //max_rmse[i_b] = max(tmpcg_rmse_square[i_b], sqrtf(pred_f_single));
//            if(adj_rmse[i_b] > tmpcg_rmse_square[i_b])
//               max_rmse[i_b] =  adj_rmse[i_b] * 0.8 + tmpcg_rmse_square[i_b] * 0.2;
//            else
//               max_rmse[i_b] = tmpcg_rmse_square[i_b];
            //max_rmse[i_b] = tmpcg_rmse_square[i_b];
            //max_rmse[i_b] = max(adj_rmse[i_b]/sqrtf(pred_f_single), tmpcg_rmse_square[i_b]);
            //max_rmse[i_b] = max(adj_rmse[i_b], sqrtf(sum_square_vt[i_b] / (double)(valid_count[i_b] - SCCD_MAX_NUM_C + 1)));
            //max_rmse = max(adj_rmse[i_b] * sqrtf(kalman_ratio), tmpcg_rmse_square[i_b]);
            //v_dif_mag[i_b][i_conse] =  (clry[i_b][cur_i + i_conse] - pred_y_single) / sqrtf(pred_f_single) / max_rmse[i_b];
            //v_dif_mag[i_b][i_conse] =  (clry[i_b][cur_i + i_conse] - pred_y_single - tmpcg_rmse[i_b])  / max_rmse[i_b];
            //v_dif_mag[i_b][i_conse] =  (clry[i_b][cur_i + i_conse] - pred_y_single - sum_vt[i_b] / (double)valid_count[i_b])  / max_rmse[i_b];
            //v_dif_mag[i_b][i_conse] =  (clry[i_b][cur_i + i_conse] - pred_y_single)  / max_rmse[i_b];
            v_dif_mag[i_b][i_conse] =  (clry[i_b][cur_i + i_conse] - pred_y_single) / max_rmse[i_b];
            //v_dif_mag[i_b][i_conse] =  (clry[i_b][cur_i + i_conse] - pred_y_single) * kalman_ratio[i_b] / max_rmse[i_b];
//            if( clrx[cur_i] > 736499 - 1)
//            {
//                //printf("%d\n", clrx[cur_i]);
//                printf("%d conse, b%d: pred_f_single = %f; pre_y_single = %f; actual = %f; v = %f\n", i_conse, i_b + 1, pred_f_single, pred_y_single, clry[i_b][cur_i + i_conse], v_dif_mag[i_b][i_conse]);
//            }
//            if(cur_i == 203)
//            {
//                printf("%d conse, b%d: pred_f_single = %f; pre_y_single = %f; actual = %f; v = %f\n", i_conse, i_b + 1, pred_f_single, pred_y_single, clry[i_b][cur_i + i_conse], v_dif_mag[i_b][i_conse]);
//            }
            for (b = 0; b < n_focus_variable; b++)
            {
                if (i_b == focus_blist[b])
                {
                    v_dif[b][i_conse] =  v_dif_mag[i_b][i_conse];
                    v_dif_mag_standard = v_dif[b][i_conse] * v_dif[b][i_conse];
                    v_dif_mag_norm[i_conse] = v_dif_mag_norm[i_conse] + v_dif_mag_standard;
                    break;
                }
            }

//            for (b = 0; b < DEFAULT_N_FOCUS_VARIABLE; b++)
//            {
//                if (i_b == DEFAULT_FOCUS_BLIST[b])
//                {
//                    v_dif_mag_norm_ini = v_dif_mag_norm_ini+ v_dif_mag[i_b][i_conse] * v_dif_mag[i_b][i_conse];
//                    break;
//                }
//            }
        }

        /* eliminate gradual change */
//        if(v_dif_mag_norm[i_conse] < t_cg_outelier)
//        {
//            abrupt_change_flag = FALSE;
//            break;
//        }
//        if( clrx[cur_i] > 730597)
//        {
//            printf("mag is %f; prob for conse %d is %f\n", v_dif_mag_norm[i_conse], i_conse, Chi_Square_Distribution(v_dif_mag_norm[i_conse], n_focus_variable));
//        }
//        if( cur_i > 203)
//        {
//            printf("mag is %f; prob for conse %d is %f\n", v_dif_mag_norm[i_conse], i_conse, Chi_Square_Distribution(v_dif_mag_norm[i_conse], n_focus_variable));
//        }
        /* if it is smaller than S_TCG, no need to continue */
        // if(v_dif_mag_norm[i_conse] < t_cg_outelier)

//        if(cur_i == 185){
//            cur_i = 185;
//        }
        if(training_type == TRAINING_TYPE_PARAMETER)
        {
            double s_tcg = X2(n_focus_variable, probability_threshold);
            //double s_tcg = 11.07;
            if (i_conse == 0)
            {
                *change_magnitude = v_dif_mag_norm[i_conse] * s_tcg/ t_cg_adjust;
            }
            else
            {
               if(v_dif_mag_norm[i_conse] * s_tcg / t_cg_adjust < *change_magnitude)
                   *change_magnitude = v_dif_mag_norm[i_conse] * s_tcg / t_cg_adjust;
            }

            if(v_dif_mag_norm[i_conse] < t_cg_adjust)
            {
                abrupt_change_flag = FALSE;
            }
        }
        else if (training_type == TRAINING_TYPE_REGULAR){
            if(v_dif_mag_norm[i_conse] < t_cg_adjust)
            {
                abrupt_change_flag = FALSE;
                 break;
            }

        }

        // prob_conse = prob_conse + Chi_Square_Distribution(v_dif_mag_norm[i_conse], n_focus_variable);


        //abrupt_change_flag = FALSE;
    }

    /* back up test: testing medium value*/
//    if(abrupt_change_flag == TRUE)
//    {
//        for(i_b = 0; i_b < n_total_variable; i_b++)
//        {
//            quick_sort_double(v_dif_mag[i_b], 0, conse-1);
//            matlab_2d_double_median(v_dif_mag, i_b, conse,
//                                   &rec_cg[*num_curve].magnitude[i_b]);
//            for (b = 0; b < n_focus_variable; b++)
//            {
//                if (i_b == focus_blist[b])
//                {
//                    medium_v_dif[b] = rec_cg[*num_curve].magnitude[i_b];
//                }
//            }
//        }

//        for(i_b = 0; i_b < n_total_variable; i_b++)
//        {
//            for(i_conse = 0; i_conse < conse; i_conse++)
//                v_dif_gradiant[i_b][i_conse] = fabsf(v_dif_mag[i_b][i_conse] - rec_cg[*num_curve].magnitude[i_b]);

//            for (b = 0; b < n_focus_variable; b++)
//            {
//                if (i_b == focus_blist[b])
//                {
//                    quick_sort_double(v_dif_gradiant[i_b], 0, conse - 1);
//                    matlab_2d_double_median(v_dif_gradiant, i_b, conse,
//                                           &medium_v_dif_gradient[b]);
//                }
//            }
//        }


//        /* Gil's spread */
////        for(i_b = 0; i_b < n_total_variable; i_b++)
////        {
////            k = 0;
////            for(i_conse = 0; i_conse < conse - 1; i_conse++)
////                for(j_conse = 0; j_conse < conse - 1; j_conse++)
////                    if(i_conse > j_conse)
////                    {
////                        v_dif_gradiant[i_b][k] = fabsf(v_dif_mag[i_b][i_conse] - v_dif_mag[i_b][j_conse]);
////                        k = k + 1;
////                    }

////            for (b = 0; b < n_focus_variable; b++)
////            {
////                if (i_b == focus_blist[b])
////                {
////                    quick_sort_double(v_dif_gradiant[i_b], 0, (conse - 1) * (conse - 2) /2 - 1);
////                    matlab_2d_double_median(v_dif_gradiant, i_b, (conse - 1) * (conse - 2) /2,
////                                           &medium_v_dif_gradient[b]);
////                }
////            }
////        }

//        double medium_minus_gradient[n_focus_variable];
//        v_dif_mag_standard = 0;
//        for(i_b = 0; i_b < n_total_variable; i_b++)
//        {
//            for (b = 0; b < n_focus_variable; b++)
//            {
//                if (i_b == focus_blist[b])
//                {
//                    if(medium_v_dif[b] < 0)
//                    {
//                        medium_minus_gradient[b] = medium_v_dif[b] + medium_v_dif_gradient[b] / 3;
//                        if(medium_minus_gradient[b] > 0)
//                            medium_minus_gradient[b] = 0;
//                    }
//                    else
//                    {
//                        medium_minus_gradient[b] = medium_v_dif[b] - medium_v_dif_gradient[b] / 3;
//                        if(medium_minus_gradient[b] < 0)
//                            medium_minus_gradient[b] = 0;
//                    }
//                    v_dif_mag_standard = v_dif_mag_standard + medium_minus_gradient[b] * medium_minus_gradient[b];
//                }
//            }
//        }

//        if (v_dif_mag_standard < t_cg_adjust)
//            abrupt_change_flag = FALSE;
//    }

     if(abrupt_change_flag == TRUE)
     {
//        quick_sort_double(v_dif_mag_norm, 0, conse - 1);
//        matlab_2d_double_median(&v_dif_mag_norm, 0, conse,
//                               &mediam_mag);
        // prob_conse_avg = Chi_Square_Distribution(mediam_mag, n_focus_variable);
        // prob_conse_avg = prob_conse / conse;
//        if(prob_conse_avg < S_TCG)
//        {
//           abrupt_change_flag = FALSE;
//        }

//        if(sum_median_tmp < t_cg_adjust)
//        {
//            abrupt_change_flag = FALSE;
//        }



        for(i_b = 0; i_b < n_total_variable; i_b++)
        {
            quick_sort_double(v_dif_mag[i_b], 0, conse-1);
            matlab_2d_double_median(v_dif_mag, i_b, conse,
                                   &rec_cg[*num_curve].magnitude[i_b]);
            for (b = 0; b < n_focus_variable; b++)
            {
                if (i_b == focus_blist[b])
                {

//                        matlab_2d_double_median(v_dif, b, conse,
//                                               &medium_v_dif[b]);  // for cold disturbance extraction
                    medium_v_dif[b] = rec_cg[*num_curve].magnitude[i_b] ;
                    // sum_median_tmp = sum_median_tmp + medium_v_dif[b] * medium_v_dif[b];
                }
            }
        }

        double mean_angle_2;
        mean_angle_2 = angl_scatter_measure(medium_v_dif, v_dif, n_focus_variable, conse);
        //mean_angle  = MeanAngl(v_dif, n_focus_variable, conse);
        if (training_type == TRAINING_TYPE_PARAMETER){
            if (mean_angle_2 > NSIGN_sccd)
            {
                *change_magnitude = -9999.0; // give a extreme value meaning that this pixel is impossible to be potential break point
                abrupt_change_flag = FALSE;
            }
            else
            {
                if(-rec_cg[*num_curve].magnitude[2]+rec_cg[*num_curve].magnitude[3]-rec_cg[*num_curve].magnitude[4] > 0) // break not related to disturbance
                {
                    abrupt_change_flag = TRUE;
                    *change_magnitude = -9999.0; /* assin it to -9999 for not related to disturbance*/
                }
                else
                {
                    /* only this condition is candidate break point */
                    abrupt_change_flag = FALSE; // abrupt_change_flag is always FALSE for training stage
                }
            }
        }else if(training_type == TRAINING_TYPE_REGULAR){
            if (mean_angle_2 > NSIGN_sccd)
            //if (mean_angle > NSIGN)
            {
                abrupt_change_flag = FALSE;
            }
        }
    }




//    /* gradual change mode*/
//    else
//    {

//        if(grad_change_flag == TRUE)
//        {
//            /* 2nd test: medium mag test */
//            double sum_median_tmp = 0;
//            for(i_b = 0; i_b < n_total_variable; i_b++)
//            {
//                quick_sort_double(v_dif_mag[i_b], 0, conse-1);
//                for (b = 0; b < n_focus_variable; b++)
//                {
//                    if (i_b == focus_blist[b])
//                    {
//                        matlab_2d_double_median(v_dif_mag, i_b, conse,
//                                               &medium_v_dif[b]);
////                        for(i_conse = 0; i_conse < conse; i_conse++)
////                        {
////                             printf("%d conse, b%d: v_dif_mag = %f\n", i_conse, i_b + 1, v_dif_mag[i_b][i_conse]);
////                         }
//                        sum_median_tmp = sum_median_tmp + medium_v_dif[b] * medium_v_dif[b];
//                    }
//                }
//            }

//            if(sum_median_tmp < t_cg_gradual)
//                 grad_change_flag = FALSE;

//            /* 3rd test: change angle test */
//            if(grad_change_flag == TRUE)
//            {
//                mean_angle = angl_scatter_measure(medium_v_dif, v_dif, n_focus_variable, conse);
//                //mean_angle = MeanAngl(v_dif, n_focus_variable, conse);

//                if (mean_angle > NSIGN_sccd)
//                {
//                    grad_change_flag = FALSE;
//                }
//            }

//        }
//    }


    if (TRUE == abrupt_change_flag)
    {
        /********************************/
        /*  changed has been detected   */
        /********************************/

        /*fitting state variables and save to rec_cg)*/
        for(i_b = 0; i_b < n_total_variable; i_b++)
        {

//            status = level_ts_fit(level_state_records, rec_cg[*num_curve].t_start,
//                    clrx[cur_i] - 1, i_b, &c0, &c1, starting_date);
//            rec_cg[*num_curve].coefs[i_b][0] = (double)c0;
//            rec_cg[*num_curve].coefs[i_b][1] = (double)c1;

//            for (k = 2; k < SCCD_MAX_NUM_C; k++) /* slope is the second element of coefs, so start from 1*/
//            {
//                rec_cg[*num_curve].coefs[i_b][k] = fit_cft1[i_b][k];
//            }

            status = auto_ts_fit(clrx, clry, i_b, i_b, i_start, cur_i - 1, SCCD_MAX_NUM_C,
                                fit_cft_average, &rmse_tmp, rec_v_dif_copy);

//            for (k = 0; k < cur_i - i_start; k++) /* slope is the second element of coefs, so start from 1*/
//            {
//                printf("%f\n", rec_v_dif_copy[i_b][k]);
//            }
            // printf("auto_ts_fit2 finished \n", i);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Calling auto_ts_fit for change detection with "
                     "enough observations\n", FUNC_NAME, FAILURE);
            }


            for (k = 0; k < SCCD_MAX_NUM_C; k++)
            {
                /**********************************/
                /*                                */
                /* Record fitted coefficients.    */
                /*                                */
                /**********************************/

                rec_cg[*num_curve].coefs[i_b][k] = fit_cft_average[i_b][k];
            }

            rec_cg[*num_curve].rmse[i_b] = sqrtf(sum_square_vt[i_b] / (valid_count[i_b]));

//            for(i_conse = 0; i_conse < conse; i_conse++)
//            {
//                printf("v_dif_mag for i_conse = %d for band = %d is %f \n", i_conse, i_b + 1,  v_dif_mag[i_b][i_conse]);
//            }

            /* record regular magnitudes */
            // quick_sort_double(v_dif_mag[i_b], 0, conse-1);

        } //for(i_b = 0; i_b < n_total_variable; i_b++)

        /* record break and end */
        rec_cg[*num_curve].t_break = clrx[cur_i];
        *prev_i_break = cur_i;
        rec_cg[*num_curve].t_end = clrx[cur_i-1];

        rec_cg[*num_curve].num_obs = cur_i -i_start;
        rec_cg[*num_curve].t_confirmed = clrx[cur_i + conse - 1];
        rec_cg[*num_curve].change_prob = 100;

        if(abrupt_change_flag == TRUE)
            rec_cg[*num_curve].category = 0;
        else
            rec_cg[*num_curve].category = 70;

//        for(i_b = 0; i_b < n_total_variable; i_b++)
//        {
//            rec_cg[*num_curve].obs_disturb[i_b] = (double)instance[i_b].H;
//            for (k = 0; k < SCCD_MAX_NUM_C - 1; k++)
//            {
//                if (k < instance[i_b].m)
//                    rec_cg[*num_curve].state_disturb[i_b][k] = (double)gsl_matrix_get(instance[i_b].Q, k, k);
//                else
//                    rec_cg[*num_curve].state_disturb[i_b][k] = 0;
//            }
//        }

        /**********************************************/
        /*                                            */
        /* Identified and move on for the next        */
        /* functional curve.                          */
        /*                                            */
        /**********************************************/

        *num_curve = *num_curve + 1;

        if (*num_curve >= NUM_FC)
        {
            /******************************************/
            /*                                        */
            /* Reallocate memory for rec_cg.          */
            /*                                        */
            /******************************************/

            rec_cg = realloc(rec_cg, (*num_curve + 1) * sizeof(Output_t_sccd));
            if (rec_cg == NULL)
            {
                RETURN_ERROR("ERROR allocating rec_cg memory",
                                     FUNC_NAME, FAILURE);
            }
        }

            /* don't update a and p to cur_i +1! leave it to the next curve initialization (ssm_initialize)*/
            /* make regular KF updates and update a and p until clrx[cur_i+1]*/
//        for(i_b = 0; i_b < n_total_variable; i_b++)
//        {
//            KF_ts_filter_regular(instance, clrx, clry[i_b], vec_next_P[i_b], vec_next_a[i_b],  cur_i,
//                                                i_b, rmse_records, level_state_records,
//                                                annual_state_records, semi_state_records,
//                                                starting_date, b_fastmode, fit_cft,
//                                                valid_count, sum_square_smooth_Q,
//                                                sum_square_smooth_H, sum_smooth_Q, sum_smooth_H);
//        }

         RETURN_VALUE = CHANGEDETECTED;
    } // if (TRUE == change_flag)
    /* outlier detected */
    //else if(first_smooth_obs_noise > S_T_MAX_CG)
    // else if(v_dif_mag_norm[0]  > S_T_MAX_CG)
    else if(v_dif_mag_norm[0]  > t_cg_outelier)
    {
        /**********************************************/
        /*                                            */
        /*    still need to update a and p,           */
        /*    but treat cur_i as missing value        */
        /*                                            */
        /**********************************************/
        // printf("%d, %d\n", cur_i, clrx[cur_i]);
        for(i_b = 0; i_b < n_total_variable; i_b++)
        {
            KF_ts_filter_falsechange(&instance[i_b], clrx, vec_next_P[i_b] , vec_next_a[i_b], cur_i,
                                 i_b, level_state_records, annual_state_records,
                                 semi_state_records, third_state_records, starting_date, b_fastmode);
        }


        /**********************************************/
        /*                                            */
        /*    Remove noise.                           */
        /*                                            */
        /**********************************************/
        for (m = cur_i; m < *end -1; m++)
        {
            clrx[m] = clrx[m+1];
            for (i_b = 0; i_b < n_total_variable; i_b++)
                clry[i_b][m] = clry[i_b][m+1];
        }

        RETURN_VALUE = FALSECHANGE;
    }
    else
    {
        /* using vec_next_P[i_b] , vec_next_a[i_b], fit_cft[i_b] because this is real update*/
        for(i_b = 0; i_b < n_total_variable; i_b++)
        {

            for(i = 0; i < SCCD_MAX_NUM_C; i++)
            {
                fit_cft4[i_b][i] = fit_cft3[i_b][i];
                fit_cft3[i_b][i] = fit_cft2[i_b][i];
                fit_cft2[i_b][i] = fit_cft1[i_b][i];
            }


            KF_ts_filter_regular(&instance[i_b], clrx, clry[i_b], vec_next_P[i_b], vec_next_a[i_b],  cur_i,
                                 i_b, rmse_records, level_state_records,
                                 annual_state_records, semi_state_records, third_state_records,
                                 starting_date, b_fastmode, fit_cft1[i_b],
                                 valid_count, sum_square_smooth_Q,
                                 sum_square_smooth_H, sum_smooth_Q, sum_smooth_H, sum_square_vt, sum_vt,
                                 sum_kalman_coef, FALSE, temporal_rmse_square_sum, temporal_rmse_sum,
                                 temporal_rmse_count);

        }

        *clrx_record4 = *clrx_record3;
        *clrx_record3 = *clrx_record2;
        *clrx_record2 = *clrx_record1;
        *clrx_record1 = clrx[cur_i];

        RETURN_VALUE = REGULAREND;
    }


    status = free_2d_array((void **)v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: v_dif\n", FUNC_NAME,
                      FAILURE);
    }


    status = free_2d_array ((void **) v_dif_mag);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: v_dif_mag\n", FUNC_NAME,
                      FAILURE);
    }

//    status = free_2d_array ((void **) v_dif_gradiant);
//    if (status != SUCCESS)
//    {
//        RETURN_ERROR ("Freeing memory: v_dif_gradiant\n", FUNC_NAME,
//                      FAILURE);
//    }

    for(i = 0; i < n_total_variable; i ++)
    {
       gsl_vector_free(copy_vec_next_a[i]);
       gsl_matrix_free(copy_vec_next_P[i]);
    }

    status = free_2d_array((void *)copy_vec_next_P);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Free memory: copy_vec_next_P memory", FUNC_NAME, FAILURE);
    }

    status = free_2d_array((void *)copy_vec_next_a);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Free memory: copy_vec_next_a memory", FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **) fit_cft_average);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fit_cft_average\n", FUNC_NAME,
                      FAILURE);
    }

    free(pred_y);
    pred_y = NULL;
    free(pred_y_f);
    pred_y_f = NULL;


    free(v_dif_mag_norm);
    v_dif_mag_norm = NULL;

    status = free_2d_array ((void **) rec_v_dif_copy);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: \n",
                      FUNC_NAME, FAILURE);
    }

    free(tmpcg_rmse_square);
    free(tmpcg_rmse);
    free(max_rmse);
//    status = free_2d_array ((void **) fit_cft);
//    if (status != SUCCESS)
//    {
//        RETURN_ERROR ("Freeing memory: fit_cft\n", FUNC_NAME,
//                      FAILURE);
//    };

    /**********************************************/
    /*                                            */
    /* Free allocated memories.                   */
    /*                                            */
    /**********************************************/
    free(d_yr);
    free(medium_v_dif);

    return (RETURN_VALUE);

}

/************************************************************************
FUNCTION: step2_chow_test

PURPOSE:
Step 2 of S-CCD: doing chow test
RETURN VALUE:
Type = int (SUCCESS OR FAILURE)

Programmer: Su Ye
**************************************************************************/

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
)
{
    int i_b, b, i;
    double rmse_first;
    double rmse_second;
    double rmse_combined;
    double rmse_first_square;
    double rmse_second_square;
    double rmse_combined_square;
    int status;
    double **fit_cft;
    double **temp_v_dif_first;
    double **temp_v_dif_second;
    double **temp_v_dif_combined;
    char FUNC_NAME[] = "step2_chow_test";
    int nums;
    int total_num = first_seg_end - first_seg_start + second_seg_end - second_seg_start + 2;
    float **clry_combined;
    int *clrx_combined;
    double f[n_focus_variable];
    int chow_pass_count = 0;
    double f_stastic;
    double f_sum = 0;
    int df = SCCD_MAX_NUM_C;

    temp_v_dif_first = (double **)allocate_2d_array(TOTAL_IMAGE_BANDS, first_seg_end - first_seg_start + 1,
                                     sizeof (double));
    if (temp_v_dif_first == NULL)
    {
        RETURN_ERROR ("Allocating temp_v_dif_first memory",FUNC_NAME, FAILURE);
    }

    temp_v_dif_second = (double **)allocate_2d_array(TOTAL_IMAGE_BANDS, second_seg_end - second_seg_start + 1,
                                     sizeof (double));
    if (temp_v_dif_second == NULL)
    {
        RETURN_ERROR ("Allocating temp_v_dif_first memory",FUNC_NAME, FAILURE);
    }

    temp_v_dif_combined = (double **)allocate_2d_array(TOTAL_IMAGE_BANDS, total_num,
                                     sizeof (double));
    if (temp_v_dif_combined == NULL)
    {
        RETURN_ERROR ("Allocating temp_v_dif_combined memory",FUNC_NAME, FAILURE);
    }

    fit_cft = (double **) allocate_2d_array (TOTAL_IMAGE_BANDS, SCCD_MAX_NUM_C, sizeof (double));
    if (fit_cft == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft memory", FUNC_NAME, FAILURE);
    }

    clry_combined = (float **) allocate_2d_array (n_focus_variable, total_num, sizeof (float));
    if (clry_combined == NULL)
    {
        RETURN_ERROR ("Allocating clry_combined memory", FUNC_NAME, FAILURE);
    }

    clrx_combined = (int *)malloc(total_num * sizeof(int));
    if (clrx == NULL)
    {
        RETURN_ERROR ("Allocating clrx memory", FUNC_NAME, FAILURE);
    }

    for (b = 0; b < n_focus_variable; b++)
    {
        for(i = first_seg_start; i < first_seg_end + 1; i++)
        {
            clry_combined[b][i - first_seg_start] = clry[focus_blist[b]][i];
            clrx_combined[i - first_seg_start] = clrx[i];
        }

        for(i = second_seg_start; i < second_seg_end + 1; i++)
        {
            clry_combined[b][first_seg_end + 1 - first_seg_start + i - second_seg_start] = clry[focus_blist[b]][i];
            clrx_combined[first_seg_end + 1 - first_seg_start + i - second_seg_start] = clrx[i];
        }


    }

    f_stastic = F(df, total_num - df * 2, f_prob);
    //f_stastic = 0;
    for (b = 0; b < n_focus_variable; b++)
    {
        status = auto_ts_fit(clrx, clry, focus_blist[b], focus_blist[b], first_seg_start, first_seg_end,
                             df, fit_cft, &rmse_first, temp_v_dif_first);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Calling auto_ts_fit for chow test\n",
                          FUNC_NAME, FAILURE);
        }

        nums = first_seg_end - first_seg_start + 1;
        rmse_first_square = rmse_first * rmse_first * ( nums - df);

        status = auto_ts_fit(clrx, clry, focus_blist[b], focus_blist[b], second_seg_start, second_seg_end,
                             df, fit_cft, &rmse_second, temp_v_dif_second);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Calling auto_ts_fit for chow test\n",
                          FUNC_NAME, FAILURE);
        }

        nums = second_seg_end - second_seg_start + 1;
        rmse_second_square = rmse_second * rmse_second * (nums - df);

        status = auto_ts_fit(clrx_combined, clry_combined, b, b, 0, total_num - 1,
                             df, fit_cft, &rmse_combined, temp_v_dif_combined);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Calling auto_ts_fit for chow test\n",
                          FUNC_NAME, FAILURE);
        }
        rmse_combined_square = rmse_combined * rmse_combined * (total_num - df);

        f[b] = ((rmse_combined_square - rmse_first_square - rmse_second_square) / (double)df)
                / ((rmse_first_square + rmse_second_square) / (double)(total_num - df * 2));

        if(f[b] > f_stastic)
             chow_pass_count = chow_pass_count + 1;
        f_sum = f_sum + f[b];
    }

    status = free_2d_array ((void **)fit_cft);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fit_cft\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **)temp_v_dif_first);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: temp_v_dif_first\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **)temp_v_dif_second);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: temp_v_dif_second\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **)temp_v_dif_combined);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: temp_v_dif_combined\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **)clry_combined);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: clry_combined\n", FUNC_NAME, FAILURE);
    }


    free(clrx_combined);

    if(chow_pass_count < NUM_LASSO_BANDS / 2 + 1)
        return FALSE;
    else
        return TRUE;

//    if(f_sum / (double)NUM_LASSO_BANDS < f_stastic)
//        return FALSE;
//    else
//        return TRUE;

}
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
)
{
    int i, j, m, k;
    int i_b;
    int id_last = cur_i;
    int b, status;
    char FUNC_NAME[] = "step3_processingend";
    double **v_dif_mag_standard;           /* normalized magnitude, =v_diff in ccd.c */
    double **v_dif_mag;
    int *bl_ids, *ids;
    int *rm_ids;
    int i_span = cur_i - i_start + 1;
    int k_new;
    int rm_ids_len;
    double *pred_y;
    double *pred_y_f;
    double tmp;
    //double rmse_sqrt_tmp;
    double tmp_d;
    double *rmse;                     /* Root Mean Squared Error array.        */
    double **temp_v_dif;
    double w, w2;
    // double tmp_q;
    double c0, c1;
    int category;
    int num_c = 8;
    int update_num_c;
    double rmse_tmp;
    int center_date;
    double weight1;
    double weight2;
    double weight3;
    double weight4;
    double **fit_cft_average;
    int n_rmse;
    int center_dyr;
    double* tmpcg_rmse_square; /* to temporarily change RMSE          */
    double* tmpcg_rmse;
    int i_conse;
    double pred_y_single;
    double pred_f_single;
    double** v_diff;
    double* medium_v_dif;
    double* v_dif_mag_norm;
    double* max_rmse;
    double min_v_dif_mag_norm;
    double t_cg_adjust;
    double** v_dif_tmp;
    double** v_dif;
    double mean_angle_2;


    w = TWO_PI / AVE_DAYS_IN_A_YEAR;
    w2 = 2.0 * w;
    int conse_end = *end - cur_i;
    if (conse_end <= 0)
        conse_end = 0;
    if(conse_end < conse)
        t_cg_adjust = X2(n_focus_variable, probability_threshold);
    else
        t_cg_adjust = X2(n_focus_variable, 1 - pow(1 - probability_threshold, (double)conse / (double)conse_end));
    //t_cg_adjust = 11.07;

    v_dif_mag_standard = (double **) allocate_2d_array (n_total_variable, conse_end, sizeof (double));
    if (v_dif_mag_standard == NULL)
    {
        RETURN_ERROR ("Allocating v_dif_mag_standard memory",
                      FUNC_NAME, FAILURE);
    }

    v_dif_mag = (double **) allocate_2d_array(n_total_variable, conse_end,
                sizeof (double));
    if (v_dif_mag == NULL)
    {
        RETURN_ERROR ("Allocating v_dif_mag memory",
                                 FUNC_NAME, FAILURE);
    }

    bl_ids = (int *)malloc((*end) * sizeof(int));
    if (bl_ids == NULL)
    {
        RETURN_ERROR("ERROR allocating bl_ids memory", FUNC_NAME, FAILURE);
    }

    ids = (int *)malloc((*end) * sizeof(int));
    if (ids == NULL)
    {
        RETURN_ERROR("ERROR allocating ids memory", FUNC_NAME, FAILURE);
    }

    rm_ids = (int *)malloc((*end) * sizeof(int));
    if (rm_ids == NULL)
    {
        RETURN_ERROR("ERROR allocating rm_ids memory", FUNC_NAME, FAILURE);
    }

    pred_y_f = (double *)malloc(conse * sizeof(double));
    if (pred_y_f == NULL)
    {
        RETURN_ERROR ("Allocating pred_y_f memory", FUNC_NAME, FAILURE);
    }

    pred_y = (double *)malloc(conse * sizeof(double));
    if (pred_y == NULL)
    {
        RETURN_ERROR ("Allocating pred_y memory", FUNC_NAME, FAILURE);
    }


    rmse = (double *)malloc(n_total_variable * sizeof(double));
    if (rmse == NULL)
    {
        RETURN_ERROR ("Allocating rmse memory", FUNC_NAME, FAILURE);
    }

    temp_v_dif = (double **)allocate_2d_array(n_total_variable, *end,
                                     sizeof (double));
    if (temp_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating temp_v_dif memory",FUNC_NAME, FAILURE);
    }


    fit_cft_average = (double **) allocate_2d_array (n_total_variable, LASSO_COEFFS, sizeof (double));
    if (fit_cft_average == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft_average memory", FUNC_NAME, FAILURE);
    }

    tmpcg_rmse_square= (double*)malloc(n_total_variable * sizeof(double)); /* to temporarily change RMSE          */
    if (tmpcg_rmse_square == NULL)
    {
        RETURN_ERROR ("Allocating tmpcg_rmse_square memory",
                     FUNC_NAME, FAILURE);
    }

    tmpcg_rmse = (double*)malloc(n_total_variable * sizeof(double));
    if (tmpcg_rmse == NULL)
    {
        RETURN_ERROR ("Allocating tmpcg_rmse memory",
                     FUNC_NAME, FAILURE);
    }


    max_rmse = (double*)malloc(n_total_variable * sizeof(double));
    if (max_rmse == NULL)
    {
        RETURN_ERROR ("Allocating max_rmse memory",
                     FUNC_NAME, FAILURE);
    }

    v_dif_mag_norm = (double *)malloc(conse_end * sizeof(double));
    if (v_dif_mag_norm == NULL)
    {
        RETURN_ERROR ("Allocating v_dif_mag_norm memory", FUNC_NAME, FAILURE);
    }

    v_diff =(double **) allocate_2d_array(n_focus_variable, conse_end, sizeof (double));
    if (v_diff == NULL)
    {
        RETURN_ERROR ("Allocating v_diff memory", FUNC_NAME, FAILURE);
    }

    if (bl_train_complete == 1)
    {

        center_date = (int)(clrx[*end - 1] + clrx[cur_i]) / 2;
        weight1 = (double)(center_date - *clrx_record1);
        weight2 = (double)(center_date - *clrx_record2);
        weight3 = (double)(center_date - *clrx_record3);
        weight4 = (double)(center_date - *clrx_record4);

        for(i_b = 0; i_b < n_total_variable; i_b++)
        {
            for(i = 0; i < SCCD_MAX_NUM_C; i++)
                fit_cft_average[i_b][i] = fit_cft1[i_b][i] * weight1 / (weight1 + weight2 + weight3 + weight4)
                        + fit_cft2[i_b][i] * weight2 / (weight1 + weight2 + weight3 + weight4)
                        + fit_cft3[i_b][i] * weight3 / (weight1 + weight2 + weight3 + weight4)
                        + fit_cft4[i_b][i] * weight4 / (weight1 + weight2 + weight3 + weight4);
        }

        n_rmse = valid_count[0];

        if(n_rmse > LASSO_COEFFS * N_TIMES)
        {

            i = 1; /* step length */
            center_dyr = center_date - (int)((int)(center_date / NUM_YEARS) * NUM_YEARS);
            center_dyr = center_dyr / TEMPORAL_RMSE_BIN;
            for(i_b = 0; i_b < n_total_variable; i_b++)
            {
                tmpcg_rmse_square[i_b] = temporal_rmse_square_sum[i_b][center_dyr];
                tmpcg_rmse[i_b] = temporal_rmse_sum[i_b][center_dyr];
            }

            n_rmse = temporal_rmse_count[center_dyr];

            //max_num = max(LASSO_COEFFS * N_TIMES, valid_count[0] / 4);
            while(n_rmse <= LASSO_COEFFS * N_TIMES)
            {
                if(center_dyr - i < 0)
                {
                    if(temporal_rmse_count[TEMPORAL_RMSE_BIN_NUM + center_dyr - i] > 0)
                    {
                        n_rmse = n_rmse + temporal_rmse_count[TEMPORAL_RMSE_BIN_NUM + center_dyr - i];
                        for(i_b = 0; i_b < n_total_variable; i_b++)
                        {
                            tmpcg_rmse_square[i_b] = tmpcg_rmse_square[i_b] + temporal_rmse_square_sum[i_b][TEMPORAL_RMSE_BIN_NUM + center_dyr - i];
                            tmpcg_rmse[i_b] = tmpcg_rmse[i_b] + temporal_rmse_sum[i_b][TEMPORAL_RMSE_BIN_NUM + center_dyr - i];
                        }
                    }
                }
                else
                {
                    if(temporal_rmse_count[center_dyr - i] > 0)
                    {
                        n_rmse = n_rmse + temporal_rmse_count[center_dyr - i];
                        for(i_b = 0; i_b < n_total_variable; i_b++)
                        {
                            tmpcg_rmse_square[i_b] = tmpcg_rmse_square[i_b] + temporal_rmse_square_sum[i_b][center_dyr - i];
                            tmpcg_rmse[i_b] = tmpcg_rmse[i_b] + temporal_rmse_sum[i_b][center_dyr - i];
                        }
                    }
                }

                if(center_dyr + i > TEMPORAL_RMSE_BIN_NUM - 1)
                {
                    if(temporal_rmse_count[center_dyr + i - TEMPORAL_RMSE_BIN_NUM] > 0)
                    {
                        n_rmse = n_rmse + temporal_rmse_count[center_dyr + i - TEMPORAL_RMSE_BIN_NUM];
                        for(i_b = 0; i_b < n_total_variable; i_b++)
                        {
                            tmpcg_rmse_square[i_b] = tmpcg_rmse_square[i_b] + temporal_rmse_square_sum[i_b][center_dyr + i - TEMPORAL_RMSE_BIN_NUM];
                            tmpcg_rmse[i_b] = tmpcg_rmse[i_b] + temporal_rmse_sum[i_b][center_dyr + i - TEMPORAL_RMSE_BIN_NUM];
                        }
                    }
                }
                else
                {
                    if(temporal_rmse_count[center_dyr + i] > 0)
                    {
                        n_rmse = n_rmse + temporal_rmse_count[center_dyr + i];
                        for(i_b = 0; i_b < n_total_variable; i_b++)
                        {
                            tmpcg_rmse_square[i_b] = tmpcg_rmse_square[i_b] + temporal_rmse_square_sum[i_b][center_dyr + i];
                            tmpcg_rmse[i_b] = tmpcg_rmse[i_b] + temporal_rmse_sum[i_b][center_dyr + i];
                        }
                    }
                }
                i++;
            }

            for (i_b = 0; i_b < n_total_variable; i_b++)
            {
                tmpcg_rmse_square[i_b] = sqrtf(tmpcg_rmse_square[i_b]) / sqrtf((double)n_rmse);
                tmpcg_rmse[i_b] = tmpcg_rmse[i_b] / (double)n_rmse;
            }

        }
        else
        {
            for (i_b = 0; i_b < n_total_variable; i_b++)
            {
                tmpcg_rmse_square[i_b] = sqrtf(sum_square_vt[i_b]) / sqrtf((double)n_rmse);
                tmpcg_rmse[i_b] = sum_vt[i_b] / (double)n_rmse;
    //            tmpcg_rmse_square[i_b] = sqrtf(sum_square_vt[i_b]) / sqrtf((double)n_rmse - SCCD_NUM_C);
    //            tmpcg_rmse[i_b] = sum_vt[i_b] / (double)(n_rmse - SCCD_NUM_C);
            }
        }


        /**********************************************************/
        /*                                                        */
        /* If no break, find at the end of the time series,       */
        /* define probability of change based on conse.           */
        /*                                                        */
        /**********************************************************/

        rec_cg[*num_curve].t_break = 0;
        rec_cg[*num_curve].num_obs = *end -1 - i_start + 1;
        rec_cg[*num_curve].t_end = clrx[*end -1];
        rec_cg[*num_curve].category = 0;
        rec_cg[*num_curve].t_start = clrx[i_start];
        rec_cg[*num_curve].t_confirmed = 0;

        for (k = 0; k < conse_end; k++)
        {
            v_dif_mag_norm[k] = 0;
        }
        /* update ssm variable for the end*/
        for(i_b = 0; i_b < n_total_variable; i_b++)
        {
            max_rmse[i_b] = max(adj_rmse[i_b], tmpcg_rmse_square[i_b]);
            for (i = 0; i < conse_end; i++)
            {
//                gsl_blas_ddot(instance->Z, vec_next_a[i_b], &yt_tmp);
//                v_dif_mag[i_b][i] =  (clry[i_b][cur_i + i] - yt_tmp) / max_rmse[i_b];
//                for (b = 0; b < n_focus_variable; b++)
//                {
//                    if (i_b == focus_blist[b])
//                    {
//                        v_diff[b][i] =  v_dif_mag[i_b][i];
//                        v_dif_mag_norm[i] = v_dif_mag_norm[i] + v_dif_mag[i_b][i] * v_dif_mag[i_b][i];
//                    }
//                }
//                KF_ts_filter_regular(&instance[i_b], clrx, clry[i_b], vec_next_P[i_b], vec_next_a[i_b],
//                                     cur_i + i, i_b, rmse_records, level_state_records,
//                                     annual_state_records, semi_state_records, third_state_records,
//                                     starting_date, b_fastmode, fit_cft_average[i_b],
//                                     valid_count, sum_square_smooth_Q,
//                                     sum_square_smooth_H, sum_smooth_Q, sum_smooth_H, sum_square_vt, sum_vt,
//                                     sum_kalman_coef, FALSE, temporal_rmse_square_sum, temporal_rmse_sum,
//                                     temporal_rmse_count);
                KF_ts_predict_conse(&instance[i_b], clrx, vec_next_P[i_b] , vec_next_a[i_b], cur_i + i,
                                    cur_i + i, &pred_y_single, &pred_f_single, b_fastmode, fit_cft_average[i_b], FALSE);
                v_dif_mag[i_b][i] =  (clry[i_b][cur_i + i] - pred_y_single) / max_rmse[i_b];
                for (b = 0; b < n_focus_variable; b++)
                {
                    if (i_b == focus_blist[b])
                    {
                        v_diff[b][i] =  v_dif_mag[i_b][i];
                        v_dif_mag_norm[i] = v_dif_mag_norm[i] + v_dif_mag[i_b][i] * v_dif_mag[i_b][i];
                    }
                }

//                KF_ts_filter_regular(&instance[i_b], clrx, clry[i_b], vec_next_P[i_b], vec_next_a[i_b],
//                                     cur_i + i, i_b, rmse_records, level_state_records,
//                                     annual_state_records, semi_state_records, third_state_records,
//                                     starting_date, b_fastmode, fit_cft_average[i_b],
//                                     valid_count, sum_square_smooth_Q,
//                                     sum_square_smooth_H, sum_smooth_Q, sum_smooth_H, sum_square_vt, sum_vt,
//                                     sum_kalman_coef, FALSE, temporal_rmse_square_sum, temporal_rmse_sum,
//                                     temporal_rmse_count);

            }

            // can't include the last obs (*end - 1)
            // needs to start from the previous break
            if (FALSE== b_fastmode)
            {

                /* for the last point clrx[*end - 1]*/
                level_state_records[i_b][clrx[*end - 1]-starting_date] = (double)gsl_vector_get(vec_next_a[i_b], 0);
                annual_state_records[i_b][clrx[*end - 1]-starting_date] = (double)gsl_vector_get(vec_next_a[i_b], 1);
                semi_state_records[i_b][clrx[*end - 1]-starting_date] = (double)gsl_vector_get(vec_next_a[i_b], 3);
                if (instance[i_b].m == 7)
                    third_state_records[i_b][clrx[*end - 1]-starting_date] = (double)gsl_vector_get(vec_next_a[i_b], 5);
                else
                    third_state_records[i_b][clrx[*end - 1]-starting_date] = -9999;
            }

//            rec_cg[*num_curve].obs_disturb[i_b] = (double)instance[i_b].H;
//            for (k = 0; k < SCCD_MAX_NUM_C - 1; k++)
//            {
//                if (k < instance[i_b].m)
//                    rec_cg[*num_curve].state_disturb[i_b][k] = (double)gsl_matrix_get(instance[i_b].Q, k, k);
//                else
//                    rec_cg[*num_curve].state_disturb[i_b][k] = 0;
//            }

            update_cft(i_span, N_TIMES, MIN_NUM_C, MID_NUM_C, SCCD_MAX_NUM_C,
                      num_c, &update_num_c);

            status = auto_ts_fit(clrx, clry, i_b, i_b, i_start, *end - 1, update_num_c,
                                fit_cft1, &rmse_tmp, temp_v_dif);
            // printf("auto_ts_fit2 finished \n", i);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Calling auto_ts_fit for change detection with "
                     "enough observations\n", FUNC_NAME, FAILURE);
            }

            for (k = 0; k < SCCD_MAX_NUM_C; k++)
            {
                /**********************************/
                /*                                */
                /* Record fitted coefficients.    */
                /*                                */
                /**********************************/

                rec_cg[*num_curve].coefs[i_b][k] = fit_cft1[i_b][k];
            }

//            array_1d_mean(rmse_records[i_b], i_start, *end - 1, &tmp);
//            rec_cg[*num_curve].rmse[i_b] = tmp;
            rec_cg[*num_curve].rmse[i_b] = sqrtf(sum_square_vt[i_b] / (valid_count[i_b]));
            rec_cg[*num_curve].magnitude[i_b] = 0;

//            if (FALSE == b_fastmode)
//            {
//                for (k = clrx[cur_i]; k < clrx[*end - 1]+1; k++)
//                {

//                    level_state_records[i_b][k-starting_date] = NA_VALUE;
//                    annual_state_records[i_b][k-starting_date] = NA_VALUE;
//                    semi_state_records[i_b][k-starting_date] = NA_VALUE;
//                    third_state_records[i_b][k-starting_date] = NA_VALUE;
//                }


//            }

        }



        for (i_conse = conse_end - 1; i_conse >= 0; i_conse--)
        {
            v_dif_tmp =(double **) allocate_2d_array(n_focus_variable, conse_end - i_conse, sizeof (double));

            v_dif =(double **) allocate_2d_array(n_focus_variable, conse_end - i_conse, sizeof (double));

            medium_v_dif = (double *)malloc(n_focus_variable * sizeof(double));


            for(j = 0; j < conse_end - i_conse; j++){
                for (b = 0; b < n_focus_variable; b++)
                {
                        v_dif_tmp[b][j] = v_diff[b][conse_end - 1 - j];
                        v_dif[b][j] = v_diff[b][conse_end - 1 - j];
                }
            }



            for (b = 0; b < n_focus_variable; b++)
            {


//                        matlab_2d_double_median(v_dif, b, conse,
//                                               &medium_v_dif[b]);  // for cold disturbance extraction
                    quick_sort_double(v_dif_tmp[b], 0, conse_end - i_conse - 1);
                    matlab_2d_double_median(v_dif_tmp, b, conse_end - i_conse, &medium_v_dif[b]);
                    // sum_median_tmp = sum_median_tmp + medium_v_dif[b] * medium_v_dif[b];

            }

            mean_angle_2 = angl_scatter_measure(medium_v_dif, v_dif, n_focus_variable, conse_end - i_conse);

            min_v_dif_mag_norm = 9999999;
            for(k = 0 ; k < conse_end - i_conse; k ++){
                if(v_dif_mag_norm[conse_end - 1 - k] < min_v_dif_mag_norm){
                    //printf("%f\n", v_dif_mag_norm[conse_end - 1 - k]);
                    min_v_dif_mag_norm = v_dif_mag_norm[conse_end - 1 - k];
                }
            }
            //double tt = vec_mag[i_conse];


            if ((min_v_dif_mag_norm  <= t_cg_adjust)||(mean_angle_2 >= NSIGN_sccd))
            {
                /**************************************************/
                /*                                                */
                /* The last stable ID.                            */
                /*                                                */
                /**************************************************/

                id_last = cur_i + i_conse + 1;
                status = free_2d_array ((void **) v_dif_tmp);
                if (status != SUCCESS)
                {
                    RETURN_ERROR ("Freeing memory: v_dif_tmp\n",
                                  FUNC_NAME, FAILURE);
                }

                status = free_2d_array ((void **) v_dif);
                if (status != SUCCESS)
                {
                    RETURN_ERROR ("Freeing memory: v_dif\n",
                                  FUNC_NAME, FAILURE);
                }

                free(medium_v_dif);
                break;
            }
            status = free_2d_array ((void **) v_dif_tmp);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: v_dif_tmp\n",
                              FUNC_NAME, FAILURE);
            }

            status = free_2d_array ((void **) v_dif);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: v_dif\n",
                              FUNC_NAME, FAILURE);
            }

            free(medium_v_dif);

        }

        double p = (double)(clrx[*end - 1] - clrx[id_last - 1]) / (double)min_days_conse;
        if (p > 1)
            p = 1;
        rec_cg[*num_curve].change_prob = (int)(p * 100);


    } // (bl_train == 1)
    else if (bl_train_complete == 0)
    {
        /**********************************************************/
        /*                                                        */
        /* If break found close to the end of the time series,    */
        /* use [adj_conse,MIN_NUM_C*N_TIMES+adj_conse) to fit curve. */
        /*                                                        */
        /* Update i_start.                                        */
        /*                                                        */
        /**********************************************************/
        if (*num_curve == 0)
        {
            /******************************************************/
            /*                                                    */
            /* First curve.                                       */
            /*                                                    */
            /******************************************************/

            i_start = 0;

        }
        else
        {

             i_start = prev_i_break + 1;
        }

        for (m = 0; m < *end - 1; m++)
        {
            bl_ids[m] = 0;
        }

        if (*end - i_start >= CONSE_END)
        {
            /******************************************************/
            /*                                                    */
            /* Multitemporal cloud mask.                          */
            /*                                                    */
            /******************************************************/

            status = auto_mask(clrx, clry, i_start, *end-1,
                           (double)(clrx[*end-1]-clrx[i_start]) / NUM_YEARS,
                           adj_rmse[1], adj_rmse[4], SCCD_T_CONST, bl_ids);
            if (status != SUCCESS)
                RETURN_ERROR("ERROR calling auto_mask at the end of time series",
                              FUNC_NAME, FAILURE);


            /******************************************************/
            /*                                                    */
            /* IDs to be removed.                                 */
            /*                                                    */
            /******************************************************/

            for (k = i_start; k < *end; k++)
            {
                ids[k-i_start] = k;
            }
            m= 0;
            i_span = 0;
            for (k = 0; k < *end-i_start; k++)
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
            for (k = 0, k_new=0; k < *end; k++)
            {
                if (m < rm_ids_len && k == rm_ids[m])
                {
                    m++;
                    continue;
                }
                clrx[k_new] = clrx[k];
                for (i_b = 0; i_b < n_total_variable; i_b++)
                    clry[i_b][k_new] = clry[i_b][k];
                k_new++;
            }
            *end = k_new;
        }

        /********************************************/
        /*    applying an inefficient procedure     */
        /********************************************/
        if (*end - i_start >= MIN_NUM_C)
        {
            rec_cg[*num_curve].t_end = clrx[*end-1];
            rec_cg[*num_curve].t_break = 0;
            rec_cg[*num_curve].t_start = clrx[i_start];
            rec_cg[*num_curve].category = 20;
            rec_cg[*num_curve].num_obs = *end - i_start;
            rec_cg[*num_curve].t_confirmed = 0;
            rec_cg[*num_curve].change_prob = 0;

            for (i_b = 0; i_b < n_total_variable; i_b++)
            {
                rec_cg[*num_curve].magnitude[i_b] = 0.0;
//                rec_cg[*num_curve].obs_disturb[i_b] = 0.0;
//                for (i = 0; i < SCCD_MAX_NUM_C - 1; i++)
//                {
//                    rec_cg[*num_curve].state_disturb[i_b][i] = 0.0;
//                }
            }
            i_span = *end - i_start;
            update_cft(i_span, N_TIMES, MIN_NUM_C, MID_NUM_C, SCCD_MAX_NUM_C,
                      num_c, &update_num_c);
            for (i_b = 0; i_b < n_total_variable; i_b++)
            {
                status = auto_ts_fit(clrx, clry, i_b, i_b, i_start, *end-1, update_num_c,
                                     fit_cft1, &rmse[i_b], temp_v_dif);
                if (status != SUCCESS)
                {
                    RETURN_ERROR ("Calling auto_ts_fit for clear persistent pixels\n",
                                  FUNC_NAME, FAILURE);
                }

                rec_cg[*num_curve].rmse[i_b] = rmse[i_b];
                /******************************************************/
                /*                                                    */
                /* Record change magnitude.                           */
                /*                                                    */
                /******************************************************/

                for (k = 0; k < SCCD_MAX_NUM_C; k++)
                {
                    /**************************************************/
                    /*                                                */
                    /* Record fitted coefficients.                    */
                    /*                                                */
                    /**************************************************/

                    rec_cg[*num_curve].coefs[i_b][k] = fit_cft1[i_b][k];
                }

            }

        } // if (*end - i_start >= CONSE_END)
        else
        {
	    //printf("delete num curve = %d \n", *num_curve);
            // delete this curve
            *num_curve = *num_curve - 1;
        }
        for (i_b = 0; i_b < n_total_variable; i_b++)
        {
            /* need to re-fit ssm to get new Q and H*/
//            if(*end - i_start > MID_NUM_C * NUM_COEFFS)
//            {

//                step1_ssm_initialize(&instance[i_b], clrx, clry[i_b], i_start, i, vec_next_a[i_b], vec_next_P[i_b],
//                                         i_b, level_state_records, annual_state_records,
//                                     semi_state_records, rmse_records, starting_date, prev_i_break,
//                                     instance[i_b].Q, instance[i_b].H);
//            }
//            /* using previous instance Q and H to update states */
//            else


//            if (FALSE == b_fastmode)
//            {
//                for (k = clrx[prev_i_break]; k < clrx[*end-1]+1; k++)
//                {

//                    level_state_records[i_b][k-starting_date] = NA_VALUE;
//                    annual_state_records[i_b][k-starting_date] = NA_VALUE;
//                    semi_state_records[i_b][k-starting_date] = NA_VALUE;
//                    third_state_records[i_b][k-starting_date] = NA_VALUE;
//                }
//            }
        }
    }else if (bl_train_complete == 2)
    {

        /**********************************************************/
        /*                                                        */
        /* If no break, find at the end of the time series,       */
        /* define probability of change based on conse.           */
        /*                                                        */
        /**********************************************************/

        rec_cg[*num_curve].t_break = 0;
        rec_cg[*num_curve].num_obs = cur_i - i_start + 1;
        rec_cg[*num_curve].change_prob = 0;
        rec_cg[*num_curve].t_confirmed = 0;
        rec_cg[*num_curve].t_end = clrx[cur_i];
        rec_cg[*num_curve].category = 0;
        rec_cg[*num_curve].t_start = clrx[i_start];
        /* update ssm variable for the end*/
        for(i_b = 0; i_b < n_total_variable; i_b++)
        {
            // can't include the last obs (*end - 1)
            // needs to start from the previous break
//            tmp_q = sum_square_smooth_Q[i_b] / valid_count[i_b] - pow((sum_smooth_Q[i_b] / valid_count[i_b]), 2);
//            rec_cg[*num_curve].obs_disturb[i_b] = 0;
//            for (k = 0; k < SCCD_MAX_NUM_C - 1; k++)
//            {
//                rec_cg[*num_curve].state_disturb[i_b][k] = 0;
//            }


            for(k = 0; k < SCCD_MAX_NUM_C; k ++)
                 rec_cg[*num_curve].coefs[i_b][k] =fit_cft1[i_b][k];



//            array_1d_mean(rmse_records[i_b], i_start, *end - 1, &tmp);
//            rec_cg[*num_curve].rmse[i_b] = tmp;
            rec_cg[*num_curve].rmse[i_b] = 0;
            rec_cg[*num_curve].magnitude[i_b] = 0;
            rec_cg[*num_curve].t_confirmed = 0;

        }

    } // (bl_train == 1) // (bl_train == 0)


    for (k = 0; k < * num_curve; k++)
    {
        category =getcategorycurve(rec_cg, k);
        //category =getcategorycurve_old(rec_cg, k);
        rec_cg[k].category = rec_cg[k].category + category;
    }

    status = free_2d_array((void **)v_dif_mag_standard);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: v_dif_mag_standard\n", FUNC_NAME,
                      FAILURE);
    }


    status = free_2d_array ((void **) v_dif_mag);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: v_dif_mag\n", FUNC_NAME,
                      FAILURE);
    }


    free(bl_ids);
    bl_ids = NULL;

    free(ids);
    ids = NULL;

    free(rm_ids);
    rm_ids = NULL;

    free(pred_y_f);
    pred_y_f = NULL;

    free(pred_y);
    pred_y = NULL;

    free(rmse);
    rmse = NULL;

    status = free_2d_array ((void **) temp_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: temp_v_dif\n",
                      FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **) fit_cft_average);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fit_cft_average\n", FUNC_NAME,
                      FAILURE);
    }

    free(tmpcg_rmse_square);
    free(tmpcg_rmse);

    free(max_rmse);

    free(v_dif_mag_norm);
    status = free_2d_array ((void **) v_diff);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: v_diff\n",
                      FUNC_NAME, FAILURE);
    }

    return SUCCESS;

}

int sccd_stand_procedure
(
    int valid_num_scenes,             /* I:  number of scenes  */
    int *valid_date_array,           /* I: valid date time series  */
    short int **buf,                 /* I:  pixel-based time series  */
    short int *fmask_buf,           /* I:  mask-based time series  */
    int *id_range,
    Output_t_sccd *rec_cg,
    int *num_curve,                 /* Intialize NUM of Functional Curves    */
    char *states_output_dir,
    bool b_fastmode,
    double user_probability_threshold,
    int min_days_conse,
    int training_type,  /* for training process*/
    int monitorwindow_lowerlim,  /* for training process*/
    int monitorwindow_upperlim,  /* for training process*/
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
)
{
    int i_b = 5;
//    gsl_vector Z;         /* A m*1 vector corresponding to the system matrix of observation equation */
//    double H;            /* A 1*1 matrix corresponding to observational disturbance epsilon      */
//    gsl_matrix T;          /* A m x m matrix corresponding to the first system matrix of state equation. */
//    gsl_matrix Q;           /* A m*m matrix for the process disturbance covariance*/
//    gsl_vector a1;         /* A m x 1 vector containing the expected values of the initial states.    */
//    gsl_matrix P1;         /* A m x m matrix containing the covariance matrix of the nondiffuse part of the initial state vector */
//    gsl_matrix P1inf;      /* A m x m matrix containing the covariance matrix  */

    int status;
    int i, j, k;
    int i_copy;
    char FUNC_NAME[] = "sccd_stand_procedure";
    //char msg_str[MAX_STR_LEN];       /* Input data scene name                 */
    int end;
    double **fit_cft1;                 /* Fitted coefficients 2-D array.        */
    double **fit_cft2;                 /* Fitted coefficients 2-D array.        */
    double **fit_cft3;
    double **fit_cft4; /* Fitted coefficients 2-D array.        */
    int clrx_record1;
    int clrx_record2;
    int clrx_record3;
    int clrx_record4;
    double **tmp_fit_cft;                 /* temporal fitted coefficients 2-D array.        */
              /* Mimimum RMSE                          */
    double *rmse_ini;                 /* Root Mean Squared Error array for initialization stages       */
    //int rec_fc;                      /* Record num. of functional curves      */
    int i_dense;
    int i_span;
    double time_span;                 /* Span of time in no. of years.         */
    int n_clr;                  /* the number of clear pixels                          */
    int *clrx;                  /* clear pixel curve in X direction (date)             */
    float **clry;               /* clear pixel curve in Y direction (spectralbands)    */
    int *clrx_copy;
    float **clry_copy;
    float max_date_difference;   /* maximum difference between two neighbor dates        */
    //double **temp_v_dif;              /*            */
    float date_vario;           /* median date                                          */
    int i_start;
    int i_start_copy;
    int bl_train;
    int bl_train_complete;
    ssmodel* instance;
    //double** initial_a1;
    double w = TWO_PI / AVE_DAYS_IN_A_YEAR;
    FILE* fptr_save;

    gsl_vector** vec_next_a;          /* a vector of a for current i for multiple band */
    gsl_matrix** vec_next_P;          /* a vector p matrix for current i for multiple band */

    //gsl_matrix *initial_Q; /*  Q matrix of the previous curve  */
   // double initial_H;

    double** level_state_records;
    double** annual_state_records;
    double** semi_state_records;
    double** third_state_records;
    double** rmse_records;
    int i_span_min;
    int prev_i_break;
    int prev_break_date;
    // int prev_date_break;
    int starting_date; /* the beginning date of s-ccd*/
    char states_output_dir_full[MAX_STR_LEN];
    double rmse_sccd;

    float* adj_rmse; /* Adjusted RMSE for all bands          */
    int* valid_count;
    double* sum_square_smooth_Q;
    double* sum_square_smooth_H;
    double* sum_smooth_Q;
    double* sum_smooth_H;
    double* sum_square_vt;
    double* sum_vt;
    double* sum_kalman_coef;
    double* var_smooth_Q;
    double* var_smooth_H;

    int i_count;
    int adj_conse;
    double t_cg_adjust;
    double t_cg_outelier;
    double t_cg_gradual;
    int b_assigninstancememory = 0;
    double** rec_v_dif;
    int update_num_c = MIN_NUM_C;

    //double Q_tmp[SCCD_MAX_NUM_C - 1];
    double** temporal_rmse_square_sum;
    double** temporal_rmse_sum;
    int* temporal_rmse_count;
    double** yearly_temporal_rmse_square_sum;
    double** yearly_temporal_rmse_sum;
    int* yearly_temporal_rmse_count;
    int i_start_record;
    int i_copy_record;
    int  chow_test_start;
    double change_magnitude_tmp; /* for training procedure*/
    double max_CM_stable = -9999.0;
    double max_CM_dist = -9999.0;

    double s_tcg = X2(user_n_focus_variable, user_probability_threshold);
    //double s_tcg = 11.07;

    double TCWET, TCGRE, TCBRI;
    short int *sensor_buf_valid;

    int n_focus_variable;
    int n_total_variable;
    int focus_blist[TOTAL_IMAGE_BANDS + TOTAL_INDICES];
    double probability_threshold;
    float rf_feature[N_FEATURE];
    //int start_year;

    valid_count = (int*)malloc(user_n_total_variable * sizeof(int));
    sum_square_smooth_Q = (double*)malloc(user_n_total_variable * sizeof(double));
    sum_square_smooth_H = (double*)malloc(user_n_total_variable * sizeof(double));
    sum_smooth_Q = (double*)malloc(user_n_total_variable * sizeof(double));
    sum_smooth_H = (double*)malloc(user_n_total_variable * sizeof(double));
    sum_square_vt = (double*)malloc(user_n_total_variable * sizeof(double));
    sum_vt = (double*)malloc(user_n_total_variable * sizeof(double));
    sum_kalman_coef = (double*)malloc(user_n_total_variable * sizeof(double));
    var_smooth_Q = (double*)malloc(user_n_total_variable * sizeof(double));
    var_smooth_H = (double*)malloc(user_n_total_variable * sizeof(double));
    adj_rmse = (float*)malloc(user_n_total_variable * sizeof(float));
    // rf_feature = (float*)malloc(N_FEATURE * sizeof(float));
    //const float* out_result = NULL;
    //int tmp_total_variable;

    for (k = 0; k < user_n_total_variable; k++)
    {
            adj_rmse[k] = 0.0;
    }

    vec_next_a = (gsl_vector**) allocate_2d_array(user_n_total_variable, 1, sizeof(gsl_vector));
    if (vec_next_a == NULL)
    {
        RETURN_ERROR ("Allocating vec_next_a memory", FUNC_NAME, FAILURE);
    }

    vec_next_P = (gsl_matrix **)allocate_2d_array(user_n_total_variable, 1, sizeof(gsl_matrix));
    if (vec_next_P == NULL)
    {
        RETURN_ERROR ("Allocating vec_next_P memory", FUNC_NAME, FAILURE);
    }


    fit_cft1 = (double **) allocate_2d_array (user_n_total_variable, LASSO_COEFFS, sizeof (double));
    if (fit_cft1 == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft memory", FUNC_NAME, FAILURE);
    }

    fit_cft2 = (double **) allocate_2d_array (user_n_total_variable, LASSO_COEFFS, sizeof (double));
    if (fit_cft2 == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft memory", FUNC_NAME, FAILURE);
    }

    fit_cft3 = (double **) allocate_2d_array (user_n_total_variable, LASSO_COEFFS, sizeof (double));
    if (fit_cft3 == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft memory", FUNC_NAME, FAILURE);
    }

    fit_cft4 = (double **) allocate_2d_array (user_n_total_variable, LASSO_COEFFS, sizeof (double));
    if (fit_cft4 == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft memory", FUNC_NAME, FAILURE);
    }
    tmp_fit_cft = (double **) allocate_2d_array (user_n_total_variable, LASSO_COEFFS, sizeof (double));
    if (tmp_fit_cft == NULL)
    {
        RETURN_ERROR ("Allocating tmp_fit_cft memory", FUNC_NAME, FAILURE);
    }


    rmse_ini = (double *)malloc(user_n_total_variable * sizeof(double));
    if (rmse_ini == NULL)
    {
        RETURN_ERROR ("Allocating rmse_ini memory", FUNC_NAME, FAILURE);
    }

    clrx = (int *)malloc(valid_num_scenes * sizeof(int));
    if (clrx == NULL)
    {
        RETURN_ERROR ("Allocating clrx memory", FUNC_NAME, FAILURE);
    }

    clry = (float **) allocate_2d_array(user_n_total_variable, valid_num_scenes,
                                         sizeof (float));
    if (clry == NULL)
    {
        RETURN_ERROR ("Allocating clry memory", FUNC_NAME, FAILURE);
    }

    clrx_copy = (int *)malloc(valid_num_scenes * sizeof(int));
    if (clrx_copy == NULL)
    {
        RETURN_ERROR ("Allocating clrx_copy memory", FUNC_NAME, FAILURE);
    }

    clry_copy = (float **) allocate_2d_array(user_n_total_variable, valid_num_scenes,
                                         sizeof (float));
    if (clry_copy == NULL)
    {
        RETURN_ERROR ("Allocating clry_copy memory", FUNC_NAME, FAILURE);
    }

    /* alloc memory for ssm matrix */
    instance = malloc(user_n_total_variable * sizeof(ssmodel));
    if(instance == NULL)
    {
       RETURN_ERROR ("Allocating instance memory", FUNC_NAME, FAILURE);
    }

    sensor_buf_valid = (short int *)malloc(valid_num_scenes * sizeof(short int));
    if (sensor_buf_valid == NULL)
    {
        RETURN_ERROR ("Allocating clrx_copy memory", FUNC_NAME, FAILURE);
    }

    n_clr = 0;
    int clry_col_count = 7;
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
                clrx_copy[n_clr] = valid_date_array[i];
                sensor_buf_valid[n_clr] = sensor_buf[i];
                for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                {

                    clry[k][n_clr] = (float)buf[k][i];

                    clry_copy[k][n_clr] = clry[k][n_clr];

                    //printf("%3.2f\n", clry[k][n_clr]);
                }
                n_clr++;
            }
        }
    }
    // printf("n_clr is %d", n_clr);

    if(user_n_total_variable > TOTAL_IMAGE_BANDS){
        clry_col_count = 7;
        if(NDVI_INCLUDED == TRUE){
            for(i = 0; i < n_clr; i++){
                if((int)(clry[3][i] + clry[2][i]) == 0)
                    clry[clry_col_count][i] = 0;
                else
                    clry[clry_col_count][i] = roundf(10000 * ((float)(clry[3][i] - clry[2][i])/(clry[3][i] + clry[2][i])));
                clry_copy[clry_col_count][i] = clry[clry_col_count][i];
            }
            clry_col_count = clry_col_count + 1;
        }

        if(NBR_INCLUDED == TRUE){
            for(i = 0; i < n_clr; i++){
                if((int)(clry[3][i] + clry[5][i]) == 0)
                    clry[clry_col_count][i] = 0;
                else
                    clry[clry_col_count][i] = roundf(10000 * ((float)(clry[3][i] - clry[5][i])/(clry[3][i] + clry[5][i])));
                clry_copy[clry_col_count][i] = clry[clry_col_count][i];
            }
            clry_col_count = clry_col_count + 1;
        }

        if(RGI_INCLUDED == TRUE){
            for(i = 0; i < n_clr; i++){
                if((int)clry[1][i] == 0)
                    clry[clry_col_count][i] = 0;
                else{
                    if (10000 * clry[2][i]/clry[1][i] > 32767)
                        clry[clry_col_count][i] = 32767.0;
                    else
                        clry[clry_col_count][i] = roundf((10000 * (float)clry[2][i])/clry[1][i]);
                }
                clry_copy[clry_col_count][i] = clry[clry_col_count][i];
            }
            clry_col_count = clry_col_count + 1;
        }

        if(TCTWETNESS_INCLUDED == TRUE){
            for(i = 0; i < n_clr; i++){
                // source: http://www.sjsu.edu/faculty/watkins/tassel.htm
                if(sensor_buf_valid[i] == LANDSAT45_TM || sensor_buf_valid[i] == LANDSAT7_ETM)
                    clry[clry_col_count][i] = roundf(0.1509 * clry[0][i] + 0.1793 * clry[1][i] + 0.3299 * clry[2][i] +
                                0.3406 * clry[3][i] - 0.7112 * clry[4][i] - 0.4572 * clry[5][i]);
                else if (sensor_buf_valid[i] == LANDSAT8_OLI)
                    //source: https://yceo.yale.edu/tasseled-cap-transform-landsat-8-oli
                    clry[clry_col_count][i] = roundf(0.1511 * clry[0][i] + 0.1973 * clry[1][i] + 0.3283 * clry[2][i] +
                                0.3407 * clry[3][i] - 0.7117 * clry[4][i] - 0.4559 * clry[5][i]);
                clry_copy[clry_col_count][i] = clry[clry_col_count][i];
            }
            clry_col_count = clry_col_count + 1;
        }

        if(TCTGREENNESS_INCLUDED == TRUE){
            for (i = 0; i < n_clr; i++)
            {
                // source: http://www.sjsu.edu/faculty/watkins/tassel.htm
                if(sensor_buf_valid[i] == LANDSAT45_TM || sensor_buf_valid[i] == LANDSAT7_ETM)
                    clry[clry_col_count][i] = roundf(-0.2848 * clry[0][i] - 0.2435 * clry[1][i] - 0.5436 * clry[2][i] +
                                0.7243 * clry[3][i] + 0.0840 * clry[4][i] - 0.1800 * clry[5][i]);
                else if (sensor_buf_valid[i] == LANDSAT8_OLI)
                    //source: https://yceo.yale.edu/tasseled-cap-transform-landsat-8-oli
                    clry[clry_col_count][i] = roundf(- 0.2941 * clry[0][i] - 0.243 * clry[1][i] - 0.5424 * clry[2][i] +
                                0.7276 * clry[3][i] + 0.0713 * clry[4][i] - 0.1608 * clry[5][i]);
                clry_copy[clry_col_count][i] = clry[clry_col_count][i];
             }
             clry_col_count = clry_col_count + 1;
        }

        if(EVI_INCLUDED == TRUE){
            for (i = 0; i < n_clr; i++)
            {
                if((int)(clry[3][i] + 6 * clry[2][i] - 7.5 * clry[0][i] + 10000) == 0)
                    clry[clry_col_count][i] = 0.0;
                else
                    clry[clry_col_count][i] = roundf(10000 * 2.5 * (clry[3][i] - clry[2][i])/(clry[3][i] + 6 * clry[2][i] - 7.5 * clry[0][i] + 10000));
                clry_copy[clry_col_count][i] = clry[clry_col_count][i];
            }
            clry_col_count = clry_col_count + 1;
        }

        if(DI_INCLUDED == TRUE){
            for (i = 0; i < n_clr; i++)
            {
                if(sensor_buf_valid[i] == LANDSAT45_TM || sensor_buf_valid[i] == LANDSAT7_ETM){
                    TCWET = roundf(0.1509 * clry[0][i] + 0.1793 * clry[1][i] + 0.3299 * clry[2][i] +
                                0.3406 * clry[3][i] - 0.7112 * clry[4][i] - 0.4572 * clry[5][i]);
                    TCGRE = roundf(-0.2848 * clry[0][i] - 0.2435 * clry[1][i] - 0.5436 * clry[2][i] +
                                0.7243 * clry[3][i] + 0.0840 * clry[4][i] - 0.1800 * clry[5][i]);
                    TCBRI = roundf(0.3037 * clry[0][i] + 0.2793 * clry[1][i] + 0.4343 * clry[2][i] +
                                0.5585 * clry[3][i] + 0.5082 * clry[4][i] + 0.1863 * clry[5][i]);
                }
                else if (sensor_buf_valid[i] == LANDSAT8_OLI){
                    //source: https://yceo.yale.edu/tasseled-cap-transform-landsat-8-oli
                    TCWET = roundf(0.1511 * clry[0][i] + 0.1973 * clry[1][i] + 0.3283 * clry[2][i] +
                                0.3407 * clry[3][i] - 0.7117 * clry[4][i] - 0.4559 * clry[5][i]);
                    TCGRE = roundf(- 0.2941 * clry[0][i] - 0.243 * clry[1][i] - 0.5424 * clry[2][i] +
                                0.7276 * clry[3][i] + 0.0713 * clry[4][i] - 0.1608 * clry[5][i]);
                    TCBRI = roundf(0.3029 * clry[0][i] + 0.2786 * clry[1][i] + 0.4733 * clry[2][i] +
                                0.5599 * clry[3][i] + 0.508 * clry[4][i] + 0.1872 * clry[5][i]);
                }
                clry[clry_col_count][i] = TCBRI - TCWET -TCGRE;
                clry_copy[clry_col_count][i] = clry[clry_col_count][i];
            }
            clry_col_count = clry_col_count + 1;
        }

        if(NDMI_INCLUDED == TRUE){
            for(i = 0; i < n_clr; i++){
                if((int)(clry[3][i] + clry[4][i]) == 0)
                    clry[clry_col_count][i] = 0;
                else
                    clry[clry_col_count][i] = roundf(10000 * ((float)(clry[3][i] - clry[4][i])/(clry[3][i] + clry[4][i])));
                clry_copy[clry_col_count][i] = clry[clry_col_count][i];
            }
            clry_col_count = clry_col_count + 1;
        }

    }


    rmse_records = (double **)allocate_2d_array (user_n_total_variable, n_clr, sizeof (double));
    if (rmse_records == NULL)
    {
        RETURN_ERROR ("Allocating rmse_records memory", FUNC_NAME, FAILURE);
    }

    temporal_rmse_square_sum = (double **) allocate_2d_array (user_n_total_variable, TEMPORAL_RMSE_BIN_NUM, sizeof (double));
    if (temporal_rmse_square_sum == NULL)
    {
        RETURN_ERROR ("Allocating temporal_rmse_square_sum memory", FUNC_NAME, FAILURE);
    }

    temporal_rmse_sum = (double **) allocate_2d_array (user_n_total_variable, TEMPORAL_RMSE_BIN_NUM, sizeof (double));
    if (temporal_rmse_sum == NULL)
    {
        RETURN_ERROR ("Allocating temporal_rmse_square_sum memory", FUNC_NAME, FAILURE);
    }

    temporal_rmse_count = (int*)malloc(TEMPORAL_RMSE_BIN_NUM * sizeof(int));
    if (temporal_rmse_count == NULL)
    {
        RETURN_ERROR ("Allocating temporal_rmse_count memory", FUNC_NAME, FAILURE);
    }

    yearly_temporal_rmse_square_sum = (double **) allocate_2d_array (user_n_total_variable, TEMPORAL_RMSE_BIN_NUM, sizeof (double));
    if (yearly_temporal_rmse_square_sum == NULL)
    {
        RETURN_ERROR ("Allocating temporal_rmse_square_sum memory", FUNC_NAME, FAILURE);
    }

    yearly_temporal_rmse_sum = (double **) allocate_2d_array (user_n_total_variable, TEMPORAL_RMSE_BIN_NUM, sizeof (double));
    if (yearly_temporal_rmse_sum == NULL)
    {
        RETURN_ERROR ("Allocating temporal_rmse_square_sum memory", FUNC_NAME, FAILURE);
    }

    yearly_temporal_rmse_count = (int*)malloc(TEMPORAL_RMSE_BIN_NUM * sizeof(int));
    if (yearly_temporal_rmse_count == NULL)
    {
        RETURN_ERROR ("Allocating temporal_rmse_count memory", FUNC_NAME, FAILURE);
    }

    /* starting date used to record state values*/
    starting_date = clrx[0];
    end = n_clr;

    level_state_records = (double**) allocate_2d_array (user_n_total_variable,
                                             clrx[end - 1] - starting_date + 1, sizeof(double));
    if (level_state_records == NULL)
        RETURN_ERROR ("Allocating level_state_records memory", FUNC_NAME, FAILURE);

    annual_state_records = (double**)allocate_2d_array (user_n_total_variable,
                                             clrx[end - 1] - starting_date + 1, sizeof(double));
    if (annual_state_records == NULL)
        RETURN_ERROR ("Allocating annual_state_records memory", FUNC_NAME, FAILURE);

    semi_state_records = (double**)allocate_2d_array (user_n_total_variable,
                                            clrx[end - 1] - starting_date + 1, sizeof(double));
    if (semi_state_records == NULL)
        RETURN_ERROR ("Allocating semi_state_records memory", FUNC_NAME, FAILURE);

    third_state_records = (double**)allocate_2d_array (user_n_total_variable,
                                            clrx[end - 1] - starting_date + 1, sizeof(double));
    if (third_state_records == NULL)
        RETURN_ERROR ("Allocating third_state_records memory", FUNC_NAME, FAILURE);

    rec_v_dif = (double **)allocate_2d_array(user_n_total_variable, n_clr,
                                     sizeof (double));
    if (rec_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating rec_v_dif memory",FUNC_NAME, FAILURE);
    }

    for(i = 0; i < TEMPORAL_RMSE_BIN_NUM; i++)
    {
        for(i_b = 0; i_b < user_n_total_variable; i_b++)
        {
            temporal_rmse_square_sum[i_b][i] = 0;
            temporal_rmse_sum[i_b][i] = 0;
            yearly_temporal_rmse_square_sum[i_b][i] = 0;
            yearly_temporal_rmse_sum[i_b][i] = 0;
        }
        temporal_rmse_count[i] = 0;
        yearly_temporal_rmse_count[i] = 0;
    }


    /**************************************************************/
    /*                                                            */
    /* Start with mininum requirement of clear obs.               */
    /*                                                            */
    /**************************************************************/

    i = N_TIMES * (MID_NUM_C)  - 1;
    if (i > end){
        i = end;
    }
    i_copy = i;

    /**************************************************************/
    /*                                                            */
    /* The first observation for TSFit.                           */
    /*                                                            */
    /**************************************************************/

    i_start = 0;

    i_start_copy = 0;

    i_dense = 0;

    /**************************************************************/
    /*                                                            */
    /* Record the start of the model initialization               */
    /*     (0=>initial;1=>done)                                   */
    /*                                                            */
    /**************************************************************/

    bl_train = 0;
    bl_train_complete = 0;

    /**************************************************************/
    /*                                                            */
    /* initialize of i_break for the first curve.                 */
    /*                                                            */
    /**************************************************************/
    prev_i_break = 0;
    prev_break_date = clrx[prev_i_break];

//    status = adjust_median_variogram(clrx, clry, n_total_variable, 0, end, &date_vario,
//                                     &max_date_difference, adj_rmse, 1);
//    if (status != SUCCESS)
//    {
//            RETURN_ERROR("ERROR calling median_variogram routine", FUNC_NAME,
//                         FAILURE);
//    }
    if(b_landspecific == TRUE && rec_cg[0].land_type != TARGETED_LABEL){
        n_focus_variable = DEFAULT_N_FOCUS_VARIABLE;
        n_total_variable = DEFAULT_TOTAL_VARIABLE;
        probability_threshold = DEFAULT_PROBABILITY;
        for(k = 0; k < n_focus_variable; k++)
            focus_blist[k] = DEFAULT_FOCUS_BLIST[k];
    }else{
        n_focus_variable = user_n_focus_variable;
        n_total_variable = user_n_total_variable;
        probability_threshold = user_probability_threshold;
        for(k = 0; k < n_focus_variable; k++)
            focus_blist[k] = user_focus_blist[k];
    }


    // initialization
    if(b_fastmode==FALSE)
    {
        for (i_b = 0; i_b < user_n_total_variable; i_b++)
        {
            for(k = 0; k < clrx[end - 1] - starting_date + 1; k++)
            {
                    level_state_records[i_b][k] = NA_VALUE;
                    annual_state_records[i_b][k] = NA_VALUE;
                    semi_state_records[i_b][k] = NA_VALUE;
                    third_state_records[i_b][k] = NA_VALUE;
            }
        }
    }

    /**************************************************************/
    /*                                                            */
    /* While loop - process til the last clear observation - conse*/
    /*                                                            */
    /**************************************************************/
    while ((clrx[end - 1] - clrx[i - 1] > min_days_conse) && (i + conse - 1 < end))
    {
        if(0 == bl_train_complete)
        {
            if(0 == bl_train)
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
                time_span = (double)(clrx[i] - clrx[i_start]) / NUM_YEARS;

                /***************************************************************/
                /*                                                             */
                /* for the first curve, we need more time points to start with */
                /*                                                             */
                /***************************************************************/
//                if(*num_curve == 0)
//                    i_span_min = N_TIMES * MID_NUM_C;
//                else
                    i_span_min = N_TIMES * MIN_NUM_C;

                status = INCOMPLETE;
                if((i_span >= i_span_min) && (time_span >= (double)MIN_YEARS))
                {
                    /**************************************************************/
                    /*                                                            */
                    /* step1: initialize a ccd model.                            */
                    /*                                                            */
                    /**************************************************************/

                    /* reset prev_i_break to be the i_start for the first curve, cuz we won't consider */
                    /*  the points before i_start*/

                    /**************************************************************/
                    /*                                                            */
                    /* calculate variogram for each band and dates.                           */
                    /*                                                            */
                    /**************************************************************/
                    status = adjust_median_variogram(clrx_copy, clry_copy, n_total_variable, i_dense, i_copy, &date_vario,
                                                     &max_date_difference, adj_rmse, 1);
                    if (status != SUCCESS)
                    {
                            RETURN_ERROR("ERROR calling median_variogram routine", FUNC_NAME,
                                         FAILURE);
                    }

                    adj_conse = round (conse * 16 / (double)date_vario);
                    if (adj_conse < conse)
                        adj_conse = conse;

                     /* adjust conse based delta days*/
                    if(adj_conse > conse)
                    {
                        // adj_TCG = chi2inv(1 - pow(1 - PROB_T_CG, (double)conse / (double)adj_conse), NUM_LASSO_BANDS);
                        t_cg_adjust = X2(n_focus_variable, 1 - pow(1 - probability_threshold, (double)conse / (double)adj_conse));
                        //t_cg_adjust = 11.07;
                    }
                    else
                    {
                        t_cg_adjust = s_tcg;
                    }

                    status = step1_ccd_initialize(adj_conse, adj_rmse, n_clr, t_cg_adjust, &i_dense, num_curve, clrx,
                                                  clry, &i, &i_start, &i_start_copy, &end, tmp_fit_cft, rec_cg, i_span_min,
                                                  &prev_i_break, rmse_ini, n_focus_variable, n_total_variable, focus_blist, min_days_conse);

                    /**************************************************************/
                    /*                                                            */
                    /* step 1: initialize a state-space model.                    */
                    /*                                                            */
                    /**************************************************************/

                }

                if(INCOMPLETE == status)
                {
                    i++;
                    i_copy++;
                    continue;
                }
                else if(SUCCESS == status)
                {
                    bl_train = 1;

                    if(training_type == TRAINING_TYPE_CLASSIFICATION)
                    {
                        for (i_b = 0; i_b < n_total_variable; i_b++)
                        {

                            status = auto_ts_fit(clrx, clry, i_b, i_b, i_start, i, SCCD_MAX_NUM_C,
                                                 fit_cft1, &rmse_ini[i_b], rec_v_dif);
                            if (status != SUCCESS)
                            {
                                RETURN_ERROR ("Calling auto_ts_fit during continuous monitoring\n",
                                              FUNC_NAME, FAILURE);
                            }

                        }
                        bl_train = 1;
                        bl_train_complete = 2; // bl_train_complete = 2 indicate this for TRAINING_TYPE_CLASSIFICATION
                        end = i + conse - 1; // set i to the max so jump to the end
                        continue;
                    }else{
                        /* the below is to determine the land type*/
//                        if (b_landspecific == TRUE){

//                            if (*num_curve == 0){
//                                if(rec_cg[0].land_type == TARGETED_LABEL){
//                                    n_focus_variable = user_n_focus_variable;
//                                    n_total_variable = user_n_total_variable;
//                                    probability_threshold = user_probability_threshold;
//                                    for(k = 0; k < n_focus_variable; k++)
//                                        focus_blist[k] = user_focus_blist[k];
//                                }else{
//                                    n_focus_variable = DEFAULT_N_FOCUS_VARIABLE;
//                                    n_total_variable = DEFAULT_TOTAL_VARIABLE;
//                                    probability_threshold = DEFAULT_PROBABILITY;
//                                    for(k = 0; k < n_focus_variable; k++)
//                                        focus_blist[k] = DEFAULT_FOCUS_BLIST[k];
//                                }
//                            } // initial curve
//                            else{


//                                for (i_b = 0; i_b < N_FEATURE - 1; i_b++)
//                                {

//                                    status = auto_ts_fit(clrx, clry, i_b, i_b, i_start, i, SCCD_MAX_NUM_C,
//                                                         tmp_fit_cft, &rmse_ini[i_b], rec_v_dif);
//                                    if (status != SUCCESS)
//                                    {
//                                        RETURN_ERROR ("Calling auto_ts_fit during continuous monitoring\n",
//                                                      FUNC_NAME, FAILURE);
//                                    }

//                                    standand_peekday = convert_standardordinal(clrx[i_start], clrx[i]);

//                                    rf_feature[i_b] = (float)(tmp_fit_cft[i_b][0] + tmp_fit_cft[i_b][1] * standand_peekday);
//                                }
//                                rf_feature[N_FEATURE - 1] = (float)auxval;
//                                DMatrixHandle out;
//                                DMatrixHandle eval_dmats[] = {};

//                                // safe_xgboost(XGBoosterCreate(eval_dmats, 0, &booster));
//                                safe_xgboost(XGDMatrixCreateFromMat(rf_feature, 1, N_FEATURE, -9999, &out));
//                                safe_xgboost(XGBoosterPredict(booster, out, 0, 0, 0, &out_len, &out_result));

//                                if(out_result[0] == TARGETED_LABEL){
//                                    n_focus_variable = user_n_focus_variable;
//                                    n_total_variable = user_n_total_variable;
//                                    probability_threshold = user_probability_threshold;
//                                    for(k = 0; k < n_focus_variable; k++)
//                                        focus_blist[k] = user_focus_blist[k];
//                                }else{
//                                    n_focus_variable = DEFAULT_N_FOCUS_VARIABLE;
//                                    n_total_variable = DEFAULT_TOTAL_VARIABLE;
//                                    probability_threshold = DEFAULT_PROBABILITY;
//                                    for(k = 0; k < n_focus_variable; k++)
//                                        focus_blist[k] = DEFAULT_FOCUS_BLIST[k];
//                                }
//                                rec_cg[*num_curve].land_type = (int)out_result[0];
//                                safe_xgboost(XGDMatrixFree(out));

//                            }


//                        }
                    }
                } /* else if(SUCCESS == status) */
            } /* if(0 == bl_train) */
            else // (1 == bl_train)
            {
                if(i - i_start + 1 < N_TIMES * SCCD_MAX_NUM_C + 1)
                {
                    if(clrx[i + conse - 1] - clrx[i - 1] < min_days_conse)
                    //if(clrx[i + conse - 1] - clrx[i] < min_days_conse)
                    {
                        adj_conse = conse;
                        while((clrx[i + adj_conse - 1] - clrx[i - 1] < min_days_conse)
                              && i + adj_conse - 1 < end
                              && adj_conse <= 3 * conse)
                        {
                            adj_conse = adj_conse + 1;
                        }
                        //t_cg_adjust = s_tcg;

                        t_cg_adjust = X2(n_focus_variable, 1 - pow(1 - probability_threshold, (double)conse / (double)adj_conse));
                        //t_cg_adjust = 11.07;
                    }
                    else
                    {
                        adj_conse = conse;
                        t_cg_adjust = s_tcg;

                    }
                    t_cg_outelier = X2(n_focus_variable, S_T_MAX_CG_PROB);
                    //t_cg_outelier = 35;
                    status =  step1_update_cft(adj_conse, adj_rmse, n_clr, t_cg_adjust, clrx, clry, i, i_start,
                                               tmp_fit_cft, rec_cg, rmse_ini, num_curve, end, &prev_i_break,
                                               n_focus_variable, n_total_variable, focus_blist, t_cg_outelier);
                    if(status == CHANGEDETECTED)
                    {
                        rec_cg[*num_curve].num_obs = i_span;
                        rec_cg[*num_curve].change_prob = 100;
                        /**********************************************/
                        /*                                            */
                        /* Start from i for the next functional     */
                        /* curve.                                     */
                        /*                                            */
                        /**********************************************/
                        i_start = i + 1;

                        i_start_copy = i_copy;

                        /**********************************************/
                        /*                                            */
                        /* Start training again.                      */
                        /*                                            */
                        /**********************************************/
                        bl_train = 0;
                        prev_break_date = clrx[prev_i_break];
                    }
                    else if(status == FALSECHANGE)
                    {
                        end--;
                        i--;
                    }
                } // if(i - i_start + 1 < N_TIMES * SCCD_MAX_NUM_C)
                else /* begin ssm initialization*/
                {
                    // initalization
                    for(i_b = 0; i_b < n_total_variable; i_b++)
                    {
                        valid_count[i_b] = 0;
                        sum_square_smooth_Q[i_b] = 0;
                        sum_square_smooth_H[i_b] = 0;
                        sum_smooth_Q[i_b] = 0;
                        sum_smooth_H[i_b] = 0;
                        sum_square_vt[i_b] = 0;
                        sum_vt[i_b] = 0;
                        sum_kalman_coef[i_b] = 0;
                    }

                    for(j = 0; j < TEMPORAL_RMSE_BIN_NUM; j++)
                    {
                        for(i_b = 0; i_b < n_total_variable; i_b++)
                            temporal_rmse_square_sum[i_b][j] = 0;
                        temporal_rmse_count[j] = 0;
                    }
                    prev_break_date = clrx[prev_i_break];

                    i--;

                    for (i_b = 0; i_b < n_total_variable; i_b++)
                    {

                        status = auto_ts_fit(clrx, clry, i_b, i_b, i_start, i, SCCD_MAX_NUM_C,
                                             tmp_fit_cft, &rmse_ini[i_b], rec_v_dif);
                        if (status != SUCCESS)
                        {
                            RETURN_ERROR ("Calling auto_ts_fit during continuous monitoring\n",
                                          FUNC_NAME, FAILURE);
                        }

                    }


                    for(i_b = 0; i_b < n_total_variable; i_b++)
                    {
                        /* eliminate slope as we only use m = 5*/
                        for(j = 0; j < SCCD_MAX_NUM_C; j++)
                        {
                            if(j == 0)
                                fit_cft1[i_b][j] = (double)tmp_fit_cft[i_b][j] + (double)tmp_fit_cft[i_b][1] * clrx[i_start];
                            else if(j == 1)
                                fit_cft1[i_b][j] = 0;
                            else
                                fit_cft1[i_b][j] = (double)tmp_fit_cft[i_b][j];
                        }

                        // determine m
                        instance[i_b].structure = 0;
                        instance[i_b].m = SCCD_MAX_NUM_C - 1;
                            /* determine num_c */

                        if (SCCD_MAX_NUM_C == 8)
                        {
                            if((fit_cft1[i_b][6] == 0) && (fit_cft1[i_b][7] == 0))
                                /* only trend and annual cycle*/
                                /*  compute a vector from fit_cft */
                                instance[i_b].m = instance[i_b].m - 2;
                            else
                                instance[i_b].structure = instance[i_b].structure + 100;

                        }

//                        if (SCCD_MAX_NUM_C > 5)
//                        {
//                            if((fit_cft1[i_b][4] == 0) && (fit_cft1[i_b][5] == 0))
//                                /* only trend and annual cycle*/
//                                instance[i_b].m = instance[i_b].m - 2;
//                            else
//                                instance[i_b].structure = instance[i_b].structure + 10;
//                        }


//                        if (SCCD_MAX_NUM_C > 3)
//                        {
//                            if((fit_cft1[i_b][2] == 0) && (fit_cft1[i_b][3] == 0))
//                                /*  compute a vector from fit_cft */
//                                instance[i_b].m = instance[i_b].m - 2;
//                            else
//                                instance[i_b].structure = instance[i_b].structure + 1;
//                        }
                        // at least 3 coefficients
                        instance[i_b].structure = instance[i_b].structure + 10;
                        instance[i_b].structure = instance[i_b].structure + 1;

                        instance[i_b].a1 = gsl_vector_alloc (instance[i_b].m);
                        instance[i_b].Z = gsl_vector_alloc (instance[i_b].m);
                        instance[i_b].T = gsl_matrix_calloc (instance[i_b].m, instance[i_b].m);
                        instance[i_b].Q = gsl_matrix_calloc (instance[i_b].m, instance[i_b].m);
                        instance[i_b].P1 = gsl_matrix_calloc (instance[i_b].m, instance[i_b].m);
                        instance[i_b].P1inf = gsl_matrix_calloc (instance[i_b].m, instance[i_b].m);
                        /* alloc memory for each element*/

                        vec_next_a[i_b] =  gsl_vector_alloc (instance[i_b].m);
                        vec_next_P[i_b] = gsl_matrix_calloc (instance[i_b].m, instance[i_b].m);

                        b_assigninstancememory = 1;

                        /**************************************************************/
                        /*                                                            */
                        /*  initialize ssm variables                                  */
                        /*                                                            */
                        /**************************************************************/
                        valid_count[i_b] = 0;
                        sum_square_smooth_Q[i_b] = 0;
                        sum_square_smooth_H[i_b] = 0;
                        sum_smooth_Q[i_b] = 0;
                        sum_smooth_H[i_b] = 0;
                        sum_square_vt[i_b] = 0;
                        sum_vt[i_b] = 0;
                        sum_kalman_coef[i_b] = 0;

                        step1_ssm_initialize(&instance[i_b], clrx, clry[i_b], i_start, i, vec_next_a[i_b], vec_next_P[i_b],
                                             i_b, level_state_records,  annual_state_records,
                                             semi_state_records, third_state_records, rmse_records, starting_date, prev_break_date,
                                             *num_curve,  rec_cg, rmse_ini[i_b], adj_rmse[i_b], b_fastmode, fit_cft1[i_b],
                                             valid_count, sum_square_smooth_Q, sum_square_smooth_H, sum_smooth_Q, sum_smooth_H, sum_square_vt,
                                             sum_vt, sum_kalman_coef, temporal_rmse_square_sum, temporal_rmse_sum, temporal_rmse_count);

                        /* record t_stable */
                        rec_cg[*num_curve].t_start = clrx[i_start];

                        /* initialization fitcft */
                        for (j = 0; j < SCCD_MAX_NUM_C; j++)
                        {
                           fit_cft2[i_b][j] = fit_cft1[i_b][j];
                           fit_cft3[i_b][j] = fit_cft1[i_b][j];
                           fit_cft4[i_b][j] = fit_cft1[i_b][j];
                        }

                        clrx_record1 = clrx[i];
                        clrx_record2 = clrx[i];
                        clrx_record3 = clrx[i];
                        clrx_record4 = clrx[i];

                    }

                    status = adjust_median_variogram(clrx_copy, clry_copy, n_total_variable, i_start_copy, i_copy, &date_vario,
                                                     &max_date_difference, adj_rmse, 1);
                    if (status != SUCCESS)
                    {
                            RETURN_ERROR("ERROR calling median_variogram routine", FUNC_NAME,
                                         FAILURE);
                    }
//                    for (k = 0; k < n_total_variable; k++)
//                    {
//                            adj_rmse[k] = 0.0;
//                    }

                    /* initialization stage stops, and continious change detection starts */
                    bl_train_complete  =  1;
                    update_num_c = 4;
                    i_count = clrx[i];
                    i++;
                    i_copy++;
                    continue;
                }
            } /* if (1 == bl_train) */
        }   /* if(0==bl_train_complete)*/
        /**************************************************************/
        /*                                                            */
        /*    step 2: kalman-filter change detection                  */
        /*                                                            */
        /**************************************************************/
        else
        {

            //int adj_conse_record;
            /* adjust conse to include at_least min_days_conse days for monitoring window */
            if(clrx[i + conse - 1] - clrx[i - 1] < min_days_conse)
            {
                adj_conse = conse;
                while((clrx[i + adj_conse - 1] - clrx[i - 1] < min_days_conse)
                      && i + adj_conse - 1 < end
                      && adj_conse <= 3 * conse)
                {
                    adj_conse = adj_conse + 1;
                }
                //adj_conse_record = adj_conse;
                //adj_conse = round((double)adj_conse / conse) * conse;
                t_cg_adjust = X2(n_focus_variable, 1 - pow(1 - probability_threshold, (double)conse / (double)adj_conse));
                //t_cg_adjust = 11.07;
            }
            else
            {
                adj_conse = conse;
                t_cg_adjust = s_tcg;
                //t_cg_outelier = S_TCG_MIN;
            }

            // t_cg_adjust = chi2inv(1 - pow(1 - S_TCG, (double)conse / (double)adj_conse), n_focus_variable);
            //t_cg_adjust = X2(n_focus_variable, 1 - pow(1 - S_TCG, (double)conse / (double)adj_conse));

            t_cg_outelier = X2(n_focus_variable, S_T_MAX_CG_PROB);
            // t_cg_outelier = 35;
//            int tt;
//            tt = clrx[i];

            i_span = i - i_start +1;

            status = step2_KF_ChangeDetection(instance, vec_next_P, vec_next_a, clrx, clry, i, level_state_records,
                                              annual_state_records, semi_state_records, third_state_records,
                                              rmse_records, rec_cg, num_curve, &end, adj_conse, starting_date,
                                              i_start, &prev_i_break, b_fastmode, fit_cft1, fit_cft2, fit_cft3, fit_cft4,
                                              &clrx_record1, &clrx_record2, &clrx_record3, &clrx_record4,
                                              valid_count, sum_square_smooth_Q, sum_square_smooth_H, sum_smooth_Q, sum_smooth_H,
                                              sum_square_vt, sum_vt, sum_kalman_coef, &i_count, adj_rmse, t_cg_adjust, t_cg_outelier,
                                              t_cg_gradual, temporal_rmse_square_sum, temporal_rmse_sum, temporal_rmse_count,
                                              yearly_temporal_rmse_square_sum, yearly_temporal_rmse_sum, yearly_temporal_rmse_count,
                                              probability_threshold, &change_magnitude_tmp, training_type, n_focus_variable, n_total_variable, focus_blist);

            if(status != CHANGEDETECTED){


                if(training_type == TRAINING_TYPE_PARAMETER)
                {
                    if (clrx[i] < monitorwindow_lowerlim){
                        if (change_magnitude_tmp > max_CM_stable)
                            max_CM_stable = change_magnitude_tmp;
                    }
                    else if((clrx[i] > monitorwindow_lowerlim) && (clrx[i] < monitorwindow_upperlim)){
                        if (change_magnitude_tmp > max_CM_dist)
                            max_CM_dist = change_magnitude_tmp;
                    }else if(clrx[i] > monitorwindow_upperlim - 1){
                        i = end - 2; // set i to the max so jump to the end
                        continue;
                    }
                }


            }

            if(status == CHANGEDETECTED)
            {
                if(b_assigninstancememory == 1)
                {
                    for(i_b = 0; i_b < n_total_variable; i_b ++)
                    {
                       gsl_vector_free(vec_next_a[i_b]);
                       gsl_matrix_free(vec_next_P[i_b]);
                    }

                    for(i_b = 0; i_b < n_total_variable; i_b++)
                    {
                        gsl_vector_free(instance[i_b].Z);
                        //(&instance[i_b])->Z = NULL;
                        //gsl_matrix_free (instance[i_b].H);
                        gsl_matrix_free(instance[i_b].T);
                        //(&instance[i_b])->T = NULL;
                        gsl_matrix_free(instance[i_b].Q);
                        //(&instance[i_b])->Q = NULL;
                        gsl_vector_free(instance[i_b].a1);
                        //(&instance[i_b])->a1 = NULL;
                        gsl_matrix_free(instance[i_b].P1);
                        //(&instance[i_b])->P1 = NULL;
                        gsl_matrix_free(instance[i_b].P1inf);
                        //(&instance[i_b])->P1inf = NULL;
                    }
                    b_assigninstancememory = 0;

                }

                rec_cg[*num_curve].num_obs = i_span;
                rec_cg[*num_curve].change_prob = 100;
                /**********************************************/
                /*                                            */
                /* Start from i for the next functional     */
                /* curve.                                     */
                /*                                            */
                /**********************************************/
                i_start_record = i_start; // make record in case that chow test fail
                i_copy_record = i_start_copy;

                i_start = i + 1;
                i_start_copy = i_copy;

                /**********************************************/
                /*                                            */
                /* Start training again.                      */
                /*                                            */
                /**********************************************/
                bl_train = 0;
                bl_train_complete = 0;

                /**********************************************/
                /*                                            */
                /* use default parameters for initialization    */
                /*                                            */
                /**********************************************/

                if(b_landspecific == TRUE){
                    n_focus_variable = DEFAULT_N_FOCUS_VARIABLE;
                    n_total_variable = DEFAULT_TOTAL_VARIABLE;
                    probability_threshold = DEFAULT_PROBABILITY;
                    for(k = 0; k < n_focus_variable; k++)
                        focus_blist[k] = DEFAULT_FOCUS_BLIST[k];
                }

                /**********************************************/
                /*                                            */
                /* INCREASE num_curve if it exceeds the limit  */
                /*                                            */
                /**********************************************/

                if (*num_curve >= NUM_FC)
                {
                    /**************************************************/
                    /*                                                */
                    /* Reallocate memory for rec_cg.                  */
                    /*                                                */
                    /**************************************************/

                    rec_cg = realloc(rec_cg, (*num_curve + 1) * sizeof(Output_t));
                    if (rec_cg == NULL)
                    {
                        RETURN_ERROR("ERROR allocating rec_cg memory",
                                     FUNC_NAME, FAILURE);
                    }
                }


            }
            else if(status == FALSECHANGE)
            {
                //printf("%d\n", clrx[i]);
                end--;
                i--;
            }

            //update coefficient for each year
            if((status == REGULAREND) && (clrx[i] - i_count > NUM_YEARS))
            {

                status = adjust_median_variogram(clrx_copy, clry_copy, n_total_variable, i_start_copy, i_copy, &date_vario,
                                                 &max_date_difference, adj_rmse, 1);
                if (status != SUCCESS)
                {
                        RETURN_ERROR("ERROR calling median_variogram routine", FUNC_NAME,
                                     FAILURE);
                }

//                for(i_b = 0; i_b < n_total_variable; i_b++)
//                {
//                    var_smooth_Q[i_b] = sum_square_smooth_Q[i_b] / (valid_count[i_b] - 1);
//                    var_smooth_H[i_b] = sum_square_smooth_H[i_b] / (valid_count[i_b] - 1);
//                    instance[i_b].H = var_smooth_H[i_b];
//                      //instance[i_b].H = sum_square_vt[i_b] / (valid_count[i_b]);
//                    for(k = 0; k< instance[i_b].m; k++)
//                        Q_tmp[k] = var_smooth_Q[i_b] *
//                            gsl_matrix_get(instance[i_b].Q, k, k) / gsl_matrix_get(instance[i_b].Q, 0, 0);

//                    for(k = 0; k< instance[i_b].m; k++)
//                        gsl_matrix_set(instance[i_b].Q, k, k, Q_tmp[k]);

//                  }

//                for(k = 0; k < TEMPORAL_RMSE_BIN_NUM; k++)
//                {
//                    for(i_b = 0; i_b < n_total_variable; i_b++)
//                    {
//                        temporal_rmse_square_sum[i_b][k] = yearly_temporal_rmse_square_sum[i_b][k] + temporal_rmse_square_sum[i_b][k];
//                        temporal_rmse_sum[i_b][k] = yearly_temporal_rmse_sum[i_b][k] + temporal_rmse_sum[i_b][k] ;
//                        yearly_temporal_rmse_square_sum[i_b][k] = 0;
//                        yearly_temporal_rmse_sum[i_b][k] = 0;
//                    }
//                    temporal_rmse_count[k] = yearly_temporal_rmse_count[k] + temporal_rmse_count[k];
//                    yearly_temporal_rmse_count[k] = 0;
//                }

                i_count = clrx[i];


            }

        }

        i++;
        i_copy++;

    } /* end for while (i < end - conse) */

    /* step 3: processing the end of time series */
    status = step3_processingend(instance, vec_next_P, vec_next_a, clrx, clry, i, level_state_records,
                                 annual_state_records, semi_state_records, third_state_records,
                                 rmse_records, rec_cg, num_curve, &end, conse, starting_date, bl_train_complete,
                                 i_start, prev_i_break, adj_rmse, b_fastmode, fit_cft1, fit_cft2, fit_cft3, fit_cft4,
                                 &clrx_record1, &clrx_record2, &clrx_record3, &clrx_record4, probability_threshold,
                                 valid_count, sum_square_smooth_Q, sum_square_smooth_H,sum_smooth_Q,
                                 sum_smooth_H, sum_square_vt, sum_vt, sum_kalman_coef, temporal_rmse_square_sum,
                                 temporal_rmse_sum, temporal_rmse_count, n_focus_variable, n_total_variable, focus_blist,
                                 min_days_conse);

    if(rec_cg[0].category < 80)
    {
        /* save out state and f info*/
        if(b_fastmode == FALSE)
        {
            for(i_b = 0; i_b < n_total_variable; i_b++)
            {
                sprintf(states_output_dir_full, "%s%d.csv", states_output_dir, i_b + 1);
                fptr_save= fopen(states_output_dir_full, "w");

                fprintf(fptr_save, "%s, %s, %s, %s, %s, %s, %s\n", "Dates", "Level",
                        "Annual", "Semiannual", "Fourmonth", "obs", "rmse");

                int obs_count = i_dense;
                const char* str_nan = "NA";
                for (i = clrx[i_dense] - starting_date; i < rec_cg[*num_curve].t_end + 1 - starting_date; i++)
                {
                    if(i  == clrx[obs_count] - starting_date){
                        fprintf(fptr_save, "%i, %f, %f, %f, %f, %f, %f\n",
                                i+starting_date, level_state_records[i_b][i],
                                annual_state_records[i_b][i], semi_state_records[i_b][i],
                                third_state_records[i_b][i], clry[i_b][obs_count],
                                rmse_records[i_b][obs_count]);
                        obs_count++;
                    }
                    else
                    {
                        fprintf(fptr_save, "%i, %f, %f, %f, %f, %s, %s\n",
                                i+starting_date, level_state_records[i_b][i],
                                annual_state_records[i_b][i], semi_state_records[i_b][i],
                                third_state_records[i_b][i],
                                str_nan, str_nan);
                    }

                }
                fclose(fptr_save);
            }

        }
    }

    if(training_type == TRAINING_TYPE_PARAMETER) /* save magnitude result for training mode*/
    {
        rec_cg[0].magnitude[0] = max_CM_stable;
        rec_cg[0].magnitude[1] = max_CM_dist;
    }

    if(b_assigninstancememory == 1)
    {
        for(i_b = 0; i_b < n_total_variable; i_b++)
        {
            gsl_vector_free(instance[i_b].Z);
            //(&instance[i_b])->Z = NULL;
            //gsl_matrix_free (instance[i_b].H);
            gsl_matrix_free(instance[i_b].T);
            //(&instance[i_b])->T = NULL;
            gsl_matrix_free(instance[i_b].Q);
            //(&instance[i_b])->Q = NULL;
            gsl_vector_free(instance[i_b].a1);
            //(&instance[i_b])->a1 = NULL;
            gsl_matrix_free(instance[i_b].P1);
            //(&instance[i_b])->P1 = NULL;
            gsl_matrix_free(instance[i_b].P1inf);
            //(&instance[i_b])->P1inf = NULL;
        }

        for(i = 0; i < n_total_variable; i ++)
        {
           gsl_vector_free(vec_next_a[i]);
           gsl_matrix_free(vec_next_P[i]);
        }

    }

    //printf("free stage 1 \n");
    free(valid_count);
    free(sum_square_smooth_Q);
    free(sum_square_smooth_H);
    free(sum_smooth_Q);
    //printf("free stage 2 \n");
    free(sum_smooth_H);
    free(sum_square_vt);
    free(sum_vt);
    free(sum_kalman_coef);
    free(var_smooth_Q);
    //printf("free stage 3 \n");
    free(var_smooth_H);
    free(adj_rmse);

    status = free_2d_array ((void **)vec_next_a);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: vec_next_a\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **)vec_next_P);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: vec_next_P\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **) fit_cft1);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fit_cft1\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **) fit_cft2);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fit_cft2\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **) fit_cft3);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fit_cft3\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **) fit_cft4);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fit_cft4\n", FUNC_NAME, FAILURE);
    }


    status = free_2d_array ((void **) tmp_fit_cft);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: tmp_fit_cft\n", FUNC_NAME, FAILURE);
    }

    //printf("free stage 4 \n");
    free(rmse_ini);
    rmse_ini = NULL;


//    status = free_2d_array((void **)instance);
//    if (status != SUCCESS)
//    {
//        RETURN_ERROR("Freeing instance: clry\n", FUNC_NAME, FAILURE);
//    }

    free(clrx);
    clrx = NULL;

    free(clrx_copy);
    clrx_copy = NULL;

    status = free_2d_array((void **)clry);
    if (status != SUCCESS)
    {
        RETURN_ERROR("Freeing memory: clry\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array((void **)clry_copy);
    if (status != SUCCESS)
    {
        RETURN_ERROR("Freeing memory: clry_copy\n", FUNC_NAME, FAILURE);
    }

    free(instance);
    instance = NULL;
    //gsl_matrix_free(initial_Q);

//    status = free_2d_array ((void **) temp_v_dif);
//    if (status != SUCCESS)
//    {
//        RETURN_ERROR ("Freeing memory: temp_v_dif\n", FUNC_NAME, FAILURE);
//    }

    free(sensor_buf_valid);

    status = free_2d_array((void **)rmse_records);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: rmse_records\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **) temporal_rmse_square_sum);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: temporal_rmse_square_sum\n",
                      FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **) temporal_rmse_sum);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: temporal_rmse_square_sum\n",
                      FUNC_NAME, FAILURE);
    }

    free(temporal_rmse_count);

    status = free_2d_array ((void **) yearly_temporal_rmse_square_sum);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: yearly_temporal_rmse_square_sum\n",
                      FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **) yearly_temporal_rmse_sum);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: yearly_temporal_rmse_square_sum\n",
                      FUNC_NAME, FAILURE);
    }

    //printf("free stage 6 \n");
    free(yearly_temporal_rmse_count);

    status = free_2d_array ((void **)level_state_records);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: level_state_records\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **)annual_state_records);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: annual_state_records\n", FUNC_NAME, FAILURE);
    }

    //printf("free stage 5 \n");
    status = free_2d_array ((void **)semi_state_records);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: semi_state_records\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **)third_state_records);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: third_state_records\n", FUNC_NAME, FAILURE);
    }


    status = free_2d_array ((void **) rec_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: rec_v_dif\n",
                      FUNC_NAME, FAILURE);
    }


    // free(rf_feature);

    return SUCCESS;

}

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
)
{
    int n_sn;
    int i, k;
    int end;
    int i_start;
    int i_span;
    int *clrx;                       /* clear pixel curve in X direction (date)    */
    float **clry;                    /* clear pixel curve in Y direction (spectralbands)    */
    double **fit_cft;                 /* Fitted coefficients 2-D array.        */
    double *rmse;                     /* Root Mean Squared Error array.        */
    double **temp_v_dif;              /* for the thermal band.......           */
    char FUNC_NAME[] = "inefficientobs_procedure";
    int n_clr = 0;
    int i_b;
    int status;

    clrx = (int*)malloc(valid_num_scenes * sizeof(int));
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

    temp_v_dif = (double **)allocate_2d_array(TOTAL_IMAGE_BANDS, valid_num_scenes,
                                     sizeof (double));
    if (temp_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating temp_v_dif memory",FUNC_NAME, FAILURE);
    }

    fit_cft = (double **) allocate_2d_array (TOTAL_IMAGE_BANDS, LASSO_COEFFS,
                                         sizeof (double));
    if (fit_cft == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft memory", FUNC_NAME, FAILURE);
    }

    rmse = (double *)malloc(TOTAL_IMAGE_BANDS * sizeof(double));
    if (rmse == NULL)
    {
        RETURN_ERROR ("Allocating rmse memory", FUNC_NAME, FAILURE);
    }

    if (sn_pct > T_SN)
    {
        n_sn = 0;

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
                     clry[k][n_sn] = (double)buf[k][i];
                }
                n_sn++;
            }
        }
        end = n_sn;

        if (n_sn < N_TIMES * MIN_NUM_C) // not enough snow pixels
        {
            RETURN_ERROR ("Not enough good snow observations\n",
                 FUNC_NAME, FAILURE);
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

                if (i_span < MIN_NUM_C * N_TIMES)
                    fit_cft[k][0] = 10000; // fixed value for saturated pixels
                else
                {
                    status = auto_ts_fit(clrx, clry, k, k, 0, i_span-1, MIN_NUM_C,
                             fit_cft, &rmse[k], temp_v_dif);

                    if (status != SUCCESS)
                        RETURN_ERROR ("Calling auto_ts_fit1\n",
                               FUNC_NAME, EXIT_FAILURE);
                }
             }
             else
             {
                for (i = 0; i < end; i++)
                {
                    //printf("thermal: %f\n", clry[k][i]);
                    //printf("clry: %f\n", (double)clry[k][i]);
                    if (clry[k][i] > -9320 && clry[k][i] < 7070)
                    {
                        clrx[i_span] = clrx[i];
                        clry[k][i_span] = clry[k][i];
                        i_span++;
                    }
                 }

                status = auto_ts_fit(clrx, clry, k, k, 0, i_span-1, MIN_NUM_C,
                         fit_cft, &rmse[k], temp_v_dif);

                if (status != SUCCESS)
                    RETURN_ERROR ("Calling auto_ts_fit1\n",
                           FUNC_NAME, EXIT_FAILURE);

              }
        }


        /**********************************************************/
        /*                                                        */
        /*                                                        */
        /**********************************************************/

        rec_cg[*num_curve].t_start = clrx[i_start-1];
        rec_cg[*num_curve].t_end = clrx[end-1];

        /**********************************************************/
        /*                                                        */
        /* No break at the moment.                                */
        /*                                                        */
        /**********************************************************/

        rec_cg[*num_curve].t_break = 0;

        /**********************************************************/
        /*                                                        */
        /* Record postion of the pixel.                           */
        /*                                                        */
        /**********************************************************/

        //rec_cg[*num_curve].pos.row = row;
        //rec_cg[*num_curve].pos.col = col;

        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
        {
            for (k = 0; k < SCCD_MAX_NUM_C; k++)
            {
                /**************************************************/
                /*                                                */
                /* Record fitted coefficients.                    */
                /*                                                */
                /**************************************************/

                rec_cg[*num_curve].coefs[i_b][k] = fit_cft[i_b][k];
            }

            /******************************************************/
            /*                                                    */
            /* Record rmse of the pixel.                          */
            /*                                                    */
            /******************************************************/

            rec_cg[*num_curve].rmse[i_b] = rmse[i_b];
        }

        /**********************************************************/
        /*                                                        */
        /* Record change probability, number of observations.     */
        /*                                                        */
        /**********************************************************/

        rec_cg[*num_curve].num_obs = n_sn;
        rec_cg[*num_curve].t_confirmed = 0;
        rec_cg[*num_curve].change_prob = 0;
        rec_cg[*num_curve].category = 50; /* snow pixel */

        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
        {
            /******************************************************/
            /*                                                    */
            /* Record change magnitude.                           */
            /*                                                    */
            /******************************************************/

            rec_cg[*num_curve].magnitude[i_b] = 0.0;
//            rec_cg[*num_curve].obs_disturb[i_b] = 0.0;
//            for (i = 0; i < SCCD_MAX_NUM_C - 1; i++)
//            {
//                rec_cg[*num_curve].state_disturb[i_b][i] = 0.0;
//            }
        }

        /**********************************************************/
        /*                                                        */
        /* NUM of Fitted Curves (*num_curve).                         */
        /*                                                        */
        /**********************************************************/

//        *num_curve = *num_curve + 1;

//        if (*num_curve >= NUM_FC)
//        {

//            /******************************************************/
//            /*                                                    */
//            /* Reallocate memory for rec_cg.                      */
//            /*                                                    */
//            /******************************************************/

//            rec_cg = realloc(rec_cg, (*num_curve + 1) * sizeof(Output_t));
//            if (rec_cg == NULL)
//            {
//                RETURN_ERROR("ERROR allocating rec_cg memory",
//                             FUNC_NAME, FAILURE);
//            }
//        }
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
                for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                    clry[k][n_clr] = (double)buf[k][i];
                n_clr++;
            }
        }
        end = n_clr;


        n_clr = 0;
        double band2_median; // probably not good practice to declare here....
        quick_sort_double(clry[1], 0, end - 1);
        matlab_2d_double_median(clry, 1, end, &band2_median);

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
            RETURN_ERROR("Not enough good clear observations\n",
                        FUNC_NAME, FAILURE);
        }
        else
        {
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
                status = auto_ts_fit(clrx, clry, i_b, i_b, 0, end-1, MIN_NUM_C,
                                     fit_cft, &rmse[i_b], temp_v_dif);
                if (status != SUCCESS)
                {
                    RETURN_ERROR ("Calling auto_ts_fit for clear persistent pixels\n",
                                  FUNC_NAME, FAILURE);
                }
            }
        }

        /**********************************************************/
        /*                                                        */
        /* Update information at each iteration.                  */
        /* Record time of curve start, time of curve end.         */
        /*                                                        */
        /**********************************************************/

        rec_cg[*num_curve].t_start = clrx[i_start-1];
        rec_cg[*num_curve].t_end = clrx[end-1];

        /**********************************************************/
        /*                                                        */
        /* No break at the moment.                                */
        /*                                                        */
        /**********************************************************/

        rec_cg[*num_curve].t_break = 0;

        /**********************************************************/
        /*                                                        */
        /* Record postion of the pixel.                           */
        /*                                                        */
        /**********************************************************/

        //rec_cg[*num_curve].pos.row = row;
        //rec_cg[*num_curve].pos.col = col;

        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
        {
            for (k = 0; k < MIN_NUM_C; k++)
            {
                /**************************************************/
                /*                                                */
                /* Record fitted coefficients.                    */
                /*                                                */
                /**************************************************/

                rec_cg[*num_curve].coefs[i_b][k] = fit_cft[i_b][k];
            }

            /******************************************************/
            /*                                                    */
            /* Record rmse of the pixel.                          */
            /*                                                    */
            /******************************************************/

            rec_cg[*num_curve].rmse[i_b] = rmse[i_b];
        }

        /**********************************************************/
        /*                                                        */
        /* Record change probability, number of observations,     */
        /* fit category.                                          */
        /*                                                        */
        /**********************************************************/
        rec_cg[*num_curve].num_obs = n_clr;
        rec_cg[*num_curve].t_confirmed = 0;
        rec_cg[*num_curve].change_prob = 0;
        rec_cg[*num_curve].category = 40;

        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
        {
            /******************************************************/
            /*                                                    */
            /* Record change magnitude.                           */
            /*                                                    */
            /******************************************************/

            rec_cg[*num_curve].magnitude[i_b] = 0.0;
//            rec_cg[*num_curve].obs_disturb[i_b] = 0.0;
//            for (i = 0; i < SCCD_MAX_NUM_C - 1; i++)
//            {
//                rec_cg[*num_curve].state_disturb[i_b][i] = 0.0;
//            }
        }

        /**********************************************************/
        /*                                                        */
        /* NUM of Fitted Curves (*num_curve).                         */
        /*                                                        */
        /**********************************************************/

//        *num_curve = *num_curve + 1;

//        if (*num_curve >= NUM_FC)
//        {
//            /******************************************************/
//            /*                                                    */
//            /* Reallocate memory for rec_cg.                      */
//            /*                                                    */
//            /******************************************************/

//            rec_cg = realloc(rec_cg, (*num_curve + 1) * sizeof(Output_t));
//            if (rec_cg == NULL)
//            {
//                RETURN_ERROR("ERROR allocating rec_cg memory",
//                             FUNC_NAME, FAILURE);
//            }
//        }
    }
    // *num_curve = *num_curve + 1;

    free(clrx);
    clrx = NULL;

    status = free_2d_array((void **) clry);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: clry\n",
                      FUNC_NAME, FAILURE);
    }

    free(rmse);
    rmse = NULL;

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
)
{
    double* H_pool;
    double* q_pool;
    int n_H_pool;
    int i, j, k;

    int n_q_pool = 3;
    double max_lik = -999999;
    char FUNC_NAME[] = "grid_searching";
    double tmp_q;
    double tmp_lik;
    int record_h;
    int record_q;
    double tmp_f_rmse = 0;
    double tmp_v_rmse = 0;
    double ini_p;

    n_H_pool = 3;
    H_pool = malloc(n_H_pool * sizeof(double));
    if (H_pool == NULL)
    {
        RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
    }

    ini_p = pow(ini_a , 2) * 0.005 / (instance ->m + 1) * 2;
    if(ini_p < rmse * rmse)
        H_pool[0] = rmse * rmse - ini_p;
    else
        H_pool[0] = rmse * rmse;

    ini_p = pow(ini_a , 2) * 0.01 / (instance ->m + 1) * 2;
    if(ini_p < rmse * rmse)
        H_pool[1] = rmse * rmse - ini_p;
    else
        H_pool[1] = rmse * rmse;

    ini_p = pow(ini_a , 2) * 0.015 / (instance ->m + 1) * 2;
    if(ini_p < rmse * rmse)
        H_pool[2] = rmse * rmse - ini_p;
    else
        H_pool[2] = rmse * rmse;

    //H_pool[3] = (mean_y * 0.125) * (mean_y * 0.125);
   //H_pool[4] = (mean_y * 0.15) * (mean_y * 0.15);
   // H_pool[3] = (mean_y * 0.125) * (mean_y * 0.125);

//    if (rmse > mean_y * 0.125)
//    {
//        n_H_pool = 4;
//        H_pool = malloc(n_H_pool * sizeof(double));
//        if (H_pool == NULL)
//        {
//            RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
//        }
//        H_pool[0] = (mean_y * 0.05) * (mean_y * 0.05);
//        H_pool[1] = (mean_y * 0.075) * (mean_y * 0.075);
//        H_pool[2] = (mean_y * 0.1) * (mean_y * 0.1);
//        H_pool[3] = (mean_y * 0.125) * (mean_y * 0.125);
//    }
//    if (rmse > mean_y * 0.1)
//    {
//        n_H_pool = 3;
//        H_pool = malloc(n_H_pool * sizeof(double));
//        if (H_pool == NULL)
//        {
//            RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
//        }
//        H_pool[0] = (mean_y * 0.05) * (mean_y * 0.05);
//        H_pool[1] = (mean_y * 0.075) * (mean_y * 0.075);
//        H_pool[2] = (mean_y * 0.1) * (mean_y * 0.1);
//    }
//    else if (rmse > mean_y * 0.075)
//    {
//        n_H_pool = 2;
//        H_pool = malloc(n_H_pool * sizeof(double));
//        if (H_pool == NULL)
//        {
//            RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
//        }
//        H_pool[0] = (mean_y * 0.05) * (mean_y * 0.05);
//        H_pool[1] = (mean_y * 0.075) * (mean_y * 0.075);
//    }
//    else if (rmse > mean_y * 0.05 )
//    {
//        n_H_pool = 1;
//        H_pool = malloc(n_H_pool * sizeof(double));
//        if (H_pool == NULL)
//        {
//            RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
//        }
//        H_pool[0] = (mean_y * 0.05 ) * (mean_y * 0.05 );
//    }
//    else // rmse < mean_y * 0.05
//    {
//        *best_h = rmse * rmse;
//        *best_q = INI_Q00;
//        return SUCCESS;
//    }
//    n_H_pool = 4;
//    H_pool = malloc(n_H_pool * sizeof(double));
//    if (H_pool == NULL)
//    {
//        RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
//    }
//    H_pool[0] = rmse * rmse * 0.8;
//    H_pool[1] = rmse * rmse * 0.9;
//    H_pool[2] = rmse * rmse * 0.99;
//    H_pool[3] = rmse * rmse * 0.999;

    n_q_pool = 1;
    q_pool = malloc(n_q_pool * sizeof(double));
    if (q_pool == NULL)
    {
        RETURN_ERROR("ERROR allocating q_pool memory", FUNC_NAME, FAILURE);
    }
    q_pool[0] = *best_q;
//    q_pool[1] = rmse * rmse * 0.01;
//    q_pool[2] = rmse * rmse * 0;

//    n_H_pool = 1;
//    H_pool = malloc(n_H_pool * sizeof(double));
//    if (H_pool == NULL)
//    {
//        RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
//    }
//    H_pool[0] = *best_h;

    for (i = 0 ; i < n_H_pool; i++)
    {
        for (j = 0; j < n_q_pool; j++)
        {
            ini_p = pow(ini_a , 2) * (i + 1) * 0.005;

            for(k = 0; k < instance->m; k++)
            {
                gsl_matrix_set(instance->P1, k, k, ini_p);
            }

            // tmp_q =  (rmse * rmse - H_pool[i]) / (pow(2.0, j+1)) / ((actual_m + 1)/2 * interval);
            // tmp_q =  (H_pool[i]) / (pow(2.0, j + 1)) / ((actual_m + 1)/2 * interval);
            //tmp_q = (rmse * rmse - H_pool[i])/ ((actual_m + 1)/2 * interval);;
            tmp_lik = ssmloglik_gridsearching(instance, H_pool[i], q_pool[j], &tmp_f_rmse, &tmp_f_rmse);
            if (tmp_lik > max_lik)
            {
                max_lik = tmp_lik;
                *best_h = H_pool[i];
                *best_q = q_pool[j];
                *best_p = ini_p;
                record_h = i;
                record_q = j;
            }
        }
    }

    //printf("best q is %d \n", record_q);

    free(H_pool);
    free(q_pool);
    return SUCCESS;


}

int Jazwinski_searching
(
     ssmodel *instance,
     double rmse,
     double mean_y,
     double *best_h,
     double *best_q,
     double interval,
     int actual_m
)
{
    double* H_pool;
    int n_H_pool;
    int i;

    double max_lik = -999999;
    char FUNC_NAME[] = "grid_searching";
    double tmp_q;
    double tmp_lik;
    double tmp_f_rmse = 0;
    double tmp_v_rmse = 0;

//    n_H_pool = 4;
//    H_pool = malloc(n_H_pool * sizeof(double));
//    if (H_pool == NULL)
//    {
//        RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
//    }
//    H_pool[0] = (mean_y * 0.05) * (mean_y * 0.05);
//    H_pool[1] = (mean_y * 0.075) * (mean_y * 0.075);
//    H_pool[2] = (mean_y * 0.1) * (mean_y * 0.1);
//    H_pool[3] = (mean_y * 0.125) * (mean_y * 0.125);
   //H_pool[4] = (mean_y * 0.15) * (mean_y * 0.15);
   // H_pool[3] = (mean_y * 0.125) * (mean_y * 0.125);

//    if (rmse > mean_y * 0.125)
//    {
//        n_H_pool = 4;
//        H_pool = malloc(n_H_pool * sizeof(double));
//        if (H_pool == NULL)
//        {
//            RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
//        }
//        H_pool[0] = (mean_y * 0.05) * (mean_y * 0.05);
//        H_pool[1] = (mean_y * 0.075) * (mean_y * 0.075);
//        H_pool[2] = (mean_y * 0.1) * (mean_y * 0.1);
//        H_pool[3] = (mean_y * 0.125) * (mean_y * 0.125);
//    }
//    if (rmse > mean_y * 0.25)
//    {
//        n_H_pool = 5;
//        H_pool = malloc(n_H_pool * sizeof(double));
//        if (H_pool == NULL)
//        {
//            RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
//        }
//        H_pool[0] = (mean_y * 0.05) * (mean_y * 0.05);
//        H_pool[1] = (mean_y * 0.1) * (mean_y * 0.1);
//        H_pool[2] = (mean_y * 0.15) * (mean_y * 0.15);
//        H_pool[3] = (mean_y * 0.2) * (mean_y * 0.2);
//        H_pool[3] = (mean_y * 0.25) * (mean_y * 0.5);


//    }
//    else if (rmse > mean_y * 0.2)
//    {
//        n_H_pool = 4;
//        H_pool = malloc(n_H_pool * sizeof(double));
//        if (H_pool == NULL)
//        {
//            RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
//        }
//        H_pool[0] = (mean_y * 0.05) * (mean_y * 0.05);
//        H_pool[1] = (mean_y * 0.1) * (mean_y * 0.1);
//        H_pool[2] = (mean_y * 0.15) * (mean_y * 0.15);
//        H_pool[3] = (mean_y * 0.2) * (mean_y * 0.2);

//    }
    if (rmse > mean_y * 0.2)
    {
        n_H_pool = 4;
        H_pool = malloc(n_H_pool * sizeof(double));
        if (H_pool == NULL)
        {
            RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
        }
        H_pool[0] = (mean_y * 0.05) * (mean_y * 0.05);
        H_pool[1] = (mean_y * 0.1) * (mean_y * 0.1);
        H_pool[2] = (mean_y * 0.15) * (mean_y * 0.15);
        H_pool[3] = (mean_y * 0.2) * (mean_y * 0.2);

    }
    else if (rmse > mean_y * 0.15)
    {
        n_H_pool = 3;
        H_pool = malloc(n_H_pool * sizeof(double));
        if (H_pool == NULL)
        {
            RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
        }
        H_pool[0] = (mean_y * 0.05) * (mean_y * 0.05);
        H_pool[1] = (mean_y * 0.1) * (mean_y * 0.1);
        H_pool[2] = (mean_y * 0.15) * (mean_y * 0.15);

    }
    else if (rmse > mean_y * 0.1)
    {
        n_H_pool = 2;
        H_pool = malloc(n_H_pool * sizeof(double));
        if (H_pool == NULL)
        {
            RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
        }
        H_pool[0] = (mean_y * 0.05) * (mean_y * 0.05);
        H_pool[1] = (mean_y * 0.1) * (mean_y * 0.1);
    }
    else if (rmse > mean_y * 0.05)
    {
        n_H_pool = 1;
        H_pool = malloc(n_H_pool * sizeof(double));
        if (H_pool == NULL)
        {
            RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
        }
        H_pool[0] = (mean_y * 0.05) * (mean_y * 0.05);
    }
    else // rmse < mean_y * 0.05
    {
        *best_h = rmse * rmse;
        *best_q = INI_Q00;
        return SUCCESS;
    }
//    n_H_pool = 1;
//    H_pool = malloc(n_H_pool * sizeof(double));
//    if (H_pool == NULL)
//    {
//        RETURN_ERROR("ERROR allocating H_pool memory", FUNC_NAME, FAILURE);
//    }
//    H_pool[0] = (mean_y * 0.05 + 50) * (mean_y * 0.05 + 50);

    for (i = 0 ; i < n_H_pool; i++)
    {
        tmp_f_rmse = 0;
        tmp_q = 0;
        tmp_lik = ssmloglik_gridsearching(instance, H_pool[i], tmp_q,  &tmp_f_rmse, &tmp_v_rmse);
        if( rmse * rmse > tmp_f_rmse)
            tmp_q = ( rmse *rmse - tmp_f_rmse) / ((actual_m + 1)/2  * interval);
        else
            tmp_q = 0;

        tmp_f_rmse = 0;
        tmp_lik = ssmloglik_gridsearching(instance, H_pool[i], tmp_q,  &tmp_f_rmse, &tmp_v_rmse);

        if (tmp_lik > max_lik)
        {
            max_lik = tmp_lik;
            *best_h = H_pool[i];
            *best_q = tmp_q;
            //record_h = i;
        }

    }

    free(H_pool);

    return SUCCESS;


}

int getcategorycurve
(
    Output_t_sccd* t,
    int num_curve
)
{
    double green_evidence = 0;
    double magnitude_mean = 0;
//    if(t[num_curve].magnitude[2] < - S_TCG_single)
//        green_evidence = green_evidence + 1;
//    if(t[num_curve].magnitude[3] > S_TCG_single)
//        green_evidence = green_evidence + 1;
//    if(t[num_curve].magnitude[4] < - S_TCG_single)
//        green_evidence = green_evidence + 1;
    int i;
//    double dist_threshold = 2;
    bool dist_flag = FALSE;
//    for(i = 0; i < n_focus_variable; i++)
//    {
//        if(t[num_curve].magnitude[focus_blist[i]] < dist_threshold)
//        {
//            dist_flag = FALSE;
//            break;
//        }
//    }

    if (dist_flag == TRUE)
    {
        return DIST_TYPE_NOISE;
    }
    else
    {
        green_evidence =  - t[num_curve].magnitude[2] +  t[num_curve].magnitude[3] - t[num_curve].magnitude[4];
        //green_evidence =  t[num_curve].magnitude[7]  - t[num_curve].magnitude[2] +  t[num_curve].magnitude[3] - t[num_curve].magnitude[4];
        // magnitude_mean = (fabsf(t[num_curve].magnitude[2]) +  fabsf(t[num_curve].magnitude[3]) + fabsf(t[num_curve].magnitude[4])) / 3.0;
        if(green_evidence > 0)
        // if((t[num_curve].magnitude[2]< - S_TCG_single) && (t[num_curve].magnitude[3] > S_TCG_single) && (t[num_curve].magnitude[4] < -S_TCG_single))
        {
            if((t[num_curve + 1].coefs[3][1] > fabsf(t[num_curve].coefs[3][1])) &&
                    (t[num_curve + 1].coefs[2][1] < -fabsf(t[num_curve].coefs[2][1])) &&
                    (t[num_curve + 1].coefs[4][1] < -fabsf(t[num_curve].coefs[4][1])))
                return DIST_TYPE_RESTOR;
            else
                return DIST_TYPE_REGROW;
        }
        else
            return DIST_TYPE_DIST;
    }

}


int getcategorycurve_old
(
    Output_t_sccd* t,
    int num_curve
)
{
    if((t[num_curve].magnitude[2]< - S_TCG_single) && (t[num_curve].magnitude[3] > S_TCG_single) && (t[num_curve].magnitude[4] < -S_TCG_single))
    {
        if((t[num_curve + 1].coefs[3][1] > fabsf(t[num_curve].coefs[3][1])) &&
                (t[num_curve + 1].coefs[2][1] < -fabsf(t[num_curve].coefs[2][1])) &&
                (t[num_curve + 1].coefs[4][1] < -fabsf(t[num_curve].coefs[4][1])))
            return DIST_TYPE_RESTOR;
        else
            return DIST_TYPE_REGROW;
    }
    else
        return DIST_TYPE_DIST;

}

