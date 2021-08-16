/* c translation for KFAS R package */
#include <math.h>
#include "KFAS.h"
#include "lbfgs.h"
#include "const.h"
#include "defines.h"

/* diffuse filtering for single time point */
void dfilter1step
(
    int cur_monitor_date,
    int next_valid_date, /*I/O*/
    int* valid_count,
    double yt,
    gsl_vector* zt,
    double ht,
    gsl_matrix* tt,
    gsl_matrix* rqr,
    gsl_vector* at,    /*I/O*/
    gsl_matrix* pt,    /*I/O*/
    double* vt,     /*I/O*/
    double* ft,    /*I/O*/
    gsl_vector* kt,    /*I/O*/
    gsl_matrix* pinf,  /*I/O*/
    double* finf,    /*I/O*/
    gsl_vector* kinf,  /*I/O*/
    int* rankp,        /*I/O*/
    double* lik,       /*I/O*/
    int m,
    bool fast_mode,
    double* fit_cft,
    int structure
)
{

    double finv, tmp;

    gsl_vector* ahelp;
    gsl_matrix* mm;
    ahelp = gsl_vector_alloc(m);
    mm = gsl_matrix_alloc(m, m);
    //memcpy(zt_sub, &zt[1], 2*sizeof(*a));

    /*kt = pt*zt*/
    gsl_blas_dsymv(CblasUpper, 1.0, pt,
                   zt, 0.0, kt);

    /*ft = kt *ztt + ht*/
    gsl_blas_ddot(zt, kt, ft);
    *ft = *ft +ht;

    if (next_valid_date == cur_monitor_date)
    {
        if(TRUE == fast_mode)
            fit_cft2vec_a(fit_cft, at, cur_monitor_date, m, structure);
        /*finf = zt*pinf*zt'*/
        gsl_blas_dsymv(CblasUpper, 1.0, pinf,
                       zt, 0.0, kinf);
        gsl_blas_ddot(zt, kinf, finf);

        /*vt = yt - zt*at*/
        gsl_blas_ddot(zt, at, &tmp);
        *vt = yt - tmp;

        if(*finf > KFAS_TOL)
        {
            finv = 1.0 / *finf;
            gsl_blas_daxpy((*vt) * finv, kinf, at);

            gsl_blas_dsyr(CblasUpper, *ft * pow(finv, 2),
                          kinf, pt);
            gsl_blas_dsyr2(CblasUpper, -finv, kt, kinf, pt);
            gsl_blas_dsyr(CblasUpper, -finv, kinf, pinf);

            *lik = *lik - 0.5 * log(*finf);

            *rankp = *rankp - 1;

        }
        else
        {
            *finf = 0.0;
            if(*ft > KFAS_TOL)
            {
                finv = 1.0 / *ft;
                gsl_blas_daxpy((*vt) * finv, kt, at);
                gsl_blas_dsyr(CblasUpper, -finv, kt, pt);
                *lik = *lik - KFAS_LIK_C - 0.5 * (log(*ft) + pow(*vt, 2) * finv);
            }
        }

        if(*ft < KFAS_TOL)
            *ft = 0.0;
        if(rankp == 0)
            return;


        (*valid_count) ++;
    }

    if(FALSE == fast_mode)
    {
        /* ahelp = 1.0 * tt* at */
        gsl_blas_dgemv(CblasNoTrans,1.0, tt, at, 0.0, ahelp);

        gsl_vector_memcpy(at, ahelp);
    }




    /* mm = pt*tt */
    gsl_blas_dsymm(CblasRight, CblasUpper, 1.0, pt, tt, 0.0, mm);

    /* pt = mm * tt^T */
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mm, tt, 0.0, pt);

    gsl_matrix_add(pt, rqr);

    /* mm = pinf*tt */
    gsl_blas_dsymm(CblasRight, CblasUpper, 1.0, pinf, tt, 0.0, mm);

    /* pinf = mm * tt^T */
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mm, tt, 0.0, pinf);

    gsl_vector_free(ahelp);
    gsl_matrix_free(mm);

    return;

}

/* non-diffuse filtering for single time point */
void filter1step
(
    int cur_monitor_date,  /*I */
    int next_valid_date,  /*I */
    int* valid_count,
    double yt,             /*I */
    gsl_vector* zt,        /*I */
    double ht,             /*I */
    gsl_matrix* tt,        /*I */
    gsl_matrix* rqr,     /*I */
    gsl_vector* at,    /*I/O*/
    gsl_matrix* pt,    /*I/O*/
    double* vt,     /*I/O*/
    double* ft,    /*O*/
    gsl_vector* kt,    /*O*/
    double* lik,       /*I/O*/
    int m,
    bool fast_mode,
    double* fit_cft,
    int structure
)
{
    double finv;
    double tmp;
    gsl_vector* ahelp;
    gsl_matrix* mm;

    ahelp = gsl_vector_alloc(m);
    mm = gsl_matrix_alloc(m, m);
    //memcpy(zt_sub, &zt[1], 2*sizeof(*a));



    if (next_valid_date == cur_monitor_date)
    {
        if(TRUE == fast_mode)
            fit_cft2vec_a(fit_cft, at, cur_monitor_date, m, structure);
        //printf("i = %d: %f, %f, %f\n", cur_monitor_date, gsl_vector_get(at, 0), gsl_vector_get(at, 1), gsl_vector_get(at, 3));
        /*kt = zt*pt*/
        gsl_blas_dsymv(CblasUpper, 1.0, pt,
                       zt, 0.0, kt);

        /*ft = kt *ztt + ht*/
        gsl_blas_ddot(zt, kt, ft);
        *ft = *ft + ht;

        /*vt = yt - zt*at*/
//        for(i = 0; i < m; i++)
//        {
//                printf("%.10g\n", gsl_vector_get(at, i));

//        }
        gsl_blas_ddot(zt, at, &tmp);
        *vt = yt - tmp;

        if(*ft > KFAS_TOL)
        {
            finv = 1.0 / *ft;

            //at = at + (*vt) * finv * kt;
            gsl_blas_daxpy((*vt) * finv, kt, at);

            //pt = pt + kt * -finv * kt;
            gsl_blas_dsyr(CblasUpper, -finv,
                          kt, pt);

            *lik = *lik - KFAS_LIK_C - 0.5 * (log(*ft) + pow(*vt, 2) * finv);

        }
        else
        {
            *ft = 0.0;
        }
        (*valid_count) ++;

        /* update fit_cft */
        if(TRUE == fast_mode)
            vec_a2fit_cft(at, fit_cft, cur_monitor_date, m, structure);

    }

    if(FALSE == fast_mode)
    {
        /* ahelp = 1.0 * tt* at */
        gsl_blas_dgemv(CblasNoTrans,1.0, tt, at, 0.0, ahelp);

        gsl_vector_memcpy(at, ahelp);
    }



    /* mm = tt*pt*/
    gsl_blas_dsymm(CblasRight, CblasUpper, 1.0, pt, tt, 0.0, mm);

//    for(i = 0; i < m; i++)
//    {
//        for( j = 0; j < m; j ++)
//        {
//            printf("%.10g ", gsl_matrix_get(pt, i, j));
//            if(j == m - 1)
//                printf("\n");
//        }
//     }



    /* pt = mm * tt^T */
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mm, tt, 0.0, pt);

//    for(i = 0; i < m; i++)
//    {
//        for( j = 0; j < m; j ++)
//        {
//            printf("%.10g ", gsl_matrix_get(pt, i, j));
//            if(j == m - 1)
//                printf("\n");
//        }
//     }

    gsl_matrix_add(pt, rqr);

//    for(i = 0; i < m; i++)
//    {
//        for( j = 0; j < m; j ++)
//        {
//            printf("%.10g ", gsl_matrix_get(pt, i, j));
//            if(j == m - 1)
//                printf("\n");
//        }
//     }

    gsl_vector_free(ahelp);
    gsl_matrix_free(mm);
    return;

}

/* non-diffuse filtering for prediction */
void filter1step_missingobs
(
    gsl_vector* zt,        /*I */
    double ht,             /*I */
    gsl_matrix* tt,        /*I */
    gsl_matrix* rqr,     /*I */
    gsl_vector* at,    /*I/O: */
    double* ft,         /*O*/
    gsl_matrix* pt,    /*I/O*/
    gsl_vector* kt,    /*I/O*/
    int m,              /*I*/
    bool fast_mode
)
{
    gsl_vector* ahelp;
    gsl_matrix* mm;

    ahelp = gsl_vector_alloc(m);
    mm = gsl_matrix_alloc(m, m);
    //memcpy(zt_sub, &zt[1], 2*sizeof(*a));

    if (fast_mode == FALSE)
    {
        /* at = 1.0 * tt* at */
        gsl_blas_dgemv(CblasNoTrans,1.0, tt, at, 0.0, ahelp);

        gsl_vector_memcpy(at, ahelp);
    }



    /* mm = tt*pt*/
    gsl_blas_dsymm(CblasRight, CblasUpper, 1.0, pt, tt, 0.0, mm);

    /* pt = mm * tt^T */
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mm, tt, 0.0, pt);

    gsl_matrix_add(pt, rqr);


    gsl_vector_free(ahelp);
    gsl_matrix_free(mm);
    return;

}

void filter1step_missingobs_noq
(
    gsl_vector* zt,        /*I */
    double ht,             /*I */
    gsl_matrix* tt,        /*I */
    gsl_matrix* rqr,     /*I */
    gsl_vector* at,    /*I/O: */
    double* ft,         /*O*/
    gsl_matrix* pt,    /*I/O*/
    gsl_vector* kt,    /*I/O*/
    int m,              /*I*/
    bool fast_mode
)
{
    gsl_vector* ahelp;
    gsl_matrix* mm;

    ahelp = gsl_vector_alloc(m);
    mm = gsl_matrix_alloc(m, m);
    //memcpy(zt_sub, &zt[1], 2*sizeof(*a));

    if (fast_mode == FALSE)
    {
        /* at = 1.0 * tt* at */
        gsl_blas_dgemv(CblasNoTrans,1.0, tt, at, 0.0, ahelp);

        gsl_vector_memcpy(at, ahelp);
    }



    /* mm = tt*pt*/
    gsl_blas_dsymm(CblasRight, CblasUpper, 1.0, pt, tt, 0.0, mm);

    /* pt = mm * tt^T */
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mm, tt, 0.0, pt);

    gsl_vector_free(ahelp);
    gsl_matrix_free(mm);
    return;

}

void filter1step_validobs
(
    double yt,             /*I */
    gsl_vector* zt,        /*I */
    double *ht,             /*I */
    gsl_matrix* tt,        /*I */
    gsl_matrix* rqr,     /*I */
    gsl_vector* at,    /*I/O*/
    gsl_matrix* pt,    /*I/O*/
    double* vt,     /*I/O*/
    double* ft,    /*I/O*/
    gsl_vector* kt,    /*I/O*/
    int m,
    bool fast_mode,
    bool update_H,
    gsl_vector* att
)
{
    double finv;
    double tmp;
    double p;
    gsl_vector* ahelp;
    gsl_matrix* mm;

    ahelp = gsl_vector_alloc(m);
    mm = gsl_matrix_alloc(m, m);

    /*kt = zt*pt*/
    gsl_blas_dsymv(CblasUpper, 1.0, pt,
                   zt, 0.0, kt);

    /*ft = kt *ztt + ht*/
    gsl_blas_ddot(zt, kt, &p);
    *ft = p + *ht;


    gsl_blas_ddot(zt, at, &tmp);
    *vt = yt - tmp;

    if(*ft > KFAS_TOL)
    {
        finv = 1.0 / *ft;
        gsl_blas_daxpy((*vt) * finv, kt, at);

        //at = at + (*vt) * finv * kt;

        // pt = pt - kt * ktt * Finv
        gsl_blas_dsyr(CblasUpper, -finv,
                      kt, pt);
        // Joseph's form covariance update
        // pt = [I - k * z] * pt * [I - k * z]' + k * HT\ * kt
    }
    else
    {
        *ft = 0.0;
    }

    gsl_vector_memcpy(att, at);
//    if(fast_mode == FALSE)
//    {
        /* ahelp = 1.0 * tt* at */
        gsl_blas_dgemv(CblasNoTrans,1.0, tt, at, 0.0, ahelp);

        gsl_vector_memcpy(at, ahelp);
//    }

    if(update_H == TRUE)
    {
        //gsl_blas_ddot(zt, at, &tmp);
        //tmp_h = yt - tmp;
        if((*vt) * (*vt) > p)
            *ht = (1 - SMOOTH_FACTOR) * (*ht) + SMOOTH_FACTOR * ((*vt) * (*vt) - p);
//        *ht = (1 - SMOOTH_FACTOR) * (*ht) + SMOOTH_FACTOR * tmp_h * tmp_h;
    }

////    /* mm = tt*pt*/
    gsl_blas_dsymm(CblasRight, CblasUpper, 1.0, pt, tt, 0.0, mm);

    /* pt = mm * tt^T */
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mm, tt, 0.0, pt);


    gsl_matrix_add(pt, rqr);


    gsl_vector_free(ahelp);
    gsl_matrix_free(mm);

    return;
}

void dfilter1step_validobs
(
    double yt,             /*I */
    gsl_vector* zt,        /*I */
    double ht,             /*I */
    gsl_matrix* tt,        /*I */
    gsl_matrix* rqr,     /*I */
    gsl_vector* at,    /*I/O*/
    gsl_matrix* pt,    /*I/O*/
    double* vt,     /*I/O*/
    double* ft,    /*I/O*/
    gsl_vector* kt,    /*I/O*/
    gsl_matrix* pinf,  /*I/O*/
    double* finf,    /*I/O*/
    gsl_vector* kinf,  /*I/O*/
    int* rankp,        /*I/O*/
    int m,
    bool fast_mode
)
{
    double finv;
    double tmp;
    double p;
    gsl_vector* ahelp;
    gsl_matrix* mm;

    ahelp = gsl_vector_alloc(m);
    mm = gsl_matrix_alloc(m, m);

    /*kt = zt*pt*/
    gsl_blas_dsymv(CblasUpper, 1.0, pt,
                   zt, 0.0, kt);

    /*ft = kt *ztt + ht*/
    gsl_blas_ddot(zt, kt, &p);
    *ft = p + ht;

    /*finf = zt*pinf*zt'*/
    gsl_blas_dsymv(CblasUpper, 1.0, pinf,
                   zt, 0.0, kinf);
    gsl_blas_ddot(zt, kinf, finf);

    gsl_blas_ddot(zt, at, &tmp);
    *vt = yt - tmp;

    if(*finf > KFAS_TOL)
    {
        finv = 1.0 / *finf;
        gsl_blas_daxpy((*vt) * finv, kinf, at);
        //at = at + (*vt) * finv * kt;
        gsl_blas_dsyr(CblasUpper, *ft * pow(finv, 2),
                      kinf, pt);
        gsl_blas_dsyr2(CblasUpper, -finv, kt, kinf, pt);
        gsl_blas_dsyr(CblasUpper, -finv, kinf, pinf);
        *rankp = *rankp - 1;
    }
    else
    {
        *finf = 0.0;
        if(*ft > KFAS_TOL)
        {
            finv = 1.0 / *ft;
            gsl_blas_daxpy((*vt) * finv, kt, at);
            gsl_blas_dsyr(CblasUpper, -finv, kt, pt);
        }
    }

    if(*ft < KFAS_TOL)
        *ft = 0.0;

    if(*rankp == 0)
    {
        gsl_vector_free(ahelp);
        gsl_matrix_free(mm);
        return;
    }

    /* ahelp = 1.0 * tt* at */
    gsl_blas_dgemv(CblasNoTrans,1.0, tt, at, 0.0, ahelp);

    gsl_vector_memcpy(at, ahelp);

    /* mm = tt*pt*/
    gsl_blas_dsymm(CblasRight, CblasUpper, 1.0, pt, tt, 0.0, mm);

    /* pt = mm * tt^T */
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mm, tt, 0.0, pt);
    gsl_matrix_add(pt, rqr);

    /* mm = pinf*tt */
    gsl_blas_dsymm(CblasRight, CblasUpper, 1.0, pinf, tt, 0.0, mm);

    /* pinf = mm * tt^T */
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mm, tt, 0.0, pinf);

    gsl_vector_free(ahelp);
    gsl_matrix_free(mm);

    return;
}

void dfilter1step_missingobs
(
    gsl_vector* zt,        /*I */
    double ht,             /*I */
    gsl_matrix* tt,        /*I */
    gsl_matrix* rqr,     /*I */
    gsl_vector* at,    /*I/O: */
    double* ft,         /*O*/
    gsl_matrix* pt,    /*I/O*/
    gsl_vector* kt,    /*I/O*/
    double* finf,    /*I/O*/
    gsl_matrix* pinf,  /*I/O*/
    gsl_vector* kinf,  /*I/O*/
    int m,              /*I*/
    bool fast_mode
)
{
    gsl_vector* ahelp;
    gsl_matrix* mm;

    ahelp = gsl_vector_alloc(m);
    mm = gsl_matrix_alloc(m, m);
    //memcpy(zt_sub, &zt[1], 2*sizeof(*a));

    if (fast_mode == FALSE)
    {
        /* at = 1.0 * tt* at */
        gsl_blas_dgemv(CblasNoTrans,1.0, tt, at, 0.0, ahelp);

        gsl_vector_memcpy(at, ahelp);
    }


    /* mm = tt*pt*/
    gsl_blas_dsymm(CblasRight, CblasUpper, 1.0, pt, tt, 0.0, mm);

    /* pt = mm * tt^T */
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mm, tt, 0.0, pt);

    gsl_matrix_add(pt, rqr);

    /* mm = pinf*tt */
    gsl_blas_dsymm(CblasRight, CblasUpper, 1.0, pinf, tt, 0.0, mm);

    /* pinf = mm * tt^T */
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mm, tt, 0.0, pinf);

    gsl_vector_free(ahelp);
    gsl_matrix_free(mm);
    return;

}

lbfgsfloatval_t evaluate
(
    ssmodel *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
)
{
    int i, j;
    lbfgsfloatval_t ndeps = 1e-3;
    lbfgsfloatval_t *x_g1; // for finite gradient 1
    lbfgsfloatval_t *x_g2; // for finite gradient 2
    lbfgsfloatval_t lik;

    lik = ssmloglik(instance, x);

    x_g1 =  lbfgs_malloc(SSM_OPTIMVAR);;

    x_g2 =  lbfgs_malloc(SSM_OPTIMVAR);;

    for (i = 0; i < SSM_OPTIMVAR; i++)
    {
        for(j = 0; j < SSM_OPTIMVAR; j++)
        {
            if(i == j)
            {
                x_g1[j] = x[j] - ndeps;
                x_g2[j] = x[j] + ndeps;
            }
            else
            {
                x_g1[j] = x[j];
                x_g2[j] = x[j];
            }

        }

        g[i] = (ssmloglik(instance, x_g2) -
                ssmloglik(instance, x_g1))/(2 * ndeps);
        //printf("x[%i]=%f, g[%i] = %f\n", i, x[i], i , g[i]);
    }



    lbfgs_free(x_g1);
    lbfgs_free(x_g2);

    return lik;

}

lbfgsfloatval_t evaluate_bounds
(
    ssmodel *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
)
{
    int i, j;
    lbfgsfloatval_t ndeps = 1e-3;
    lbfgsfloatval_t *x_g1; // for finite gradient 1
    lbfgsfloatval_t *x_g2; // for finite gradient 2
    lbfgsfloatval_t lik;

    lik = ssmloglik(instance, x);

    x_g1 =  lbfgs_malloc(SSM_OPTIMVAR);;

    x_g2 =  lbfgs_malloc(SSM_OPTIMVAR);;

    for (i = 0; i < SSM_OPTIMVAR; i++)
    {
        for(j = 0; j < SSM_OPTIMVAR; j++)
        {
            if(i == j)
            {
                x_g1[j] = x[j] - ndeps;
                x_g2[j] = x[j] + ndeps;
            }
            else
            {
                x_g1[j] = x[j];
                x_g2[j] = x[j];
            }

        }

        g[i] = (ssmloglik_bounds(instance, x_g2, instance->q_lower, instance->q_upper) -
                ssmloglik_bounds(instance, x_g1, instance->q_lower, instance->q_upper))/(2 * ndeps);
        //printf("x[%i]=%f, g[%i] = %f\n", i, x[i], i , g[i]);
    }



    lbfgs_free(x_g1);
    lbfgs_free(x_g2);

    return lik;

}


lbfgsfloatval_t evaluate_2bounds
(
    ssmodel *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
)
{
    int i, j;
    lbfgsfloatval_t ndeps = 1e-3;
    lbfgsfloatval_t *x_g1; // for finite gradient 1
    lbfgsfloatval_t *x_g2; // for finite gradient 2
    lbfgsfloatval_t lik;

    lik = ssmloglik(instance, x);

    x_g1 =  lbfgs_malloc(SSM_OPTIMVAR);;

    x_g2 =  lbfgs_malloc(SSM_OPTIMVAR);;

    for (i = 0; i < SSM_OPTIMVAR; i++)
    {
        for(j = 0; j < SSM_OPTIMVAR; j++)
        {
            if(i == j)
            {
                x_g1[j] = x[j] - ndeps * fabs(x[j]);
                x_g2[j] = x[j] + ndeps * fabs(x[j]);
            }
            else
            {
                x_g1[j] = x[j];
                x_g2[j] = x[j];
            }

        }

        g[i] = (ssmloglik_2bounds(instance, x_g2, instance->q_lower, instance->q_upper, instance->h_lower, instance->h_upper) -
                ssmloglik_2bounds(instance, x_g1, instance->q_lower, instance->q_upper, instance->h_lower, instance->h_upper))/(2 * ndeps);
        //printf("x[%i]=%f, g[%i] = %f\n", i, x[i], i , g[i]);
    }



    lbfgs_free(x_g1);
    lbfgs_free(x_g2);

    return lik;

}

lbfgsfloatval_t evaluate_b6
(
    ssmodel *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
)
{
    int i, j;
    lbfgsfloatval_t ndeps = 1e-3;
    lbfgsfloatval_t *x_g1; // for finite gradient 1
    lbfgsfloatval_t *x_g2; // for finite gradient 2
    lbfgsfloatval_t lik;

    lik = ssmloglik(instance, x);

    x_g1 =  lbfgs_malloc(SSM_OPTIMVAR);;

    x_g2 =  lbfgs_malloc(SSM_OPTIMVAR);;

    for (i = 0; i < SSM_OPTIMVAR; i++)
    {
        for(j = 0; j < SSM_OPTIMVAR; j++)
        {
            if(i == j)
            {
                x_g1[j] = x[j] - ndeps;
                x_g2[j] = x[j] + ndeps;
            }
            else
            {
                x_g1[j] = x[j];
                x_g2[j] = x[j];
            }

        }

        g[i] = (ssmloglik(instance, x_g2) -
                ssmloglik(instance, x_g1))/(2 * ndeps);
        //printf("x[%i]=%f, g[%i] = %f\n", i, x[i], i , g[i]);
    }



    lbfgs_free(x_g1);
    lbfgs_free(x_g2);

    return lik;

}


lbfgsfloatval_t ssmloglik(
    ssmodel *instance,
    const lbfgsfloatval_t *x
)
{
    int i, k;
    lbfgsfloatval_t lik = 0.0;
    int d = 0;
    int rankp = 0;
    int valid_date_count = 0;

    double* vt;
    double* ft;
    double* finf;
    int initial_date;
    gsl_vector* at;
    gsl_matrix* pinf;
    gsl_matrix* pt;
    gsl_vector* kinf;
    gsl_vector* kt;
    gsl_matrix* Q;
    double* fit_cft;
    fit_cft = (double*) malloc((instance->m + 1) * sizeof(float));
    /* initialize fit_cft*/
    vec_a2fit_cft(instance->a1, fit_cft, instance->clear_date_array[0], instance->m, instance->structure);

    /* alloc memory  */
    vt = malloc(sizeof(double));
    ft = malloc(sizeof(double));
    finf =  malloc(sizeof(double));
    at = gsl_vector_alloc(instance->m);
    kinf = gsl_vector_alloc(instance->m);
    kt = gsl_vector_alloc(instance->m);
    pinf = gsl_matrix_alloc(instance->m, instance->m);
    pt = gsl_matrix_alloc(instance->m, instance->m);
    Q = gsl_matrix_calloc(instance->m, instance->m);

    /* initialize parameter  */
    gsl_matrix_memcpy(pt, instance->P1);
    gsl_vector_memcpy(at, instance->a1);
    gsl_matrix_memcpy(pinf, instance->P1inf);

    for (i = 0; i < (int)instance->P1inf->size1; i++)
        for(k = 0; k < (int)instance->P1inf->size2; k++)
            rankp = rankp + (int)(*gsl_matrix_ptr(instance->P1inf, i, k));


    /*option 1*/
//    gsl_matrix_set(Q, 1, 1, pow(exp(x[0]/2),2));
//    gsl_matrix_set(Q, 2, 2, pow(exp(x[1]/2),2));
//    gsl_matrix_set(Q, 3, 3, pow(exp(x[2]/2),2));
//    gsl_matrix_set(Q, 4, 4, pow(exp(x[3]/2),2));
//    gsl_matrix_set(Q, 5, 5, pow(exp(x[4]/2),2));

    /*option 3*/
    gsl_matrix_set(Q, 0, 0, (double)(Q00_UPPER + Q00_LOWER)/2 + x[0] * (double)(Q00_UPPER - Q00_LOWER)*0.5/pow(1 + x[0]*x[0],0.5));
    gsl_matrix_set(Q, 1, 1, (double)(Q11_UPPER + Q11_LOWER)/2 + x[1] * (double)(Q11_UPPER - Q11_LOWER)*0.5/pow(1 + x[1]*x[1],0.5));
    gsl_matrix_set(Q, 2, 2, (double)(Q11_UPPER + Q11_LOWER)/2 + x[1] * (double)(Q11_UPPER - Q11_LOWER)*0.5/pow(1 + x[1]*x[1],0.5));
    gsl_matrix_set(Q, 3, 3, (double)(Q22_UPPER + Q22_LOWER)/2 + x[2] * (double)(Q22_UPPER - Q22_LOWER)*0.5/pow(1 + x[2]*x[2],0.5));
    gsl_matrix_set(Q, 4, 4, (double)(Q22_UPPER + Q22_LOWER)/2 + x[2] * (double)(Q22_UPPER - Q22_LOWER)*0.5/pow(1 + x[2]*x[2],0.5));




    //H = pow(exp(x[5]/2),2);
    // H = (q_upper + q_lower)/2 + x[3] * (double)(q_upper - q_lower)*0.5/pow(1 + x[3]*x[3],0.5);


    initial_date = instance->clear_date_array[0];
    /* diffuse initialization */
    if(rankp > 0)
    {
         while((d < instance->n)&&(rankp > 0))
         {
             dfilter1step(d + initial_date, instance->clear_date_array[valid_date_count], &valid_date_count,
                          (double)instance->yt[valid_date_count], instance->Z,
                          instance->H, instance->T, Q,
                          at, pt, vt, ft, kt, pinf, finf, kinf,
                          &rankp, &lik, instance->m, TRUE, fit_cft, instance->structure);
             d = d + 1;

         }

    }

   int t = 0;
    /* non-diffuse filtering continues from i=d */
    for( i = d; i < instance->n; i++)
    {
        if(i == 32)
        {
           t = t + 1;
        }
        filter1step(i + initial_date, instance->clear_date_array[valid_date_count], &valid_date_count,
                    (double)instance->yt[valid_date_count], instance->Z,
                     instance->H, instance->T, Q, at, pt, vt, ft, kt, &lik, instance->m, TRUE, fit_cft, instance->structure);
    }


    free(vt);
    free(ft);
    free(finf);
    free(fit_cft);
    gsl_vector_free(at);
    gsl_vector_free(kinf);
    gsl_vector_free(kt);
    gsl_matrix_free(pinf);
    gsl_matrix_free(pt);
    gsl_matrix_free(Q);

    return -lik;
}


lbfgsfloatval_t ssmloglik_bounds(
    ssmodel *instance,
    const lbfgsfloatval_t *x,
    double q_lower,
    double q_upper
)
{
    int i, k;
    lbfgsfloatval_t lik = 0.0;
    int d = 0;
    int rankp = 0;
    int valid_date_count = 0;
    double tmp;
    double* fit_cft;
    double* vt;
    double* ft;
    double* finf;
    int initial_date;
    gsl_vector* at;
    gsl_matrix* pinf;
    gsl_matrix* pt;
    gsl_vector* kinf;
    gsl_vector* kt;
    gsl_matrix* Q;


    fit_cft = (double*) malloc((instance->m + 1) * sizeof(double));
    /* initialize fit_cft*/
    vec_a2fit_cft(instance->a1, fit_cft, instance->clear_date_array[0], instance->m, instance->structure);

    /* alloc memory  */
    vt = malloc(sizeof(double));
    ft = malloc(sizeof(double));
    finf =  malloc(sizeof(double));
    at = gsl_vector_alloc(instance->m);
    kinf = gsl_vector_alloc(instance->m);
    kt = gsl_vector_alloc(instance->m);
    pinf = gsl_matrix_alloc(instance->m, instance->m);
    pt = gsl_matrix_alloc(instance->m, instance->m);
    Q = gsl_matrix_calloc(instance->m, instance->m);

    /* initialize parameter  */
    gsl_matrix_memcpy(pt, instance->P1);
    gsl_vector_memcpy(at, instance->a1);
    gsl_matrix_memcpy(pinf, instance->P1inf);

    for (i = 0; i < (int)instance->P1inf->size1; i++)
        for(k = 0; k < (int)instance->P1inf->size2; k++)
            rankp = rankp + (int)(*gsl_matrix_ptr(instance->P1inf, i, k));


    /*option 1*/
//    gsl_matrix_set(Q, 1, 1, pow(exp(x[0]/2),2));
//    gsl_matrix_set(Q, 2, 2, pow(exp(x[1]/2),2));
//    gsl_matrix_set(Q, 3, 3, pow(exp(x[2]/2),2));
//    gsl_matrix_set(Q, 4, 4, pow(exp(x[3]/2),2));
//    gsl_matrix_set(Q, 5, 5, pow(exp(x[4]/2),2));

    /*option 3*/
    tmp = (double)(q_upper + q_lower)/2 + x[0] * (double)(q_upper - q_lower)*0.5/pow(1 + x[0]*x[0], 0.5);
    gsl_matrix_set(instance->Q, 0, 0, tmp);
    gsl_matrix_set(instance->Q, 1, 1, tmp);
    gsl_matrix_set(instance->Q, 2, 2, tmp);
    gsl_matrix_set(instance->Q, 3, 3, tmp);
    gsl_matrix_set(instance->Q, 4, 4, tmp);




    //H = pow(exp(x[5]/2),2);
    // H = (q_upper + q_lower)/2 + x[3] * (double)(q_upper - q_lower)*0.5/pow(1 + x[3]*x[3],0.5);


    initial_date = instance->clear_date_array[0];
    /* diffuse initialization */
    if(rankp > 0)
    {
         while((d < instance->n)&&(rankp > 0))
         {
             dfilter1step(d + initial_date, instance->clear_date_array[valid_date_count], &valid_date_count,
                          (double)instance->yt[valid_date_count], instance->Z,
                          instance->H, instance->T, Q,
                          at, pt, vt, ft, kt, pinf, finf, kinf,
                          &rankp, &lik, instance->m, TRUE, fit_cft, instance->structure);
             d = d + 1;

         }

    }

   int t = 0;
    /* non-diffuse filtering continues from i=d */
    for( i = d; i < instance->n; i++)
    {
        if(i == 32)
        {
           t = t + 1;
        }
        filter1step(i + initial_date, instance->clear_date_array[valid_date_count], &valid_date_count,
                    (double)instance->yt[valid_date_count], instance->Z,
                     instance->H, instance->T, Q, at, pt, vt, ft, kt, &lik, instance->m, TRUE, fit_cft, instance->structure);
    }


    free(vt);
    free(ft);
    free(finf);
    free(fit_cft);
    gsl_vector_free(at);
    gsl_vector_free(kinf);
    gsl_vector_free(kt);
    gsl_matrix_free(pinf);
    gsl_matrix_free(pt);
    gsl_matrix_free(Q);

    return -lik;
}


lbfgsfloatval_t ssmloglik_2bounds(
    ssmodel *instance,
    const lbfgsfloatval_t *x,
    double q_lower,
    double q_upper,
    double h_lower,
    double h_upper
)
{
    int i, k;
    // double c = 0.5 * log(8.0 * atan(1.0));
    lbfgsfloatval_t lik = 0.0;
    int d = 0;
    int valid_date_count = 0;
    double H;
    double tmp;
    double* fit_cft;
    double* vt;
    double* ft;
    int initial_date;
    gsl_vector* at;
    gsl_matrix* pt;
    gsl_vector* kt;
    gsl_matrix* Q;

    fit_cft = (double*) malloc((instance->m + 1) * sizeof(double));
    /* initialize fit_cft*/
    vec_a2fit_cft(instance->a1, fit_cft, instance->clear_date_array[0], instance->m, instance->structure);

    /* alloc memory  */
    vt = malloc(sizeof(double));
    ft = malloc(sizeof(double));
    at = gsl_vector_alloc(instance->m);
    kt = gsl_vector_alloc(instance->m);;
    pt = gsl_matrix_alloc(instance->m, instance->m);
    Q = gsl_matrix_calloc(instance->m, instance->m);

    /* initialize parameter  */
    gsl_matrix_memcpy(pt, instance->P1);
    gsl_vector_memcpy(at, instance->a1);


//    for (i = 0; i < (int)instance->P1inf->size1; i++)
//        for(k = 0; k < (int)instance->P1inf->size2; k++)
//            rankp = rankp + (int)(*gsl_matrix_ptr(instance->P1inf, i, k));


    /*option 1*/
//    gsl_matrix_set(Q, 1, 1, pow(exp(x[0]/2),2));
//    gsl_matrix_set(Q, 2, 2, pow(exp(x[1]/2),2));
//    gsl_matrix_set(Q, 3, 3, pow(exp(x[2]/2),2));
//    gsl_matrix_set(Q, 4, 4, pow(exp(x[3]/2),2));
//    gsl_matrix_set(Q, 5, 5, pow(exp(x[4]/2),2));

    /*option 3*/
    tmp = (double)(q_upper + q_lower)/2 + x[0] * (double)(q_upper - q_lower)*0.5/pow(1 + x[0]*x[0], 0.5);
    gsl_matrix_set(Q, 0, 0, tmp);
    gsl_matrix_set(Q, 1, 1, tmp);
    gsl_matrix_set(Q, 2, 2, tmp);
    gsl_matrix_set(Q, 3, 3, tmp);
    gsl_matrix_set(Q, 4, 4, tmp);

    H = (h_upper + h_lower)/2 + x[1] * (double)(h_upper - h_lower)*0.5/pow(1 + x[1]*x[1], 0.5);


    //H = pow(exp(x[5]/2),2);
    // H = (q_upper + q_lower)/2 + x[3] * (double)(q_upper - q_lower)*0.5/pow(1 + x[3]*x[3],0.5);


    initial_date = instance->clear_date_array[0];
    /* diffuse initialization */
//    if(rankp > 0)
//    {
//         while((d < instance->n)&&(rankp > 0))
//         {
//             dfilter1step(d + initial_date, instance->clear_date_array[valid_date_count], &valid_date_count,
//                          (double)instance->yt[valid_date_count], instance->Z,
//                          H, instance->T, Q,
//                          at, pt, vt, ft, kt, pinf, finf, kinf,
//                          &rankp, &lik, DEFAULT_M);
//             d = d + 1;

//         }

//    }

    /* non-diffuse filtering continues from i=d */
    for( i = d; i < instance->n; i++)
    {
        filter1step(i + initial_date, instance->clear_date_array[valid_date_count], &valid_date_count,
                    (double)instance->yt[valid_date_count], instance->Z,
                     H, instance->T, Q, at, pt, vt, ft, kt, &lik, instance->m, TRUE, fit_cft, instance->structure);
    }


    free(vt);
    free(ft);
    free(fit_cft);
    gsl_vector_free(at);
    gsl_vector_free(kt);
    gsl_matrix_free(pt);
    gsl_matrix_free(Q);

    return -lik;
}

double ssmloglik_gridsearching
(
    ssmodel *instance,
    double H,
    double q,
    double *f_rmse,
    double *v_rmse
)
{
    int i, k;
    double c = 0.5 * log(8.0 * atan(1.0));
    double lik = 0.0;
    int d = 0;
    // int rankp = 0;
    int valid_date_count = 0;
    // double tmp;
    double* fit_cft;
    double* vt;
    double* ft;
    int initial_date;
    gsl_vector* at;
    gsl_matrix* pt;
    gsl_vector* kt;
    gsl_matrix* Q;
    int valid_count = 0;

    fit_cft = (double*) malloc(MAX_NUM_C * sizeof(double));
    /* initialize fit_cft*/
    vec_a2fit_cft(instance->a1, fit_cft, instance->clear_date_array[0], instance->m, instance->structure);

    /* alloc memory  */
    vt = (double*)malloc(sizeof(double));
    ft = (double*)malloc(sizeof(double));
    at = gsl_vector_alloc(instance->m);
    kt = gsl_vector_alloc(instance->m);
    pt = gsl_matrix_alloc(instance->m, instance->m);
    Q = gsl_matrix_calloc(instance->m, instance->m);

    /* initialize parameter  */
    gsl_matrix_memcpy(pt, instance->P1);
    gsl_vector_memcpy(at, instance->a1);

    for(i = 0; i < instance->m; i++)
       gsl_matrix_set(Q, i, i, q);



    //H = pow(exp(x[5]/2),2);
    // H = (q_upper + q_lower)/2 + x[3] * (double)(q_upper - q_lower)*0.5/pow(1 + x[3]*x[3],0.5);


    initial_date = instance->clear_date_array[0];

    /* non-diffuse filtering continues from i=d */
    for( i = d; i < instance->n; i++)
    {
        filter1step(i + initial_date, instance->clear_date_array[valid_date_count], &valid_date_count,
                    (double)instance->yt[valid_date_count], instance->Z,
                     H, instance->T, Q, at, pt, vt, ft, kt, &lik, instance->m, TRUE, fit_cft, instance->structure);

        /* first obs shouldn't be counted as arbitary initial P*/
        if((i + initial_date) == instance->clear_date_array[valid_date_count - 1])
        {

            if(i > d)
            //if(valid_date_count > SCCD_MAX_NUM_C - 1)
            {
                if((*vt) * (*vt) < 1.96 * 1.96 * H)
                {
                    *v_rmse = *v_rmse + (*vt) * (*vt);
                    *f_rmse = *f_rmse + (*ft);
                    valid_count++;
                }

            }
        }

    }

    *f_rmse = *f_rmse / valid_count;
    *v_rmse = *v_rmse / valid_count;

    free(vt);
    free(ft);
    free(fit_cft);
    gsl_vector_free(at);
    gsl_vector_free(kt);
    gsl_matrix_free(pt);
    gsl_matrix_free(Q);

    return lik;
}

int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
)
{
//    printf("Iteration %d:\n", k);
//    printf("  fx = %f, x[0] = %f, x[1] = %f",
//           fx, x[0], x[1]);
//    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
//    printf("\n");
    return 0;
}

/* the function that optimized state-space model to get initial SSM_OPTIMVAR */
int fitSSM
(
    ssmodel *instance
)
{
    int result = 0;
    lbfgs_parameter_t param;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x;
    int n_param;

    double tmp;
    /* if curve number is larger than 1, H is excluded*/
    n_param = SSM_OPTIMVAR;
    x = lbfgs_malloc(n_param);

     /* Option 1: convert them for initialize the variables. */
//    x[0] = log(sqrt(gsl_matrix_get(instance->Q, 1, 1))) * 2;
//    //printf("x[0]=%f\n",x[0]);
//    x[1] = log(sqrt(gsl_matrix_get(instance->Q, 2, 2))) * 2;
//    //printf("x[1]=%f\n",x[1]);
//    x[2] = log(sqrt(gsl_matrix_get(instance->Q, 3, 3))) * 2;
//    x[3] = log(sqrt(gsl_matrix_get(instance->Q, 4, 4))) * 2;
//    x[4] = log(sqrt(gsl_matrix_get(instance->Q, 5, 5))) * 2;


    /* option2: don't convert  */
//   x[0] = gsl_matrix_get(instance->Q, 1, 1);
//   //printf("x[0]=%f\n",x[0]);
//   x[1] = gsl_matrix_get(instance->Q, 2, 2);
//   //printf("x[1]=%f\n",x[1]);
//   x[2] = gsl_matrix_get(instance->Q, 3, 3);
//   x[3] = gsl_matrix_get(instance->Q, 4, 4);
//   x[4] = gsl_matrix_get(instance->Q, 5, 5);
//   x[5] = instance->H;

     /* option3: using bounded update functions  */

   // printf("gsl_matrix_get(instance->Q, 1, 1)= %f\n", gsl_matrix_get(instance->Q, 1, 1));
    tmp = pow(gsl_matrix_get(instance->Q, 0, 0)-(double)(Q00_UPPER + Q00_LOWER)/2, 2);
    if (gsl_matrix_get(instance->Q, 0, 0) > (double)(Q00_UPPER + Q00_LOWER)/2)
        x[0] =pow(tmp/(pow((double)(Q00_UPPER -Q00_LOWER)/2, 2) - tmp), 0.5);
    else
        x[0] = -pow(tmp/(pow((double)(Q00_UPPER -Q00_LOWER)/2, 2) - tmp), 0.5);

    tmp = pow(gsl_matrix_get(instance->Q, 1, 1)-(double)(Q11_UPPER + Q11_LOWER)/2, 2);
    if (gsl_matrix_get(instance->Q, 1, 1) > (double)(Q11_UPPER + Q11_LOWER)/2)
        x[1] =pow(tmp/(pow((double)(Q11_UPPER -Q11_LOWER)/2, 2) - tmp), 0.5);
    else
        x[1] = -pow(tmp/(pow((double)(Q11_UPPER -Q11_LOWER)/2, 2) - tmp), 0.5);

    tmp = pow(gsl_matrix_get(instance->Q, 2, 2)-(double)(Q22_UPPER + Q22_LOWER)/2, 2);
    if (gsl_matrix_get(instance->Q, 2, 2) > (double)(Q22_UPPER + Q22_LOWER)/2)
        x[2] = pow(tmp/(pow((double)(Q22_UPPER -Q22_LOWER)/2, 2) - tmp), 0.5);
    else
        x[2] = -pow(tmp/(pow((double)(Q22_UPPER -Q22_LOWER)/2, 2) - tmp), 0.5);



    /* if it is for temperature band*/
//    if (b_BT == TRUE)
//    {
//        tmp = pow(instance->H - (double)(q_upper_B6 + q_lower_B6)/2, 2);
//        if (instance->H > (double)(q_upper_B6 + q_lower_B6)/2)
//            x[3] = pow(tmp/(pow((double)(q_upper_B6 -q_lower_B6)/2, 2) - tmp), 0.5);
//        else
//            x[3] =-pow(tmp/(pow((double)(q_upper_B6 -q_lower_B6)/2, 2) - tmp), 0.5);
//    }
//    else
//    {
//        tmp = pow(instance->H - (double)(q_upper + q_lower)/2, 2);
//        if (instance->H > (double)(q_upper + q_lower)/2)
//            x[3] = pow(tmp/(pow((double)(q_upper -q_lower)/2, 2) - tmp), 0.5);
//        else
//            x[3] =-pow(tmp/(pow((double)(q_upper -q_lower)/2, 2) - tmp), 0.5);
//    }


    /* Initialize the parameters for the L-BFGS optimization. */
    lbfgs_parameter_init(&param);
    param.max_iterations = 50;
    //param.past = 5;
    //param.epsilon = 1e-7;
    param.delta = 1e-2;
    //param.linesearch =  LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
    param.linesearch = LBFGS_LINESEARCH_DEFAULT ;
    //param.wolfe = 0.01;
    param.min_step = 1e-7;
//    param.max_step = 1e+20;
//   param.gtol = 0.9;
//    param.max_linesearch = 20;

    /* Start the L-BFGS optimization; this will invoke the callback functions
            evaluate() and progress() when necessary.*/


   result = lbfgs(n_param, x, &fx, evaluate, progress, instance,
               &param);



    /* set the Q and H values in instance */

    /*option 1*/
//    gsl_matrix_set(instance->Q, 1, 1, pow(exp(x[0]/2),2));
//    gsl_matrix_set(instance->Q, 2, 2, pow(exp(x[1]/2),2));
//    gsl_matrix_set(instance->Q, 3, 3, pow(exp(x[2]/2),2));
//    gsl_matrix_set(instance->Q, 4, 4, pow(exp(x[3]/2),2));
//    gsl_matrix_set(instance->Q, 5, 5, pow(exp(x[4]/2),2));

    /*option 3*/
    gsl_matrix_set(instance->Q, 0, 0, (double)(Q00_UPPER + Q00_LOWER)/2 + x[0] * (double)(Q00_UPPER - Q00_LOWER)*0.5/pow(1 + x[0]*x[0],0.5));
    gsl_matrix_set(instance->Q, 1, 1, (double)(Q11_UPPER + Q11_LOWER)/2 + x[1] * (double)(Q11_UPPER - Q11_LOWER)*0.5/pow(1 + x[1]*x[1],0.5));
    gsl_matrix_set(instance->Q, 2, 2, (double)(Q11_UPPER + Q11_LOWER)/2 + x[1] * (double)(Q11_UPPER - Q11_LOWER)*0.5/pow(1 + x[1]*x[1],0.5));
    gsl_matrix_set(instance->Q, 3, 3, (double)(Q22_UPPER + Q22_LOWER)/2 + x[2] * (double)(Q22_UPPER - Q22_LOWER)*0.5/pow(1 + x[2]*x[2],0.5));
    gsl_matrix_set(instance->Q, 4, 4, (double)(Q22_UPPER + Q22_LOWER)/2 + x[2] * (double)(Q22_UPPER - Q22_LOWER)*0.5/pow(1 + x[2]*x[2],0.5));

//    printf(" x[5] = %f\n", x[5]);

//    if(b_BT == TRUE)
//    {
//        instance->H = (q_upper_B6 + q_lower_B6)/2 + x[3] * (double)(q_upper_B6 - q_lower_B6)*0.5/pow(1 + x[3]*x[3],0.5);

//    }
//    else
//    {

//        instance->H = (q_upper + q_lower)/2 + x[3] * (double)(q_upper - q_lower)*0.5/pow(1 + x[3]*x[3], 0.5);

//    }


//    printf("instance->Q[1] = %f\n", gsl_matrix_get(instance->Q, 1, 1));
//    printf("instance->Q[2] = %f\n", gsl_matrix_get(instance->Q, 2, 2));
//    printf("instance->Q[3] = %f\n", gsl_matrix_get(instance->Q, 3, 3));
//    printf("instance->Q[4] = %f\n", gsl_matrix_get(instance->Q, 4, 4));
//    printf("instance->Q[5] = %f\n", gsl_matrix_get(instance->Q, 5, 5));

//    gsl_matrix_set(instance->Q, 1, 1, x[0]);
//    printf("x[0]= %f\n", x[0]);
//    gsl_matrix_set(instance->Q, 2, 2, x[1]);
//    printf("x[1]= %f\n", x[1]);
//    gsl_matrix_set(instance->Q, 3, 3, x[2]);
//    printf("x[2]= %f\n", x[2]);
//    gsl_matrix_set(instance->Q, 4, 4, x[3]);
//    printf("x[3]= %f\n", x[3]);
//    gsl_matrix_set(instance->Q, 5, 5, x[4]);
//    printf("x[4]= %f\n", x[4]);
//    instance->H = x[5];
//    printf("x[5]= %f\n", x[5]);

    lbfgs_free(x);

    return result;


}

/* the function that optimized state-space model to get initial SSM_OPTIMVAR */
int fitSSM_bounds
(
    ssmodel *instance,
    double q_lower,
    double q_upper
)
{
    int result = 0;
    lbfgs_parameter_t param;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x;
    int n_param;

    double tmp;
    double tmp2;
    double tmp3;
    /* if curve number is larger than 1, H is excluded*/
    n_param = SSM_OPTIMVAR;
    x = lbfgs_malloc(n_param);


     /* option3: using bounded update functions  */

   // printf("gsl_matrix_get(instance->Q, 1, 1)= %f\n", gsl_matrix_get(instance->Q, 1, 1));
    tmp = pow(gsl_matrix_get(instance->Q, 0, 0)-(double)(q_lower + q_upper)/2, 2);
    if (gsl_matrix_get(instance->Q, 0, 0) > (double)(q_lower + q_upper)/2)
        x[0] =pow(tmp/(pow((double)(q_upper - q_lower)/2, 2) - tmp), 0.5);
    else
        x[0] = -pow(tmp/(pow((double)(q_upper - q_lower)/2, 2) - tmp), 0.5);


    /* Initialize the parameters for the L-BFGS optimization. */
    lbfgs_parameter_init(&param);
    param.max_iterations = 50;
    //param.past = 5;
    //param.epsilon = 1e-7;
    param.delta = 1e-2;
    //param.linesearch =  LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
    param.linesearch = LBFGS_LINESEARCH_DEFAULT ;
    //param.wolfe = 0.01;
    param.min_step = 1e-7;
//    param.max_step = 1e+20;
//   param.gtol = 0.9;
//    param.max_linesearch = 20;

    /* Start the L-BFGS optimization; this will invoke the callback functions
            evaluate() and progress() when necessary.*/

    instance->q_lower = q_lower;
    instance->q_upper = q_upper;
    result = lbfgs(n_param, x, &fx, evaluate_bounds, progress, instance,
               &param);



    /* set the Q and H values in instance */

    /*option 1*/
//    gsl_matrix_set(instance->Q, 1, 1, pow(exp(x[0]/2),2));
//    gsl_matrix_set(instance->Q, 2, 2, pow(exp(x[1]/2),2));
//    gsl_matrix_set(instance->Q, 3, 3, pow(exp(x[2]/2),2));
//    gsl_matrix_set(instance->Q, 4, 4, pow(exp(x[3]/2),2));
//    gsl_matrix_set(instance->Q, 5, 5, pow(exp(x[4]/2),2));

    /*option 3*/
    tmp2 = (double)(q_upper + q_lower)/2 + x[0] * (double)(q_upper - q_lower)*0.5/pow(1 + x[0]*x[0], 0.5);
    gsl_matrix_set(instance->Q, 0, 0, tmp2);
    gsl_matrix_set(instance->Q, 1, 1, tmp2);
    gsl_matrix_set(instance->Q, 2, 2, tmp2);
    gsl_matrix_set(instance->Q, 3, 3, tmp2);
    gsl_matrix_set(instance->Q, 4, 4, tmp2);


    lbfgs_free(x);

    return result;

}

/* the function that optimized state-space model to get initial SSM_OPTIMVAR */
int fitSSM_2bounds
(
    ssmodel *instance,
    double q_lower,
    double q_upper,
    double h_lower,
    double h_upper
)
{
    int result = 0;
    lbfgs_parameter_t param;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x;


    double tmp;
    double tmp2;
    /* if curve number is larger than 1, H is excluded*/
    x = lbfgs_malloc(SSM_OPTIMVAR);


     /* option3: using bounded update functions  */

   // printf("gsl_matrix_get(instance->Q, 1, 1)= %f\n", gsl_matrix_get(instance->Q, 1, 1));
    tmp = pow(gsl_matrix_get(instance->Q, 0, 0)-(double)(q_lower + q_upper)/2, 2);
    if (gsl_matrix_get(instance->Q, 0, 0) > (double)(q_lower + q_upper)/2)
        x[0] =pow(tmp/(pow((double)(q_upper - q_lower)/2, 2) - tmp), 0.5);
    else
        x[0] = -pow(tmp/(pow((double)(q_upper - q_lower)/2, 2) - tmp), 0.5);


    tmp = pow(instance->H - (double)(h_upper + h_lower)/2, 2);
    if (instance->H > (double)(h_upper + h_lower)/2)
        x[1] = pow(tmp/(pow((double)(h_upper - h_lower)/2, 2) - tmp), 0.5);
    else
        x[1] =-pow(tmp/(pow((double)(h_upper - h_lower)/2, 2) - tmp), 0.5);


    /* Initialize the parameters for the L-BFGS optimization. */
    lbfgs_parameter_init(&param);
    param.max_iterations = 50;
    //param.past = 5;
    //param.epsilon = 1e-7;
    param.delta = 1e-2;
    //param.linesearch =  LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
    param.linesearch = LBFGS_LINESEARCH_DEFAULT ;
    //param.wolfe = 0.01;
    //param.min_step = 1e-7;
    //param.epsilon = 1e-2;
//    param.max_step = 1e+20;
//   param.gtol = 0.9;
    param.max_linesearch = 5;

    /* Start the L-BFGS optimization; this will invoke the callback functions
            evaluate() and progress() when necessary.*/

    instance->q_lower = q_lower;
    instance->q_upper = q_upper;
    instance->h_lower = h_lower;
    instance->h_upper = h_upper;
    result = lbfgs(SSM_OPTIMVAR, x, &fx, evaluate_2bounds,  progress, instance,&param);



    /* set the Q and H values in instance */

    /*option 1*/
//    gsl_matrix_set(instance->Q, 1, 1, pow(exp(x[0]/2),2));
//    gsl_matrix_set(instance->Q, 2, 2, pow(exp(x[1]/2),2));
//    gsl_matrix_set(instance->Q, 3, 3, pow(exp(x[2]/2),2));
//    gsl_matrix_set(instance->Q, 4, 4, pow(exp(x[3]/2),2));
//    gsl_matrix_set(instance->Q, 5, 5, pow(exp(x[4]/2),2));

    /*option 3*/
    tmp2 = (double)(q_upper + q_lower)/2 + x[0] * (double)(q_upper - q_lower)*0.5/pow(1 + x[0]*x[0], 0.5);
    gsl_matrix_set(instance->Q, 0, 0, tmp2);
    gsl_matrix_set(instance->Q, 1, 1, tmp2);
    gsl_matrix_set(instance->Q, 2, 2, tmp2);
    gsl_matrix_set(instance->Q, 3, 3, tmp2);
    gsl_matrix_set(instance->Q, 4, 4, tmp2);

    instance->H = (h_upper + h_lower)/2 + x[1] * (double)(h_upper - h_lower) * 0.5/pow(1 + x[1]*x[1],0.5);
    lbfgs_free(x);

    return result;


}


void fit_cft2vec_a
(
     double *fit_cft,
     gsl_vector* vec_next_a,
     int cur_date,
     int m,
     int structure
)
{
    double w = TWO_PI / AVE_DAYS_IN_A_YEAR;
    int count_m = 1;
    gsl_vector_set(vec_next_a, 0, (double)(fit_cft[0]));

    if(count_m < m){
        if(structure % 10 == 1){
            gsl_vector_set(vec_next_a, count_m , (double)(fit_cft[2] * cos((double)cur_date * w)
                           + fit_cft[3] * sin((double)cur_date * w)));
            gsl_vector_set(vec_next_a, count_m + 1, (double)(- fit_cft[2] * sin((double)cur_date * w)
                           + fit_cft[3] * cos((double)cur_date * w)));
            count_m = count_m + 2;
        }
    }


    if(count_m < m){
        if(structure / 10 == 1){
            gsl_vector_set(vec_next_a, count_m, (double)(fit_cft[4] * cos((double)cur_date * 2 * w)
                           + fit_cft[5] * sin((double)cur_date * 2 * w)));
            gsl_vector_set(vec_next_a, count_m + 1, (double)(- fit_cft[4] * sin((double)cur_date * 2 * w)
                           + fit_cft[5] * cos((double)cur_date * 2 * w)));
            count_m = count_m + 2;
        }
    }

    if(count_m < m){
        if(structure / 100 == 1){
            gsl_vector_set(vec_next_a, count_m, (double)(fit_cft[6] * cos((double)cur_date * 3 * w)
                           + fit_cft[7] * sin((double)cur_date * 3 * w)));
            gsl_vector_set(vec_next_a, count_m + 1, (double)(- fit_cft[6] * sin((double)cur_date * 3 * w)
                           + fit_cft[7] * cos((double)cur_date * 3 * w)));
            count_m = count_m + 2;
        }
    }

//    else if(m = 5)
//    {
//        gsl_vector_set(vec_next_a, 5, 0);
//        gsl_vector_set(vec_next_a, 6, 0);
//    }

    //printf("cos((double)cur_date * 2 * w) = %f\n", cos((double)cur_date * 2 * w));

}


void vec_a2fit_cft
(
     gsl_vector *vec_next_a,
     double *fit_cft,
     int cur_date,
     int m,
     int structure
)
{
    double w = TWO_PI / AVE_DAYS_IN_A_YEAR;
    int count_m = 1;
    fit_cft[0] = gsl_vector_get(vec_next_a, 0);
    fit_cft[1] = 0;
    int i;
    for (i = 1; i < LASSO_COEFFS;i++)
        fit_cft[i] = 0;

    if(count_m < m){
        if(structure % 10 == 1){
            fit_cft[2] = cos((double)cur_date * w) * gsl_vector_get(vec_next_a, count_m) -
                                sin((double)cur_date * w) * gsl_vector_get(vec_next_a, count_m + 1);
            fit_cft[3] = cos((double)cur_date * w) * gsl_vector_get(vec_next_a, count_m + 1) +
                                sin((double)cur_date * w) * gsl_vector_get(vec_next_a, count_m);
            count_m = count_m + 2;
        }
    }

    if(count_m < m){
        if(structure / 10 == 1){
            fit_cft[4] = cos((double)cur_date * 2 * w) * gsl_vector_get(vec_next_a, count_m) -
                                sin((double)cur_date * 2 * w) * gsl_vector_get(vec_next_a, count_m + 1);
            fit_cft[5] = cos((double)cur_date * 2 * w) * gsl_vector_get(vec_next_a, count_m + 1) +
                                sin((double)cur_date * 2 * w) * gsl_vector_get(vec_next_a, count_m);
            count_m = count_m + 2;
        }
    }

    if(count_m < m){
        if(structure / 100 == 1){
            fit_cft[6] = cos((double)cur_date * 3 * w) * gsl_vector_get(vec_next_a, count_m) -
                                sin((double)cur_date * 3 * w) * gsl_vector_get(vec_next_a, count_m + 1);
            fit_cft[7] = cos((double)cur_date * 3 * w) * gsl_vector_get(vec_next_a, count_m + 1) +
                                sin((double)cur_date * 3 * w) * gsl_vector_get(vec_next_a, count_m);
            count_m = count_m + 2;
        }
    }


    //printf("cos((double)cur_date * 2 * w) = %f\n", cos((double)cur_date * 2 * w));



}
//void vmmin(int n0, double *b, double *Fmin, optimfn fminfn, optimgr fmingr,
//      int maxit, int trace, int *mask,
//      double abstol, double reltol, int nREPORT, void *ex,
//      int *fncount, int *grcount, int *fail)
//{
//    bool accpoint, enough;
//    double *g, *t, *X, *c, **B;
//    int   count, funcount, gradcount;
//    double f, gradproj;
//    int   i, j, ilast, iter = 0;
//    double s, steplength;
//    double D1, D2;
//    int   n, *l;

//    if (maxit <= 0) {
//    *fail = 0;
//    *Fmin = fminfn(n0, b, ex);
//    *fncount = *grcount = 0;
//    return;
//    }

//    if (nREPORT <= 0)
//    error(_("REPORT must be > 0 (method = \"BFGS\")"));
//    l = (int *) R_alloc(n0, sizeof(int));
//    n = 0;
//    for (i = 0; i < n0; i++) if (mask[i]) l[n++] = i;
//    g = vect(n0);
//    t = vect(n);
//    X = vect(n);
//    c = vect(n);
//    B = Lmatrix(n);
//    f = fminfn(n0, b, ex);
//    if (!R_FINITE(f))
//    error(_("initial value in 'vmmin' is not finite"));
//    if (trace) Rprintf("initial  value %f \n", f);
//    *Fmin = f;
//    funcount = gradcount = 1;
//    fmingr(n0, b, g, ex);
//    iter++;
//    ilast = gradcount;

//    do {
//    if (ilast == gradcount) {
//        for (i = 0; i < n; i++) {
//        for (j = 0; j < i; j++) B[i][j] = 0.0;
//        B[i][i] = 1.0;
//        }
//    }
//    for (i = 0; i < n; i++) {
//        X[i] = b[l[i]];
//        c[i] = g[l[i]];
//    }
//    gradproj = 0.0;
//    for (i = 0; i < n; i++) {
//        s = 0.0;
//        for (j = 0; j <= i; j++) s -= B[i][j] * g[l[j]];
//        for (j = i + 1; j < n; j++) s -= B[j][i] * g[l[j]];
//        t[i] = s;
//        gradproj += s * g[l[i]];
//    }

//    if (gradproj < 0.0) {	/* search direction is downhill */
//        steplength = 1.0;
//        accpoint = FALSE;
//        do {
//        count = 0;
//        for (i = 0; i < n; i++) {
//            b[l[i]] = X[i] + steplength * t[i];
//            if (reltest + X[i] == reltest + b[l[i]]) /* no change */
//            count++;
//        }
//        if (count < n) {
//            f = fminfn(n0, b, ex);
//            funcount++;
//            accpoint = R_FINITE(f) &&
//            (f <= *Fmin + gradproj * steplength * acctol);
//            if (!accpoint) {
//            steplength *= stepredn;
//            }
//        }
//        } while (!(count == n || accpoint));
//        enough = (f > abstol) &&
//        fabs(f - *Fmin) > reltol * (fabs(*Fmin) + reltol);
//        /* stop if value if small or if relative change is low */
//        if (!enough) {
//        count = n;
//        *Fmin = f;
//        }
//        if (count < n) {/* making progress */
//        *Fmin = f;
//        fmingr(n0, b, g, ex);
//        gradcount++;
//        iter++;
//        D1 = 0.0;
//        for (i = 0; i < n; i++) {
//            t[i] = steplength * t[i];
//            c[i] = g[l[i]] - c[i];
//            D1 += t[i] * c[i];
//        }
//        if (D1 > 0) {
//            D2 = 0.0;
//            for (i = 0; i < n; i++) {
//            s = 0.0;
//            for (j = 0; j <= i; j++)
//                s += B[i][j] * c[j];
//            for (j = i + 1; j < n; j++)
//                s += B[j][i] * c[j];
//            X[i] = s;
//            D2 += s * c[i];
//            }
//            D2 = 1.0 + D2 / D1;
//            for (i = 0; i < n; i++) {
//            for (j = 0; j <= i; j++)
//                B[i][j] += (D2 * t[i] * t[j]
//                    - X[i] * t[j] - t[i] * X[j]) / D1;
//            }
//        } else {	/* D1 < 0 */
//            ilast = gradcount;
//        }
//        } else {	/* no progress */
//        if (ilast < gradcount) {
//            count = 0;
//            ilast = gradcount;
//        }
//        }
//    } else {		/* uphill search */
//        count = 0;
//        if (ilast == gradcount) count = n;
//        else ilast = gradcount;
//        /* Resets unless has just been reset */
//    }
//    if (trace && (iter % nREPORT == 0))
//        Rprintf("iter%4d value %f\n", iter, f);
//    if (iter >= maxit) break;
//    if (gradcount - ilast > 2 * n)
//        ilast = gradcount;	/* periodic restart */
//    } while (count != n || ilast != gradcount);
//    if (trace) {
//    Rprintf("final  value %f \n", *Fmin);
//    if (iter < maxit) Rprintf("converged\n");
//    else Rprintf("stopped after %i iterations\n", iter);
//    }
//    *fail = (iter < maxit) ? 0 : 1;
//    *fncount = funcount;
//    *grcount = gradcount;
//}



