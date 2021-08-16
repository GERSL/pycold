#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdbool.h>
#include "lbfgs.h"
#include "const.h"

typedef struct
{
    int clear_date_array[MAXIMUM_INI_LEN];    /* I: valid date time series  */
    float yt[MAXIMUM_INI_LEN];           /* I: valid observations  */
    gsl_vector* Z;         /* A m matrix corresponding to the system matrix of observation equation */
    double H;            /* A 1*1 matrix corresponding to observational disturbances epsilon      */
    gsl_matrix* T;          /* A m x m matrix corresponding to the first system matrix of state equation. */
    gsl_matrix* Q;           /* A m*m matrix for the state disturbance covariance*/
    gsl_vector* a1;         /* A m x 1 vector containing the expected values of the initial states.    */
    gsl_matrix* P1;         /* A m x m matrix containing the covariance matrix of the nondiffuse part of the initial state vector */
    gsl_matrix* P1inf;
    double q_lower;
    double q_upper;
    double h_lower;
    double h_upper;
    int m;
    int n;                /* number of time points */
    int structure; /* 100: 4 month, 10: semi-annual; 1: annual */
} ssmodel;

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
);

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
);

/* non-diffuse filtering for single time point */
void filter1step
(
    int cur_monitor_date,  /*I */
    int next_valid_date,  /*I/O */
    int* valid_count,
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
    double* lik,       /*I/O*/
    int m,
    bool fast_mode,
    double* fit_cft,
    int structure
);

/* non-diffuse filtering for missing obs */
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
);


/* diffuse filtering for missing obs */
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
);

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
);

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
);

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
);

int fitSSM
(
    ssmodel *instance
);

int fitSSM_bounds
(
    ssmodel *instance,
    double q_lower,
    double q_upper
);

int fitSSM_2bounds
(
    ssmodel *instance,
    double q_lower,
    double q_upper,
    double h_lower,
    double h_upper
);

lbfgsfloatval_t ssmloglik(
    ssmodel *instance,
    const lbfgsfloatval_t *x
);

double ssmloglik_gridsearching
(
    ssmodel *instance,
    double H,
    double q,
    double *f_rmse,
    double *v_rmse
);

lbfgsfloatval_t ssmloglik_bounds(
    ssmodel *instance,
    const lbfgsfloatval_t *x,
    double q_lower,
    double q_upper
);

lbfgsfloatval_t ssmloglik_2bounds(
    ssmodel *instance,
    const lbfgsfloatval_t *x,
    double q_lower,
    double q_upper,
    double h_lower,
    double h_upper
);

lbfgsfloatval_t evaluate
(
    ssmodel *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
);

lbfgsfloatval_t evaluate_bounds
(
    ssmodel *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
);

lbfgsfloatval_t evaluate_2bounds
(
    ssmodel *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
);


lbfgsfloatval_t evaluate_b6
(
    ssmodel *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
);



void filter1step_pred
(
    int pred_i,  /*I */
    gsl_vector* zt,        /*I */
    double ht,             /*I */
    gsl_matrix* tt,        /*I */
    gsl_matrix* rqr,     /*I */
    gsl_vector* at,    /*I: */
    double ft,
    gsl_matrix* pt,    /*I/O*/
    double* ft_records,    /*I/O: ft records */
    double** at_records,   /*I/O: at records */
    gsl_vector* kt,    /*I/O*/
    int m
);

void fit_cft2vec_a
(
     double *fit_cft,
     gsl_vector *vec_next_a,
     int cur_date,
     int m,
     int structure
);

void vec_a2fit_cft
(
     gsl_vector *vec_next_a,
     double *fit_cft,
     int cur_date,
     int m,
     int structure
);

//typedef double optimfn(int n, double *par, void *ex);
//typedef void optimgr(int n, double *par, double *gr, void *ex);

//void vmmin(int n, double *x, double *Fmin, optimfn fn, optimgr gr, int maxit, int trace,
//int *mask, double abstol, double reltol, int nREPORT, void *ex, int *fncount, int *grcount, int *fail);
