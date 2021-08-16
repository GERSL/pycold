/* this is an effort to consolidate all defines.  Previously, they    */
/* scattered throughout .c and/or .h files, and not always included   */
/* and/or avilable everywhere or where needed.  Also, some reduncancy */
/* and conflicts existed.                                             */

/* from ccdc.c */
#define TOTAL_IMAGE_BANDS 7 /* Number of image bands, for loops.   HLS   */
#define TOTAL_INDICES 0 /* Number of image bands, for loops. */
#define TOTAL_BANDS 8     /* Total image plus mask bands, for loops.  HLS */
#define MIN_NUM_C 4       /* Minimum number of coefficients           */
#define MID_NUM_C 6       /* Mid-point number of coefficients         */
#define SCCD_MAX_NUM_C 6   /* Maximum number of coefficients           */
#define SCCD_NUM_C 5       /* Mid-point number of coefficients         */
#define MAX_NUM_C 8       /* Maximum number of coefficients           */
#define N_TIMES 3         /* number of clear observations/coefficients*/
#define NUM_YEARS 365.25  /* average number of days per year          */
#define NUM_FC 40       /* Values change with number of pixels run  */
#define T_CONST 4.42      /* Threshold for cloud, shadow, and snow detection */
                          /* Old T_CONST = 4.89 SY 11/14/2018         */
#define MIN_YEARS 1       /* minimum year for model intialization     */
#define T_SN 0.75         /* no change detection for permanent snow pixels */
#define T_CLR 0.25        /* Fmask fails threshold                    */
#define T_MAX_CG 30.8562  /* chi-square inversed T_max_cg (1e-5) for    */
                          /* last step noise removal                    */
                          /*  SY 11/19/2018: change from 35.8882 to 30.8562)                  */

/* threshold (degree) of mean included angle                          */
#define NSIGN 45
#define NSIGN_sccd 30

/* from 2darray.c */
/* Define a unique (i.e. random) value that can be used to verify a pointer
   points to an LSRD_2D_ARRAY. This is used to verify the operation succeeds to
   get an LSRD_2D_ARRAY pointer from a row pointer. */
#define SIGNATURE 0x326589ab

/* Given an address returned by the allocate routine, get a pointer to the
   entire structure. */
#define GET_ARRAY_STRUCTURE_FROM_PTR(ptr) \
    ((LSRD_2D_ARRAY *)((char *)(ptr) - offsetof(LSRD_2D_ARRAY, memory_block)))


/* from input.c */
//#define TOTAL_IMAGE_BANDS 7

/* from misc.c */
/* 12-31-1972 is 720258 in julian day since year 0001 */
/* matlab use 720624 since year 0000 */
#define JULIAN_DATE_LAST_DAY_1972 720624
#define LANDSAT_START_YEAR 1973
#define LEAP_YEAR_DAYS 366
#define NON_LEAP_YEAR_DAYS 365
#define AVE_DAYS_IN_A_YEAR 365.25
#define ROBUST_COEFFS 5
#define LASSO_COEFFS 8
//#define SCCD_COEFFS 6
//#define TOTAL_IMAGE_BANDS 7

/* from input.h */
/* possible cfmask values */
#define CFMASK_CLEAR   0
#define CFMASK_WATER   1
#define CFMASK_SHADOW  2
#define CFMASK_SNOW    3
#define CFMASK_CLOUD   4
#define CFMASK_FILL  255
#define IMAGE_FILL -9999


#define CFMASK_BAND    7

/* from output.h */
#define FILL_VALUE 255
#define NUM_COEFFS 8
#define NUM_BANDS 7

/* from s_ccd.c */
#define ANNUAL_STATE 1
#define SEMIANNUAL_STATE 2

/* step3_KF_ChangeDetection RETURN_VALUE*/
#define REGULAREND 0
#define CHANGEDETECTED 1
#define FALSECHANGE 2

/* for ini_mode in sccd */
#define SIMPLE 1
#define CONSE_END 6

#define UPDATE_FREQ 3
// #define CONSE 6           /* No. of CONSEquential pixels 4 bldg. model*/
#define LASSO_MIN 6           /* No. of CONSEquential pixels 4 bldg. model*/

#define INITIAL_P_RATIO 0.05
#define SCCD_T_CONST 4.42

#define DEFAULT_S_TCG 11.0705      /* chi-square inversed T_cg (0.95)  */
#define NUM_LASSO_BANDS 5 /* Number of bands for Least Absolute Shrinkage */
                          /* and Selection Operator LASSO regressions */
#define DEFAULT_COLD_TCG 15.0863

//#define DEFAULT_M 5
#define SMOOTH_FACTOR 0
/* variable number to be optimized in State space model */
#define SSM_OPTIMVAR 2

#define S_TCG_single 0


#define INI_Q00 0

#define INI_Q 1
#define INI_Q11 5
#define INI_Q22 5
#define INI_P 1000
#define INI_H 0


#define Q00_LOWER 1E-2
#define Q11_LOWER 1E-2
#define Q22_LOWER 1E-2
#define H_LOWER 900
#define H_LOWER_B6 1000

#define Q00_UPPER 400
#define Q11_UPPER 400
#define Q22_UPPER 400
#define H_UPPER 90000
#define H_UPPER_B6 120000

#define KFAS_LIK_C 0.5 * log(8.0 * atan(1.0))
#define KFAS_TOL 1E-14

#define T_MAX_CG 30.8562  // prob = 1-1e-5
#define S_T_MAX_CG_PROB 0.9999


/* from ccd.c detection method*/
//#define CCD 1
#define COLD 1
#define SCCD 2
// #define EXTRACTION 3 // PIDS legacy
#define OBCOLD 3 // COLD with outputting CM magnitudes
#define OBCOLD_RECONSTRUCT 4


#define MIN_DAYS_CONSE 80
#define SKIP_PERCENTAGE 0.03

// #define WINDOW_SIZE
#define WEIGHT1 0.1
#define WEIGHT2 0.2
#define WEIGHT3 0.3
#define WEIGHT4 0.4

#define LANDSAT45_TM 1
#define LANDSAT7_ETM 2
#define LANDSAT8_OLI 3

# define DIST_TYPE_DIST 1 // abrupt disturbance
# define DIST_TYPE_RESTOR 2 // reforestation/afforestation
# define DIST_TYPE_REGROW 3 // regrowth
# define DIST_TYPE_NOISE 4 // noise

# define TEMPORAL_RMSE_BIN 6
# define TEMPORAL_RMSE_BIN_NUM 61 // i.e, 366/6

#define NDVI_INDEX        7
#define NBR_INDEX       8
#define RGI_INDEX       9
#define TCTWETNESS_INDEX      10
#define TCTGREENNESS_INDEX        11
#define EVI_INDEX       12
#define DI_INDEX       13
#define NDMI_INDEX 14

#define TIFF_FORMAT 1
#define ENVI_FORMAT 2

#define TRAINING_TYPE_REGULAR 0
#define TRAINING_TYPE_PARAMETER 1
#define TRAINING_TYPE_CLASSIFICATION 2

#define NA_VALUE -9999

#define DEFAULT_N_FOCUS_VARIABLE 5
#define DEFAULT_TOTAL_VARIABLE 8

#define TARGETED_LABEL 2

#define DEFAULT_PROBABILITY 0.99
#define DEFAULT_CONSE 6
#define N_FEATURE 7

#define ORDINALDAY_19710101 719528
#define JULY1ST_DOY 183

#define SPATIAL_OUTPUT_INTERVAL 30


