#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <omp.h>
#include <stdbool.h>
#include <unistd.h>
#include "defines.h"
#include "cold.h"
#include "const.h"
#include "utilities.h"
#include "2d_array.h"
#include "input.h"
#include "output.h"
#include "misc.h"
#include "s_ccd.h"

/******************************************************************************
MODULE:  get_args
PURPOSE:  Gets the command-line arguments and validates that the required
arguments were specified.
RETURN VALUE:
Type = int
Value           Description
-----           -----------
FAILURE         Error getting the command-line arguments or a command-line
                argument and associated value were not specified
SUCCESS         No errors encountered
HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
11/15/2018  Su Ye            orginal develop
******************************************************************************/
int get_variables_hpc
(
    int argc,              /* I: number of cmd-line args                    */
    char *argv[],          /* I: string of cmd-line args                    */
    char *in_path,         /* O: directory locaiton for input data          */
    char *out_path,        /* O: directory location for output files        */
    int *method,            /* O: 1 - COLD, 2- S-CCD          */
    char *user_mask_path,   /* O: directory location for user_mask           */
    double *probability_threshold,  /* O: probability threshold between 0 and 1  */
    int *min_days_conse         /* O: minimum peek days, S-CCD only  */
)
{
    char cwd[MAX_STR_LEN]; // current directory path
    //char cwd[] = "/Documents/QTProjects/build-CCDC_C-Desktop-Debug";
    char var_path[MAX_STR_LEN];
    FILE *var_fp;
    char line1[MAX_STR_LEN], line2[MAX_STR_LEN], line3[MAX_STR_LEN],
            line4[MAX_STR_LEN], line5[MAX_STR_LEN], line6[MAX_STR_LEN];
    char FUNC_NAME[] = "get_variables";

    // when there is no variable command-line argument,
    // use the default variable text path


    getcwd(cwd, sizeof(cwd));
    //printf("getvariable");
    sprintf(var_path, "%s/%s", cwd, "variables");

    var_fp = fopen(var_path, "r");

    if(var_fp == NULL)
    {
        RETURN_ERROR("There is no variable file in the bin folder", FUNC_NAME, ERROR);
    }

    fscanf(var_fp, "%s\n", line1);
    //in_path = strchr(line2, '=') + 1;getvariable
    strcpy(in_path, strchr(line1, '=') + 1);

    fscanf(var_fp, "%s\n", line2);
    strcpy(out_path, strchr(line2, '=') + 1);

//    fscanf(var_fp, "%s\n", line4);
//    *n_cores = atoi(strchr(line4, '=') + 1);

//    fscanf(var_fp, "%s\n", line5);
//    *n_nodes = atoi(strchr(line5, '=') + 1);
    fscanf(var_fp, "%s\n", line3);
    *method = atoi(strchr(line3, '=') + 1);

    fscanf(var_fp, "%s\n", line4);
    strcpy(user_mask_path, strchr(line4, '=') + 1);


//    fscanf(var_fp, "%s\n", line7);
//    *row = atoi(strchr(line7, '=') + 1);

//    fscanf(var_fp, "%s\n", line8);
//    *col = atoi(strchr(line8, '=') + 1);

    fscanf(var_fp, "%s\n", line5);
    *probability_threshold = atof(strchr(line5, '=') + 1);

    fscanf(var_fp, "%s\n", line6);
    *min_days_conse = atoi(strchr(line6, '=') + 1);
    //printf("debug: step4\n");


    // check if variable inputs are complete
    if (*method > 2 || *method < 0)
    {
        RETURN_ERROR("Variables error: method variable must be either 1 or 2", FUNC_NAME, ERROR);

    }
    if (*in_path == '\0')
    {
        RETURN_ERROR("Variables error: 'in_path' "
                     "cannot be empty in 'Variables' file!", FUNC_NAME, ERROR);

    }
    if (*out_path == '\0')
    {
        RETURN_ERROR("Variables error: 'out_path' "
                     "cannot be empty in 'Variables' file!", FUNC_NAME, ERROR);
    }

//    fscanf(var_fp, "%s\n", line);
//    fscanf(var_fp, "%s\n", line);
//    fscanf(var_fp, "%s\n", line);
//    fscanf(var_fp, "in_path=%s\n", in_path);
//    fscanf(var_fp, "out_path=%s\n", out_path);
//    fscanf(var_fp, "n_cores=%d\n", n_cores);
//    fscanf(var_fp, "user_mask_path=%s\n", user_mask_path);
//    fscanf(var_fp, "row=%d\n", row);
//    fscanf(var_fp, "col=%d\n", col);

    fclose(var_fp);
    return SUCCESS;

}


int sccd_hpc
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
    double probability_threshold,
    int min_days_conse,
    short int *sensor_buf,
    int conse
)
{
    /* hpc version temporally set these parameter as fixed */
    bool b_fastmode = TRUE;
    char* states_output_dir;
    int training_type = 0; /* '0' is normal */
    int monitorwindow_lowerlin = 0;  /* for training process*/
    int monitorwindow_upperlim = 0;   /* for training process*/
    int n_focus_variable = DEFAULT_N_FOCUS_VARIABLE;
    int n_total_variable = TOTAL_IMAGE_BANDS;
    int focus_blist[TOTAL_IMAGE_BANDS] = {1, 2, 3, 4, 5};
    bool NDVI_INCLUDED = FALSE;
    bool NBR_INCLUDED = FALSE;
    bool RGI_INCLUDED = FALSE;
    bool TCTWETNESS_INCLUDED = FALSE;
    bool TCTGREENNESS_INCLUDED = FALSE;
    bool EVI_INCLUDED = FALSE;
    bool DI_INCLUDED = FALSE;
    bool NDMI_INCLUDED = FALSE;
    bool b_landspecific  = FALSE;
    short int auxval = 0;
    /* hpc version temporally set these parameter as fixed */

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

    id_range = (int*)malloc(valid_num_scenes * sizeof(int));
    // printf("valid_num_scenes is %d\n", valid_num_scenes);

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
        rec_cg[i].category = NA_VALUE;
        rec_cg[i].land_type = NA_VALUE;

        for (j = 0; j < TOTAL_IMAGE_BANDS+TOTAL_INDICES; j++){
            // rec_cg[i].obs_disturb[j] = NA_VALUE;
            rec_cg[i].rmse[j] = NA_VALUE;
            rec_cg[i].magnitude[j] = NA_VALUE;
            for(k = 0; k < SCCD_MAX_NUM_C; k++){ 
                // rec_cg[i].state_disturb[j][k] = NA_VALUE;
                rec_cg[i].coefs[j][k] = NA_VALUE;
            }
            
        }
    }

    //printf("num_fc_initial is %d\n", *num_fc);
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

       result = sccd_stand_procedure(valid_num_scenes, valid_date_array, buf,
                                     fmask_buf, id_range, rec_cg, num_fc,
                                     states_output_dir, b_fastmode, probability_threshold,
                                     min_days_conse, training_type, monitorwindow_lowerlin,
                                     monitorwindow_upperlim, sensor_buf, n_focus_variable,
                                     n_total_variable, focus_blist, NDVI_INCLUDED,
                                     NBR_INCLUDED, RGI_INCLUDED,TCTWETNESS_INCLUDED,TCTGREENNESS_INCLUDED,
                                     EVI_INCLUDED, DI_INCLUDED, NDMI_INCLUDED, b_landspecific, auxval, conse);


    }

    //printf("free stage 7 \n");
    //printf("num_fc is %d\n", *num_fc);
    free(id_range);
    id_range = NULL;
    //printf("free stage 8 \n");
    if (result == SUCCESS)
    {
        return (SUCCESS);
    }
    else
    {
        return (FAILURE);
    }
}


int main(int argc, char *argv[])
{

//    char in_path[] = "/media/su/DataFactory/Dissertation_Meta/Landsat/ENVI_LandsatARD";
//    char mask_path[] = "/home/su/Documents/RProjects/Dissertation/Results/hostspecies_mask.envi";
//    char out_path[] = "/media/su/DataFactory/Dissertation_Meta/Landsat/ccdc_result/RUN0";

//    int n_cores = 8;
//    int row = 310;
//    int col = 10;
//    //int row = 1;
//    //int col =
//    int mode = 2;                    /* CCD detection mode
//                                        3: whole images; 1: pixel-based;
//                                        2: scanline-based*/
    /* inputted argument */                /* CCD detection mode
                                        3: whole images; 1: pixel-based;
                                        2: scanline-based*/
    int method;
    double probability_threshold;
    int output_mode;
    //bool verbose = TRUE;

    char scene_list_filename[] = "scene_list.txt"; /* file name containing list of input sceneIDs */
    char breakdatemap_list_filename[] = "breakdatemap_list.txt"; /* file name containing list of breakdatemap file names */
    char msg_str[MAX_STR_LEN];       /* Input data scene name                 */
    int i, j;                           /* Loop counters                         */
    char scene_list_directory[MAX_STR_LEN]; /* full directory of scene list*/
    int status;                      /* Return value from function call       */
    char FUNC_NAME[] = "main";       /* For printing error messages           */

    int *sdate;                      /* Pointer to list of acquisition dates  */
    char **scene_list;                /* 2-D array for list of scene IDs       */

    int num_scenes;                  /* Number of input scenes defined        */
    Input_meta_t *meta;              /* Structure for ENVI metadata hdr info  */
    char tmpstr[MAX_STR_LEN];        /* char string for text manipulation      */

    time_t now;                  /* For logging the start, stop, and some     */
    int result;

    int num_fc;                        /* the number of functional curve        */
    char out_filename[MAX_STR_LEN];
    char errmsg[MAX_STR_LEN];   /* for printing error text to the log.  */
    //block_num = (int)meta->lines / threads;

    /**************************************************************************/
    /*                   Parallel scanline processing                         */
    /**************************************************************************/

    int n_process;
    int n_block;
    int process_id;
    int n_remain_line;

    int **valid_date_array_scanline;
    short int **fmask_buf_scanline;        /* fmask buf, valid pixels only*/
    int* valid_scene_count_scanline;
    short int **buf;                       /* This is the image bands buffer, valid pixel only*/
    short int *tmp_buf;                   /* This is the image bands buffer, valid pixel only*/
    bool isUserMaskExist;     /*  if user mask exists or not        */
    char *user_mask_scanline;
    bool user_mask_hasvalidcount = FALSE;
    // int n_working_process; /* total processes - 1 (master process)*/
    int i_col;
    int* num_fc_scanline;
    FILE *fh[MAX_SCENE_LIST];
    FILE *fusermask_bip;
    FILE *fhoutput;
    FILE *fhoutput_cm;
    FILE *fhoutput_cm_date;
    FILE* fd;
    char img_filename[MAX_STR_LEN]; /* file name constructed from sceneID   */
    int k, b;
    short int **tmp_buf_2d;                   /* This is the image bands buffer, valid pixel only*/
    //short int *tmp_fmask_buf;              /* fmask buf, valid pixels only          */
    //int *tmp_valid_date_array;             /* Sdate array after cfmask filtering    */
    int n_rows;
    int n_cols;
    Output_t* rec_cg;
    Output_t_sccd* s_rec_cg;
    int new_n_block; // the row number for each block
    int starting_row; // the starting row for each block
    short int tmp_sensor;
    int interval;
    short int *sensor_buf;
    double tcg;
    int starting_date; // the initial date of the whole dataset
    int n_CM_maps;
    short int* CM_outputs;
    char* CM_outputs_date;
    // bool b_singleline = FALSE;
    FILE *fh_breakdatemap[MAX_YEAR_RANGE];
    int *tmp_buf_breakdatemap;
    char **breakdatemap_list;
    char breakdatemap_list_directory[MAX_STR_LEN];
    int num_breakdatemaps;
    int **breakdates_scanline;
    int n_probability_maps;
    int sample_row = 0;
    int sample_col = 0;
    int conse;

    int mode = 3; // 1 - pixel-based; 2 - scanline-based; 3 - area-based

    char CM_filename[MAX_STR_LEN];
    char CM_date_filename[MAX_STR_LEN];
    char CM_fullpath[MAX_STR_LEN];
    char CM_date_fullpath[MAX_STR_LEN];
    char cold_config_fullname[MAX_STR_LEN];
    char cold_config_fullpath[MAX_STR_LEN];
    char out_fullpath[MAX_STR_LEN];
    char out_fullpath_tmp[MAX_STR_LEN];
    bool b_obcold_reconstruction = FALSE; // indicate if it is used to reconstruct rec_cg in obcold
    bool b_outputCM = FALSE;

    int nvals;
    if (argv[1][0] == 'r'){
        b_obcold_reconstruction = TRUE;
        method = COLD;
    }else if(argv[1][0] == 's'){
        method = SCCD;
    }
    else if (argv[1][0] =='c'){
        method = COLD;
    }
    else if (argv[1][0] == 'o'){
        b_outputCM = TRUE;
        method = COLD;
    }
    else
        RETURN_ERROR("The second input parameter has to be r, s, o or c", FUNC_NAME, FAILURE);

    process_id =  (int)strtol(argv[2], NULL, 10) - 1;
    n_process = (int)strtol(argv[3], NULL, 10);
    //argv[4] and argv[5] are in_path and out_path
    // optional parameters
    if (argc == 7){ // debug use for single line
        sample_row = (int)strtol(argv[6], NULL, 10);
        mode = 2;
    }

    //probability_threshold = (double)strtol(argv[6], NULL, 10);
    probability_threshold = DEFAULT_PROBABILITY;
    conse = DEFAULT_CONSE;

    meta = (Input_meta_t *)malloc(sizeof(Input_meta_t));
    if (meta == NULL)
    {
        RETURN_ERROR ("ERROR allocating meta",
                      FUNC_NAME, FAILURE);
    }

    scene_list = (char **) allocate_2d_array (MAX_SCENE_LIST, ARD_STR_LEN,
                                                             sizeof (char));
    if (scene_list == NULL)
    {
        RETURN_ERROR("ERROR allocating scene_list memory", FUNC_NAME, FAILURE);
    }


    /**************************************************************/
    /*                                                            */
    /*   record the start time of just the CDCD         */
    /*     algorithm.  Up until here, it has all just been        */
    /*     setting it up......                                    */
    /*                                                            */
    /**************************************************************/
    time (&now);                 /*     intermediate times.                   */
    snprintf (msg_str, sizeof(msg_str), "CCDC start_time\n", ctime (&now));

    //printf("debug: step1\n");
    LOG_MESSAGE (msg_str, FUNC_NAME);

    // printf("The actual assigned process number is %d\n", n_process);



    //printf("debug: step2\n");
    /**************************************************************/
    /*                                                            */
    /*   read CCD variable                                        */
    /*                                                            */
    /**************************************************************/
//    result = get_variables_hpc(argc, argv, in_path, out_path,
//                            &method, mask_path, &probability_threshold,
//                            &min_days_conse);
//    if(result == ERROR)
//    {
//        RETURN_ERROR("CCDC procedure fails. The program stops!", FUNC_NAME, FAILURE);
//    }

    /* if out_path not exist, create it*/
    struct stat st_tmp = {0};
    if (stat(argv[5], &st_tmp) == -1) {
        mkdir(argv[5], 0700);
    }


    sprintf(scene_list_directory, "%s/%s", argv[4], scene_list_filename);

    if (access(scene_list_directory, F_OK) != 0) /* File does not exist */
    {
        status = create_scene_list(argv[4], &num_scenes, scene_list_filename);
        if(status != SUCCESS){
            sprintf(errmsg,  "Running create_scene_list file: %s\n", scene_list_directory);
            RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
        }
    }
    else
    {
        num_scenes = MAX_SCENE_LIST;
    }

    /**************************************************************/
    /*                                                            */
    /* Fill the scene list array with full path names.            */
    /*                                                            */
    /**************************************************************/
    fd = fopen(scene_list_directory, "r");
    if (fd == NULL)
    {
        sprintf(errmsg,  "Opening scene_list file fails: %s\n", scene_list_directory);
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }

    for (i = 0; i < num_scenes; i++)
    {
        if (fscanf(fd, "%s", tmpstr) == EOF)
            break;
        strcpy(scene_list[i], tmpstr);
    }
    num_scenes = i;

    fclose(fd);

    sdate = (int *)malloc(num_scenes * sizeof(int));

    if (sdate == NULL)
    {
        RETURN_ERROR("ERROR allocating sdate memory", FUNC_NAME, FAILURE);
    }


    /**************************************************************/
    /*                                                            */
    /* Sort scene_list based on year & julian_day, then do the    */
    /* swath filter, but read it above first.                     */
    /*                                                            */
    /**************************************************************/

    status = sort_scene_based_on_year_doy_row(scene_list, num_scenes, sdate, ENVI_FORMAT);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling sort_scene_based_on_year_jday",
                  FUNC_NAME, FAILURE);
    }


    status = read_envi_header(argv[4], scene_list[0], meta);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling read_envi_header",
                      FUNC_NAME, FAILURE);
    }

    n_cols = meta->samples;
    n_rows = meta->lines;
    if(n_process > 1){
        n_block = (int)(meta->lines  / n_process);
    }
    else{
        n_block = meta->lines;
    }

    n_remain_line = meta->lines -  n_block * n_process;
    interval = n_process;

    /* save basic info file used for obia info*/
    if(process_id == 0){
        if(b_obcold_reconstruction == FALSE){
            sprintf(cold_config_fullname, "cold_config.yaml");
            sprintf(cold_config_fullpath, "%s/%s", argv[5], cold_config_fullname);
            fd = fopen(cold_config_fullpath, "w");
            if (fd == NULL)
            {
                RETURN_ERROR("Opening cold_config file", FUNC_NAME, ERROR);
            }
            starting_date = sdate[0];
            n_probability_maps = (sdate[num_scenes - 1] - sdate[0]) / SPATIAL_OUTPUT_INTERVAL + 1;
            fprintf(fd, "starting_date: %d\n", starting_date);
            fprintf(fd, "n_cm_maps: %d\n", n_probability_maps);
            fprintf(fd, "stack_path: %s\n", argv[4]);
            fprintf(fd, "in_cm_path: %s\n", argv[5]);
            fprintf(fd, "cm_interval: %d\n", SPATIAL_OUTPUT_INTERVAL);
            fprintf(fd, "change probability: %f\n", probability_threshold);
            fprintf(fd, "conse: %d\n", conse);
            fclose(fd);
        }else{
            struct stat st = {0};
            sprintf(out_fullpath, "%s/obcold", argv[5]);
            if (stat(out_fullpath, &st) == -1) {
                mkdir(out_fullpath, 0700);
            }
        }
    }

    buf = (short int **) allocate_2d_array (TOTAL_IMAGE_BANDS, num_scenes * n_cols, sizeof (short int));
    if(buf == NULL)
    {
        RETURN_ERROR ("Allocating buf", FUNC_NAME, FAILURE);
    }

    tmp_buf_2d = (short int **) allocate_2d_array (TOTAL_IMAGE_BANDS, num_scenes, sizeof (short int));
    if (tmp_buf_2d == NULL)
    {
        RETURN_ERROR ("ERROR allocating tmp_buf_2d",
                      FUNC_NAME, FAILURE);
    }

    fmask_buf_scanline = (short int **)allocate_2d_array (n_cols, num_scenes, sizeof (short int));
    if(fmask_buf_scanline == NULL)
    {
        RETURN_ERROR ("Allocating fmask_buf_scanline", FUNC_NAME, FAILURE);
    }

    valid_scene_count_scanline = (int*)malloc(n_cols * sizeof(int));
    if(valid_scene_count_scanline == NULL)
    {
        RETURN_ERROR ("Allocating valid_scene_count_scanline", FUNC_NAME, FAILURE);
    }

    valid_date_array_scanline = (int**)allocate_2d_array (n_cols, num_scenes, sizeof (int));
    if(valid_date_array_scanline == NULL)
    {
        RETURN_ERROR ("Allocating valid_date_array_scanline", FUNC_NAME, FAILURE);
    }

    sensor_buf = (short int *) malloc(num_scenes * sizeof(short int));
    if(sensor_buf == NULL)
    {
        RETURN_ERROR ("Allocating sensor_buf fails", FUNC_NAME, FAILURE);
    }
//    scene_list = (char **) allocate_2d_array (MAX_SCENE_LIST, ARD_STR_LEN,
//                                                       sizeof (char));
//    if(scene_list == NULL)
//    {
//        RETURN_ERROR ("Allocating scene_list", FUNC_NAME, FAILURE);
//    }


    /**************************************************************/
    /*                                                            */
    /* Now that we know the actual number of scenes, allocate     */
    /* memory for date array.                                     */
    /*                                                            */
    /**************************************************************/

    num_fc_scanline = (int *)malloc(n_cols * sizeof(short int));
    if (num_fc_scanline == NULL)
    {
        RETURN_ERROR ("ERROR allocating num_fc_scanline",
                      FUNC_NAME, FAILURE);
    }


    user_mask_scanline = (char* )malloc(n_cols * sizeof(char));
    if (user_mask_scanline == NULL)
    {
        RETURN_ERROR ("ERROR allocating user_mask_scanline",
                      FUNC_NAME, FAILURE);
    }


    tmp_buf = malloc(sizeof(short int) * n_cols * TOTAL_BANDS);
    if (tmp_buf == NULL)
    {
        RETURN_ERROR ("ERROR allocating tmp_buf",
                      FUNC_NAME, FAILURE);
    }


    tmp_buf_breakdatemap = malloc(sizeof(int) * n_cols);
    if (tmp_buf_breakdatemap == NULL)
    {
        RETURN_ERROR ("ERROR allocating tmp_buf_breakdatemap",
                      FUNC_NAME, FAILURE);
    }


    starting_date = sdate[0];
    n_CM_maps = (sdate[num_scenes - 1] - sdate[0]) / SPATIAL_OUTPUT_INTERVAL + 1;
    CM_outputs = malloc(sizeof (short int) * n_CM_maps);

    if(CM_outputs == NULL){
         RETURN_ERROR("ERROR allocating CM_outputs", FUNC_NAME, FAILURE);
    }

    CM_outputs_date =  malloc(sizeof (char) * n_CM_maps);
    if(CM_outputs_date == NULL){
         RETURN_ERROR("ERROR allocating CM_outputs_date", FUNC_NAME, FAILURE);
    }



    breakdatemap_list = (char **) allocate_2d_array (MAX_YEAR_RANGE, ARD_STR_LEN, sizeof (char));
    if (breakdatemap_list == NULL)
    {
        RETURN_ERROR("ERROR allocating breakdatemap_list memory", FUNC_NAME, FAILURE);
    }

    breakdates_scanline = (int**)allocate_2d_array (n_cols, MAX_YEAR_RANGE, sizeof (int));
    if(breakdates_scanline == NULL)
    {
        RETURN_ERROR ("Allocating breakdates_scanline", FUNC_NAME, FAILURE);
    }

    if (mode == 3){
        if(process_id < n_remain_line)
        {
            new_n_block = n_block + 1;
            //starting_row = new_n_block * (process_id);
            //prev_n_block = process_id * n_block;
        }
        else
        {
            new_n_block = n_block;
            //starting_row =  (n_block + 1) * n_remain_line + (process_id - n_remain_line) * n_block;
        }
    }else{
        new_n_block = 1;
    }

//    if (b_singleline == TRUE){
//        new_n_block = 1;
//    }

//    n_block = 0;
//    process_id = 0;
//    prev_n_block = process_id * n_block;

    /*********************************************/
    // debug modes
    // printf("process_id = %d; new_n_block = %d; num_scenes = %d \n", process_id, new_n_block, num_scenes);
    //prev_n_block = 119;

    /* for reconstructing rec_cg in sp-cold*/
    if(b_obcold_reconstruction == TRUE)
    {
        sprintf(breakdatemap_list_directory, "%s/breakdate_maps/%s", argv[5], breakdatemap_list_filename);
        if (access(breakdatemap_list_directory, F_OK) != 0)
            RETURN_ERROR("Can't locate breakdate_map_list file", FUNC_NAME, FAILURE);

        fd = fopen(breakdatemap_list_directory, "r");
        if (fd == NULL)
        {
            RETURN_ERROR("Opening breakdatemap_list_directory file", FUNC_NAME, FAILURE);
        }

        for (i = 0; i < MAX_YEAR_RANGE; i++)
        {
            if (fscanf(fd, "%s", tmpstr) == EOF)
                break;
            strcpy(breakdatemap_list[i], tmpstr);
        }
        num_breakdatemaps = i;
        fclose(fd);

        for(i = 0; i < num_breakdatemaps; i++){
            sprintf(img_filename, "%s/breakdate_maps/%s", argv[5], breakdatemap_list[i]);
            fh_breakdatemap[i] = open_raw_binary(img_filename,"rb");
        }
    }


    for(i = 0; i < num_scenes; i++)
    {
        // printf("%s:%d\n", scene_list[i], i);
        sprintf(img_filename, "%s/%s/%s_MTLstack", argv[4], scene_list[i], scene_list[i]);
        fh[i] = open_raw_binary(img_filename,"rb");
    }

    int current_row;
    // for(current_row = starting_row; current_row < starting_row + new_n_block; current_row++)
    tcg = X2(NUM_LASSO_BANDS, probability_threshold);

    for(j = 0; j < new_n_block; j++)
    {

        // for debug
//        current_row = 58;
        //j = new_n_block - 1;
        if (mode == 2)
            current_row = sample_row - 1;
        else
            current_row = process_id + j * interval;

        for (i = 0 ; i < n_cols; i++)
        {
            valid_scene_count_scanline[i] = 0;
        }


        snprintf (msg_str, sizeof(msg_str), "Row %d processing started\n", current_row + 1);
        LOG_MESSAGE (msg_str, FUNC_NAME)

        if(b_obcold_reconstruction == TRUE)
        {
            sprintf(out_filename, "record_change_row%d_obcold.dat", current_row + 1);
            sprintf(out_fullpath, "%s/obcold/%s", argv[5], out_filename);
            sprintf(out_fullpath_tmp, "%s/obcold/%s.part", argv[5], out_filename);
        }
        else
        {
            if(method == COLD){
                sprintf(out_filename, "record_change_row%d_cold.dat", current_row + 1);
                sprintf(CM_filename, "CM_row%d.dat", current_row + 1);
                sprintf(CM_date_filename, "CM_date_row%d.dat", current_row + 1);
                sprintf(out_fullpath, "%s/%s", argv[5], out_filename);
                sprintf(out_fullpath_tmp, "%s/%s.part", argv[5], out_filename);
                sprintf(CM_fullpath, "%s/%s", argv[5], CM_filename);
                sprintf(CM_date_fullpath, "%s/%s", argv[5], CM_date_filename);
                if(b_outputCM == TRUE){
                    fhoutput_cm = fopen(CM_fullpath,"w");
                    fhoutput_cm_date = fopen(CM_date_fullpath,"w");
                }
            }else if(method == SCCD)
            {
                sprintf(out_filename, "record_change_row%d_sccd.dat", current_row + 1);
                sprintf(out_fullpath, "%s/%s", argv[5], out_filename);
                sprintf(out_fullpath, "%s/%s.part", argv[5], out_filename);
            }

        }


        /******************************************************************/
        /*                                                                */
        /*       check if row has be processed or not                     */
        /*                                                                */
        /******************************************************************/
        if( access(out_fullpath, F_OK ) == 0 ) {
            // jump to the next row
            continue;
        }

        fhoutput = fopen(out_fullpath_tmp,"w");


        /******************************************************************/
        /*                                                                */
        /* Read the image bands for this scene.                           */
        /*                                                                */
        /******************************************************************/

         for (i = 0; i < num_scenes; i++)
         {
             //printf("num_scenes is %d with proccess id %d\n", num_scenes, process_id);
             if(scene_list[i][2] == '4' || scene_list[i][2] == '5')
                 tmp_sensor = LANDSAT45_TM;
             else if(scene_list[i][2] == '7')
                 tmp_sensor = LANDSAT7_ETM;
             else if(scene_list[i][2] == '8')
                 tmp_sensor = LANDSAT8_OLI;

             fseek(fh[i], (current_row * n_cols * TOTAL_BANDS) * sizeof(short int), SEEK_SET);
             fread(tmp_buf, sizeof(short int), TOTAL_BANDS * n_cols, fh[i]);
             if (tmp_buf == NULL)
             {
                 sprintf(errmsg, "error reading %d scene, %d row\n", i, current_row + 1);
                 RETURN_ERROR(errmsg, FUNC_NAME, ERROR);
             }

             for(k = 0; k < n_cols; k++)
             {
                 if (tmp_buf[k * TOTAL_BANDS + TOTAL_BANDS - 1] < CFMASK_FILL)
                 {
                     for(b = 0; b < TOTAL_IMAGE_BANDS; b++)
                     {
                         //printf("%d\n", valid_scene_count[k]);
                         buf[b][k * num_scenes + valid_scene_count_scanline[k]] = tmp_buf[k * TOTAL_BANDS + b];
                     }
                     fmask_buf_scanline[k][valid_scene_count_scanline[k]] = (short int)tmp_buf[k * TOTAL_BANDS + TOTAL_BANDS - 1];
                     valid_date_array_scanline[k][valid_scene_count_scanline[k]] = sdate[i];
                     sensor_buf[valid_scene_count_scanline[k]] = tmp_sensor;
                     valid_scene_count_scanline[k] = valid_scene_count_scanline[k] + 1;
                 }
             }
        }

         if(b_obcold_reconstruction == TRUE)
         {
             for (i = 0; i < num_breakdatemaps; i++)
             {
                 fseek(fh_breakdatemap[i], (current_row * n_cols) * sizeof(int), SEEK_SET);
                 fread(tmp_buf_breakdatemap, sizeof(int), n_cols, fh_breakdatemap[i]);
                 for(k = 0; k < n_cols; k++){
                     breakdates_scanline[k][i] = tmp_buf_breakdatemap[k];
                 }
             }
         }

        /******************************************************************/
        /*                                                                */
        /*           Begin CCDC scanline processing                       */
        /*                                                                */
        /******************************************************************/
        for(i_col = 0; i_col < n_cols; i_col++)
        {
            num_fc = 0;
            for(b = 0; b < TOTAL_IMAGE_BANDS; b++)
            {
               tmp_buf_2d[b]  = buf[b] + i_col * num_scenes;
            }

            if(b_obcold_reconstruction == TRUE){
                rec_cg = malloc(NUM_FC * sizeof(Output_t));
                if (rec_cg == NULL)
                {
                    RETURN_ERROR ("ERROR allocating rec_cg",
                                  FUNC_NAME, FAILURE);
                }
                result = obcold_reconstruction_procedure(tmp_buf_2d, fmask_buf_scanline[i_col], valid_date_array_scanline[i_col],
                                                         valid_scene_count_scanline[i_col], rec_cg, &num_fc, num_breakdatemaps,
                                                         breakdates_scanline[i_col],n_cols, i_col + 1, current_row + 1, conse);
                for(i = 0; i < num_fc; i++)
                {
                   result = write_output_binary(fhoutput, rec_cg[i]);
                }
                free(rec_cg);

            }
            else{
                if (method == COLD)
                {

                    rec_cg = malloc(NUM_FC * sizeof(Output_t));
                    if (rec_cg == NULL)
                    {
                        RETURN_ERROR ("ERROR allocating rec_cg",
                                      FUNC_NAME, FAILURE);
                    }

                    for(i = 0; i < n_CM_maps; i++){
                        CM_outputs[i] = NA_VALUE;
                        CM_outputs_date[i] = 255;
                    }

    //                printf("temporal success2 with processid %d\n", process_id);
    //                printf("valid_scene_count_scanline[i_col] is %d;i_col is %d; current_row is %d; probability threshold is %f\n",
    //                       valid_scene_count_scanline[i_col], i_col, current_row, probability_threshold);
                    result = cold(tmp_buf_2d, fmask_buf_scanline[i_col], valid_date_array_scanline[i_col],
                                  valid_scene_count_scanline[i_col], n_cols, i_col + 1, current_row + 1,
                                  tcg, conse, b_outputCM, starting_date, rec_cg, &num_fc, CM_outputs, CM_outputs_date);

                    for(i = 0; i < num_fc; i++)
                    {
                        result = write_output_binary(fhoutput, rec_cg[i]);
                    }

                    if(b_outputCM == TRUE){
                        nvals = fwrite (CM_outputs, sizeof(int16), n_CM_maps, fhoutput_cm);
                        if (nvals != n_CM_maps)
                        {
                            RETURN_ERROR("Incorrect amount of data written", FUNC_NAME, ERROR);
                        }
                        nvals = fwrite (CM_outputs_date, sizeof(int8), n_CM_maps, fhoutput_cm_date);
                        if (nvals != n_CM_maps)
                        {
                            RETURN_ERROR("Incorrect amount of data written", FUNC_NAME, ERROR);
                        }
                    }


                    free(rec_cg);
                    // for debug
                    // printf("free stage 7 \n");
                }
                else
                {
                    s_rec_cg = malloc(NUM_FC * sizeof(Output_t_sccd));
                    if (s_rec_cg == NULL)
                    {
                        RETURN_ERROR ("ERROR allocating s_rec_cg",
                                      FUNC_NAME, FAILURE);
                    }

                    result = sccd_hpc(tmp_buf_2d, fmask_buf_scanline[i_col], valid_date_array_scanline[i_col], valid_scene_count_scanline[i_col],
                                      s_rec_cg, &num_fc,n_cols, i_col + 1, current_row + 1, probability_threshold, MIN_DAYS_CONSE, sensor_buf,conse);
                    //printf("free stage 9 \n");
                    for(i = 0; i < num_fc; i++)
                    {
                        result = write_output_binary_sccd(fhoutput, s_rec_cg[i]);
                    }

                    //printf("free stage 10 \n");
                    free(s_rec_cg);

                }
            }


            if (result != SUCCESS)
            {
                RETURN_ERROR("ccd procedure fails at row %d and col %d \n", current_row + 1, i_col + 1);

            }
//            int count;


        } // end for(i_col = 0; i_col < num_samples; i_col++)


       fclose(fhoutput);

       if(b_outputCM == TRUE){
           fclose(fhoutput_cm);
           fclose(fhoutput_cm_date);
       }

       if (rename(out_fullpath_tmp, out_fullpath) == 0)
       {
           snprintf (msg_str, sizeof(msg_str), "Row %d finished\n", current_row + 1);
           LOG_MESSAGE (msg_str, FUNC_NAME);
       }
       else
       {
           snprintf (msg_str, sizeof(msg_str), "Row %d renaming failed\n", current_row + 1);
           LOG_MESSAGE (msg_str, FUNC_NAME);
       }

    } // for(j = 0; j < n_block; j++)


    for(i = 0; i < num_scenes; i++){
       fclose(fh[i]);
    }

    if (b_obcold_reconstruction == TRUE){
        for(i = 0; i < num_breakdatemaps; i++){
            fclose(fh_breakdatemap[i]);
        }
    }



    /********************************************************/
    /*                   free memory                        */
    /********************************************************/

    free(meta);

    status = free_2d_array((void **)scene_list);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: scene_list\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array((void **)buf);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: buf\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array((void **)tmp_buf_2d);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: tmp_buf_2d\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array((void **) fmask_buf_scanline);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fmask_buf_scanline\n", FUNC_NAME, FAILURE);
    }

    free(valid_scene_count_scanline);

    status = free_2d_array((void **)valid_date_array_scanline);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: valid_date_array_scanline\n", FUNC_NAME, FAILURE);
    }

    free(num_fc_scanline);

    free(sensor_buf);


    free(user_mask_scanline);
    free(tmp_buf);
    free(sdate);
    free(CM_outputs);
    free(CM_outputs_date);

    status = free_2d_array((void **)breakdatemap_list);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: breakdatemap_list\n", FUNC_NAME, FAILURE);
    }

    status = free_2d_array((void **)breakdates_scanline);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: breakdates_scanline\n", FUNC_NAME, FAILURE);
    }

    free(tmp_buf_breakdatemap);

    // free(fh);
    // return SUCCESS;

}
