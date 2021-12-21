#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <stdio.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <stdlib.h>
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

int main(int argc, char *argv[])
{
    int result;
    /* inputted argument, exe mode */
    int mode;                           /* CCD detection mode
                                        3: whole images; 1: pixel-based;
                                        2: scanline-based*/
    char in_path[MAX_STR_LEN];
    float probability_threshold;
    int conse;


    //printf("nonono");

    /**************************************************************/
    /*                                                            */
    /*   read CCD variable                                        */
    /*                                                            */
    /**************************************************************/
    /* need to recover for exe */
//    result = get_variables(argc, argv, &mode, in_path, out_path, &n_cores,
//                           &row, &col, &task, mask_path, &b_fastmode, &b_outputCSV);

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
    //bool verbose = TRUE;

    char msg_str[MAX_STR_LEN];       /* Input data scene name                 */
    int i;                           /* Loop counters                         */
    int starting_date;
    int last_date;
    int status;                      /* Return value from function call       */
    char FUNC_NAME[] = "main";       /* For printing error messages           */

    int *sdate;                      /* Pointer to list of acquisition dates  */
    char **scene_list;                /* 2-D array for list of scene IDs       */

    int num_scenes;                  /* Number of input scenes defined        */
    char tmpstr[MAX_STR_LEN];        /* char string for text manipulation      */

    time_t now;                  /* For logging the start, stop, and some     */

    int num_fc;                        /* the number of functional curve        */
    char out_filename[MAX_STR_LEN];
    char errmsg[MAX_STR_LEN];   /* for printing error text to the log.  */
    //block_num = (int)meta->lines / threads;

    /**************************************************************************/
    /*                   Parallel scanline processing                         */
    /**************************************************************************/

    int n_process;
    int process_id;

    short int *fmask_buf;        /* fmask buf, valid pixels only*/
    short int **buf;                       /* This is the image bands buffer, valid pixel only*/
    short int *tmp_buf;                   /* This is the image bands buffer, valid pixel only*/
    bool isUserMaskExist;     /*  if user mask exists or not        */
    bool user_mask_hasvalidcount = FALSE;
    // int n_working_process; /* total processes - 1 (master process)*/
    int i_col;
    FILE *fh;
    FILE *fusermask_bip;
    FILE *fhoutput;
    FILE *fhoutput_cm;
    FILE *fhoutput_cm_date;
    FILE *fhoutput_cm_direction;
    FILE *fd;
    FILE *fd2;
    char img_filename[MAX_STR_LEN]; /* file name constructed from sceneID   */
    int k, b;
    short int **tmp_buf_2d;                   /* This is the image bands buffer, valid pixel only*/
    //short int *tmp_fmask_buf;              /* fmask buf, valid pixels only          */
    //int *tmp_valid_date_array;             /* Sdate array after cfmask filtering    */
    Output_t* rec_cg;
    Output_t_sccd* s_rec_cg;
    int starting_row; // the starting row for each block
    short int tmp_sensor;
    int interval;
    short int *sensor_buf;
    double tcg;
    int n_cm_maps = 0;
    short int* CM_outputs;
    unsigned char* CM_outputs_date;
    unsigned char* CMdirection_outputs;
    // bool b_singleline = FALSE;
    FILE *fh_breakdatemap;
    int *tmp_buf_breakdatemap;
    char **breakdatemap_list;
    char breakdatemap_list_directory[MAX_STR_LEN];
    int num_breakdatemaps;
    int *breakdates_block;
    // int sample_row = 0;
    // int sample_col = 0;

    char CM_filename[MAX_STR_LEN];
    char CM_direction_filename[MAX_STR_LEN];
    char CM_date_filename[MAX_STR_LEN];
    char CM_fullpath[MAX_STR_LEN];
    char CM_direction_fullpath[MAX_STR_LEN];
    char CM_date_fullpath[MAX_STR_LEN];
    char cold_config_fullname[MAX_STR_LEN];
    char cold_config_fullpath[MAX_STR_LEN];
    char out_fullpath[MAX_STR_LEN];
    char out_fullpath_tmp[MAX_STR_LEN];
    bool b_obcold_reconstruction = FALSE; // indicate if it is used to reconstruct rec_cg in obcold
    bool b_outputCM = FALSE;
    bool b_partition = TRUE; //
    int n_rows = 0;
    int n_cols = 0;
    // int ROW_STEP = 0;
    // int PARTITION = 0;
    int CM_OUTPUT_INTERVAL;
    int nvals;
    int original_row;   // the current processing row
    int original_col;
    Input_meta_t *meta;
    int n_block_v;
    int n_block_h;
    int current_block_x = 0;
    int current_block_y = 0;

    int sample_process;
    int block_width = 1;
    int block_height = 1;
    FILE *sampleFile;
    char * csv_row;
    int valid_scene_count = 0;
    int pixel_qa;
    num_scenes = MAX_SCENE_LIST;


    /**************************************************************/
    /*                                                            */
    /*   preprocessing inputted parameters                        */
    /*                                                            */
    /**************************************************************/
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

    /* stack is complete image, not partitions */
    if(argv[1][1] == 'c'){
        b_partition = FALSE;
    }

    //argv[2] and argv[3] are in_path and out_path
    // optional parameters


    get_coldparameters(&n_rows, &n_cols, &n_block_h, &n_block_v, &CM_OUTPUT_INTERVAL, &probability_threshold, &conse);


    //probability_threshold = (double)strtol(argv[6], NULL, 10);
//    probability_threshold = DEFAULT_PROBABILITY;
//    conse = DEFAULT_CONSE;


    scene_list = (char **) allocate_2d_array (MAX_SCENE_LIST, ARD_STR_LEN,
                                                             sizeof (char));
    if (scene_list == NULL)
    {
        RETURN_ERROR("ERROR allocating scene_list memory", FUNC_NAME, FAILURE);
    }

    meta = (Input_meta_t *)malloc(sizeof(Input_meta_t));
    if (meta == NULL)
    {
        RETURN_ERROR ("ERROR allocating meta",
                      FUNC_NAME, FAILURE);
    }


    /**************************************************************/
    /*                                                            */
    /*   record the start time of just the COLD                   */
    /*     algorithm.  Up until here, it has all just been        */
    /*     setting it up......                                    */
    /*                                                            */
    /**************************************************************/
    time (&now);                 /*     intermediate times.                   */
    snprintf (msg_str, sizeof(msg_str), "COLD starts \n", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);

    /* if out_path not exist, create it*/
    struct stat st_tmp = {0};
    if (stat(argv[3], &st_tmp) == -1) {
        mkdir(argv[3], 0700);
    }

    // chi-square probability
    tcg = X2(NUM_LASSO_BANDS, probability_threshold);


    /**************************************************************/
    /*                                                            */
    /* Now that we know the actual number of scenes, allocate     */
    /* memory for date array.                                     */
    /*                                                            */
    /**************************************************************/
    sdate = (int *)malloc(num_scenes * sizeof(int));

    if (sdate == NULL)
    {
        RETURN_ERROR("ERROR allocating sdate memory", FUNC_NAME, FAILURE);
    }

    buf = (short int **) allocate_2d_array (TOTAL_IMAGE_BANDS, MAX_SCENE_LIST, sizeof (short int));
    if(buf == NULL)
    {
        RETURN_ERROR ("Allocating buf", FUNC_NAME, FAILURE);
    }


    fmask_buf = (short int *) malloc(num_scenes * sizeof(short int));
    if(fmask_buf == NULL)
    {
        RETURN_ERROR ("Allocating fmask_buf", FUNC_NAME, FAILURE);
    }

    sensor_buf = (short int *) malloc(MAX_SCENE_LIST * sizeof (short int));
    if(sensor_buf  == NULL)
    {
        RETURN_ERROR ("Allocating sensor_buf ", FUNC_NAME, FAILURE);
    }

    csv_row = malloc(MAX_STR_LEN * sizeof(char));
    if(sensor_buf  == NULL)
    {
        RETURN_ERROR ("Allocating  csv_row", FUNC_NAME, FAILURE);
    }


    // printf("%d \n", num_scenes);
    sampleFile = fopen(argv[2], "r");
    int row_count = 0;
    while (fgets(csv_row, 255, sampleFile) != NULL)
    {
        if(row_count != 0) // we skip first line because it is a header
        {
            sdate[valid_scene_count] = atoi(strtok(csv_row, ","));
            buf[0][valid_scene_count] = (short)atoi(strtok(NULL, ","));
            buf[1][valid_scene_count] = (short)atoi(strtok(NULL, ","));
            buf[2][valid_scene_count] = (short)atoi(strtok(NULL, ","));
            buf[3][valid_scene_count] = (short)atoi(strtok(NULL, ","));
            buf[4][valid_scene_count] = (short)atoi(strtok(NULL, ","));
            buf[5][valid_scene_count] = (short)atoi(strtok(NULL, ","));
            buf[6][valid_scene_count] = (short)atoi(strtok(NULL, ","));
            pixel_qa = atoi(strtok(NULL, ","));
    //                if (training_type == 1)
    //                {
    //                   fmask_buf[valid_scene_count] = (short)pixel_qa;
    //                   sensor_buf[valid_scene_count] = (short)atoi(strtok(NULL, ","));
    //                }
    //                else
    //                   fmask_buf[valid_scene_count] = (short)qabitval(pixel_qa);
            //fmask_buf[valid_scene_count] = (short)qabitval(pixel_qa);
            fmask_buf[valid_scene_count] = (short)pixel_qa;
            sensor_buf[valid_scene_count] = (short)atoi(strtok(NULL, ","));

            valid_scene_count++;
        }
        row_count ++;
    }
    fclose(sampleFile);

    n_cm_maps = (sdate[valid_scene_count - 1] - sdate[0]) / CM_OUTPUT_INTERVAL + 1;
    starting_date = sdate[0];
    /* save basic info file used for obia info*/

    if(b_obcold_reconstruction == FALSE){
        sprintf(cold_config_fullname, "cold_config.yaml");
        sprintf(cold_config_fullpath, "%s/%s", argv[3], cold_config_fullname);
        fd = fopen(cold_config_fullpath, "w");
        if (fd == NULL)
        {
            RETURN_ERROR("Opening cold_config file", FUNC_NAME, ERROR);
        }
        // starting_date = ORDINALDAY_19710101;

        fprintf(fd, "n_cm_maps: %d\n", n_cm_maps);
        fprintf(fd, "stack_path: %s\n", argv[2]);
        fprintf(fd, "coldresults_path: %s\n", argv[3]);
        fprintf(fd, "starting_date: %d\n", starting_date);
        fprintf(fd, "change_probability: %f\n", probability_threshold);
        fprintf(fd, "conse: %d\n", conse);
        fclose(fd);
    }else{
        struct stat st2 = {0};
        sprintf(out_fullpath, "%s/obcold", argv[3]);
        if (stat(out_fullpath, &st2) == -1) {
            mkdir(out_fullpath, 0700);
        }
    }



//    sensor_buf = (short int *) malloc(num_scenes * sizeof(short int));
//    if(sensor_buf == NULL)
//    {
//        RETURN_ERROR ("Allocating sensor_buf fails", FUNC_NAME, FAILURE);
//    }

//    scene_list = (char **) allocate_2d_array (MAX_SCENE_LIST, ARD_STR_LEN,
//                                                       sizeof (char));
//    if(scene_list == NULL)
//    {
//        RETURN_ERROR ("Allocating scene_list", FUNC_NAME, FAILURE);
//    }

    CM_outputs = malloc(sizeof (short int) * n_cm_maps);
    if(CM_outputs == NULL){
        printf("%d h%d v%d", n_cm_maps, current_block_x, current_block_y);
         RETURN_ERROR("ERROR allocating CM_outputs", FUNC_NAME, FAILURE);
    }

    CMdirection_outputs =  malloc(sizeof (unsigned char) * n_cm_maps);
    if(CMdirection_outputs == NULL){
         RETURN_ERROR("ERROR allocating CM_outputs_date", FUNC_NAME, FAILURE);
    }

    CM_outputs_date =  malloc(sizeof (unsigned char) * n_cm_maps);
    if(CM_outputs_date == NULL){
         RETURN_ERROR("ERROR allocating CM_outputs_date", FUNC_NAME, FAILURE);
    }



    /******************************************************************/
    /*                                                                */
    /* Start reading dataset into the memory                          */
    /*                                                                */
    /******************************************************************/
    if(b_obcold_reconstruction == TRUE)
    {

        breakdates_block = malloc(MAX_YEAR_RANGE * sizeof (int));
        if(breakdates_block == NULL)
        {
            RETURN_ERROR ("Allocating breakdates_block", FUNC_NAME, FAILURE);
        }
        snprintf (msg_str, sizeof(msg_str), "breakdatemap reading finished\n");
        LOG_MESSAGE (msg_str, FUNC_NAME)

        sampleFile = fopen(argv[4], "r");
        num_breakdatemaps = 0;
        row_count = 0;
        while (fgets(csv_row, 255, sampleFile) != NULL)
        {
            if(row_count != 0) // we skip first line because it is a header
            {
                breakdates_block[num_breakdatemaps] = atoi(strtok(csv_row, ","));
                num_breakdatemaps++;
            }
            row_count++;
        }
        fclose(sampleFile);
    }


    /******************************************************************/
    /*                                                                */
    /*       Reading dataset finished                                 */
    /*                                                                */
    /******************************************************************/

    snprintf (msg_str, sizeof(msg_str), "fopen and fread finished \n");
    LOG_MESSAGE (msg_str, FUNC_NAME)


    // for(original_row = starting_row; original_row < starting_row + n_block; original_row++)


    /******************************************************************/
    /*                                                                */
    /*        Set up outputting dataset                                */
    /*                                                                */
    /******************************************************************/
    if(b_obcold_reconstruction == TRUE)
    {
        sprintf(out_filename, "record_change_h%d_v%d_obcold.dat", current_block_x, current_block_y);
        sprintf(out_fullpath, "%s/obcold/%s", argv[3], out_filename);
        sprintf(out_fullpath_tmp, "%s/obcold/%s.part", argv[3], out_filename);
    }
    else
    {
        if(method == COLD){
            sprintf(out_filename, "record_change_h%d_v%d_cold.dat", current_block_x, current_block_y);
            sprintf(CM_filename, "CM_h%d_v%d.dat", current_block_x, current_block_y);
            sprintf(CM_date_filename, "CM_date_h%d_v%d.dat", current_block_x, current_block_y);
            sprintf(CM_direction_filename, "CM_direction_h%d_v%d.dat", current_block_x, current_block_y);
            sprintf(out_fullpath, "%s/%s", argv[3], out_filename);
            sprintf(out_fullpath_tmp, "%s/%s.part", argv[3], out_filename);
            sprintf(CM_fullpath, "%s/%s", argv[3], CM_filename);
            sprintf(CM_direction_fullpath, "%s/%s", argv[3], CM_direction_filename);
            sprintf(CM_date_fullpath, "%s/%s", argv[3], CM_date_filename);
        }else if(method == SCCD)
        {
            sprintf(out_filename, "record_change_h%d_v%d_sccd.dat", current_block_x, current_block_y);
            sprintf(out_fullpath, "%s/%s", argv[3], out_filename);
            sprintf(out_fullpath_tmp, "%s/%s.part", argv[3], out_filename);
        }
    }


    /******************************************************************/
    /*                                                                */
    /*       check if row has be processed or not                     */
    /*                                                                */
    /******************************************************************/

    if( access(out_fullpath, F_OK ) == 0 ) {
        snprintf (msg_str, sizeof(msg_str), "Block_h%d_v%d finished\n", current_block_x, current_block_y);
        LOG_MESSAGE (msg_str, FUNC_NAME);
        // jump to the next row
        return SUCCESS;
    }

    fhoutput = fopen(out_fullpath_tmp,"w");

    if(b_outputCM == TRUE){
        fhoutput_cm = fopen(CM_fullpath,"w");
        fhoutput_cm_direction = fopen(CM_direction_fullpath,"w");
        fhoutput_cm_date = fopen(CM_date_fullpath,"w");
    }

    snprintf (msg_str, sizeof(msg_str), "COLD processing for Block_h%d_v%d started\n",
              current_block_x, current_block_y);
    LOG_MESSAGE (msg_str, FUNC_NAME)

    /******************************************************************/
    /*                                                                */
    /*           Begin CCDC scanline processing                       */
    /*                                                                */
    /******************************************************************/
    for(i_col = 0; i_col < block_width * block_height; i_col++)
    {
        num_fc = 0;
        if(b_obcold_reconstruction == TRUE){
            rec_cg = malloc(NUM_FC * sizeof(Output_t));
            if (rec_cg == NULL)
            {
                RETURN_ERROR ("ERROR allocating rec_cg",
                              FUNC_NAME, FAILURE);
            }
            result = obcold_reconstruction_procedure(buf, fmask_buf, sdate,
                                                     num_scenes, rec_cg, &num_fc, num_breakdatemaps,
                                                     breakdates_block,n_cols, original_col, original_row, conse);
            // snprintf (msg_str, sizeof(msg_str), "col %d reconstruction calculation finished\n", i_col+1);
            // LOG_MESSAGE (msg_str, FUNC_NAME)
            for(i = 0; i < num_fc; i++)
            {
               result = write_output_binary(fhoutput, rec_cg[i]);
            }
            free(rec_cg);

        }
        else
        {
               if (method == COLD)
               {

                    rec_cg = malloc(NUM_FC * sizeof(Output_t));
                    if (rec_cg == NULL)
                    {
                        RETURN_ERROR ("ERROR allocating rec_cg",
                                      FUNC_NAME, FAILURE);
                    }
                    for(i = 0; i < n_cm_maps; i++){
                        CM_outputs[i] = NA_VALUE;
                        CM_outputs_date[i] = 255;
                        CMdirection_outputs[i] = 255;
                    }

        //                printf("temporal success2 with processid %d\n", process_id);
        //                printf("valid_scene_count_scanline[i_col] is %d;i_col is %d; original_row is %d; probability threshold is %f\n",
        //                       valid_scene_count_scanline[i_col], i_col, original_row, probability_threshold);
                    result = cold(buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], fmask_buf, sdate, valid_scene_count,
                                  n_cols, original_col, original_row, tcg, conse, b_outputCM,
                                  starting_date, rec_cg, &num_fc, CM_OUTPUT_INTERVAL, CM_outputs, CMdirection_outputs, CM_outputs_date);

                    // snprintf (msg_str, sizeof(msg_str), "pixel %d COLD calculation finished\n", i_col+1);
                    // LOG_MESSAGE (msg_str, FUNC_NAME)

                    for(i = 0; i < num_fc; i++)
                    {
                        result = write_output_binary(fhoutput, rec_cg[i]);
                    }

                    if(b_outputCM == TRUE){
                        nvals = fwrite (CM_outputs, sizeof(int16), n_cm_maps, fhoutput_cm);
                        if (nvals != n_cm_maps)
                        {
                            RETURN_ERROR("Incorrect amount of data written", FUNC_NAME, ERROR);
                        }
                        nvals = fwrite (CMdirection_outputs, sizeof(int8), n_cm_maps, fhoutput_cm_direction);
                        if (nvals != n_cm_maps)
                        {
                            RETURN_ERROR("Incorrect amount of data written", FUNC_NAME, ERROR);
                        }
                        nvals = fwrite (CM_outputs_date, sizeof(int8), n_cm_maps, fhoutput_cm_date);
                        if (nvals != n_cm_maps)
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

        //                result = sccd_hpc(tmp_buf_2d, fmask_buf_scanline[i_col], valid_date_array_scanline[i_col], valid_scene_count_scanline[i_col],
        //                                  s_rec_cg, &num_fc,n_cols, i_col + 1, original_row, probability_threshold, MIN_DAYS_CONSE, sensor_buf,conse);
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
            RETURN_ERROR("ccd procedure fails at Block_h%d_v%d \n", current_block_x, current_block_y);

        }
    //            int count;
     } // for(i_col = 0; i_col < block_width * block_height; i_col++)


    fclose(fhoutput);

    if(b_outputCM == TRUE){
       fclose(fhoutput_cm);
       fclose(fhoutput_cm_direction);
       fclose(fhoutput_cm_date);
    }

    if (rename(out_fullpath_tmp, out_fullpath) == 0)
    {
       snprintf (msg_str, sizeof(msg_str), "Block_h%d_v%d finished\n", current_block_x, current_block_y);
       LOG_MESSAGE (msg_str, FUNC_NAME);
    }
    else
    {
       snprintf (msg_str, sizeof(msg_str), "Block_h%d_v%d renaming failed\n", current_block_x, current_block_y);
       LOG_MESSAGE (msg_str, FUNC_NAME);
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

    free(fmask_buf);


    // free(sensor_buf);

    free(sdate);
    free(CM_outputs);
    free(CM_outputs_date);
    free(CMdirection_outputs);

    if(b_obcold_reconstruction == TRUE)
    {
        free(breakdates_block);
    }


    // free(fh);
    return SUCCESS;

}
