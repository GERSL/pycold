#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <stdio.h>
#include <sys/timeb.h>
#include <sys/time.h>
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
//int sccd_executor(
//    int mode,
//    char* in_path,
//    char* out_path,
//    int n_cores,
//    int row,
//    int col,
//    int METHOD,
//    char* mask_path,
//    double probability_threshold,
//    int min_days_conse,
//    int output_mode,
//    int verbose,
//    int training_type, /* for training process*/
//    int monitorwindow_lowerlin, /* for training process*/
//    int monitorwindow_upperlim, /* for training process*/
//    int bandselection_bit,
//    char* classification_config /* classification config file */
//)
{
    int result;
    /* inputted argument, exe mode */
    int mode;                           /* CCD detection mode
                                        3: whole images; 1: pixel-based;
                                        2: scanline-based*/
    char in_path[MAX_STR_LEN];
    char out_path[MAX_STR_LEN];
    int n_cores;
    int row;
    int col;
    int task;
    char mask_path[MAX_STR_LEN];
    int min_days_conse;
    double probability_threshold;
    int conse;
    int output_mode;
    bool verbose = TRUE;
    bool b_fastmode;
    bool b_outputCSV;           /* output point-based time series csv, only for debug   */


    //printf("nonono");

    /**************************************************************/
    /*                                                            */
    /*   read CCD variable                                        */
    /*                                                            */
    /**************************************************************/
    /* need to recover for exe */
    result = get_variables(argc, argv, &mode, in_path, out_path, &n_cores,
                           &row, &col, &task, mask_path, &probability_threshold, conse,
                           &min_days_conse, &b_fastmode, &b_outputCSV);

    tsalgorithm_executor(mode, in_path, out_path, n_cores, row, col, task, mask_path, probability_threshold,
                         conse, min_days_conse,output_mode,verbose);

    return SUCCESS;
}
