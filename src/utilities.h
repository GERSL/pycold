#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdbool.h>
#include <stdio.h>


#define LOG_MESSAGE(message, module) \
            write_message((message), (module), "INFO", \
                          __FILE__, __LINE__, stdout);


#define WARNING_MESSAGE(message, module) \
            write_message((message), (module), "WARNING", \
                          __FILE__, __LINE__, stdout);


#define ERROR_MESSAGE(message, module) \
            write_message((message), (module), "ERROR", \
                          __FILE__, __LINE__, stdout);


#define RETURN_ERROR(message, module, status) \
           {write_message((message), (module), "ERROR", \
                          __FILE__, __LINE__, stdout); \
            return (status);}


void write_message
(
    const char *message, /* I: message to write to the log */
    const char *module,  /* I: module the message is from */
    const char *type,    /* I: type of the error */
    char *file,          /* I: file the message was generated in */
    int line,            /* I: line number in the file where the message was
                               generated */
    FILE * fd            /* I: where to write the log message */
);


char *sub_string         /* explicit control of a substring function  */
(
    const char *source,  /* I: input string                           */
    size_t start,        /* I: index for start of sub string          */
    size_t length        /* I: number of characters to grab           */
);

void quick_sort(int arr[], char *brr[], int crr[], int left, int right);
int partition(int arr[], char *brr[], int crr[], int left, int right);
void quick_sort_float(float arr[], int left, int right);
void quick_sort_double(double arr[], int left, int right);
int partition_float (float arr[], int left, int right);
int partition_double (double arr[], int left, int right);
int partition_index (int arr[],  int *index, int left, int right);
void quick_sort_buf (int arr[], short int **brr, short int *fmask_buf, int left, int right);
void quick_sort_index (int arr[], int *index,  int left, int right);

int get_variables
(
    int argc,              /* I: number of cmd-line args                    */
    char *argv[],          /* I: string of cmd-line args                    */
    int *mode,              /* O: 1 - pixel-based; 2 - scanline-based;       */
                           /*    3 - wholescene                             */
    char *in_path,         /* O: directory locaiton for input data          */
    char *out_path,        /* O: directory location for output files        */
    int *n_cores,            /* O: number of cores                            */
    int *row,              /* O: (required for mode 1 and mode 2)           */
                           /*    row number for the pixel                   */
    int *col,              /* O: (required for mode 1)                      */
                           /* O: col number for the pixel                   */
    int *task,
    char *user_mask_path,   /* O: directory location for user_mask           */
    double *probability_threshold,
    int *conse,
    int *min_days_conse,
    bool *b_fastmode,
    bool *b_outputCSV
);

int get_args
(
    int argc,              /* I: number of cmd-line args                    */
    char *argv[],          /* I: string of cmd-line args                    */
    int *mode,             /* O: 1 - pixel-based; 2 - scanline-based; 3 - wholescene; 4-validation */
    char *in_path,         /* O: directory locaiton for input data          */
    char *out_path,        /* O: directory location for output files        */
    int *n_cores,          /* O: the number of cores                        */
    int *row,              /* O: row number for the pixel                   */
    int *col,              /* O: col number for the pixel                   */
    int *method,
    char *user_mask_path,
    double *probability_threshold,
    int *min_days_conse,
    bool *b_fastmode,
    bool *b_outputCSV
);

int get_classificationconfig
(
    char* var_path,
    char *xgboost_model_path,
    int *specific_label,
    char *auxiliary_var1
);

int find_index_clrx(
    int *clrx,
    int end,
    int input_ordinal_date
);

int usage_message ();

#endif /* UTILITIES_H */
