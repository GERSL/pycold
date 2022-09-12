#include "defines.h"
#include "utilities.h"
#include "const.h"
#include "output.h"


int convert_jday_to_year_doy_from_0001
(
    int jday,      /* I: julian date since year 0000 */
    int *year,      /* O: year */
    int *doy       /* O: day of the year */
)
{
    char FUNC_NAME[] = "convert_jday_to_year_doy_from_0001";
    int status;
    int day_count = ORDINAL_DATE_LAST_DAY_1972;
    int year_count = LANDSAT_START_YEAR - 1;

    if (jday < ORDINAL_DATE_LAST_DAY_1972)
    {
        RETURN_ERROR ("Landsat data starts from 1973", FUNC_NAME, ERROR);
    }

    while(day_count < jday)
    {
        year_count = year_count + 1;
        status = is_leap_year(year_count);
        if (status == TRUE)
        {
            day_count += LEAP_YEAR_DAYS;
        }
        else
        {
            day_count += NON_LEAP_YEAR_DAYS;
        }

    }

    if (is_leap_year(year_count))
    {
        day_count -= LEAP_YEAR_DAYS;
    }
    else
    {
        day_count -= NON_LEAP_YEAR_DAYS;
    }

    *doy = jday - day_count;
    *year = year_count;

    return (SUCCESS);
}

int firstDegradationYear(Output_t* t, int num_fc, int n_lassoband, bool updirection, int *dyear)
{
    int i;
    int degradation_doy;
    int degradation_year = 0;
    if (num_fc > 1)
    {
        for(i = 0; i < num_fc; i++)
        {
            if(updirection == TRUE)
            {
                if(t[i].magnitude[n_lassoband] > 0)
                {
                    if(convert_jday_to_year_doy_from_0001(t[i].t_break, &degradation_year, &degradation_doy) == SUCCESS)
                    {
                        *dyear = degradation_year;
                        return SUCCESS;
                    }
                    else
                    {
                        return FAILURE;
                    }

                }
            }
            else
            {
                if(t[i].magnitude[n_lassoband] < 0)
                {
                    if(convert_jday_to_year_doy_from_0001(t[i].t_break, &degradation_year, &degradation_doy) == SUCCESS)
                    {
                        *dyear = degradation_year;
                        return SUCCESS;
                    }
                    else
                    {
                        return FAILURE;
                    }
                }
            }

        }
    }

    *dyear = degradation_year;
    return SUCCESS;

}

int write_raw_binary
(
    FILE *rb_fptr,      /* I: pointer to the raw binary file */
    int nlines,         /* I: number of lines to write to the file */
    int nsamps,         /* I: number of samples to write to the file */
    int size,           /* I: number of bytes per pixel (ex. sizeof(uint8)) */
    int *img_array     /* I: array of nlines * nsamps * size to be written
                              to the raw binary file */
)
{
    int nvals;               /* number of values written to the file */
    char FUNC_NAME[] = "write_raw_binary"; /* function name */

    /* Write the data to the raw binary file */
    nvals = fwrite (img_array, size, nlines * nsamps, rb_fptr);
//    for(i = 0 ; i < nsamps; i++)
//    {
//        printf("%d\n", img_array[i]);
//    }
    if (nvals != nlines * nsamps)
    {
        RETURN_ERROR("Incorrect amount of data written", FUNC_NAME, ERROR);
    }

    return (SUCCESS);
}

void write_envi_header
(
    FILE *fptr,      /* I: pointer to the raw binary file */
    Input_meta_t *meta     /* I: saved header file info */
)
{
    fprintf(fptr, "ENVI\n");
    fprintf(fptr, "description = {CCDC outputs}\n");
    fprintf(fptr, "samples = %d\n", meta->samples);
    fprintf(fptr, "lines = %d\n", meta->lines);
    fprintf(fptr, "bands = %d\n", 1);
    fprintf(fptr, "head offset = 0\n");
    fprintf(fptr, "file type = ENVI Standard\n");
    fprintf(fptr, "data type = %d\n", 3);
    fprintf(fptr, "interleave = bip\n");
    fprintf(fptr, "map info = {UTM, 1.000, 1.000, %d, %d, %d, %d, 13, North, WGS-84, units=Meters}\n",
            meta->upper_left_x, meta->upper_left_y, meta->pixel_size, meta->pixel_size);
}


int write_output_binary
(
    FILE *fptr,      /* I: pointer to the binary file */
    Output_t t          /* I: outputted structure     */
)
{
    int nvals;               /* number of values written to the file */
    char FUNC_NAME[] = "write_output_binary"; /* function name */
    //int size = sizeof(Output_t);

    nvals = fwrite(&t, sizeof(Output_t), 1, fptr);

    if (nvals != 1)
    {
        RETURN_ERROR("Incorrect amount of data written", FUNC_NAME, ERROR);
    }

    return (SUCCESS);
}

int write_output_binary_sccd
(
    FILE *fptr,      /* I: pointer to the binary file */
    Output_sccd t          /* I: outputted structure     */
)
{
    int nvals;               /* number of values written to the file */
    char FUNC_NAME[] = "write_output_binary_sccd"; /* function name */

    nvals = fwrite(&t, sizeof(Output_sccd), 1, fptr);

    if (nvals != 1)
    {
        RETURN_ERROR("Incorrect amount of data written", FUNC_NAME, ERROR);
    }

    return (SUCCESS);
}
//int write_output_mat
//(
//    FILE *csv_fptr,      /* I: pointer to the csv file */
//    Output_t t          /* I: outputted structure     */
//)
//{
//    MATFile *pmat;
//    mxArray *coefs, *rmse, *magnitude;

//    coefs = mxCreateDoubleMatrix(8,7,mxREAL);
//    rmse = mxCreateDoubleMatrix(1,7,mxREAL);
//    magnitude = mxCreateDoubleMatrix(1,7,mxREAL);


//}
