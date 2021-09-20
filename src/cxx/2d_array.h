#ifndef MISC_2D_ARRAY_H
#define MISC_2D_ARRAY_H
#include <stdio.h>


void **allocate_2d_array
(
    int rows,             /* I: Number of rows for the 2D array */
    int columns,          /* I: Number of columns for the 2D array */
    size_t member_size    /* I: Size of the 2D array element */
);


int get_2d_array_size
(
    void **array_ptr,     /* I: Pointer returned by the alloc routine */
    int *rows,            /* O: Pointer to number of rows */
    int *columns          /* O: Pointer to number of columns */
);


int free_2d_array
(
    void **array_ptr     /* I: Pointer returned by the alloc routine */
);


#endif
