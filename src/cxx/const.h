
#ifndef CONST_H
#define CONST_H

#include <math.h>

typedef signed short int16;
typedef unsigned char uint8;
typedef signed char int8;

#ifndef min
    #define min(a,b) (((a) < (b)) ? (a) : (b))
#endif


#ifndef max
    #define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef max
    #define max(a,b) (((a) > (b)) ? (a) : (b))
#endif



#define TOTAL_BANDS 8

#define PI 3.1415926535897935
#define TWO_PI (2.0 * PI)
#define HALF_PI (PI / 2.0)


#define DEG (180.0 / PI)
#define RAD (PI / 180.0)


#ifndef SUCCESS
    #define SUCCESS  0
#endif

#ifndef ERROR
    #define ERROR -1
#endif


#ifndef FAILURE
    #define FAILURE 1
#endif

#ifndef TRUE
    #define TRUE 1
#endif

#ifndef FALSE
    #define FALSE 0
#endif

#ifndef INCOMPLETE
    #define INCOMPLETE 2
#endif


#define MINSIGMA 1e-5

#define MAX_STR_LEN 512
#define MAX_SCENE_LIST 3922
#define MAX_YEAR_RANGE 50
#define ARD_STR_LEN 50

#define MAXIMUM_INI_LEN 400
#endif
