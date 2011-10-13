#ifndef COMMON_H
#define COMMON_H

typedef enum {
    False=0,
    True,
} bool;

#ifdef SINGLEPRECISION
#define real float
#else
#define real double
#endif


#endif
