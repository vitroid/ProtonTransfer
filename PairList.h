#ifndef PAIRLIST_H
#define PAIRLIST_H
#include "varray.h"
#include "common.h"
typedef struct sPairList
{
    struct sPairList* next;
    int i,j;
    real sqrdistance;
} sPairList;


extern sPairList* PairList(int n, real* x, real* y, real* z, real r);

#endif
