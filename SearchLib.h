#ifndef PATHSEARCH_H
#define PATHSEARCH_H

#include "LattIce.h"

typedef struct{
    int     maxdepth;
    int*    positions;
    int*    lastO;
    real* solvation;
    sLattIce* ice;
} sSearch;

extern sSearch* new_Search( sLattIce* ice, int maxdepth );

extern void SearchRecursively( sSearch* search, int depth, real maxradius, real pot, bool dryrun, bool bidir );

extern real SearchMove( sSearch* search, int next, int depth, bool dryrun );

#endif
