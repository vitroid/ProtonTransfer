#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "IntHash.h"
#include "PairList.h"

typedef struct sNodeList {
    int node;
    struct sNodeList* next;
} sNodeList;

typedef struct sGridInfo{
    sNodeList* nodeList;
    int x,y,z;
    struct sGridInfo* next;
} sGridInfo;

#define Serialize(x,y,z) ( ( ( (z) + 500 )*1000 + (y) + 500 )*1000 + (x) + 500 )
/*空疎3次元配列の値を参照する。*/
sGridInfo* Get(sIntHash* ih, int x, int y, int z)
{
    return (sGridInfo*) IntHash_QueryValue( ih, Serialize( x,y,z ) );
}

/*空疎3次元配列の値を代入する。*/
void Set(sIntHash* ih, int x, int y, int z, sGridInfo* value )
{
    IntHash_RegisterValue( ih, Serialize( x,y,z ), (int)value );
}

/*周期境界を用いない場合。*/
sPairList* PairList(int n, real* x, real* y, real* z, real r)
{
    int ix, iy, iz;
    sGridInfo* grid = NULL;
    sGridInfo* g;
    sIntHash*  ih;
    int i;
    sPairList* pairList = NULL;

    ih = IntHash_Init( 22 );
    for( i=0; i<n; i++ ){
        sGridInfo* gridInfo;
        sNodeList* nodeList=malloc(sizeof(sNodeList));
        nodeList->node = i;
        //サイズrのグリッドに割り振る
        ix = rint( x[i] / r );
        iy = rint( y[i] / r );
        iz = rint( z[i] / r );
        //登録する
        gridInfo = Get( ih, ix,iy,iz );
        if ( gridInfo == NULL ){
            gridInfo    = malloc( sizeof( sGridInfo ) );
            gridInfo->x = ix;
            gridInfo->y = iy;
            gridInfo->z = iz;
            gridInfo->nodeList = nodeList;
            gridInfo->next = grid;
            grid = gridInfo;
            nodeList->next = NULL;
            Set( ih, ix,iy,iz, gridInfo );
        }
        else{
            nodeList->next = gridInfo->nodeList;
            gridInfo->nodeList = nodeList;
        }
    }
    //中実な要素について、その隣接グリッド内の全サイトとペアリストを生成する。
    g = grid;
    while ( g != NULL ){
        int xx,yy,zz;
        sNodeList* n1;
        sNodeList* n2;

        for( xx=-1; xx<=1; xx++ ){
            for( yy=-1; yy<=1; yy++ ){
                for( zz=-1; zz<=1; zz++ ){
                    if ( xx==0 && yy==0 && zz==0 ) goto ext;
                    sGridInfo* g2 = Get( ih, xx + g->x, yy + g->y, zz + g->z );
                    if ( g2 ){
                        for( n1=g->nodeList; n1 != NULL; n1 = n1->next ){
                            for( n2=g2->nodeList; n2 != NULL; n2 = n2->next ){
                                int i,j;
                                real dx,dy,dz,dd;
                                i = n1->node;
                                j = n2->node;
                                dx = x[i] - x[j];
                                dy = y[i] - y[j];
                                dz = z[i] - z[j];
                                dd = dx*dx+dy*dy+dz*dz;
                                if ( dd < r*r ){
                                    sPairList* pl = malloc(sizeof(sPairList));
                                    pl->i = i;
                                    pl->j = j;
                                    pl->sqrdistance = dd;
                                    pl->next = pairList;
                                    pairList = pl;
                                }
                            }
                        }
                    }
                }
            }
        }
    ext:
        for( n1=g->nodeList; n1 != NULL; n1 = n1->next ){
            for( n2 = n1->next; n2 != NULL; n2 = n2->next ){
                int i,j;
                real dx,dy,dz,dd;
                i = n1->node;
                j = n2->node;
                dx = x[i] - x[j];
                dy = y[i] - y[j];
                dz = z[i] - z[j];
                dd = dx*dx+dy*dy+dz*dz;
                if ( dd < r*r ){
                    sPairList* pl = malloc(sizeof(sPairList));
                    pl->i = i;
                    pl->j = j;
                    pl->sqrdistance = dd;
                    pl->next = pairList;
                    pairList = pl;
                }
            }
        }
        g = g->next;
    }
    //メモリーの開放。これにかなり時間を喰われている。
    g = grid;
    while ( g != NULL ){
        sGridInfo* next=g->next;
        sNodeList* n=g->nodeList;
        while ( n!=NULL ){
            sNodeList* next=n->next;
            free(n);
            n=next;
        }
        free(g);
        g=next;
    }
    IntHash_Done( ih );
    return pairList;
}

void PairListCheck( int n, real* x, real* y, real* z, real r )
{
    int i,j;
    for( i=0; i<n; i++ ){
        for( j=i+1; j<n; j++ ){
            real dx,dy,dz,dd;
            dx = x[i] - x[j];
            dy = y[i] - y[j];
            dz = z[i] - z[j];
            dd = dx*dx+dy*dy+dz*dz;
            if ( dd < r*r ){
                printf(":%d %d %f\n", i,j, dd );
            }
        }
    }
}

#define N 1000000
int test_main( int argc, char* argv[] )
{
    int i;
    real x[N],y[N],z[N];
    for( i=0; i<N; i++ ){
        x[i] = drand48() * 1000;
        y[i] = drand48() * 1000;
        z[i] = drand48() * 1000;
        //printf( "%d %f %f %f\n", i, x[i], y[i], z[i] ); 
    }
    
    PairList( N, x,y,z, 10.0 );
    //PairListCheck( N, x,y,z, 10.0 );
    exit(0);
}
