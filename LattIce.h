#ifndef PROTONXFER_H
#define PROTONXFER_H
#include "common.h"

#define MAXNEIBOR 6

typedef struct 
{
    /*charge is 0 for virtual site (on the surface of the cluster)*/
    /*charge$B$,(B0$B$N%5%$%H$K$O?eJ,;R$OCV$+$J$$!#C1$K!"?eAG7k9g$r:n$k$?$a$KCV$/%@%_!<%5%$%H(B*/
    /*charge$B$NCM<+BN$K$O0UL#$O$J$$!#(B*/
    int charge;
    /*number of adjacent sites*/
    int nAdj;
    /*ID of neighboring Oxygen atoms and spanning Hydrogen bonds*/
    int oxygen[MAXNEIBOR], hydrogenBond[MAXNEIBOR];
    //number of outgoing edges from the oxygen (should be 2 for water)
    int outgo;
    real x,y,z;
}
sOxygen;

typedef struct
{
    /*IDs of two oxygen atoms*/
    int oxygen[2];
    /*direction of the bond*/
    /*oxygen[0]->oxygen[1]$B8~$-$N?eAG7k9g(B($B?eAG$O(Boxygen[0]$B$K6a$$B&$K$/$k(B)
     *$B$r(B+1$B$H$7!"5U8~$-$r(B-1$B$H$9$k!#Cf1{$K0LCV$9$k>l9g(B($B%W%m%H%s(B)$B$O(B0$B$H$9(B
     *$B$k!#I,$:(Boxygen[0] < oxygen[1]$B$,@.$jN)$D$h$&$K$9$k(B*/
    int direction;
    /*Coordinate for each position.*/
    real px[3],py[3],pz[3];
}
sHydrogenBond;

typedef struct
{
    sOxygen*       oxygens;
    sHydrogenBond* bonds;
    int currentBond;
    real currentPotential;
    //total number of oxygen atoms
    int nOxygen;
    int nBond;
    //center-of-mass
    real comx,comy,comz;
    //box size
    bool periodic;
    real bx,by,bz;
    /*$B%(%M%k%.!<7W;;$r>JN,$9$k;~$O(B1*/
    //bool dryrun;
} sLattIce;

//

/*$B?eAG7k9g$N8~$-$K$h$k?eAG$N0LCV$N0c$$!#(B*/
/* 2007-7-12 user$B$,;XDj$9$k$h$&$KJQ99!#(B($B%G%U%)%k%HCM$O0J2<$NCM$N$^$^(B)
#define IN     0.7L
#define MID    0.5L
#define OUT    0.3L
*/

#define QH    +1.0L
//#define QO    (-2.0L * QH)
#define QP    +1.0L

#define OUTGOING +1
#define INCOMING -1
#define PROTON   0
#define FORWARD  +1
#define BACKWARD -1

//$B;@AG$rG[CV$7!"7k9g$r:n$j!"%W%m%H%s$rCV$$$F!"=i4|G[CV$r7A@.$9$k!#(B
extern sLattIce* new_ProtonatedIceIc( int n, real hposition, bool dryrun );

extern void SnapShot( sLattIce* ice, int aux );

extern real TotalEnergy( sLattIce* ice );

extern int BondDirection( sLattIce* ice, int o1, int i );

extern real ProtonPotentialEnergy( sLattIce* ice, int proton, real radius );

extern real ProtonInPCIPotentialEnergy( sLattIce* ice, int proton, real radius );
extern real ProtonInPDIPCIPotentialEnergy( sLattIce* ice, int proton, real radius );
extern real ProtonInPDIPotentialEnergy( sLattIce* ice, int proton, real radius );

extern real TotalEnergy( sLattIce* ice );

extern void CheckConsistency( sLattIce* ice, real radius );

extern real Radius( sLattIce* ice, int proton );

extern sLattIce* LoadLatticeAndProtonate( FILE* file, bool randomize, int grind, real hposition, bool dryrun );

extern sLattIce* load_LattIce( FILE* file, real hposition );

extern void DisorderProton( sLattIce* ice );

extern void CheckConvergence( sLattIce* ice );

extern real Move( sLattIce* ice, int next, bool dryrun );

extern real SharedOxygen( sLattIce* ice, int next );

/*radius > 0 : save edges inside the radius only; save all otherwise.*/
extern void saveNGPH( sLattIce* ice, real radius, FILE* file );

extern void saveProtonPosition( sLattIce* ice, FILE* file );

extern void saveSite( sLattIce* ice, FILE* file );

#endif
