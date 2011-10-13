#ifndef PROTONXFER_H
#define PROTONXFER_H
#include "common.h"

#define MAXNEIBOR 6

typedef struct 
{
    /*charge is 0 for virtual site (on the surface of the cluster)*/
    /*chargeが0のサイトには水分子は置かない。単に、水素結合を作るために置くダミーサイト*/
    /*chargeの値自体には意味はない。*/
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
    /*oxygen[0]->oxygen[1]向きの水素結合(水素はoxygen[0]に近い側にくる)
     *を+1とし、逆向きを-1とする。中央に位置する場合(プロトン)は0とす
     *る。必ずoxygen[0] < oxygen[1]が成り立つようにする*/
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
    /*エネルギー計算を省略する時は1*/
    //bool dryrun;
} sLattIce;

//

/*水素結合の向きによる水素の位置の違い。*/
/* 2007-7-12 userが指定するように変更。(デフォルト値は以下の値のまま)
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

//酸素を配置し、結合を作り、プロトンを置いて、初期配置を形成する。
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
