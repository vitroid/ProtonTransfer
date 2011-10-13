#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SearchLib.h"
#include "Utility.h"

sSearch* new_Search( sLattIce* ice, int maxdepth )
{
    sSearch* search = malloc(sizeof(sSearch) );
    
    search->maxdepth = maxdepth;
    search->ice      = ice;
    
    search->positions = calloc( maxdepth+2, sizeof(int) );
    search->lastO     = calloc( maxdepth+2, sizeof(int) );
    search->solvation = calloc( maxdepth+2, sizeof(real) );

    search->positions[maxdepth+1] = 0;
    search->positions[maxdepth]   = search->ice->currentBond;
    search->lastO[maxdepth]       = -1;

    return search;
}


#ifdef DEBUG
static real etotal;
#endif

void SearchRecursively( sSearch* search, int depth, real maxradius, real pot, bool dryrun, bool bidir )
{
    int newpos,j,i;
    real radius0;
    
    search->positions[ depth ] = search->ice->currentBond;
    search->solvation[ depth ] = search->ice->currentPotential;

#ifdef DEBUG
    printf("(%f)", TotalEnergy( search->ice ) - etotal );
#endif
    printf( "%d %22.18f ", depth, pot );
    radius0 = Radius( search->ice, search->positions[ depth ] );
    printf( "%f ", radius0 );
    for( i=search->maxdepth; depth<=i; i-- ){
        printf("%d ", search->positions[i] );
    }
    //show distance from the center of mass
    printf( "\n" );
#ifdef DEBUG
    CheckConsistency( search->ice );
#endif
    
    if ( depth == 0 ){
        printf("\n");
        return;
    }
    //重心からある程度以上プロトンが遠ざかったら探索を中止する。
    if ( maxradius < radius0 ){
        printf("\n");
        return;
    }

    for( newpos=0; newpos<2; newpos++ ){
        int o1;
        
        o1 = search->ice->bonds[ search->ice->currentBond ].oxygen[ newpos ];
        /*もと来た方向に戻らないようにlastOと照合する。*/
        //printf("%d %d %d\n", search->ice->oxygens[ o1 ].charge,search->lastO[ depth ], o1 );
        
        if ( bidir || (search->ice->oxygens[ o1 ].charge != 0 && search->lastO[ depth ] != o1) ){
            int first=1;
            for( j=0; j<search->ice->oxygens[o1].nAdj; j++ ){
                if ( BondDirection( search->ice, o1, j ) == OUTGOING ) {
                    int destBond = search->ice->oxygens[ o1 ].hydrogenBond[ j ];
                    int oo1 = search->ice->bonds[ destBond ].oxygen[0];
                    int oo2 = search->ice->bonds[ destBond ].oxygen[1];
                    oo1 = search->ice->oxygens[ oo1 ].charge;
                    oo2 = search->ice->oxygens[ oo2 ].charge;
                    /*元来た道に戻るのを禁止する*/
                    /*さらに、dangling bondにprotonが乗らないようにする。ダングリングボンドは、片端がダミーサイト。平成16年2月18日(水)追加*/
                    if ( ( bidir || (search->positions[ depth+1 ] != destBond && destBond != search->ice->currentBond ) ) && ! ( oo1==0 || oo2==0 )){
                        real deltae,oldEnergy, newEnergy;
                        int newDirection;
                        int lastDirection;
                        
                        /*プロトン移動する先の結合の向き(復旧に使用)*/
                        lastDirection        = search->ice->bonds[ destBond ].direction;
                        /*プロトン移動する先を記録*/
                        search->positions[ depth-1 ] = destBond;
                        search->lastO[ depth-1 ]     = o1;
                        /*プロトン移動する先が、現在の結合のどちら側か*/
                        newDirection = ( newpos == 0 ) ? FORWARD : BACKWARD;
                        
                        // 前者は現在のプロトンの溶媒和エネルギー、後者は移動先の水素のエネルギー
                        if ( ! dryrun )
			  oldEnergy = search->solvation[ depth ] + ProtonPotentialEnergy( search->ice, destBond, 0.0L );
                        //move protons
                        search->ice->bonds[ search->ice->currentBond ].direction = newDirection;
                        search->ice->bonds[ destBond ].direction    = 0;
                        
                        // 移動先をプロトンにする(結合の中央に位置)場合のエネルギー
                        if ( ! dryrun )
			  newEnergy = ProtonPotentialEnergy( search->ice, search->ice->currentBond, 0.0L );
                        search->ice->currentBond = destBond;
                        // 移動元を通常のHBにした場合の溶媒和エネルギー
                        if ( ! dryrun ){
			  search->ice->currentPotential = ProtonPotentialEnergy( search->ice, search->ice->currentBond, 0.0L );
			  newEnergy += search->ice->currentPotential;
                        
			  deltae = ( newEnergy - oldEnergy );
			}

                        //gnuplotでtreeを描くための細工。
                        if ( first ){
                            first=0;
                        } else{
                            printf( "%d %22.18f %f\n", depth, pot, radius0 );
                        }
			
#ifdef DEBUG
			if ( ! dryrun )
			  printf( "DeltaSigma=%f, Absolute=%f\n", pot + deltae, TotalEnergy( search->ice ) );
#endif
                        
                        /*再帰呼びだし*/
			SearchRecursively( search, depth-1, maxradius, pot + deltae, dryrun, bidir );


                        /*もとに戻す。*/
                        destBond = search->positions[ depth-1 ];
                        search->ice->currentBond = search->positions[ depth ];
                        
                        search->ice->bonds[ search->ice->currentBond ].direction = 0;
                        search->ice->bonds[ destBond ].direction    = lastDirection;

                        search->ice->currentPotential = search->solvation[ depth ];
                    }
                }
            }
        }
    }
}

                
//プロトンを動かす時に、酸素番号を記録する。
real SearchMove( sSearch* search, int next, int depth, bool dryrun )
{
    int o1;
    
    o1     = SharedOxygen( search->ice, next );
    search->lastO[ depth ]     = o1;
    return Move( search->ice, next, dryrun );
}
