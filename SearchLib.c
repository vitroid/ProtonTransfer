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
    //$B=E?4$+$i$"$kDxEY0J>e%W%m%H%s$,1s$6$+$C$?$iC5:w$rCf;_$9$k!#(B
    if ( maxradius < radius0 ){
        printf("\n");
        return;
    }

    for( newpos=0; newpos<2; newpos++ ){
        int o1;
        
        o1 = search->ice->bonds[ search->ice->currentBond ].oxygen[ newpos ];
        /*$B$b$HMh$?J}8~$KLa$i$J$$$h$&$K(BlastO$B$H>H9g$9$k!#(B*/
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
                    /*$B85Mh$?F;$KLa$k$N$r6X;_$9$k(B*/
                    /*$B$5$i$K!"(Bdangling bond$B$K(Bproton$B$,>h$i$J$$$h$&$K$9$k!#%@%s%0%j%s%0%\%s%I$O!"JRC<$,%@%_!<%5%$%H!#J?@.(B16$BG/(B2$B7n(B18$BF|(B($B?e(B)$BDI2C(B*/
                    if ( ( bidir || (search->positions[ depth+1 ] != destBond && destBond != search->ice->currentBond ) ) && ! ( oo1==0 || oo2==0 )){
                        real deltae,oldEnergy, newEnergy;
                        int newDirection;
                        int lastDirection;
                        
                        /*$B%W%m%H%s0\F0$9$k@h$N7k9g$N8~$-(B($BI|5l$K;HMQ(B)*/
                        lastDirection        = search->ice->bonds[ destBond ].direction;
                        /*$B%W%m%H%s0\F0$9$k@h$r5-O?(B*/
                        search->positions[ depth-1 ] = destBond;
                        search->lastO[ depth-1 ]     = o1;
                        /*$B%W%m%H%s0\F0$9$k@h$,!"8=:_$N7k9g$N$I$A$iB&$+(B*/
                        newDirection = ( newpos == 0 ) ? FORWARD : BACKWARD;
                        
                        // $BA0<T$O8=:_$N%W%m%H%s$NMOG^OB%(%M%k%.!<!"8e<T$O0\F0@h$N?eAG$N%(%M%k%.!<(B
                        if ( ! dryrun )
			  oldEnergy = search->solvation[ depth ] + ProtonPotentialEnergy( search->ice, destBond, 0.0L );
                        //move protons
                        search->ice->bonds[ search->ice->currentBond ].direction = newDirection;
                        search->ice->bonds[ destBond ].direction    = 0;
                        
                        // $B0\F0@h$r%W%m%H%s$K$9$k(B($B7k9g$NCf1{$K0LCV(B)$B>l9g$N%(%M%k%.!<(B
                        if ( ! dryrun )
			  newEnergy = ProtonPotentialEnergy( search->ice, search->ice->currentBond, 0.0L );
                        search->ice->currentBond = destBond;
                        // $B0\F085$rDL>o$N(BHB$B$K$7$?>l9g$NMOG^OB%(%M%k%.!<(B
                        if ( ! dryrun ){
			  search->ice->currentPotential = ProtonPotentialEnergy( search->ice, search->ice->currentBond, 0.0L );
			  newEnergy += search->ice->currentPotential;
                        
			  deltae = ( newEnergy - oldEnergy );
			}

                        //gnuplot$B$G(Btree$B$rIA$/$?$a$N:Y9)!#(B
                        if ( first ){
                            first=0;
                        } else{
                            printf( "%d %22.18f %f\n", depth, pot, radius0 );
                        }
			
#ifdef DEBUG
			if ( ! dryrun )
			  printf( "DeltaSigma=%f, Absolute=%f\n", pot + deltae, TotalEnergy( search->ice ) );
#endif
                        
                        /*$B:F5"8F$S$@$7(B*/
			SearchRecursively( search, depth-1, maxradius, pot + deltae, dryrun, bidir );


                        /*$B$b$H$KLa$9!#(B*/
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

                
//$B%W%m%H%s$rF0$+$9;~$K!";@AGHV9f$r5-O?$9$k!#(B
real SearchMove( sSearch* search, int next, int depth, bool dryrun )
{
    int o1;
    
    o1     = SharedOxygen( search->ice, next );
    search->lastO[ depth ]     = o1;
    return Move( search->ice, next, dryrun );
}
