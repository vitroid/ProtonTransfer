/*
  2008-2-2 LattIce.c --> LattIce2.c
  cluster$B30<~$N?eAG7k9g$r!"(BLattIce.c$B$G$O!"$^$8$a$KJ,;RC10L$G07$C$F$$$k!#(B
  $B$D$^$j!"308~$-$N7k9g$O3J;R$+$iFM$-=P$F$$$k$,!"5U8~$-$N7k9g$O4"$j$H$C(B
  $B$F$"$k!#$3$l$r(BPCI$B$K%U%#%C%H$5$;$h$&$H$9$k$H!"(BPCI$B$b$3$N%/%i%9%?!<$N30(B
  $B<~$N7A$K$"$o$;$F4"$j<h$i$6$k$r$($:!"Hs>o$K$H$j$"$D$+$$$,$d$d$3$7$$!#(B
  $B$3$N$h$&$J07$$$K$7$?M}M3$O!"%W%m%H%s$,%/%i%9%?!<I=LL$^$GMh$?$H$-$KIT(B
  $B@09g$,$*$3$i$J$$$h$&$K$9$k$?$a$@$C$?$N$@$,!"<B:]$K$O8=:_$N%W%m%0%i%`(B
  $B$G$O%W%m%H%s$OI=LL$KMh$J$$$h$&M^;_$5$l$F$$$k$N$G!"LdBj$K$J$i$J$$!#(B

  LattIce2.c$B$G$O!"%@%_!<$NFb8~$-?eAG7k9g$rDI2C$9$k$3$H$G!"(BPDI-PCI$B6-3&(B
  $B$N$G$3$\$3$r$J$/$7!"<h$j07$$$r4JC1$K$9$k!#(B
 */
/*
 * $B%W%m%H%s0\F0$N%b%s%F%+%k%m!#%/%i%9%?$r:n$j!"Cf?4$K%W%m%H%s$rCV$/!#(B
 $B%(%M%k%.!<$r7W;;$7!"%W%m%H%s0\F0$r<!!9$K9T$&!#3J;R$NBg$-$5$O$G$-$k$@(B
 $B$1Bg$-$/$G$-$k$h$&$K$9$k!#(B

 MonteCarlo$B$N4pK\E*$J%9%F%C%W$O(B

 1. $B8=:_$NG[CV$+$i!"?7$7$$G[CV$r@8@.(B
 2. $B?7$7$$G[CV$N%(%M%k%.!<$r7W;;(B
 3. accept/reject$B$NI>2A(B(energy$B$r4p=`$K(B)
 4a. accept$B$J$i8E$$G[CV$r4~5Q$7$F?7$7$$G[CV$r:NMQ(B
 4b. reject$B$J$i?7$7$$G[CV$r4~5Q!#(B

 $B%W%m%H%s$N0LCV$,JQ$o$k$H!";@AG$N0LCV$bHyL/$KJQ$o$k$,!":#2s$OLO7?7W;;(B
 $B$H$$$&$3$H$G!";@AG$O>o$K3J;RE@$K!"?eAG$O>o$K(BO-O$B@~>e$K$"$k$b$N$H$9$k!#(B

 $B?eAG$N0LCV$OFsCM$GI=$;$k$+$i!"J,;R?t$,(BN$B$@$H!"?eAG7k9g$NK\?t$O(B2N$B!#(B2N$B%S%C(B
 $B%H$N%a%b%j$G$H$j$"$($:F0$/$O$:!#(B

 $B0l2s$N%W%m%H%s0\F0$G(B2$B8D$N%W%m%H%s$N0LCV$,JQ$o$k!#$=$l$K$h$C$F(B2 x 3N$BBP(B
 $B$NAj8_:nMQ$N:F7W;;$,I,MW!#$?$@$7!";@AG(B-$B?eAG4V%(%M%k%.!<$O$"$i$+$8$a7W(B
 $B;;$7$F$*$1$k!#(BH-H$B4V%(%M%k%.!<$O$"$i$+$8$a7W;;$7$F$?$a$F$*$/$K$O?t$,B?(B
 $B$9$.$k$N$G!"Kh2s:F7W;;$9$k$b$N$H$9$l$P!"Aj8_:nMQ$N7W;;$OA4It$G(B2 x 2 x 
 2N=8N$BBP$H$J$k!#(B

 $B;@AG$N0LCV$KDL$7HV9f$rM?$($k!#?eAG7k9g$r$J$9(B2$B$D$N;@AG$N4V$N?eAG$N0LCV(B
 $B$O!"DL$7HV9f$,>.$5$$;@AG$K6a$$0LCV$r(B-1$B!"Bg$-$$;@AGB&$r(B+1$B$H$7!"Cf1{(B($B%W(B
 $B%m%H%s(B)$B$r(B0$B$H$9$k!#(B

 $B?eAG$N0LCV$O!"(B2$B$D$N;@AG$NHV9f$H!"0LCVCM$NAH$_$"$o$;$G:BI8$,7h(B
 $B$^$k!#?eAG$K$bDL$7HV9f$rM?$($k!#(B

 $B$"$k;@AG$KNY@\$9$k;@AG(B4$B$D$r0lMw$K$7$F$*$/!#$^$?$=$N4V$N?eAG$NHV9f$b0l(B
 $BMw$K$7$F$*$/!#(B

 $B%/%i%9%?I=LL$K0LCV$9$k?eJ,;R$O!"7k9g$r(B4$B$D$b$?$J$$$N$G!"(Bdummy$B%5%$%H$H(B
 $B$7$F07$$!"?eAG7k9g$NDj5A$K$@$1;HMQ$9$k!#(B

 $BAj8_:nMQ%(%M%k%.!<$NC10L$O!"3J;R>e$N5wN%(B1$BN%$l$?(B2$B$D$N(B($B?eJ,;R$N(B)$B?eAG$N(B
 $B4V$N%(%M%k%.!<$r(B1$B$H$9$k!#;@AG$NEE2Y$O?eAG$N(B2$BG\$H$9$k!#%W%m%H%s$@$1$O(B
 $B0[$J$kEE2Y(B(QP)$B$r;}$D$b$N$H$9$k!#(B

 $B?eAG$N0LCV$O!"$=$l$>$lF10lJ,;RFb!"%W%m%H%s!"NY@\J,;RFb$K$"$k$H$-!";@(B
 $BAG;@AG4V$ND>@~>e$N!"(B0.3, 0.5, 0.7$B$N0LCV$K$*$/!#<B:]$NI9$O!"(BO-O$B4V5w(B
 $BN%(B2.75A$B!"(BO-H$B5wN%(B1.00A$B!#HfN($H$7$F$O(B0.36$B$0$i$$$K$J$k!#9b05$NI9$J$I(B
 $B$r%7%_%e%l!<%7%g%s$9$k2DG=@-$b$"$k$N$G!"$3$N?t;z$O;XDj$G$-$k$[$&$,$$$$!#(B

 $B%W%m%H%s$,IUCe$7$F$$$k?eJ,;R$G$bDL>o$N?eJ,;R$G$b!"%W%m%H%s$r=|$/?eAG(B
 $B$N?t$O(B2$B$GJQ$o$j$J$$$+$i!"0lJ,;R$"$?$j$N=P7k9g?t(B(outgo)$B$O%W%m%H%s0\F0(B
 $B$D$l$F99?7$9$kI,MW$O$J$$!#$?$@!"D9$/%7%_%e%l!<%7%g%s$7$F$$$k$&$A$K!"(B
 $B%W%m%H%s$,%/%i%9%?$NI=LL$K$-$F$7$^$&$H$^$:$$!#(BILA($BL58B3J;R6a;w(B)$B$G(B
 $BBP=h$9$k!#(B
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LattIce.h"
#include "PairList.h"
#include "Utility.h"

/*infinite lattice approximation*/
//ILA$B$r;H$&$+$I$&$+$O(BMakefile$B$G;XDj$9$k!#(B
//#define ILA

/**/

/*ice Ic$B$G$N:BI8$H(BID$B$N4X78(B*/

/*cubic ice$B$N9=B$$O!"LL?4N)J}3J;R$r(B2$B$D$:$i$;$F$+$5$M$?9=B$$K$J$C$F$$$k!#(B
  $BLL?4N)J}3J;R$b$^$?!"C1=c3J;R$r(B3$BJ}8~$K$:$i$;$F=E$M$?9=B$$K$J$C$F$$$k!#(B
  10x10x10=1000$BJ,;R$N(Bcubic ice$B9=B$$r:n$k>l9g!"$^$:(B5x5x5$B$NC1=c3J;R$r$D(B
  $B$/$j!"3J;RE@$K(B0..124$BHV$r3d$j$U$k!#<!$K!"(Byz$BJ?LL>e$G$:$i$7$?C1=c3J;R$K(B
  $B$O!"$b$H$NC1=c3J;R$K(B125$B$rB-$7$?HV9f$r!"(Bxz$BJ?LL>e$G$:$i$7$?C1=c3J;R$K(B
  $B$O(B250$B$rB-$7$?HV9f$r!"(Byz$BJ?LL$G$:$i$7$?C1=c3J;R$K$O(B375$B$r2C;;$9$k!#$5$i(B
  $B$K!"LL?4N)J}3J;RA4BN$r(B(1,1,1)$BJ}8~$K$:$i$7$?I{3J;R$K$O(B500$B$r2C;;$9$k!#(B

  $B$d$d$3$7$$$,!"3J;RHV9f$N?6$j$+$?$O$3$N(Broutine$B$K$N$_0MB8$7$F$$$k!#(B
  $B3F3J;RE@$K0l0U$KHV9f$,$D$/$J$i!"B>$NJ}K!$G$b9=$o$J$$!#(B

*/
int Coordinate2OxygenID( int n, int x, int y, int z )
{


    int xh,yh,zh, xq,yq,zq;
    int id;
    
    /*$B<g3J;R(B($BC1=c3J;R(B)$B$NHV9f(B*/
    xq = x/4;
    yq = y/4;
    zq = z/4;
    id = (zq * (n/2) + yq) * (n/2) + xq;
    
    xh = x/2;
    yh = y/2;
    zh = z/2;
    
    if ( xq*2 == xh && yq*2 == yh && zq*2 == zh ){
    }
    else
        if ( xq * 2 == xh )
            id += n*n*n*1/8;
        else if ( yq * 2 == yh )
            id += n*n*n*2/8;
        else if ( zq * 2 == zh )
            id += n*n*n*3/8;

    if ( xh * 2 != x ) id += n*n*n/2;
    
#ifdef DEBUG
    if ( id == 150 ){
        printf("%d %d %d %d %d\n", x,y,z,id,n);
    }
#endif
    
    //if ( 999 < id )
    //    printf("ID! %d\n", id );
    
    return id;
}


/*o1$BHVL\$N;@AG$N(Bi$BHVL\$N%\%s%I$,=P7k9g$J$i(B1$B!"F~7k9g$J$i(B-1, $B<+M3M[;R$J$i(B0*/
int BondDirection( sLattIce* const ice, int o1, int i )
{
    int bond = ice->oxygens[o1].hydrogenBond[ i ];
    int dir = ( ice->bonds[ bond ].oxygen[0] == o1 ) ? +1 : -1;
    dir *= ice->bonds[ bond ].direction;
    return dir;
}



/*o1$BHVL\$N;@AG$+$i=P$F$$$k7k9g$NK\?t$r?t$($k!#(B2$B$J$iDL>o$N?eJ,;R!"(B3$B$J$i%W%m%H%s$,IUCe$7$??eJ,;R!#(B*/
int OutBonds( sLattIce* const ice, int o1 )
{
    int i;
    int det;
    
    det = 0;
    
    for ( i=0; i<ice->oxygens[o1].nAdj; i++ ){
        det += ( BondDirection( ice, o1, i ) == OUTGOING );
    }
    
    return det;
}


void PurgeDefects( sLattIce* ice )
{
    int ndefect;
    int* defect;
    int o1;
    
    defect = malloc( sizeof( int ) * ice->nOxygen );
    
    /*initialize*/
    for ( o1=0; o1<ice->nOxygen; o1++ ){
        defect[o1] = o1;
    }
    ndefect = ice->nOxygen;

    /*Move them randomly*/
    while ( 0 < ndefect ) {
        //choose one
        int i = lrand48() % ndefect;
        o1 = defect[i];
#ifdef DEBUG
        printf("ndefect=%d\n", ndefect );
#endif
        //charge==0$B$N%5%$%H$O%@%_!<%5%$%H(B($B%/%i%9%?!<I=LL$r=$>~$9$k$?$a$KIU2C$7$?%5%$%H(B)
        if ( 2 == ice->oxygens[o1].outgo || 0 == ice->oxygens[o1].charge ){
            // $B7g4Y$G$J$$%5%$%H$O8+IU$1$7$@$$%j%9%H$+$i>C$9!#(B
            ndefect --;
            defect[i] = defect[ndefect];
        }
        //o1$B$,(Bdummy site$B$G$O$J$/!"=P7k9g$,B?$$>l9g(B
        else if ( 2 < ice->oxygens[o1].outgo ){
            int bond;
            int o2;
            //$B%i%s%@%`$K(Bo1$B$+$i$N=P7k9g$r$5$,$9(B
            while ( 1 ){
                bond = lrand48() % ice->oxygens[o1].nAdj;
                if ( BondDirection( ice, o1, bond ) == OUTGOING )
                    break;
            }
            o2 = ice->oxygens[o1].oxygen[bond];
            ice->oxygens[o1].outgo --;
            //$B7k9g$NH?E>(B
            ice->bonds[ ice->oxygens[o1].hydrogenBond[bond] ].direction = 
                - ice->bonds[ ice->oxygens[o1].hydrogenBond[bond] ].direction;
            // if o2 is not dummy site
            if ( ice->oxygens[o2].charge != 0 ){
                // if o2 is originally a non-defect site,
                if ( ice->oxygens[o2].outgo == 2 ){
                    //newly add to defect list
                    defect[ndefect++] = o2;
                }
                //ice->oxygens[o2].outgo ++;
            }
            //$B$3$A$i$K0\$7$?!#(Bdummy site$B$G$"$C$F$b!"0l1~=P7k9g?t$O?t$($F$*$/!#(B
            //$B$?$@$7(Bdummy site$B$O(Bdefect$B$K$O$J$i$J$$!#(B
            ice->oxygens[o2].outgo ++;
        }
        //o1$B$,(Bdummy site$B$G$O$J$/!"=P7k9g$,B-$j$J$$>l9g(B
        else {
            int bond;
            int o2;
            //$B%i%s%@%`$K(Bo1$B$+$i$NF~7k9g$r$5$,$9(B
            while ( 1 ){
                bond = lrand48() % ice->oxygens[o1].nAdj;
                if ( BondDirection( ice, o1, bond ) == INCOMING )
                    break;
            }
            o2 = ice->oxygens[o1].oxygen[bond];
            ice->oxygens[o1].outgo ++;
            //$B7k9g$NH?E>(B
            ice->bonds[ ice->oxygens[o1].hydrogenBond[bond] ].direction = 
                - ice->bonds[ ice->oxygens[o1].hydrogenBond[bond] ].direction;
            // if o2 is not dummy site
            if ( ice->oxygens[o2].charge != 0 ){
                //if o2 is originally a non-defect site
                if ( ice->oxygens[o2].outgo == 2 ){
                    //newly add to defect list
                    defect[ndefect++] = o2;
                }
                //ice->oxygens[o2].outgo --;
            }
            ice->oxygens[o2].outgo --;
        }
    }
    free( defect );
}



/*$B?eAG$N5o>l=j(BH-Bin$B$r=`Hw$9$k(B*/
void SetHBin( sLattIce* ice, int h, real hposition )
{
    int o1,o2;
    real dx,dy,dz;
    o1 = ice->bonds[h].oxygen[0];
    o2 = ice->bonds[h].oxygen[1];
    
    dx = ice->oxygens[o2].x - ice->oxygens[o1].x;
    dy = ice->oxygens[o2].y - ice->oxygens[o1].y;
    dz = ice->oxygens[o2].z - ice->oxygens[o1].z;
    if ( ice->periodic ){
        dx -= ice->bx * rint( dx / ice->bx );
        dy -= ice->by * rint( dy / ice->by );
        dz -= ice->bz * rint( dz / ice->bz );
    }
    
    ice->bonds[h].px[0] = ice->oxygens[o1].x + (1-hposition) * dx;
    ice->bonds[h].py[0] = ice->oxygens[o1].y + (1-hposition) * dy;
    ice->bonds[h].pz[0] = ice->oxygens[o1].z + (1-hposition) * dz;
    ice->bonds[h].px[1] = ice->oxygens[o1].x + 0.5 * dx;
    ice->bonds[h].py[1] = ice->oxygens[o1].y + 0.5 * dy;
    ice->bonds[h].pz[1] = ice->oxygens[o1].z + 0.5 * dz;
    ice->bonds[h].px[2] = ice->oxygens[o1].x + hposition * dx;
    ice->bonds[h].py[2] = ice->oxygens[o1].y + hposition * dy;
    ice->bonds[h].pz[2] = ice->oxygens[o1].z + hposition * dz;
}



//-1, 0, +1$B$N>l=j$K$"$k;~$N:BI8$r$"$i$+$8$a7W;;$7$F$*$/!#(B
void SetHBins( sLattIce* ice, real hposition )
{
    int h;
    for( h=0; h<ice->nBond; h++ ){
      SetHBin( ice, h, hposition );
    }
}

    

//$B;@AG$N=E?4$r5a$a$k!#(B
void CenterOfMass( sLattIce* ice )
{
    int o1;
    ice->comx=0;
    ice->comy=0;
    ice->comz=0;
    if ( ! ice->periodic ) {
        for( o1=0; o1<ice->nOxygen; o1++ ){
            ice->comx += ice->oxygens[o1].x;
            ice->comy += ice->oxygens[o1].y;
            ice->comz += ice->oxygens[o1].z;
        }
        ice->comx /= ice->nOxygen;
        ice->comy /= ice->nOxygen;
        ice->comz /= ice->nOxygen;
    }
}



/*create N x N x N cubic ice */
sLattIce* CubicIce( int n, real hposition )
{
    int x,y,z;
    int nhb;
    int o1;
    
    sLattIce* ice;
    
    ice = malloc( sizeof(sLattIce) );

    ice->oxygens = calloc( n*n*n,   sizeof( sOxygen ) );
    ice->bonds   = calloc( n*n*n*2, sizeof( sHydrogenBond ) );
    ice->nOxygen = n*n*n;
    ice->periodic = False;
    //ice->dryrun  = False;
    
    for ( x=0; x<2*n; x+=2 ){
        for ( y=0; y<2*n; y+=2 ){
            for ( z=0; z<2*n; z+=2 ){
                if ( ( ( x+y+z ) % 4 ) == 0 ){
                    int xi,yi,zi;
                    //$B0LCV$+$i;@AG(BID$B$rF@$k(B
                    int self = Coordinate2OxygenID( n, x, y, z );
                    
                    if ( ice->oxygens[self].x != 0 || ice->oxygens[self].y != 0
                         || ice->oxygens[self].z != 0 ){
                        printf("%d %d %d %d is used\n", self,x,y,z);
                    }
                    //$B0LCV$rJ]B8$9$k(B
                    ice->oxygens[self].x = x;
                    ice->oxygens[self].y = y;
                    ice->oxygens[self].z = z;

                    //$BN)J}BNI=LL$K$"$k$J$i(B
                    if ( x == 0 || x == 2*n-2 ||
                         y == 0 || y == 2*n-2 ||
                         z == 0 || z == 2*n-2 )
                      //$B%@%_!<%5%$%H$G$"$k!#(B
		      ice->oxygens[self].charge = 0;
                    else{
                      //$B?eAG$NEE2Y$N(B-2$BG\!#(B
		      ice->oxygens[self].charge = -2;
                      ice->oxygens[self].oxygen[0] =
                        Coordinate2OxygenID( n, x+1, y+1, z+1 );
                      ice->oxygens[self].oxygen[1] =
                        Coordinate2OxygenID( n, x+1, y-1, z-1 );
                      ice->oxygens[self].oxygen[2] =
                        Coordinate2OxygenID( n, x-1, y+1, z-1 );
                      ice->oxygens[self].oxygen[3] =
                        Coordinate2OxygenID( n, x-1, y-1, z+1 );
                      ice->oxygens[self].nAdj = 4;
                    }
                    
                    xi = x + 1;
                    yi = y + 1;
                    zi = z + 1;
                    
                    self = Coordinate2OxygenID( n, xi, yi, zi );
                    ice->oxygens[self].x = xi;
                    ice->oxygens[self].y = yi;
                    ice->oxygens[self].z = zi;

                    if ( xi == 2*n-1 ||
                         yi == 2*n-1 ||
                         zi == 2*n-1 )
                        ice->oxygens[self].charge = 0;
                    else{
                        ice->oxygens[self].charge = -2;
                        ice->oxygens[self].oxygen[0] =
                            Coordinate2OxygenID( n, xi-1, yi-1, zi-1 );
                        ice->oxygens[self].oxygen[1] =
                            Coordinate2OxygenID( n, xi-1, yi+1, zi+1 );
                        ice->oxygens[self].oxygen[2] =
                            Coordinate2OxygenID( n, xi+1, yi-1, zi+1 );
                        ice->oxygens[self].oxygen[3] =
                            Coordinate2OxygenID( n, xi+1, yi+1, zi-1 );
                        ice->oxygens[self].nAdj = 4;
                    }
                }
            }
        }
    }
    
    CenterOfMass( ice );
    /*$B$D$.$K?eAG7k9g$rDj5A$9$k!#(B*/
    nhb = 0;
    /*$B$9$Y$F$N;@AG86;R$K$D$$$F(B*/
    for ( o1=0; o1<n*n*n; o1++ ){
        /*$B$b$7%@%_!<%5%$%H(B($BI=LL%5%$%H(B)$B$G$J$1$l$P(B*/
        if ( ice->oxygens[o1].charge != 0 ){
            int j;
            /*4$B$D$NNY@\;@AG$K$D$$$F(B*/
            for( j=0; j<ice->oxygens[o1].nAdj; j++ ){
                int o2;
                
                o2 = ice->oxygens[o1].oxygen[j];
                /*o1<o2$B$J$i!"%\%s%I$rEPO?$9$k(B($B=EJ#EPO?$r$5$1$k$?$a(B)*/
                if ( o1 < o2 ) {
                    /*bond$B$N(B2$B$D$NC<E@$rEPO?$9$k(B*/
                    ice->bonds[nhb].oxygen[0] = o1;
                    ice->bonds[nhb].oxygen[1] = o2;
                    /*$B$3$N7k9g$r(Bo2$B$X$N7k9g$H$7$FEPO?$9$k!#(B*/
                    ice->oxygens[o1].hydrogenBond[j] = nhb;
                    /*$BNY@\;@AG$,%@%_!<%5%$%H$G$J$$$J$i(B*/
                    if ( ice->oxygens[o2].charge != 0 ){
                        int k;
                    
                        /*$BNY@\;@AG$NNY@\;@AG$K$D$$$F(B*/
                        for ( k=0; k<ice->oxygens[o2].nAdj; k++ )
                            /*$B$=$l$,(Bo1$B$J$i(B*/
                            if ( ice->oxygens[o2].oxygen[k] == o1 )
                                /*$B$3$N7k9g$r(Bo2$B$+$i(Bo1$B$X$N7k9g$H$7$FEPO?$9$k(B*/
                                ice->oxygens[o2].hydrogenBond[k] = nhb;
                    }
                    
                    /*$BEPO?7k9g?t$rA}$d$9(B*/
                    nhb ++;
                }
                /*$BNY@\;@AG$,%@%_!<%5%$%H$J$i(B*/
                else if ( ice->oxygens[o2].charge == 0 ) {
                    /*o2 < o1$B$G$b%\%s%I$rEPO?$9$k!#(B*/
                    ice->bonds[nhb].oxygen[1] = o1;
                    ice->bonds[nhb].oxygen[0] = o2;
                    ice->oxygens[o1].hydrogenBond[j] = nhb;
                    nhb ++;
                }
            }
        }
    }

    //$B7k9gAm?t$K$O!"%@%_!<%\%s%I$X$N7k9g$,4^$^$l$k$N$G!";@AG(B($BEE2Y$"$j(B)$B$NAm?t$N(B2$BG\$h$j$bBg$-$$$O$:!#(B
    ice->nBond = nhb;
    SetHBins( ice, hposition );
    
    return ice;
    
}




/*$B%W%m%H%s(B($B>l=j$O(Bproton)$B$NMOG^OB%(%M%k%.!<$r7W;;$9$k!#(B*/
/*radius$B$O!"Aj8_:nMQ7W;;$rBG$A@Z$kH>7B$r;XDj$9$k!#(B0$B$r;XDj$9$k$HBG$A@Z(B
  $B$i$J$$!#DL>o$N7W;;$G$OBG$A@Z$j$OI,MW$J$$$,!"<}B+%F%9%H7W;;$G;H$&$N$G(B
  $BDI2C$7$?!#(B
*/
real ProtonPotentialEnergy( sLattIce* const ice, int proton, real radius )
{
  //real ep = ProtonInPDIPotentialEnergy( ice, proton, radius );
#ifdef ILA
  //$B%H%]%m%8!<7g4Y$r4^$`7O$N>l9g$b(BPCI$B$K$O7g4Y$N$J$$3J;R$rM?$($kI,MW$,(B
  //$B$"$k!#(B
  //ep -= ProtonInPCIPotentialEnergy( ice, proton, radius );
#endif
  //for debug
  sprintf(stderr,"ProtonPotentialEnergy\n");
  real ep = ProtonInPDIPCIPotentialEnergy( ice, proton, radius );
  return ep;
}

        

/*$B%W%m%H%s(B($B>l=j$O(Bproton)$B$N(BPDI(Proton Disordered Ice)$B$X$NMOG^OB%(%M%k%.!<$r7W;;$9$k!#(B*/
/*radius$B$O!"Aj8_:nMQ7W;;$rBG$A@Z$kH>7B$r;XDj$9$k!#(B0$B$r;XDj$9$k$HBG$A@Z(B
  $B$i$J$$!#DL>o$N7W;;$G$OBG$A@Z$j$OI,MW$J$$$,!"<}B+%F%9%H7W;;$G;H$&$N$G(B
  $BDI2C$7$?!#(B
*/
real ProtonInPDIPotentialEnergy( sLattIce* ice, int proton, real radius )
{
  //LattIce.c$B$G$OJ,;RC10L$GAmOB$r7W;;$7$F$$$?$,!"(BLattIce2.c$B$G$O%@%_!<%\%s%I$r4^$`A47k9g$H$NAj8_:nMQ$r7W;;$9$k!#>pJsMn$A$,?4G[$@$,!"$R$H$^$:;@AG$OL5;k$9$k!#(B

    real px,py,pz,ch,esum;
    int position;
    int h;

    /*dryrun$B;~$O>o$K(B0*/
    //if ( ice->dryrun )
    //    return 0;

    position = ice->bonds[ proton ].direction + 1;
    px = ice->bonds[ proton ].px[position];
    py = ice->bonds[ proton ].py[position];
    pz = ice->bonds[ proton ].pz[position];
    ch = ice->bonds[proton].direction==0 ? QP : QH;
    
    esum = 0.0L;
    
    for( h=0; h<ice->nBond; h++ ){
      if ( h != proton ){
        int hposition = ice->bonds[ h ].direction + 1;
        real hpx = ice->bonds[ h ].px[hposition];
        real hpy = ice->bonds[ h ].py[hposition];
        real hpz = ice->bonds[ h ].pz[hposition];
	
        //Interaction with hydrogen
        real dx = px - hpx;
        real dy = py - hpy;
        real dz = pz - hpz;
	
        real r = sqrt( dx*dx + dy*dy + dz*dz );
	if( radius == 0 || r < radius ){
          real e = ch * QH / r;
          esum += e;
          //printf("PDI %f %f\n", r, e);
        }
      }
    }
    //fprintf( stderr, "%f PDI\n", count);
    return esum;
}



real ProtonInPDIPCIPotentialEnergy( sLattIce* const ice, int proton, real radius )
{
  //$B7k9g$4$H$K%-%c%s%;%k$5$;$?$[$&$,@:EY$O>e$,$k!#(B

    real px,py,pz,ch,esum;
    int position;
    int h;

    /*dryrun$B;~$O>o$K(B0*/
    //if ( ice->dryrun )
    //    return 0;

    position = ice->bonds[ proton ].direction + 1;
    px = ice->bonds[ proton ].px[position];
    py = ice->bonds[ proton ].py[position];
    pz = ice->bonds[ proton ].pz[position];
    ch = ice->bonds[proton].direction==0 ? QP : QH;
    
    esum = 0.0L;
    
    for( h=0; h<ice->nBond; h++ ){
      if ( h != proton ){
        real hpx = ice->bonds[ h ].px[1];
        real hpy = ice->bonds[ h ].py[1];
        real hpz = ice->bonds[ h ].pz[1];
	
        //Interaction with hydrogen
        real dx = px - hpx;
        real dy = py - hpy;
        real dz = pz - hpz;
	
        real r = sqrt( dx*dx + dy*dy + dz*dz );
	if( radius == 0 || r < radius ){
          real e = ch * QH / r;
          esum -= e;
          //printf("PCIPDI %f %f\n", r, -e);

          int hposition = ice->bonds[ h ].direction + 1;
          real hpx = ice->bonds[ h ].px[hposition];
          real hpy = ice->bonds[ h ].py[hposition];
          real hpz = ice->bonds[ h ].pz[hposition];
          
          //Interaction with hydrogen
          real dx = px - hpx;
          real dy = py - hpy;
          real dz = pz - hpz;
	
          real r = sqrt( dx*dx + dy*dy + dz*dz );
          e = ch * QH / r;
          esum += e;

          //printf("PDIPCI %f %f\n", r, e);
        }
      }
    }
    //fprintf( stderr, "%f PDI\n", count);
    return esum;
}



real ProtonInPCIPotentialEnergy( sLattIce* ice, int proton, real radius )
{
  //LattIce.c$B$G$OJ,;RC10L$GAmOB$r7W;;$7$F$$$?$,!"(BLattIce2.c$B$G$O%@%_!<%\%s%I$r4^$`A47k9g$H$NAj8_:nMQ$r7W;;$9$k!#>pJsMn$A$,?4G[$@$,!"$R$H$^$:;@AG$OL5;k$9$k!#(B

    real px,py,pz,ch,esum;
    int position;
    int h;

    /*dryrun$B;~$O>o$K(B0*/
    //if ( ice->dryrun )
    //    return 0;

    position = ice->bonds[ proton ].direction + 1;
    px = ice->bonds[ proton ].px[position];
    py = ice->bonds[ proton ].py[position];
    pz = ice->bonds[ proton ].pz[position];
    ch = ice->bonds[proton].direction==0 ? QP : QH;
    
    esum = 0.0L;
    
    for( h=0; h<ice->nBond; h++ ){
      if ( h != proton ){
        //PCI$B$G$O?eAG$O>o$KCf1{0LCV$K$"$k!#(B
        real hpx = ice->bonds[ h ].px[1];
        real hpy = ice->bonds[ h ].py[1];
        real hpz = ice->bonds[ h ].pz[1];
	
        //Interaction with hydrogen
        real dx = px - hpx;
        real dy = py - hpy;
        real dz = pz - hpz;
	
        real r = sqrt( dx*dx + dy*dy + dz*dz );
	if( radius == 0 || r < radius ){
          real e = ch * QH / r;
          esum += e;
          //printf("PCI %f %f\n", r, e);
        }
      }
    }
    //fprintf( stderr, "%f PDI\n", count);
    return esum;
}

                    

real Yaplot( sLattIce* ice )
{
    /*$B3NG'$N$?$a(Byaplot$BMQ$N=PNO(B*/

    real px,py,pz,ch,esum;
    int o1;
    int position;
    int proton;
    
    proton = ice->currentBond;

    position = ice->bonds[ proton ].direction;
    px = ice->bonds[ proton ].px[position+1];
    py = ice->bonds[ proton ].py[position+1];
    pz = ice->bonds[ proton ].pz[position+1];
    ch = ice->bonds[proton].direction==0 ? QP : QH;

    printf("@ 3\nr 0.2\no %f %f %f\n", px, py, pz );
    
    esum = 0.0L;
    
    for( o1=0; o1<ice->nOxygen; o1++ ){
        if ( ice->oxygens[o1].charge != 0.0 ){
            real emol=0.0L;
            real dx,dy,dz;
            int i;
            
            //Interaction with oxygen
            dx = px - ice->oxygens[o1].x;
            dy = py - ice->oxygens[o1].y;
            dz = pz - ice->oxygens[o1].z;
            printf("@ 3\nr 0.1\no %f %f %f\n", ice->oxygens[o1].x,
                   ice->oxygens[o1].y, ice->oxygens[o1].z );

            emol = ch * (-2)*QH / sqrt( dx*dx + dy*dy + dz*dz );

            printf("@ 4\n");
            //Interaction with hydrogen
            for( i=0; i<ice->oxygens[o1].nAdj; i++ ){
                if ( BondDirection( ice, o1, i ) == OUTGOING ){
                    int hb;
                    int pos;
                    hb = ice->oxygens[ o1 ].hydrogenBond[i];
                    pos = ice->bonds[ hb ].direction + 1;
                    
                    dx = px - ice->bonds[ hb ].px[pos];
                    dy = py - ice->bonds[ hb ].py[pos];
                    dz = pz - ice->bonds[ hb ].pz[pos];
                    
                    printf("l %f %f %f %f %f %f\n", ice->oxygens[o1].x, ice->oxygens[o1].y, ice->oxygens[o1].z, ice->bonds[ hb ].px[pos], ice->bonds[ hb ].py[pos], ice->bonds[ hb ].pz[pos] );
                }
            }
            esum += emol;
        }
    }
    return esum;
}

                    


void CheckNetwork( sLattIce* ice )
{
    int o1;
    
    for( o1=0; o1<ice->nOxygen; o1++ ){
        int x,y,z;
        
        x = ice->oxygens[o1].x;
        y = ice->oxygens[o1].y;
        z = ice->oxygens[o1].z;
        if ( ice->oxygens[o1].charge == 0 ){
            printf("%4d ( %4d,%4d,%4d ) dummy\n", o1, x,y,z);
        }
        else {
            int i;
            int n[3];
            for ( i=0; i<3; i++ )
                n[i] = 0;
            for ( i=0; i<ice->oxygens[o1].nAdj; i++ ) {
                n[ BondDirection( ice,o1, i ) + 1 ] ++;
            }
            printf("%4d ( %4d,%4d,%4d ) in%d out%d pro%d\n", o1, x,y,z, n[0],n[2],n[1]);
        }
    }
}

            
void RandomizeHBDirections( sLattIce* ice )
{
    int i;
    
    /*$B?eAG$N0LCV$r%i%s%@%`$K7h$a$k!#(B*/
    for ( i=0; i< ice->nBond; i++ ){
        ice->bonds[i].direction = ( lrand48() % 2 ) ? +1 : -1;
    }
    
    // count outgoing bonds of each oxygen atom
    for ( i=0; i< ice->nOxygen; i++ ){
        // if not dummy site
        if ( ice->oxygens[i].charge != 0 )
            ice->oxygens[i].outgo = OutBonds( ice, i );
        else
            ice->oxygens[i].outgo = 0;
    }
}



void Protonate( sLattIce* ice, bool dryrun )
{
    int i,o1,o2,bond;

    //$BCf1{IU6a$K(Bproton$B$rCV$/(B
    //$BCf?4$O:BI8$G(B(n,n,n)$BIU6a$+$J!#(B
    for( ;; ){
        int dx,dy,dz;
        //choose one
        o1 = lrand48() % ice->nOxygen;
        dx = ice->oxygens[o1].x - ice->comx;
        dy = ice->oxygens[o1].y - ice->comy;
        dz = ice->oxygens[o1].z - ice->comz;
        if ( ice->oxygens[o1].charge != 0 && dx*dx + dy*dy + dz*dz < 18 )
            break;
    }
    //choose a bond
    i = lrand48() % ice->oxygens[o1].nAdj;
    o2 = ice->oxygens[o1].oxygen[ i ];
    bond = ice->oxygens[o1].hydrogenBond[ i ];
    if ( OUTGOING == BondDirection( ice, o1, i ) ) {
        ice->oxygens[o1].outgo--;
    }
    else {
        ice->oxygens[o2].outgo--;
    }
    // set proton in the middle
    ice->bonds[bond].direction = PROTON;
    fprintf(stderr, "%d Protonated\n", bond);
    //$B8=:_%W%m%H%s$,$"$k>l=j$r5-21$9$k!#(B
    ice->currentBond = bond;

    //defect$B$r2r>C$9$k!#<~4|6-3&>r7o$G$O$J$$$N$G!"6-3&$^$G%W%m%H%s$r2!$7$@$9$H(Bdefect$B$O2r>C$9$k!#(B
    PurgeDefects( ice );

#ifdef DEBUG
    //Check
    //CheckNetwork( );
#endif

    //proton$B$H<~0O$NAj8_:nMQ$r;;=P$9$k!#(B
    if ( ! dryrun )
      ice->currentPotential = ProtonPotentialEnergy( ice, ice->currentBond, 0.0L );

#ifdef YAPLOT
    //Check
    Yaplot( ice );
#endif
    
}

//$B;@AG$rG[CV$7!"7k9g$r:n$j!"%W%m%H%s$rCV$$$F!"=i4|G[CV$r7A@.$9$k!#(B
//hposition$B$O!"(BO-O$B4V5wN%$N$I$N0LCV$K?eAG$rCV$/$+(B($BDL>o$N?eAG7k9g$N>l9g$K(B)
sLattIce* new_ProtonatedIceIc( int n, real hposition, bool dryrun )
{
    sLattIce* ice;
    
    //$B;@AG$N0LCV$H7k9g$r:n$k(B
    ice = CubicIce( n, hposition );
    RandomizeHBDirections( ice );
    Protonate( ice, dryrun );
    return ice;
}



//$B;@AG0LCV$r$b$H$K!"7k9g$rDj5A$9$k!#(B
void DetermineBonds( sLattIce* ice, real rhb )
{
    int nhb = 0;

    //$BC1=c$K(B2$B=E%k!<%W$K$7$F$7$^$&$H!"$b$N$9$4$/;~4V$,$+$+$k$N$G!"%0%j%C%IJ,3d$r9T$&!#(B
    if ( ice->periodic ){
        int o1;
        for ( o1=0; o1< ice->nOxygen; o1++ ){
            int o2;
            for( o2=o1+1; o2< ice->nOxygen; o2++ ){
                real dx,dy,dz;
                dx = ice->oxygens[o1].x - ice->oxygens[o2].x;
                dy = ice->oxygens[o1].y - ice->oxygens[o2].y;
                dz = ice->oxygens[o1].z - ice->oxygens[o2].z;
                if ( ice->periodic ){
                    dx -= ice->bx * rint( dx / ice->bx );
                    dy -= ice->by * rint( dy / ice->by );
                    dz -= ice->bz * rint( dz / ice->bz );
                }
                if ( dx*dx + dy*dy + dz*dz < rhb*rhb ){
                    int i;
                    
                    ice->bonds[nhb].oxygen[0] = o1;
                    ice->bonds[nhb].oxygen[1] = o2;
                    i = ice->oxygens[o1].nAdj;
                    if ( i == MAXNEIBOR ){
                        fprintf( stderr, "Error: number of neighborhood exceeds the limit at oxygen #%d.\n", o1 );
                        exit(1);
                    }
                    ice->oxygens[o1].oxygen[i] = o2;
                    ice->oxygens[o1].hydrogenBond[i] = nhb;
                    ice->oxygens[o1].nAdj ++;
                    
                    i = ice->oxygens[o2].nAdj;
                    if ( i == MAXNEIBOR ){
                        fprintf( stderr, "Error: number of neighborhood exceeds the limit at oxygen #%d.\n", o2 );
                        exit(1);
                    }
                    ice->oxygens[o2].oxygen[i] = o1;
                    ice->oxygens[o2].hydrogenBond[i] = nhb;
                    ice->oxygens[o2].nAdj ++;
                    
                    nhb ++;
                }
            }
        }
    }
    else {
        real* x;
        real* y;
        real* z;
        int i;
        sPairList* pairList;
        
        x = malloc( sizeof(real) * ice->nOxygen );
        y = malloc( sizeof(real) * ice->nOxygen );
        z = malloc( sizeof(real) * ice->nOxygen );
        
        for( i=0; i<ice->nOxygen; i++ ){
            x[i] = ice->oxygens[i].x;
            y[i] = ice->oxygens[i].y;
            z[i] = ice->oxygens[i].z;
        }
        pairList = PairList( ice->nOxygen, x,y,z, rhb );
        
        while ( pairList != NULL ){
            sPairList* next;
            int i;
            int o1,o2;
            
            o1 = pairList->i;
            o2 = pairList->j;
            if ( o2 < o1 ){
                o1 = pairList->j;
                o2 = pairList->i;
            }
            
            ice->bonds[nhb].oxygen[0] = o1;
            ice->bonds[nhb].oxygen[1] = o2;
            i = ice->oxygens[o1].nAdj;
            if ( i == MAXNEIBOR ){
                fprintf( stderr, "Error: number of neighborhood exceeds the limit at oxygen #%d.\n", o1 );
                exit(1);
            }
            ice->oxygens[o1].oxygen[i] = o2;
            ice->oxygens[o1].hydrogenBond[i] = nhb;
            ice->oxygens[o1].nAdj ++;
            
            i = ice->oxygens[o2].nAdj;
            if ( i == MAXNEIBOR ){
                fprintf( stderr, "Error: number of neighborhood exceeds the limit at oxygen #%d.\n", o2 );
                exit(1);
            }
            ice->oxygens[o2].oxygen[i] = o1;
            ice->oxygens[o2].hydrogenBond[i] = nhb;
            ice->oxygens[o2].nAdj ++;
                    
            nhb ++;

            next = pairList->next;
            free( pairList );
            pairList = next;
        }
        free( x );
        free( y );
        free( z );
    
    }
    ice->nBond = nhb;
}
    
//$B;@AG$N0LCV$rFI$_$3$`!#(B
sLattIce* load_LattIce( FILE* file, real hposition )
{
    char buf[1000];
    sLattIce* ice;
    bool coord=False, network=False;
    real rhb=0;

    ice = calloc( 1, sizeof( sLattIce ) );
    ice->periodic = False;
    //ice->dryrun   = False;
    /*$B%W%m%H%s$N>l=j$OITL@(B*/
    ice->currentBond = -1;
    ice->currentPotential = 0;
    while( NULL != fgets( buf, sizeof( buf ), file ) ){
        if ( strncmp( buf, "@BOX3",5 ) == 0 ){
            fgets( buf, sizeof( buf ), file );
            ice->periodic = True;
#ifdef SINGLEPRECISION
            sscanf( buf, "%f %f %f", &ice->bx, &ice->by, &ice->bz );
#else
            sscanf( buf, "%lf %lf %lf", &ice->bx, &ice->by, &ice->bz );
#endif
        }
        else if ( strncmp( buf, "@BXLA",5 ) == 0 ){
            fgets( buf, sizeof( buf ), file );
            ice->periodic = True;
#ifdef SINGLEPRECISION
            sscanf( buf, "%f", &ice->bx );
#else
            sscanf( buf, "%lf", &ice->bx );
#endif
            ice->by = ice->bz = ice->bx;
        }
        else if ( strncmp( buf, "@NGPH", 5 ) == 0 ){
            int o1,o2;
            int no;
            int nhb=0;

            network = True;
            fgets( buf, sizeof( buf ), file );
            no = atoi( buf );
            if ( ice->nOxygen == 0 ){
                ice->nOxygen = no;
                ice->comx = ice->comy = ice->comz = 0;
                ice->oxygens = calloc( ice->nOxygen,     sizeof( sOxygen ) );
                ice->bonds   = calloc( ice->nOxygen * 2, sizeof( sHydrogenBond ) );
            }
            else if ( ice->nOxygen != no ){
                fprintf( stderr, "Error: number of sites differ: %d != %d\n",
                         ice->nOxygen, no );
                exit(1);
            }
            /*$B%M%C%H%o!<%/$rFI$_$3$`!#0l1~M-8~%0%i%U$H9M$($k!#(B*/
            while ( NULL != fgets( buf, sizeof( buf ), file ) ){
                sOxygen*       oxy1;
                sOxygen*       oxy2;
                sHydrogenBond* hb;
                sscanf( buf, "%d %d", &o1, &o2 );
                if ( o1 < 0 ) break;
                oxy1 = &ice->oxygens[ o1 ];
                oxy2 = &ice->oxygens[ o2 ];
                hb   = &ice->bonds[ nhb ];
                oxy1->hydrogenBond[ oxy1->nAdj ] = nhb;
                oxy1->oxygen[ oxy1->nAdj ]       = o2;
                oxy2->hydrogenBond[ oxy2->nAdj ] = nhb;
                oxy2->oxygen[ oxy2->nAdj ]       = o1;
                if ( o1 < o2 ){
                    hb->oxygen[0] = o1;
                    hb->oxygen[1] = o2;
                    hb->direction = +1;
                }
                else {
                    hb->oxygen[0] = o2;
                    hb->oxygen[1] = o1;
                    hb->direction = -1;
                }
                oxy1->outgo ++;
                oxy1->nAdj ++;
                oxy2->nAdj ++;
                nhb ++;
            }
            ice->nBond = nhb;
        }            
        else if ( strncmp( buf, "@PPOS", 5 ) == 0 ){
            int o1,o2;
            int no;
            int nhb;
            sOxygen*       oxy1;
            sOxygen*       oxy2;
            sHydrogenBond* hb;

            /*@PPOS$B$h$jA0$K(B@NGPH$B$rFI$_$3$s$G$$$kI,MW$,$"$k!#(B*/
            if ( ! network ){
                fprintf( stderr, "Network topology must be read before proton position.\n" );
                exit(1);
            }
            fgets( buf, sizeof( buf ), file );
            no = atoi( buf );
            /*$B%M%C%H%o!<%/$rFI$_$3$`!#0l1~M-8~%0%i%U$H9M$($k!#(B*/
            nhb = ice->nBond;
            /*$B:G=i$N%W%m%H%s$7$+FI$^$J$$!#(B*/
            fgets( buf, sizeof( buf ), file );
            sscanf( buf, "%d %d", &o1, &o2 );
            oxy1 = &ice->oxygens[ o1 ];
            oxy2 = &ice->oxygens[ o2 ];
            hb   = &ice->bonds[ nhb ];
            oxy1->hydrogenBond[ oxy1->nAdj ] = nhb;
            oxy1->oxygen[ oxy1->nAdj ]       = o2;
            oxy2->hydrogenBond[ oxy2->nAdj ] = nhb;
            oxy2->oxygen[ oxy2->nAdj ]       = o1;
            hb->oxygen[0] = o1;
            hb->oxygen[1] = o2;
            hb->direction = 0;
            oxy1->nAdj ++;
            oxy2->nAdj ++;
            ice->currentBond = nhb;
            nhb ++;
            ice->nBond = nhb;
        }            
        else if ( strncmp( buf, "@RCOO", 5 ) == 0 ){
            /*$B7k9g$7$F$$$k$H$_$J$9(BO-O$B4V5wN%$NogCM(B*/
            fgets( buf, sizeof( buf ), file );
            rhb = atof( buf );
        }
        else if ( strncmp( buf, "@NX4A", 5 ) == 0 ||
                  strncmp( buf, "@AR3A", 5 ) == 0 ||
                  strncmp( buf, "@NX3A", 5 ) == 0 ){
            int o1;
            int no;

            coord = True;
            fgets( buf, sizeof( buf ), file );
            no = atoi( buf );
            if ( ice->nOxygen == 0 ){
                ice->nOxygen = no;
                ice->comx = ice->comy = ice->comz = 0;
                ice->oxygens = calloc( ice->nOxygen,     sizeof( sOxygen ) );
                ice->bonds   = calloc( ice->nOxygen * 2, sizeof( sHydrogenBond ) );
            }
            else if ( ice->nOxygen != no ){
                fprintf( stderr, "Error: number of sites differ: %d != %d\n",
                         ice->nOxygen, no );
                exit(1);
            }
            /*$B:BI8$rFI$_$3$`!#(B*/
            for( o1=0; o1< ice->nOxygen; o1++ ){
                fgets( buf, sizeof( buf ), file );
#ifdef SINGLEPRECISION
                sscanf( buf, "%f %f %f",
                        &ice->oxygens[o1].x, 
                        &ice->oxygens[o1].y, 
                        &ice->oxygens[o1].z );
#elif real == double
                sscanf( buf, "%lf %lf %lf",
                        &ice->oxygens[o1].x, 
                        &ice->oxygens[o1].y, 
                        &ice->oxygens[o1].z );
#else
#error
#endif
                /*$B%U%!%$%k$+$iFI$`>l9g$O$H$j$"$($:%@%_!<%5%$%H$G$O$J$$$H9M$($k!#"*J?@.(B16$BG/(B2$B7n(B18$BF|(B($B?e(B)$BJQ99!#7k9g?t$K$h$C$F%@%_!<$+$I$&$+$r7h$a$k!#(B*/
                ice->oxygens[o1].charge = -2;
            }
#ifdef DEBUG
            fprintf( stderr, "CoM: %f %f %f\n", ice->comx, ice->comy, ice->comz );
#endif
        }
        //if ( coord && network ) break;
    }

    if ( ! network && 0.0 < rhb ) {
        DetermineBonds( ice, rhb );
    }

    CenterOfMass( ice );
    /*O-O$B4V$N%W%m%H%s$N$H$j$&$k0LCV(B3$B$D$N:BI8$r7W;;$9$k!#(B*/
    SetHBins( ice, hposition );
    return ice;
    //fprintf( stderr, "Error: No inputs.\n" );
    //exit(1);
}


/*$B!V8&Ka!W3J;R$NI=LL$K0LCV$7!"NY@\%5%$%H?t$,ITB-$7$F$$$k%5%$%H$r!"%@%_!<(B
  $B%5%$%H$H$7$F$7$^$&!#(B

  $B%@%_!<%5%$%H$O!"?eAG7k9g$N%?!<%2%C%H$K$O$J$l$k$,!"$=$l<+?H$O;@AG$r$b(B
  $B$?$J$$!#?eAG$O;}$D(B(lLattIce.c$B$H$N0c$$(B)$B!#I=LL$N%@%s%0%j%s%07k9g$N8~$-(B
  $B$r$F$-$H$&$K<}$a$k$?$a$KF3F~$9$k!#(B

*/
void Grind( sLattIce* ice, int grind, real hposition )
{
    int o1;
    for( o1=0; o1<ice->nOxygen; o1++ ){
        if ( ice->oxygens[o1].nAdj < grind ){
            int i;
            //dummy site$B2=!#(Bdummy site$B$rDL>o$N;@AG$H6hJL$9$kJ}K!$OEE(B
            //$B2Y$,(B0$B$G$"$k$+H]$+!#(B
            ice->oxygens[o1].charge = 0;
            //dummy site$B$O!"(Bbond$B$r9=@.$9$k(B2$B%5%$%H$N$&$A$N(B[0]$BHVL\$N%5%$%H$K$O$J$l$J$$$N$G!"(B
            //$BA4ItF~$l$+$($k!#(B
            for( i=0; i<ice->oxygens[o1].nAdj; i++ ){
                int hb  = ice->oxygens[o1].hydrogenBond[i];
                int tmp = ice->bonds[hb].oxygen[0];
                if ( tmp == o1 ){
                    ice->bonds[hb].oxygen[0] = ice->bonds[hb].oxygen[1];
                    ice->bonds[hb].oxygen[1] = tmp;
                    ice->bonds[hb].direction = - ice->bonds[hb].direction;
                    SetHBin( ice, hb, hposition );
                }
            }
        }
    }
}


//$B9=B$$rFI$_$3$_!"%W%m%H%s$rCV$$$F!"=i4|G[CV$r7A@.$9$k!#(B
sLattIce* LoadLatticeAndProtonate( FILE* file, bool randomize, int grind, real hposition, bool dryrun )
{
    sLattIce* ice;
    
    //$B;@AG$N0LCV$rFI$_$3$`!#(B
    ice = load_LattIce( file, hposition );
    if ( ice->periodic ){
        fprintf( stderr, "Error: protonated ice should not be in periodic boundary condition.\n" );
        exit(1);
    }
    //$BNY@\%5%$%H?t$,B-$j$J$$%5%$%H$r%@%_!<%5%$%H$K$7$F$7$^$&!#(B
    Grind( ice, grind, hposition );
    /*$B%W%m%H%s$N>l=j$,FI$_$3$^$l$J$+$C$?>l9g$O!"(Bnetwork$BA4BN$r%i%s%@%^%$%:$7!"=E?4IU6a$K%W%m%H%s$N>l=j$r%i%s%@%`$KA*$V!#(B*/
    if ( ice->currentBond < 0 ){
        if ( randomize )
            RandomizeHBDirections( ice );
        Protonate( ice, dryrun );
    }
    else {
      if ( ! dryrun )
	ice->currentPotential = ProtonPotentialEnergy( ice, ice->currentBond, 0.0L );
    }
    return ice;
}




//$B?eAG7k9g$N8~$-$r%i%s%@%`$K$9$k!#(B
void DisorderProton( sLattIce* ice )
{
    RandomizeHBDirections( ice );
    
    //defect$B$r2r>C$9$k!#(B
    PurgeDefects( ice );

#ifdef DEBUG
    //Check
    //CheckNetwork( );
#endif

    //$B8=:_%W%m%H%s$,$"$k>l=j$r5-21$9$k!#(B
    ice->currentBond = -1;

    //proton$B$H<~0O$NAj8_:nMQ$r;;=P$9$k!#(B
    ice->currentPotential = 0;

#ifdef YAPLOT
    //Check
    Yaplot( ice );
#endif
}



real PotentialEnergy( sLattIce* ice )
{
    return ice->currentPotential;
}



void SnapShot( sLattIce* ice, int aux )
{
    //
    //proton$B$N0LCV$r=PNO$9$k!#(B
    //
    real x,y,z;
    x = ice->bonds[ ice->currentBond ].px[1] - ice->comx;
    y = ice->bonds[ ice->currentBond ].py[1] - ice->comy;
    z = ice->bonds[ ice->currentBond ].pz[1] - ice->comz;
    
    printf("t %f %f %f %d(%f)\n", x,y,z,aux,ice->currentPotential);
}



//$B%G%P%C%0MQ!#A4%(%M%k%.!<$r7W;;$9$k!#(B
real TotalEnergy( sLattIce* ice )
{

    real px,py,pz,ch,esum;
    int h1;
    int position;

    esum = 0.0L;

    for( h1=0; h1<ice->nBond; h1++ ){
        int h2,o1;
        
        position = ice->bonds[ h1 ].direction + 1;
        px = ice->bonds[ h1 ].px[position];
        py = ice->bonds[ h1 ].py[position];
        pz = ice->bonds[ h1 ].pz[position];
        ch = ice->bonds[ h1 ].direction==0 ? QP : QH;
    
        //H-O
        for( o1=0; o1<ice->nOxygen; o1++ ){
            if ( ice->oxygens[o1].charge != 0 ){
#ifdef ILA
#else
                real dx,dy,dz;
                //Interaction with oxygen
                dx = px - ice->oxygens[o1].x;
                dy = py - ice->oxygens[o1].y;
                dz = pz - ice->oxygens[o1].z;
                
                esum += ch * (-2)*QH / sqrt( dx*dx + dy*dy + dz*dz );
#endif
            }
        }
        //H-H
        for( h2=0; h2<h1; h2++ ){
            int pos;
            real ch2;
            real dx,dy,dz;
                
            pos = ice->bonds[ h2 ].direction + 1;
            ch2 = (ice->bonds[ h2 ].direction == 0) ? QP : QH;
            //Interaction with oxygen
            dx = px - ice->bonds[h2].px[pos];
            dy = py - ice->bonds[h2].py[pos];
            dz = pz - ice->bonds[h2].pz[pos];
            
            esum += ch * ch2 / sqrt( dx*dx + dy*dy + dz*dz );
#ifdef ILA
            dx = px - ice->bonds[h2].px[1];
            dy = py - ice->bonds[h2].py[1];
            dz = pz - ice->bonds[h2].pz[1];
            
            esum -= ch * ch2 / sqrt( dx*dx + dy*dy + dz*dz );
#endif
        }
    }
    return esum;
}



void CheckConsistency( sLattIce* ice, real radius )
{
    //debug$BMQ!#%W%m%H%s$,0l$D$@$1$+$I$&$+!"(Boutgo$B$N?t$,@5$7$$$+!"$J$I$rD4$Y$k!#(B
    int o1;
    int nproton=0;
    int ndummy=0;
    int statNAdj[10];
    int i;
    int volume=0;
    int numWater = 0; // number of H-O-H
    int numOH    = 0; // number of O-H
    int numH3O   = 0; // number of H3O
    int numOthers= 0; // other molecule (what?)
    
    for( i=0; i<10; i++)
        statNAdj[i] = 0;
    
    for( o1=0; o1<ice->nOxygen; o1++ ){
#ifdef DEBUG
        //fprintf( stderr, "O %d nAdj %d charge %d\n", o1, ice->oxygens[o1].nAdj, ice->oxygens[o1].charge );
#endif
        if ( ice->oxygens[o1].charge == 0 ){
            ndummy ++;
        }
        else{
            int j;
            int noutgo=0;
            statNAdj[ ice->oxygens[o1].nAdj ]++;
            for( j=0; j< ice->oxygens[o1].nAdj; j++ ){
                int o2  = ice->oxygens[o1].oxygen[j];
                int dir = BondDirection( ice, o1, j );
                if ( dir == OUTGOING )
                    noutgo ++;
                if ( ice->oxygens[o2].charge != 0 )
                    if ( dir == PROTON ){
                        //fprintf( stderr, "Proton@%d\n", ice->oxygens[o1].hydrogenBond[ j ] );
                        nproton ++;
                    }
            }
	    if ( noutgo == 2 ){
	      numWater++;
	    }
	    else if ( noutgo == 3 ){
	      numH3O++;
	    }
	    else if ( noutgo == 1 ){
	      numOH++;
	    }
	    else{
	      numOthers++;
	    }
            if ( noutgo != ice->oxygens[o1].outgo )
                printf( "Error: Outgo %d should be %d at %d\n",ice->oxygens[o1].outgo, noutgo, o1 );
            real dx,dy,dz;
            dx = ice->oxygens[o1].x - ice->comx;
            dy = ice->oxygens[o1].y - ice->comy;
            dz = ice->oxygens[o1].z - ice->comz;
            if ( dx*dx + dy*dy + dz*dz < radius*radius )
                volume ++;
        }
    }
    //two water molecules are in Zundel form around the proton.
    if ( nproton != 2 )
        printf( "Error: Number of proton %d should be 1\n", nproton/2 );
    printf( "Number of nodes in radius %f: %d\n", radius, volume );
    printf( "Number of dummy sites: %d\n", ndummy );
    printf( "Number of water molecules: %d\n", numWater );
    printf( "Number of H3O molecules: %d\n", numH3O );
    printf( "Number of OH molecules: %d\n", numOH );
    printf( "Number of other molecules: %d\n", numOthers );
    printf( "Histogram of the neighboring sites\n" );
    for( i=0; i<10; i++ ){
        printf("%d : %d\n", i, statNAdj[i] );
    }
}



//$B%/%i%9%?=E?4$+$i%W%m%H%s$^$G$N5wN%!#(B
real Radius( sLattIce* ice, int proton )
{
    real px,py,pz;
    int position;

    position = ice->bonds[ proton ].direction;
    px = ice->bonds[ proton ].px[position+1] - ice->comx;
    py = ice->bonds[ proton ].py[position+1] - ice->comy;
    pz = ice->bonds[ proton ].pz[position+1] - ice->comz;

    return sqrt( px*px + py*py + pz*pz );
}



/*$BL58B3J;R6a;w$N<}B+@:EY$r%A%'%C%/$9$k!#(B*/
void CheckConvergence( sLattIce* ice )
{
  real epdi, epci;
  real r=20;
  int proton = ice->currentBond;

  printf("# convergence\n");
  //while ( r < 180.0 ){
  while ( r < 20.05 ){
    //fprintf( stderr, "%f r\n", r);
    epdi = ProtonInPDIPotentialEnergy( ice, proton, r );
    epci = ProtonInPCIPotentialEnergy( ice, proton, r );
    //$B>pJsMn$A$rKI$07W;;=g=x$r;HMQ$7$?>l9g(B
    real epi = ProtonInPDIPCIPotentialEnergy( ice, proton, r );
    printf( "%f %24.17f %24.17f %24.17f %24.17f\n", r, epdi-epci, epdi, epci, epi );
    r += 0.1;
  }
}


//$B8=:_$N%W%m%H%s$N0LCV$+$i!"(Bnext$B$KF0$1$k$+$I$&$+$r%A%'%C%/$7!"%W%m%H%s$rF0$+$9!#(B
//$B%(%M%k%.!<:9$rJV$9(B
real Move( sLattIce* ice, int next, bool dryrun  )
{
    int newpos;

    for( newpos=0; newpos<2; newpos++ ){
        int o1;
        
        o1 = ice->bonds[ ice->currentBond ].oxygen[ newpos ];
        if ( ice->oxygens[ o1 ].charge != 0 ){
            int j;
            for( j=0; j<4; j++ ){
                if ( BondDirection( ice, o1, j ) == OUTGOING ) {
                    int destBond = ice->oxygens[ o1 ].hydrogenBond[ j ];
                    if ( destBond == next ) {
                        real oldEnergy, newEnergy;
                        int newDirection;
                        
                        /*$B%W%m%H%s0\F0$9$k@h$,!"8=:_$N7k9g$N$I$A$iB&$+(B*/
                        newDirection = ( newpos == 0 ) ? FORWARD : BACKWARD;
                        
                        // calculate energy for the destination hydrogen
                        if ( ! dryrun )
                          oldEnergy = ProtonPotentialEnergy( ice, ice->currentBond, 0.0L ) + ProtonPotentialEnergy( ice, destBond, 0.0L );
                        
                        //move protons
                        ice->bonds[ ice->currentBond ].direction = newDirection;
                        ice->bonds[ destBond ].direction    = 0;
                        
                        // calculate energy for new configuration
                        if ( ! dryrun )
                          newEnergy = ProtonPotentialEnergy( ice, ice->currentBond, 0.0L );
                        ice->currentBond = destBond;
                        if ( ! dryrun ){
                          ice->currentPotential = ProtonPotentialEnergy( ice, ice->currentBond, 0.0L );
                          newEnergy += ice->currentPotential;
                        
                          return ( newEnergy - oldEnergy );
                        }
                        else{
                          return 0;
                        }
                    }
                }
            }
        }
    }
    die( "Illegal path: %d", next );
}



//$B8=:_$N%W%m%H%s$N0LCV$H!"<!$N%W%m%H%s$N0LCV$N4V$G6&M-$7$F$$$k;@AG$NHV9f$rJV$9!#(B
real SharedOxygen( sLattIce* ice, int next )
{
    int newpos;

    for( newpos=0; newpos<2; newpos++ ){
        int o1;
        
        o1 = ice->bonds[ ice->currentBond ].oxygen[ newpos ];
        if ( ice->oxygens[ o1 ].charge != 0 ){
            int j;
            for( j=0; j<4; j++ ){
                if ( BondDirection( ice, o1, j ) == OUTGOING ) {
                    int destBond = ice->oxygens[ o1 ].hydrogenBond[ j ];
                    if ( destBond == next ) {
                        return o1;
                    }
                }
            }
        }
    }
    die( "Illegal path: %d", next );
}




void saveNGPH( sLattIce* ice, real radius, FILE* file )
{
    int hb;
    
    fprintf( file, "@NGPH\n%d\n", ice->nOxygen );
    for( hb=0; hb<ice->nBond; hb++ ){
        int o1,o2;
        o1 = ice->bonds[hb].oxygen[0];
        o2 = ice->bonds[hb].oxygen[1];
	if ( 0.0 < radius ){
	  if ( radius < Radius( ice, hb ) ){
	    continue;
	  }
	}
	//third value is dummy.
        if ( ice->bonds[hb].direction == +1 )
	  fprintf( file, "%d %d %d\n", o1,o2,hb );
        else if ( ice->bonds[hb].direction == -1 )
	  fprintf( file, "%d %d %d\n", o2,o1,hb );
    }
    fprintf( file, "-1 -1\n" );
}



void saveProtonPosition( sLattIce* ice, FILE* file )
{
    int hb;
    
    /*$B%W%m%H%s$N>l=j$O(BPPOS$B7A<0$GJ]B8!#(BPPOS$B$N=q<0$O(BNGPH$B$HF10l!#(B($B<B:]$K$O%W%m%H%s$,(B2$B8D0J>e$"$k$H8mF0:n$9$k!#(B)*/
    fprintf( file, "@PPOS\n%d\n", ice->nOxygen );
    for( hb=0; hb<ice->nBond; hb++ ){
        int o1,o2;
        
        o1 = ice->bonds[hb].oxygen[0];
        o2 = ice->bonds[hb].oxygen[1];
        if ( ice->bonds[hb].direction == 0 )
	  //third value is dummy
	  fprintf( file, "%d %d %d\n", o1,o2,hb );
    }
    fprintf( file, "-1 -1\n" );
}



void saveSite( sLattIce* ice, FILE* file )
{
    int o1;
    /*$B<~4|6-3&>r7o$G$"$l$PH"$NBg$-$5$r=PNO$9$k!#(B*/
    if ( ice->periodic )
        fprintf( file, "@BOX3\n%f %f %f\n", ice->bx, ice->by, ice->bz );
    /*$B;@AG$N0LCV$r(BAR3A$B7A<0$G=PNO$9$k!#(B*/
    fprintf( file, "@AR3A\n%d\n", ice->nOxygen );
    for( o1=0; o1<ice->nOxygen; o1++ ){
        real dx,dy,dz;
        dx = ice->oxygens[o1].x;
        dy = ice->oxygens[o1].y;
        dz = ice->oxygens[o1].z;
        fprintf( file, "%f %f %f\n", dx,dy,dz );
    }
}
