/*
  2008-2-2 LattIce.c --> LattIce2.c
  cluster外周の水素結合を、LattIce.cでは、まじめに分子単位で扱っている。
  つまり、外向きの結合は格子から突き出ているが、逆向きの結合は刈りとっ
  てある。これをPCIにフィットさせようとすると、PCIもこのクラスターの外
  周の形にあわせて刈り取らざるをえず、非常にとりあつかいがややこしい。
  このような扱いにした理由は、プロトンがクラスター表面まで来たときに不
  整合がおこらないようにするためだったのだが、実際には現在のプログラム
  ではプロトンは表面に来ないよう抑止されているので、問題にならない。

  LattIce2.cでは、ダミーの内向き水素結合を追加することで、PDI-PCI境界
  のでこぼこをなくし、取り扱いを簡単にする。
 */
/*
 * プロトン移動のモンテカルロ。クラスタを作り、中心にプロトンを置く。
 エネルギーを計算し、プロトン移動を次々に行う。格子の大きさはできるだ
 け大きくできるようにする。

 MonteCarloの基本的なステップは

 1. 現在の配置から、新しい配置を生成
 2. 新しい配置のエネルギーを計算
 3. accept/rejectの評価(energyを基準に)
 4a. acceptなら古い配置を棄却して新しい配置を採用
 4b. rejectなら新しい配置を棄却。

 プロトンの位置が変わると、酸素の位置も微妙に変わるが、今回は模型計算
 ということで、酸素は常に格子点に、水素は常にO-O線上にあるものとする。

 水素の位置は二値で表せるから、分子数がNだと、水素結合の本数は2N。2Nビッ
 トのメモリでとりあえず動くはず。

 一回のプロトン移動で2個のプロトンの位置が変わる。それによって2 x 3N対
 の相互作用の再計算が必要。ただし、酸素-水素間エネルギーはあらかじめ計
 算しておける。H-H間エネルギーはあらかじめ計算してためておくには数が多
 すぎるので、毎回再計算するものとすれば、相互作用の計算は全部で2 x 2 x 
 2N=8N対となる。

 酸素の位置に通し番号を与える。水素結合をなす2つの酸素の間の水素の位置
 は、通し番号が小さい酸素に近い位置を-1、大きい酸素側を+1とし、中央(プ
 ロトン)を0とする。

 水素の位置は、2つの酸素の番号と、位置値の組みあわせで座標が決
 まる。水素にも通し番号を与える。

 ある酸素に隣接する酸素4つを一覧にしておく。またその間の水素の番号も一
 覧にしておく。

 クラスタ表面に位置する水分子は、結合を4つもたないので、dummyサイトと
 して扱い、水素結合の定義にだけ使用する。

 相互作用エネルギーの単位は、格子上の距離1離れた2つの(水分子の)水素の
 間のエネルギーを1とする。酸素の電荷は水素の2倍とする。プロトンだけは
 異なる電荷(QP)を持つものとする。

 水素の位置は、それぞれ同一分子内、プロトン、隣接分子内にあるとき、酸
 素酸素間の直線上の、0.3, 0.5, 0.7の位置におく。実際の氷は、O-O間距
 離2.75A、O-H距離1.00A。比率としては0.36ぐらいになる。高圧の氷など
 をシミュレーションする可能性もあるので、この数字は指定できるほうがいい。

 プロトンが付着している水分子でも通常の水分子でも、プロトンを除く水素
 の数は2で変わりないから、一分子あたりの出結合数(outgo)はプロトン移動
 つれて更新する必要はない。ただ、長くシミュレーションしているうちに、
 プロトンがクラスタの表面にきてしまうとまずい。ILA(無限格子近似)で
 対処する。
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LattIce.h"
#include "PairList.h"
#include "Utility.h"

/*infinite lattice approximation*/
//ILAを使うかどうかはMakefileで指定する。
//#define ILA

/**/

/*ice Icでの座標とIDの関係*/

/*cubic iceの構造は、面心立方格子を2つずらせてかさねた構造になっている。
  面心立方格子もまた、単純格子を3方向にずらせて重ねた構造になっている。
  10x10x10=1000分子のcubic ice構造を作る場合、まず5x5x5の単純格子をつ
  くり、格子点に0..124番を割りふる。次に、yz平面上でずらした単純格子に
  は、もとの単純格子に125を足した番号を、xz平面上でずらした単純格子に
  は250を足した番号を、yz平面でずらした単純格子には375を加算する。さら
  に、面心立方格子全体を(1,1,1)方向にずらした副格子には500を加算する。

  ややこしいが、格子番号の振りかたはこのroutineにのみ依存している。
  各格子点に一意に番号がつくなら、他の方法でも構わない。

*/
int Coordinate2OxygenID( int n, int x, int y, int z )
{


    int xh,yh,zh, xq,yq,zq;
    int id;
    
    /*主格子(単純格子)の番号*/
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


/*o1番目の酸素のi番目のボンドが出結合なら1、入結合なら-1, 自由陽子なら0*/
int BondDirection( sLattIce* const ice, int o1, int i )
{
    int bond = ice->oxygens[o1].hydrogenBond[ i ];
    int dir = ( ice->bonds[ bond ].oxygen[0] == o1 ) ? +1 : -1;
    dir *= ice->bonds[ bond ].direction;
    return dir;
}



/*o1番目の酸素から出ている結合の本数を数える。2なら通常の水分子、3ならプロトンが付着した水分子。*/
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
        //charge==0のサイトはダミーサイト(クラスター表面を修飾するために付加したサイト)
        if ( 2 == ice->oxygens[o1].outgo || 0 == ice->oxygens[o1].charge ){
            // 欠陥でないサイトは見付けしだいリストから消す。
            ndefect --;
            defect[i] = defect[ndefect];
        }
        //o1がdummy siteではなく、出結合が多い場合
        else if ( 2 < ice->oxygens[o1].outgo ){
            int bond;
            int o2;
            //ランダムにo1からの出結合をさがす
            while ( 1 ){
                bond = lrand48() % ice->oxygens[o1].nAdj;
                if ( BondDirection( ice, o1, bond ) == OUTGOING )
                    break;
            }
            o2 = ice->oxygens[o1].oxygen[bond];
            ice->oxygens[o1].outgo --;
            //結合の反転
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
            //こちらに移した。dummy siteであっても、一応出結合数は数えておく。
            //ただしdummy siteはdefectにはならない。
            ice->oxygens[o2].outgo ++;
        }
        //o1がdummy siteではなく、出結合が足りない場合
        else {
            int bond;
            int o2;
            //ランダムにo1からの入結合をさがす
            while ( 1 ){
                bond = lrand48() % ice->oxygens[o1].nAdj;
                if ( BondDirection( ice, o1, bond ) == INCOMING )
                    break;
            }
            o2 = ice->oxygens[o1].oxygen[bond];
            ice->oxygens[o1].outgo ++;
            //結合の反転
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



/*水素の居場所H-Binを準備する*/
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



//-1, 0, +1の場所にある時の座標をあらかじめ計算しておく。
void SetHBins( sLattIce* ice, real hposition )
{
    int h;
    for( h=0; h<ice->nBond; h++ ){
      SetHBin( ice, h, hposition );
    }
}

    

//酸素の重心を求める。
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
                    //位置から酸素IDを得る
                    int self = Coordinate2OxygenID( n, x, y, z );
                    
                    if ( ice->oxygens[self].x != 0 || ice->oxygens[self].y != 0
                         || ice->oxygens[self].z != 0 ){
                        printf("%d %d %d %d is used\n", self,x,y,z);
                    }
                    //位置を保存する
                    ice->oxygens[self].x = x;
                    ice->oxygens[self].y = y;
                    ice->oxygens[self].z = z;

                    //立方体表面にあるなら
                    if ( x == 0 || x == 2*n-2 ||
                         y == 0 || y == 2*n-2 ||
                         z == 0 || z == 2*n-2 )
                      //ダミーサイトである。
		      ice->oxygens[self].charge = 0;
                    else{
                      //水素の電荷の-2倍。
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
    /*つぎに水素結合を定義する。*/
    nhb = 0;
    /*すべての酸素原子について*/
    for ( o1=0; o1<n*n*n; o1++ ){
        /*もしダミーサイト(表面サイト)でなければ*/
        if ( ice->oxygens[o1].charge != 0 ){
            int j;
            /*4つの隣接酸素について*/
            for( j=0; j<ice->oxygens[o1].nAdj; j++ ){
                int o2;
                
                o2 = ice->oxygens[o1].oxygen[j];
                /*o1<o2なら、ボンドを登録する(重複登録をさけるため)*/
                if ( o1 < o2 ) {
                    /*bondの2つの端点を登録する*/
                    ice->bonds[nhb].oxygen[0] = o1;
                    ice->bonds[nhb].oxygen[1] = o2;
                    /*この結合をo2への結合として登録する。*/
                    ice->oxygens[o1].hydrogenBond[j] = nhb;
                    /*隣接酸素がダミーサイトでないなら*/
                    if ( ice->oxygens[o2].charge != 0 ){
                        int k;
                    
                        /*隣接酸素の隣接酸素について*/
                        for ( k=0; k<ice->oxygens[o2].nAdj; k++ )
                            /*それがo1なら*/
                            if ( ice->oxygens[o2].oxygen[k] == o1 )
                                /*この結合をo2からo1への結合として登録する*/
                                ice->oxygens[o2].hydrogenBond[k] = nhb;
                    }
                    
                    /*登録結合数を増やす*/
                    nhb ++;
                }
                /*隣接酸素がダミーサイトなら*/
                else if ( ice->oxygens[o2].charge == 0 ) {
                    /*o2 < o1でもボンドを登録する。*/
                    ice->bonds[nhb].oxygen[1] = o1;
                    ice->bonds[nhb].oxygen[0] = o2;
                    ice->oxygens[o1].hydrogenBond[j] = nhb;
                    nhb ++;
                }
            }
        }
    }

    //結合総数には、ダミーボンドへの結合が含まれるので、酸素(電荷あり)の総数の2倍よりも大きいはず。
    ice->nBond = nhb;
    SetHBins( ice, hposition );
    
    return ice;
    
}




/*プロトン(場所はproton)の溶媒和エネルギーを計算する。*/
/*radiusは、相互作用計算を打ち切る半径を指定する。0を指定すると打ち切
  らない。通常の計算では打ち切りは必要ないが、収束テスト計算で使うので
  追加した。
*/
real ProtonPotentialEnergy( sLattIce* const ice, int proton, real radius )
{
  //real ep = ProtonInPDIPotentialEnergy( ice, proton, radius );
#ifdef ILA
  //トポロジー欠陥を含む系の場合もPCIには欠陥のない格子を与える必要が
  //ある。
  //ep -= ProtonInPCIPotentialEnergy( ice, proton, radius );
#endif
  //for debug
  sprintf(stderr,"ProtonPotentialEnergy\n");
  real ep = ProtonInPDIPCIPotentialEnergy( ice, proton, radius );
  return ep;
}

        

/*プロトン(場所はproton)のPDI(Proton Disordered Ice)への溶媒和エネルギーを計算する。*/
/*radiusは、相互作用計算を打ち切る半径を指定する。0を指定すると打ち切
  らない。通常の計算では打ち切りは必要ないが、収束テスト計算で使うので
  追加した。
*/
real ProtonInPDIPotentialEnergy( sLattIce* ice, int proton, real radius )
{
  //LattIce.cでは分子単位で総和を計算していたが、LattIce2.cではダミーボンドを含む全結合との相互作用を計算する。情報落ちが心配だが、ひとまず酸素は無視する。

    real px,py,pz,ch,esum;
    int position;
    int h;

    /*dryrun時は常に0*/
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
  //結合ごとにキャンセルさせたほうが精度は上がる。

    real px,py,pz,ch,esum;
    int position;
    int h;

    /*dryrun時は常に0*/
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
  //LattIce.cでは分子単位で総和を計算していたが、LattIce2.cではダミーボンドを含む全結合との相互作用を計算する。情報落ちが心配だが、ひとまず酸素は無視する。

    real px,py,pz,ch,esum;
    int position;
    int h;

    /*dryrun時は常に0*/
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
        //PCIでは水素は常に中央位置にある。
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
    /*確認のためyaplot用の出力*/

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
    
    /*水素の位置をランダムに決める。*/
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

    //中央付近にprotonを置く
    //中心は座標で(n,n,n)付近かな。
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
    //現在プロトンがある場所を記憶する。
    ice->currentBond = bond;

    //defectを解消する。周期境界条件ではないので、境界までプロトンを押しだすとdefectは解消する。
    PurgeDefects( ice );

#ifdef DEBUG
    //Check
    //CheckNetwork( );
#endif

    //protonと周囲の相互作用を算出する。
    if ( ! dryrun )
      ice->currentPotential = ProtonPotentialEnergy( ice, ice->currentBond, 0.0L );

#ifdef YAPLOT
    //Check
    Yaplot( ice );
#endif
    
}

//酸素を配置し、結合を作り、プロトンを置いて、初期配置を形成する。
//hpositionは、O-O間距離のどの位置に水素を置くか(通常の水素結合の場合に)
sLattIce* new_ProtonatedIceIc( int n, real hposition, bool dryrun )
{
    sLattIce* ice;
    
    //酸素の位置と結合を作る
    ice = CubicIce( n, hposition );
    RandomizeHBDirections( ice );
    Protonate( ice, dryrun );
    return ice;
}



//酸素位置をもとに、結合を定義する。
void DetermineBonds( sLattIce* ice, real rhb )
{
    int nhb = 0;

    //単純に2重ループにしてしまうと、ものすごく時間がかかるので、グリッド分割を行う。
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
    
//酸素の位置を読みこむ。
sLattIce* load_LattIce( FILE* file, real hposition )
{
    char buf[1000];
    sLattIce* ice;
    bool coord=False, network=False;
    real rhb=0;

    ice = calloc( 1, sizeof( sLattIce ) );
    ice->periodic = False;
    //ice->dryrun   = False;
    /*プロトンの場所は不明*/
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
            /*ネットワークを読みこむ。一応有向グラフと考える。*/
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

            /*@PPOSより前に@NGPHを読みこんでいる必要がある。*/
            if ( ! network ){
                fprintf( stderr, "Network topology must be read before proton position.\n" );
                exit(1);
            }
            fgets( buf, sizeof( buf ), file );
            no = atoi( buf );
            /*ネットワークを読みこむ。一応有向グラフと考える。*/
            nhb = ice->nBond;
            /*最初のプロトンしか読まない。*/
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
            /*結合しているとみなすO-O間距離の閾値*/
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
            /*座標を読みこむ。*/
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
                /*ファイルから読む場合はとりあえずダミーサイトではないと考える。→平成16年2月18日(水)変更。結合数によってダミーかどうかを決める。*/
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
    /*O-O間のプロトンのとりうる位置3つの座標を計算する。*/
    SetHBins( ice, hposition );
    return ice;
    //fprintf( stderr, "Error: No inputs.\n" );
    //exit(1);
}


/*「研磨」格子の表面に位置し、隣接サイト数が不足しているサイトを、ダミー
  サイトとしてしまう。

  ダミーサイトは、水素結合のターゲットにはなれるが、それ自身は酸素をも
  たない。水素は持つ(lLattIce.cとの違い)。表面のダングリング結合の向き
  をてきとうに収めるために導入する。

*/
void Grind( sLattIce* ice, int grind, real hposition )
{
    int o1;
    for( o1=0; o1<ice->nOxygen; o1++ ){
        if ( ice->oxygens[o1].nAdj < grind ){
            int i;
            //dummy site化。dummy siteを通常の酸素と区別する方法は電
            //荷が0であるか否か。
            ice->oxygens[o1].charge = 0;
            //dummy siteは、bondを構成する2サイトのうちの[0]番目のサイトにはなれないので、
            //全部入れかえる。
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


//構造を読みこみ、プロトンを置いて、初期配置を形成する。
sLattIce* LoadLatticeAndProtonate( FILE* file, bool randomize, int grind, real hposition, bool dryrun )
{
    sLattIce* ice;
    
    //酸素の位置を読みこむ。
    ice = load_LattIce( file, hposition );
    if ( ice->periodic ){
        fprintf( stderr, "Error: protonated ice should not be in periodic boundary condition.\n" );
        exit(1);
    }
    //隣接サイト数が足りないサイトをダミーサイトにしてしまう。
    Grind( ice, grind, hposition );
    /*プロトンの場所が読みこまれなかった場合は、network全体をランダマイズし、重心付近にプロトンの場所をランダムに選ぶ。*/
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




//水素結合の向きをランダムにする。
void DisorderProton( sLattIce* ice )
{
    RandomizeHBDirections( ice );
    
    //defectを解消する。
    PurgeDefects( ice );

#ifdef DEBUG
    //Check
    //CheckNetwork( );
#endif

    //現在プロトンがある場所を記憶する。
    ice->currentBond = -1;

    //protonと周囲の相互作用を算出する。
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
    //protonの位置を出力する。
    //
    real x,y,z;
    x = ice->bonds[ ice->currentBond ].px[1] - ice->comx;
    y = ice->bonds[ ice->currentBond ].py[1] - ice->comy;
    z = ice->bonds[ ice->currentBond ].pz[1] - ice->comz;
    
    printf("t %f %f %f %d(%f)\n", x,y,z,aux,ice->currentPotential);
}



//デバッグ用。全エネルギーを計算する。
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
    //debug用。プロトンが一つだけかどうか、outgoの数が正しいか、などを調べる。
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



//クラスタ重心からプロトンまでの距離。
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



/*無限格子近似の収束精度をチェックする。*/
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
    //情報落ちを防ぐ計算順序を使用した場合
    real epi = ProtonInPDIPCIPotentialEnergy( ice, proton, r );
    printf( "%f %24.17f %24.17f %24.17f %24.17f\n", r, epdi-epci, epdi, epci, epi );
    r += 0.1;
  }
}


//現在のプロトンの位置から、nextに動けるかどうかをチェックし、プロトンを動かす。
//エネルギー差を返す
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
                        
                        /*プロトン移動する先が、現在の結合のどちら側か*/
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



//現在のプロトンの位置と、次のプロトンの位置の間で共有している酸素の番号を返す。
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
    
    /*プロトンの場所はPPOS形式で保存。PPOSの書式はNGPHと同一。(実際にはプロトンが2個以上あると誤動作する。)*/
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
    /*周期境界条件であれば箱の大きさを出力する。*/
    if ( ice->periodic )
        fprintf( file, "@BOX3\n%f %f %f\n", ice->bx, ice->by, ice->bz );
    /*酸素の位置をAR3A形式で出力する。*/
    fprintf( file, "@AR3A\n%d\n", ice->nOxygen );
    for( o1=0; o1<ice->nOxygen; o1++ ){
        real dx,dy,dz;
        dx = ice->oxygens[o1].x;
        dy = ice->oxygens[o1].y;
        dz = ice->oxygens[o1].z;
        fprintf( file, "%f %f %f\n", dx,dy,dz );
    }
}
