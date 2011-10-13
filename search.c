#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include "LattIce.h"
#include "SearchLib.h"
#include "Utility.h"
#include "common.h"

typedef struct sNode
{
    struct sNode* next;
    int           value;
}
sNode;


sNode* LookUp( sNode* node, int value )
{
    while ( node != NULL ){
        if ( node->value == value )
            return node;
        node = node->next;
    }
    return NULL;
}


sNode* Insert( sNode* node, int value )
{
    if ( NULL == LookUp( node, value ) ){
        sNode* newnode;
        newnode = malloc(sizeof(sNode));
        newnode->next=node;
        newnode->value = value;
        return newnode;
    }
    else
        return node;
}



//同じノードを何度も出力しない工夫が必要。Insert時に突きあわせることにする。
void ShowProtonInYaplot( sLattIce* ice, int bond, FILE* file )
{
    int    position,o1,o2;
    real px,py,pz,o1x,o1y,o1z,o2x,o2y,o2z;

    position = ice->bonds[ bond ].direction + 1;
    px = ice->bonds[ bond ].px[position] - ice->comx;
    py = ice->bonds[ bond ].py[position] - ice->comy;
    pz = ice->bonds[ bond ].pz[position] - ice->comz;

    o1 = ice->bonds[ bond ].oxygen[0];
    o1x = ice->oxygens[ o1 ].x - ice->comx;
    o1y = ice->oxygens[ o1 ].y - ice->comy;
    o1z = ice->oxygens[ o1 ].z - ice->comz;
    
    o2 = ice->bonds[ bond ].oxygen[1];
    o2x = ice->oxygens[ o2 ].x - ice->comx;
    o2y = ice->oxygens[ o2 ].y - ice->comy;
    o2z = ice->oxygens[ o2 ].z - ice->comz;
    
    //green
    fprintf( file, "@ 3\n");
    //draw line
    fprintf( file, "l %f %f %f %f %f %f\n", o1x,o1y,o1z,o2x,o2y,o2z );
    //white
    if ( position == 1 )
        fprintf( file, "@ 5\n");
    else
        fprintf( file, "@ 2\n");
    //radius=0.1
    fprintf( file, "r 0.1\n");
    //put circle
    fprintf( file, "o %f %f %f\n", px,py,pz );

}


void ShowWaterInYaplot( sLattIce* ice, int o1, FILE* file )
{
    real px,py,pz,o1x,o1y,o1z;
    int    j;
    
    o1x = ice->oxygens[ o1 ].x - ice->comx;
    o1y = ice->oxygens[ o1 ].y - ice->comy;
    o1z = ice->oxygens[ o1 ].z - ice->comz;
    //white
    fprintf( file, "@ 2\n");
    //radius=0.3
    fprintf( file, "r 0.3\n");
    fprintf( file, "o %f %f %f\n", o1x,o1y,o1z );
    
    for( j=0; j< ice->oxygens[o1].nAdj; j++ ){
        int dir = BondDirection( ice, o1, j );
        if ( dir == OUTGOING ){
            int position;
            int bond = ice->oxygens[o1].hydrogenBond[j];
            position = ice->bonds[ bond ].direction + 1;
            px = ice->bonds[ bond ].px[position] - ice->comx;
            py = ice->bonds[ bond ].py[position] - ice->comy;
            pz = ice->bonds[ bond ].pz[position] - ice->comz;
            //green
            fprintf( file, "@ 3\n");
            //radius=0.2
            fprintf( file, "r 0.2\n");
            fprintf( file, "o %f %f %f\n", px,py,pz );
            //draw line : OH-bond
            fprintf( file, "l %f %f %f %f %f %f\n", o1x, o1y, o1z, px,py,pz );
        }
    }
}


void ShowPathsInYaplot( sLattIce* ice, sNode* node, FILE* file )
{
    while ( node != NULL ){
        ShowProtonInYaplot( ice, node->value, file );
        node = node->next;
    }
    fprintf( file, "\n");
}

void ShowRimInYaplot( sLattIce* ice, FILE* file )
{
    int o1;
    
    //すべての酸素について
    for( o1=0; o1<ice->nOxygen; o1++ ){
        //ダミーでなければ
        if ( ice->oxygens[o1].charge != 0.0 ){
            //隣接数が4でなければ
            if ( ice->oxygens[o1].nAdj != 4 ){
                int j;
                //水分子と隣接分子を表示
                fprintf( file, "y %d\n",ice->oxygens[o1].nAdj );
                ShowWaterInYaplot( ice, o1, file );
                for( j=0; j< ice->oxygens[o1].nAdj; j++ ){
                    int o2 = ice->oxygens[o1].oxygen[j];
                    ShowWaterInYaplot( ice, o2, file );
                }
            }
        }
    }
    fprintf( file, "y 6\n");
    for( o1=0; o1<ice->nOxygen; o1++ ){
        //ダミーなら
        if ( ice->oxygens[o1].charge == 0.0 ){
            real o1x,o1y,o1z;
            
            //水分子と隣接分子を表示
            o1x = ice->oxygens[ o1 ].x - ice->comx;
            o1y = ice->oxygens[ o1 ].y - ice->comy;
            o1z = ice->oxygens[ o1 ].z - ice->comz;
            //white
            fprintf( file, "@ 4\n");
            //radius=0.3
            fprintf( file, "r 0.3\n");
            fprintf( file, "o %f %f %f\n", o1x,o1y,o1z );
        }
    }
    fprintf( file, "\n");
}





typedef struct {
  int seed;
  int latticeSize;

  bool bidir;
  int    noenergy;
  //重心から距離15までのサイトのみを探索する。
  //座標をファイルから読んだ場合は座標の距離と同じ単位で指定する
  real radius;

  //int optind;
  int depth;
  //char buf[100];
  //隣接サイト数が4に満たないサイトには水分子を置かないでダミーサイトとして扱う。
  int cleanse;
  FILE*  coordFile;
  FILE*  pathFile;
  //check mode. default=0
  int check;
  /*通常は、プロトンの場所を読みこまなかった場合は、ネットワーク全体をランダム化してからprotonateする。しかし、shuffleがFalseの場合は、ランダム化せずにProtonateする。*/
  bool   shuffle;

  bool   visualize;

  FILE*  outputFile;
  FILE*  traceFile;
  int    traceEvery;
  bool   graph;
  bool   site;
  real hposition;
} sSettings;



void usage( int argc, char* argv[], sSettings* settings )
{
    fprintf( stderr, "Usage: %s [options] [randomseed]\n", argv[0] );
    fprintf( stderr, "Options:\n" );
    fprintf( stderr, "\t--lattice=file Read lattice structure from the file.\n" ); 
    fprintf( stderr, "\t--pathfile=file Read premotion path from the file.\n" ); 
    fprintf( stderr, "\t--compose=N    Compose ice Ic lattice with size=NxNxN ( N is an even number ).\n" ); 
    fprintf( stderr, "\t--cleanse=N    Treat the nodes having less than N neighborhoods as dummy nodes ( default=4 ).\n" ); 
    fprintf( stderr, "\t--hposition=x  Place the H atom at fraction x between two neighboring oxygens. ( default=0.3 ).\n" ); 
    fprintf( stderr, "\t--noshuffle    Use read network as is and do not shuffle before protonation ( default=do shuffle ).\n" ); 
    fprintf( stderr, "\t--depth=N      Search depth ( default=8 ).\n" ); 
    fprintf( stderr, "\t--check=N      Check mode. 0: none, 1:ILA convergence ( default=0 ).\n" ); 
    fprintf( stderr, "\t--noenergy     Do not calculate energy change during proton xfer.\n" ); 
    fprintf( stderr, "\t--bidir        Allow to go backward during proton xfer.\n" ); 

    if ( settings->visualize ){
      fprintf( stderr, "\t--trace=filename Save trace of proton in yaplot format.\n" ); 
      fprintf( stderr, "\t--trace-every=x  Save trace of proton every x frames.\n" ); 
      fprintf( stderr, "\t--graph         Network topology information is saved to the file ( in @NGPH format ).\n" ); 
      fprintf( stderr, "\t--site         Oxygen site positions are saved to to the file ( in @AR3A format ).\n" ); 
      fprintf( stderr, "\t--output=filename Output file. ( default=stdout ).\n" ); 
      fprintf( stderr, "\t--radius=x     Set outer limit radius to visualize ( default=0,\n\t\t\ti.e. all the edges will be shown. ).\n" ); 
    }
    else{
      fprintf( stderr, "\t--radius=x     Set outer limit radius for search ( default=1e30, i.e. unlimited ).\n" ); 
    }
    exit( 1 );
}



sSettings* Options( int argc, char *argv[] )
{
  sSettings* settings = malloc( sizeof( sSettings ) );

  settings->seed=1234;
  settings->latticeSize=0;
  settings->noenergy=0;
  settings->bidir = False;
  //重心から距離15までのサイトのみを探索する。
  //座標をファイルから読んだ場合は座標の距離と同じ単位で指定する
  settings->radius = 1e30;
  //int optind;
  settings->depth = 8;
  //隣接サイト数が4に満たないサイトには水分子を置かないでダミーサイトとして扱う。
  settings->cleanse=4;
  settings->coordFile = NULL;
  settings->pathFile = stdin;
  //check mode. default=0
  settings->check=0;
  /*通常は、プロトンの場所を読みこまなかった場合は、ネットワーク全体をランダム化してからprotonateする。しかし、shuffleがFalseの場合は、ランダム化せずにProtonateする。*/
  settings->shuffle=True;
  /*O-O間のどのあたりにHを置くか。互換性のためデフォルトは30%とするが、実際の氷の場合は1/2.75=0.3636...*/
  settings->hposition=0.3;
    
  /*for getopt*/
  int c;

  settings->outputFile = stdout;
  settings->traceFile = NULL;
  settings->traceEvery = 1;
  settings->graph=False;
  settings->site=False;

  settings->visualize = False;

  //通常はsearchだがコマンド名がvisualizeを含む場合はvisualizeとして機能する。
  if ( strstr( argv[0], "visualize" ) ){
    settings->visualize = True;
    settings->radius = 0;
  }

  static struct option long_options[] = {
    {"lattice", 1, 0, 0},
    {"cleanse", 1, 0, 0},
    {"hposition", 1, 0, 0},
    {"compose", 1, 0, 0},
    {"pathfile", 1, 0, 0},
    {"depth"  , 1, 0, 0},
    {"check"  , 1, 0, 0},
    {"noshuffle"  , 0, 0, 0},
    {"radius",  1, 0, 0},
    {"trace" , 1, 0, 0},
    {"trace-every" , 1, 0, 0},
    {"output" , 1, 0, 0},
    {"graph"  , 0, 0, 0},
    {"site"  , 0, 0, 0},
    {"noenergy", 0, 0, 0},
    {"bidir", 0, 0, 0},
    {0, 0, 0, 0}
  };

  while (1) {
    int option_index = 0;
    
    c = getopt_long (argc, argv, "",
                       long_options, &option_index);

    if (c == -1)
      break;
    
    bool done=False;
    switch (c) {
    case 0:
      done = True;
      if ( strcmp( long_options[option_index].name, "lattice" ) == 0 ){
        settings->coordFile = fopen( optarg, "r" );
      }
      else if ( strcmp( long_options[option_index].name, "pathfile" ) == 0 ){
        settings->pathFile = fopen( optarg, "r" );
      }
      else if ( strcmp( long_options[option_index].name, "compose" ) == 0 ){
        settings->latticeSize = atoi( optarg );
        //乱数の種を設定してからでないとまずい。
        //ice = new_ProtonatedIceIc( latticeSize );
      }
      else if ( strcmp( long_options[option_index].name, "hposition" ) == 0 ){
        settings->hposition = atof( optarg );
      }
      else if ( strcmp( long_options[option_index].name, "cleanse" ) == 0 ){
        settings->cleanse = atoi( optarg );
      }
      else if ( strcmp( long_options[option_index].name, "noshuffle" ) == 0 ){
        settings->shuffle=False;
      }
      else if ( strcmp( long_options[option_index].name, "depth" ) == 0 ){
        settings->depth = atoi( optarg );
      }
      else if ( strcmp( long_options[option_index].name, "check" ) == 0 ){
        settings->check = atoi( optarg );
      }
      else if ( strcmp( long_options[option_index].name, "radius" ) == 0 ){
	settings->radius = atof( optarg );
      }
      else if ( strcmp( long_options[option_index].name, "noenergy" ) == 0 ){
	settings->noenergy ++;
      }
      else if ( strcmp( long_options[option_index].name, "bidir" ) == 0 ){
	settings->bidir = True;
      }
      else {
        if ( settings->visualize ){
          if ( strcmp( long_options[option_index].name, "graph" ) == 0 ){
            settings->graph ++;
          }
          else if ( strcmp( long_options[option_index].name, "site" ) == 0 ){
            settings->site ++;
          }
          else if ( strcmp( long_options[option_index].name, "output" ) == 0 ){
            settings->outputFile = fopen( optarg, "w" );
          }
          else if ( strcmp( long_options[option_index].name, "trace" ) == 0 ){
            settings->traceFile = fopen( optarg, "w" );
          }
          else if ( strcmp( long_options[option_index].name, "trace-every" ) == 0 ){
            settings->traceEvery = atoi( optarg );
          }
          else
            done = False;
        }
        else{
          //else
	  done = False;
        }
      }
      if ( ! done ){
        printf ("option %s", long_options[option_index].name);
        if (optarg)
          printf (" with arg %s", optarg);
        printf ("\n");
        usage( argc, argv, settings );
        exit(1);
      }
      break;
            
    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
      usage( argc, argv, settings );
    }
  }
    
  //1つの場合はそれを乱数の種とみなす。
  //(同じ種を与えれば必ず同じ初期配置が得られる。)
  if ( optind < argc ){
    settings->seed =  atoi( argv[optind++] );
    if ( optind < argc )
      usage( argc, argv, settings );
  }
  return settings;
}





/*全探索を行う。とりあえず、プログラムが簡単な深さ優先探索を行う。*/
int main( int argc, char *argv[] )
{
  sSearch*  search = NULL;
  sLattIce* ice    = NULL;
  //etotalは、プロトンの溶媒和エネルギー+周辺の構造変化による全エネルギー変化の和。
  real    etotal = 0;
  sNode*    node   = NULL;
  int traceLen;
  char buf[100];

  sSettings* settings = Options( argc, argv );
  /*dryrunを1にしておくと、エネルギー計算が省かれる。*/
  bool dryrun = settings->noenergy;

  fprintf( stderr, "Seed: %d\n", settings->seed );
  srand48( settings->seed );
  
  if ( settings->latticeSize ){
    ice = new_ProtonatedIceIc( settings->latticeSize, settings->hposition, dryrun );
  }
  
  if ( settings->coordFile != NULL ) {
    if ( ice != NULL ){
      fprintf( stderr, "Error: Lattice is specified twice.\n" );
      usage( argc, argv, settings );
    }
    //水の位置と結合情報を読みこむ。
    ice = LoadLatticeAndProtonate( settings->coordFile, settings->shuffle, settings->cleanse, settings->hposition, dryrun );
    fclose( settings->coordFile );
  }
  
  if ( settings->check == 2 ){
    CheckConsistency( ice, settings->radius );
    exit(0);
  }
#ifdef ILA
  else if ( settings->check == 1 ){
    CheckConvergence( ice );
    exit(0);
  }
#endif
  search = new_Search( ice, settings->depth );
  
#ifdef DEBUG
  //just for debug
  if ( settings->visualize )
    ShowRimInYaplot( search->ice, settings->traceFile );
#endif
  if ( ! dryrun )
    etotal = ProtonPotentialEnergy( ice, ice->currentBond, 0.0L );
  //etotal = 0;
#ifdef DEBUG
  etotal = TotalEnergy( ice );
#endif
  //   fprintf( stderr, "[0: %f]\n", etotal );
  traceLen = 0;
  //標準入力から、事前にプロトンを動かしておく経路を読みこむ。
  //2008-2-10加筆。path fileの第2カラムにエネルギーが書いてあれば、そちらを採用し、エネルギー計算を省く。これでかなり速くできる。
  while( NULL != fgets( buf, sizeof(buf), settings->pathFile ) ){
    int bond;
    real absenergy;
#ifdef SINGLEPRECISION
    int col = sscanf( buf, "%d %lf", &bond, &absenergy );
#else
    int col = sscanf( buf, "%d %lf", &bond, &absenergy );
#endif
    if ( settings->visualize )
      node = Insert( node, bond );
    if ( settings->traceFile )
      ShowPathsInYaplot( search->ice, node, settings->traceFile );
    if ( search->ice->currentBond != bond ){
      die( "First bond of the proton should be %d\n", search->ice->currentBond );
    }
    while ( NULL != fgets( buf, sizeof(buf), settings->pathFile ) ){
#ifdef SINGLEPRECISION
      col = sscanf( buf, "%d %lf", &bond, &absenergy );
#else
      col = sscanf( buf, "%d %lf", &bond, &absenergy );
#endif
      if ( bond != search->ice->currentBond ){
        if ( col == 1 ){
          etotal += SearchMove( search, bond, settings->depth, False/*dryrun*/ );
          //fprintf( stderr, "%f %f\n", etotal, absenergy );
          }
        else{
          SearchMove( search, bond, settings->depth, True/*dryrun*/ );
          etotal = absenergy;
        }
      }
      fprintf( stderr, "[%d: %f]\r", ++traceLen, etotal );
      
#ifdef DEBUG
      printf( "DeltaSigma=%f, Absolute=%f\n", etotal, TotalEnergy( search->ice ) );
#endif
      if ( settings->visualize ){
        node = Insert( node, bond );
        if ( settings->traceFile && (traceLen % settings->traceEvery == 0 ) )
          ShowPathsInYaplot( search->ice, node, settings->traceFile );
      }
    }
    if ( col == 2 ){
      //途中のエネルギー計算を省いたので、最後に再計算して、確認する。
      if ( ! dryrun )
	ice->currentPotential = ProtonPotentialEnergy( ice, ice->currentBond, 0.0L );
      //fprintf("%f %f\n", etotal, ice->currentPotential);
      //exit(1);
    }
  }
  if ( settings->visualize ){
    if ( settings->graph ){
      //水素結合ネットワークトポロジーを有向グラフ形式で出力
      //半径が0でない場合は、半径内のみを出力。
      saveNGPH( search->ice, settings->radius, settings->outputFile );
      //saveNGPHはプロトンの位置を保存しないので、他の方法で保存する必要がある。
      saveProtonPosition( search->ice, settings->outputFile );
    }
    if ( settings->site ){
      /*酸素原子の位置をセーブする。NGPHの場合と違い、半径指定は無効。*/
      saveSite( search->ice, settings->outputFile );
    }
    fclose( settings->outputFile );
    if ( settings->traceFile )
      fclose( settings->traceFile );
  }
  else{
    /*深さ8まで全探索*/
    SearchRecursively( search, settings->depth, settings->radius, etotal, dryrun, settings->bidir );
  }
  exit(0);
}
