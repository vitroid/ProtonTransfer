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



//$BF1$8%N!<%I$r2?EY$b=PNO$7$J$$9)IW$,I,MW!#(BInsert$B;~$KFM$-$"$o$;$k$3$H$K$9$k!#(B
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
    
    //$B$9$Y$F$N;@AG$K$D$$$F(B
    for( o1=0; o1<ice->nOxygen; o1++ ){
        //$B%@%_!<$G$J$1$l$P(B
        if ( ice->oxygens[o1].charge != 0.0 ){
            //$BNY@\?t$,(B4$B$G$J$1$l$P(B
            if ( ice->oxygens[o1].nAdj != 4 ){
                int j;
                //$B?eJ,;R$HNY@\J,;R$rI=<((B
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
        //$B%@%_!<$J$i(B
        if ( ice->oxygens[o1].charge == 0.0 ){
            real o1x,o1y,o1z;
            
            //$B?eJ,;R$HNY@\J,;R$rI=<((B
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
  //$B=E?4$+$i5wN%(B15$B$^$G$N%5%$%H$N$_$rC5:w$9$k!#(B
  //$B:BI8$r%U%!%$%k$+$iFI$s$@>l9g$O:BI8$N5wN%$HF1$8C10L$G;XDj$9$k(B
  real radius;

  //int optind;
  int depth;
  //char buf[100];
  //$BNY@\%5%$%H?t$,(B4$B$KK~$?$J$$%5%$%H$K$O?eJ,;R$rCV$+$J$$$G%@%_!<%5%$%H$H$7$F07$&!#(B
  int cleanse;
  FILE*  coordFile;
  FILE*  pathFile;
  //check mode. default=0
  int check;
  /*$BDL>o$O!"%W%m%H%s$N>l=j$rFI$_$3$^$J$+$C$?>l9g$O!"%M%C%H%o!<%/A4BN$r%i%s%@%`2=$7$F$+$i(Bprotonate$B$9$k!#$7$+$7!"(Bshuffle$B$,(BFalse$B$N>l9g$O!"%i%s%@%`2=$;$:$K(BProtonate$B$9$k!#(B*/
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
  //$B=E?4$+$i5wN%(B15$B$^$G$N%5%$%H$N$_$rC5:w$9$k!#(B
  //$B:BI8$r%U%!%$%k$+$iFI$s$@>l9g$O:BI8$N5wN%$HF1$8C10L$G;XDj$9$k(B
  settings->radius = 1e30;
  //int optind;
  settings->depth = 8;
  //$BNY@\%5%$%H?t$,(B4$B$KK~$?$J$$%5%$%H$K$O?eJ,;R$rCV$+$J$$$G%@%_!<%5%$%H$H$7$F07$&!#(B
  settings->cleanse=4;
  settings->coordFile = NULL;
  settings->pathFile = stdin;
  //check mode. default=0
  settings->check=0;
  /*$BDL>o$O!"%W%m%H%s$N>l=j$rFI$_$3$^$J$+$C$?>l9g$O!"%M%C%H%o!<%/A4BN$r%i%s%@%`2=$7$F$+$i(Bprotonate$B$9$k!#$7$+$7!"(Bshuffle$B$,(BFalse$B$N>l9g$O!"%i%s%@%`2=$;$:$K(BProtonate$B$9$k!#(B*/
  settings->shuffle=True;
  /*O-O$B4V$N$I$N$"$?$j$K(BH$B$rCV$/$+!#8_49@-$N$?$a%G%U%)%k%H$O(B30%$B$H$9$k$,!"<B:]$NI9$N>l9g$O(B1/2.75=0.3636...*/
  settings->hposition=0.3;
    
  /*for getopt*/
  int c;

  settings->outputFile = stdout;
  settings->traceFile = NULL;
  settings->traceEvery = 1;
  settings->graph=False;
  settings->site=False;

  settings->visualize = False;

  //$BDL>o$O(Bsearch$B$@$,%3%^%s%IL>$,(Bvisualize$B$r4^$`>l9g$O(Bvisualize$B$H$7$F5!G=$9$k!#(B
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
        //$BMp?t$N<o$r@_Dj$7$F$+$i$G$J$$$H$^$:$$!#(B
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
    
  //1$B$D$N>l9g$O$=$l$rMp?t$N<o$H$_$J$9!#(B
  //($BF1$8<o$rM?$($l$PI,$:F1$8=i4|G[CV$,F@$i$l$k!#(B)
  if ( optind < argc ){
    settings->seed =  atoi( argv[optind++] );
    if ( optind < argc )
      usage( argc, argv, settings );
  }
  return settings;
}





/*$BA4C5:w$r9T$&!#$H$j$"$($:!"%W%m%0%i%`$,4JC1$J?<$5M%@hC5:w$r9T$&!#(B*/
int main( int argc, char *argv[] )
{
  sSearch*  search = NULL;
  sLattIce* ice    = NULL;
  //etotal$B$O!"%W%m%H%s$NMOG^OB%(%M%k%.!<(B+$B<~JU$N9=B$JQ2=$K$h$kA4%(%M%k%.!<JQ2=$NOB!#(B
  real    etotal = 0;
  sNode*    node   = NULL;
  int traceLen;
  char buf[100];

  sSettings* settings = Options( argc, argv );
  /*dryrun$B$r(B1$B$K$7$F$*$/$H!"%(%M%k%.!<7W;;$,>J$+$l$k!#(B*/
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
    //$B?e$N0LCV$H7k9g>pJs$rFI$_$3$`!#(B
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
  //$BI8=`F~NO$+$i!";vA0$K%W%m%H%s$rF0$+$7$F$*$/7PO)$rFI$_$3$`!#(B
  //2008-2-10$B2CI.!#(Bpath file$B$NBh(B2$B%+%i%`$K%(%M%k%.!<$,=q$$$F$"$l$P!"$=$A$i$r:NMQ$7!"%(%M%k%.!<7W;;$r>J$/!#$3$l$G$+$J$jB.$/$G$-$k!#(B
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
      //$BESCf$N%(%M%k%.!<7W;;$r>J$$$?$N$G!":G8e$K:F7W;;$7$F!"3NG'$9$k!#(B
      if ( ! dryrun )
	ice->currentPotential = ProtonPotentialEnergy( ice, ice->currentBond, 0.0L );
      //fprintf("%f %f\n", etotal, ice->currentPotential);
      //exit(1);
    }
  }
  if ( settings->visualize ){
    if ( settings->graph ){
      //$B?eAG7k9g%M%C%H%o!<%/%H%]%m%8!<$rM-8~%0%i%U7A<0$G=PNO(B
      //$BH>7B$,(B0$B$G$J$$>l9g$O!"H>7BFb$N$_$r=PNO!#(B
      saveNGPH( search->ice, settings->radius, settings->outputFile );
      //saveNGPH$B$O%W%m%H%s$N0LCV$rJ]B8$7$J$$$N$G!"B>$NJ}K!$GJ]B8$9$kI,MW$,$"$k!#(B
      saveProtonPosition( search->ice, settings->outputFile );
    }
    if ( settings->site ){
      /*$B;@AG86;R$N0LCV$r%;!<%V$9$k!#(BNGPH$B$N>l9g$H0c$$!"H>7B;XDj$OL58z!#(B*/
      saveSite( search->ice, settings->outputFile );
    }
    fclose( settings->outputFile );
    if ( settings->traceFile )
      fclose( settings->traceFile );
  }
  else{
    /*$B?<$5(B8$B$^$GA4C5:w(B*/
    SearchRecursively( search, settings->depth, settings->radius, etotal, dryrun, settings->bidir );
  }
  exit(0);
}
