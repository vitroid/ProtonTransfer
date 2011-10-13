/*$B@0?t(B(0$B$r4^$^$J$$(B)$B$rMWAG$H$9$k(BHash$B!#(B*/

#include <stdio.h>
#include <stdlib.h>
#include "IntHash.h"
#define EMPTY -1

/*shift$B$H(Bxor$B$K$h$k4JC1$J(Bhash key$B$N@8@.(B*/
unsigned int _IntHash_Encode(sIntHash *s,unsigned int key)
{
    unsigned int mod=0;
    while(key)
      {
          mod ^= (key & (s->hashsize-1));
          key >>= s->shift;
      }
    return mod;
}

/*hash$BMWAG$NHV9f$rJV$9!#$b$7B8:_$7$J$$MWAG$J$i!"6u$-MWAG$rJV$9!#MWAG$ODI2C$5$l$J$$!#(B*/
int _IntHash_QueryElement(sIntHash *ih,unsigned int key)
{
  int e=_IntHash_Encode(ih,key);
  while(1){
    if(ih->key[e]==EMPTY){
      ih->value[e]=0;
      return e;
    }
    if(ih->key[e]==key){
      return e;
    }
    e+=13;
    if(e>=ih->hashsize)
      e-=ih->hashsize;
  }
}

/*$BCM$rEPO?$9$k!#99?7$J$i(B0$B!"DI2C$J$i(B1$B$rJV$9!#(B*/
int IntHash_RegisterValue(sIntHash *ih,unsigned int key,int value)
{
  int e=_IntHash_QueryElement(ih,key);
  int v=ih->value[e];
  if(v==0){
    ih->key[e]=key;
    ih->nentry++;
    if(ih->nentry > ih->hashsize/2){
      fprintf(stderr,"Warning: hash size seems too small.\n");
    }
    if(ih->nentry >= ih->hashsize){
      fprintf(stderr,"Error: hash overflow.\n");
      exit(1);
    }
  }
  ih->value[e]=value;
  return (v==0);
}

/*$BCM$r;2>H$9$k!#(B*/
int IntHash_QueryValue(sIntHash *ih,unsigned int key)
{
  return ih->value[_IntHash_QueryElement(ih,key)];
}

/*$BCM$rKu>C$9$k!#(B*/
void IntHash_EraseOne(sIntHash *ih,unsigned int key)
{
  ih->key[_IntHash_QueryElement(ih,key)]=EMPTY;
}

sIntHash *IntHash_Init(int size)
{
  sIntHash *ih=malloc(sizeof(sIntHash));
  int m,i;
  m=1<<size;
  ih->shift=size;
  ih->hashsize=m;
  ih->key=malloc(m*sizeof(int));
  ih->value=calloc(m,sizeof(int));
  ih->nentry=0;
  for(i=0;i<m;i++)
    ih->key[i]=EMPTY;
  return ih;
}

void IntHash_Done(sIntHash *ih)
{
  free(ih->key);
  free(ih->value);
  free(ih);
}
