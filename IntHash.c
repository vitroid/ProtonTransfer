/*整数(0を含まない)を要素とするHash。*/

#include <stdio.h>
#include <stdlib.h>
#include "IntHash.h"
#define EMPTY -1

/*shiftとxorによる簡単なhash keyの生成*/
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

/*hash要素の番号を返す。もし存在しない要素なら、空き要素を返す。要素は追加されない。*/
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

/*値を登録する。更新なら0、追加なら1を返す。*/
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

/*値を参照する。*/
int IntHash_QueryValue(sIntHash *ih,unsigned int key)
{
  return ih->value[_IntHash_QueryElement(ih,key)];
}

/*値を抹消する。*/
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
