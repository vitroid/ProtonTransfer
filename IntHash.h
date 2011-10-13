#ifndef _INTHASH_H
#define _INTHASH_H

typedef struct
{
  int nentry;
  int hashsize,shift;
  int *key;
  int *value;
}
sIntHash;

int IntHash_RegisterValue(sIntHash *ih,unsigned int key,int value);
int IntHash_QueryValue(sIntHash *ih,unsigned int key);
void IntHash_EraseOne(sIntHash *ih,unsigned int key);
sIntHash *IntHash_Init(int size);
void IntHash_Done(sIntHash *ih);

#endif
