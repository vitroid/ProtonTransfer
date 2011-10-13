#ifndef _VARRAY_H
#define _VARRAY_H
/*可変長配列*/

typedef struct
{
  int initval;
  size_t size;  /*size of an entity*/
  int nmemb;    /*number of members == last member + 1*/
  int capacity; /*allocated array size*/
  char *a;
} varray;

extern varray * varray_init(size_t nmemb,size_t size,int initval);
extern void *varray_ptr(varray *v,int i);
extern int  varray_size(varray *v);
extern void varray_done(varray *v);


#endif
