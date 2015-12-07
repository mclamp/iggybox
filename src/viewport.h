#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "datamodel.h"

#define _ (-1)
#define B(n) ((1ULL<<(n*2))-1)

extern const signed char fasta_encoding[256];

typedef struct hashpos {
  unsigned int pos;	   
  struct hashpos *next;
  int count;
} Hash;

typedef struct _Lis {
  int query_start;
  int query_end;
  int target_start;
  int target_end;
  int count;
  struct _Lis *next;
} Lis;

typedef struct _BitBlock {
  int length;
  int leftshift;
  int rightshift;
  struct _BitBlock *next;
} BitBlock;


typedef struct _Comb {
  int n;
  int m;
  BitBlock *blocks;
  int numblocks;
} Comb;


typedef struct _Viewport {
  int char_width;
  int char_height;

  int old_hvalue;
  int old_vvalue;

  int width;
  int height;

  int color_by_consensus;
  int color_by_blocks;
  int first_sequence_color;
  int color_by_hash;

  Sequence *consensus;
  Block    **blocks;
  Hash     **hash;
  Comb     *comb;

} Viewport, *Viewport_ptr;

extern void  find_hash(Sequence *head,int n,Viewport *viewport);
extern int  *generate_bitstring(int n, int m);
extern Comb *bits2comb(int *bits, int n,int m);
extern int   nmer2pow(int n);
extern unsigned int map_comb(Comb *comb,unsigned long long i, int print);


