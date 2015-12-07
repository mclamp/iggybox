
#include "viewport.h"

const signed char fasta_encoding[] = {
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, 0, _, 1, _, _, _, 2, _, _, _, _, _, _, _, _,
 _, _, _, _, 3, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
};


int nmer2pow(int n) {
  int i=0;
  int val = 1;

  while (i < n*2) {
    val = val * 2;
    i++;
  }
  return val;
}

int *generate_bitstring(int n, int m) {
  /* Generates a random n of m bitstring */

  int rand_num;
  int i = 0;
  int *bits;

  bits = (int *)malloc(m*sizeof(int));
  
  for (i=0; i < m;i++) {
    bits[i] = 1;
  }

  i=0;

  while (i < (m-n)) {
    rand_num = (int)((double)m*rand()*1.0/(1.0*RAND_MAX));
    
    if (bits[rand_num] != 0) {
      bits[rand_num] = 0;
      i++;
    }
  }

  for(i=0;i< m;i++) {
    printf("%d",bits[i]);
  }

  printf("\n\n");

  return bits;
}

Comb *bits2comb(int *bits, int n,int m) {
  int i      = 0;
  int length = 0;
  int start  = 0;

  int lshift;

  BitBlock *prev;
  Comb     *comb;

  int numblocks = 0;

  lshift = 0;

  comb = (Comb *)malloc(sizeof(Comb));

  comb->n = n;
  comb->m = m;

  for (i=0; i < m; i++) {
    printf("%d",bits[i]);
  }
  
  printf("\n");
  
  i=0;

  while (i < m) {

    if (i > 0) {
      if (bits[i] == 1 && bits[i-1] == 0) {
	start = i;
	length = 1;
	lshift++;
	numblocks++;
      } else if (bits[i] == 0 && bits[i-1] == 1) {
	if (numblocks == 1) {
	  comb->blocks = (BitBlock *)malloc(sizeof(BitBlock));

	  prev = comb->blocks;
	  prev->next = NULL;
	} else {
	  prev->next = (BitBlock *)malloc(sizeof(BitBlock));

	  prev = prev->next;
	  prev->next = NULL;
	}

	prev->length      = length;

	prev->leftshift   = (m - start - length);
	prev->rightshift  = (prev->leftshift  - n + lshift); 

      } else if (bits[i] == 1 && bits[i-1] == 1) {
	length++;
	lshift++;
      } 
    } else {
      if (bits[i] == 1) {
	numblocks++;
	start = i;
	length = 1;
	lshift++;
      }
    }
      
    i++;
  }

  if (bits[m-1] == 1) {
    if (numblocks == 1) {
      comb->blocks = (BitBlock *)malloc(sizeof(BitBlock));
      
      prev = comb->blocks;
      prev->next = NULL;
    } else {
      prev->next = (BitBlock *)malloc(sizeof(BitBlock));
      
      prev = prev->next;
      prev->next = NULL;
    }
    
    prev->length      = length;
    prev->leftshift   = (m - start - length);
    prev->rightshift  = (prev->leftshift  - n + lshift); 
  }
  
  comb->numblocks = numblocks;

  return comb;
}

unsigned int map_comb(Comb *comb,unsigned long long i, int print) {

  unsigned int j = 0;
  
  BitBlock *block;
  
  block = comb->blocks;
  
  while (block != NULL) {
    j = j | ((i &(B(block->length) << (2*block->leftshift))) >> (2*block->rightshift));

    block = block->next;
  }
  return j;
}

void find_hash(Sequence *head,int n,Viewport *viewport) {
  Hash ** hash;
  Hash  * tmphash;
  Comb  * comb;
  int   * bits;
  int     i;

  static unsigned long long ei = 0;
  static unsigned int ecode;

  int e;

  Sequence *seq;

  hash = (Hash **)malloc(nmer2pow(n) * sizeof(Hash *));

  bits = (int *)malloc(sizeof(int)*n);
  bits = generate_bitstring(n,n);

  comb = bits2comb(bits,n,n);

  memset(hash,0,nmer2pow(n)*sizeof(Hash *));

  seq = head;

  while (seq != NULL) {
    i = 0;

    while (i < n-1) {
      
      e = fasta_encoding[(int)seq->sequence[i]];
      ei = (ei << 2) | e;
      i++;
    }
    
    i = n-1;
    
    while (i < seq->length) {
      e = fasta_encoding[(int)seq->sequence[i]];
      
      if (e >=0) {
	ei = (ei << 2) | e;
	ecode = map_comb(comb,ei,0);
	
	tmphash = (Hash *)malloc(sizeof(Hash));
	tmphash->pos = i;
	
	if (hash[ecode] != NULL) {
	  tmphash->count = hash[ecode]->count + 1;
	} else {
	  tmphash->count = 1;
	}
	
	tmphash->next = hash[ecode];
	hash[ecode] = tmphash;
	//print_ll(stdout,ei,m);
	//printf(" Ecode %d ",ecode);
	//print_l(stdout,ecode,n);
	//printf("\n");
      }
      i++;
    }
    seq = seq->next;
  }

  viewport->hash = hash;
  viewport->comb = comb;

  printf("Comb %d\n",comb->m);

}
