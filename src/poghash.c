#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "datamodel.h"
#include "pog_utils.h"

#define _ (-1)

extern Sequence *read_sequence(FILE *file);
extern const signed char fasta_encoding[256];

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

typedef struct hashpos {
  int pos;		/* position where word starts in sequence 1 */
  struct hashpos *next;
} Hash;

static void print_ll(FILE *fp, unsigned long long x)
{
    int i;
    for (i = 36; i >= 0; i -= 2)
	fprintf(fp, "%c", "ACGT"[(x>>i)&3]);
}

static void print_l(FILE *fp, unsigned long x)
{
    int i;
    for (i = 22; i >= 0; i -= 2)
	fprintf(fp, "%c", "ACGT"[(x>>i)&3]);
}


#define B(n) ((1ULL<<(n*2))-1)

//              8765432109876543210
// template     1110100110010101111 (12 of 19)

static unsigned int map_19to12(unsigned long long i)
{
        unsigned int j =
                ((i & (B(3) << (2*16))) >> (2*7)) |
                ((i & (B(1) << (2*14))) >> (2*6)) |
                ((i & (B(2) << (2*10))) >> (2*4)) |
                ((i & (B(1) << (2* 7))) >> (2*2)) |
                ((i & (B(1) << (2* 5))) >> (2*1)) |
                ((i & (B(4) << (2* 0))) >> (2*0)) ;
	// DEBUG
	fprintf(stdout, "%016llx %016x\n", i, j);
	print_ll(stdout, i); 
	fprintf(stdout, " ");
	print_l(stdout, j); 
	fprintf(stdout, "\n");
        return j;
}

int main(int argc, char **argv) {
  
  /*char *seq = "ATGCATGCATGCATGCATG";*/


  char *infilename = argv[1];


  FILE *infile = fopen(infilename,"r");

  if (infile == NULL) {
    printf("No file [%s] found\n",infilename);
    exit(0);
  }

  
  Sequence *seqptr;

  while ((seqptr = read_sequence(infile)) != NULL) {

    int i;
    int e;
    unsigned long long ei = 0;
    int ecode;
    
    char *seq = seqptr->sequence;

    for (i=0; i < seqptr->length;i++) {

      printf("%d %c\n",fasta_encoding[seq[i]],seq[i]);
      
      e = fasta_encoding[seq[i]];
      
      printf("%d %lld\n",e,ei);
      
      ei = (ei << 2) | e;
      
      printf("%lld\n",ei);
      
      
    }
    ecode = map_19to12(ei);
  }
  return 1;
}
