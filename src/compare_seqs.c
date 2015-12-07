#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "datamodel.h"

extern Sequence *read_fasta(char *filename);

int main(int argc, char * argv[]) {
    FILE *file;

    Sequence *seq;
    Sequence *seqptr;
    Sequence *newseqptr;
    Sequence *newseq;

    int numseqs = 0;
    int i = 0;
    int count = 0;
    int num = 0;

    char c;


    seqptr = read_fasta(argv[1]);

    seq = seqptr;

    while (seq != NULL) {
      numseqs++;
      seq = seq->next;

    }

    if (numseqs != 3) {
      printf("Number of sequences for comparison needs to be 3 - currently %d\n",numseqs);
      exit(0);
    }

    int j = 0;

    seq = seqptr;

    int cov2 = 0;
    int cov3 = 0;

    int len1 = 0;

    int bases_same = 0;
    int bases_differ = 0;

    Sequence *seq1;
    Sequence *seq2;
    Sequence *seq3;

    seq1 = seq;

    seq2 = seq1->next;
    seq3 = seq2->next;
    
    int len = strlen(seq1->sequence);

    while (j < len) {
      char c1 = seq1->sequence[j];
      char c2 = seq2->sequence[j];
      char c3 = seq3->sequence[j];

      if (c1 != '-') {
	len1++;

	if (c2 != '-' &&
	    c3 != '-') {
	  
	  if (c2 == c3) {
	    bases_same++;
	  } else {
	    //	    printf("Differing bases %c %c\n",c2,c3);
	    bases_differ++;
	  }
	}
	if (c2 != '-') {
	  cov2++;
	}
	if (c3 != '-') {
	  cov3++;
	}
	
      }

      j++;
    }

    printf("Length %d cov1 %d %d cov2 %d %d percent tot coverage %d same %d differ %d fraction different %d\n",len1,cov2,(int)(100*cov2/len1),cov3,(int)(100*cov3/len1),(int)(100*cov3/cov2),bases_same,bases_differ,(int)(100*bases_differ/(bases_same+bases_differ)));
}
