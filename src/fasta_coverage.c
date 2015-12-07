#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "datamodel.h"

extern Sequence *read_fasta(char *filename);

int main(int argc, char * argv[]) {

    Sequence *seq;
    Sequence *seqptr;

    int i = 0;
    int numseqs = 0;

    seqptr = read_fasta(argv[1]);

    seq = seqptr;
    
    int cov[200];

    int len = seqptr->length;

    while (seq != NULL) {
      numseqs++;
      seq = seq->next;
    }

    while (i < 200) {
      cov[i] = 0;
      i++;
    }

    i = 0;

    int count = 0;

    while (i < len) {
      seq = seqptr;
      seq = seqptr->next;

      int nongap = 0;

      while (seq != NULL) {
	if (seq->sequence[i] != '-') {
	  nongap++;
	}
	seq = seq->next;
      }
      cov[nongap]++;

      if (nongap > 0) {
	count++;
      }

      i++;
    }

    i = 0;

    while (i < numseqs) {
      printf("%d\t%d\t%d\n",i,cov[i],count);
      i++;
    }
}

