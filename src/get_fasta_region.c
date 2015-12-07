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

    int start = atoi(argv[2]);
    int end   = atoi(argv[3]);

    seqptr = read_fasta(argv[1]);

    seq = seqptr;

    while (seq != NULL) {
     
      int i   = start-1;

      printf(">%s\n",seq->id);

      while (i < end) {
	printf("%c",seq->sequence[i]);

	if (i > start && (i-start)%72 == 0) {
	  printf("\n");
	}
	i++;
      }
      printf("\n");
      seq = seq->next;
    }
}

