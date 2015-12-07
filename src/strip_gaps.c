#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "datamodel.h"

extern Sequence *read_fasta(char *filename);

int main(int argc, char * argv[]) {

    Sequence *seq;
    Sequence *seqptr;
    Sequence *tmpseq;

    int i = 0;
    int numseqs = 0;

    seqptr = read_fasta(argv[1]);

    seq = seqptr;

    int len = seq->length;

    while (seq != NULL) {
      numseqs++;
      seq = seq->next;
    }

    int pos = -1;

    while (i < len) {

      tmpseq = seqptr;

      int found = 0;

      while (tmpseq != NULL) {
	if (tmpseq->sequence[i] != '-') {
	  found = 1;
	  tmpseq = NULL;
	} else {
	  tmpseq = tmpseq->next;
	}
      }

      if (found == 1) {
	tmpseq = seqptr;
	pos++;
	
	
	while (tmpseq != NULL) {
	  tmpseq->sequence[pos] = tmpseq->sequence[i];
	  tmpseq = tmpseq->next;
	}
      }

      i++;
    }
    


    tmpseq = seqptr;

    while (tmpseq != NULL) {
      tmpseq->sequence[pos] = '\0';
      tmpseq = tmpseq->next;
    }

    tmpseq = seqptr;

    int count = 0;
    


    while (count < numseqs) {
      i = 0;
      int len = strlen(tmpseq->sequence);

      printf(">%s\n",tmpseq->id);
      while (i < len) {
        printf("%c",tmpseq->sequence[i]);

        if (i > 0 && i%72 == 0) {
          printf("\n");
        }
        i++;
      }

      printf("\n");
      tmpseq = tmpseq->next;
      count++;
    }

}


