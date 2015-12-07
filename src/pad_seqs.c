#include "datamodel.h"
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>

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
    int maxlen = 0;

    while (seq != NULL) {
      numseqs++;

      if (strlen(seq->sequence) > maxlen) {
	maxlen = strlen(seq->sequence);
      }
      seq = seq->next;

    }

    i = 0;
    seq = seqptr;

    while (i < numseqs) {
      
      int len = strlen(seq->sequence);

      if (len < maxlen) {

	seq->sequence = (char *)realloc(seq->sequence,(maxlen+1)*sizeof(char));

	int j = strlen(seq->sequence);

	while (j < maxlen) {
	  seq->sequence[j] = '-';
	  j++;
	}
	seq->sequence[maxlen] = '\0';
	
      }
      seq = seq->next;
      i++;
    }

    count = 0;

    newseq = seqptr;

    while (count < numseqs) {

      printf (">%s\n",newseq->id);
      
      i = 0;
      
      int len = strlen(newseq->sequence);
      
      while (i < len) {
        printf("%c",newseq->sequence[i]);
	
        if (i > 0 && i%72 == 0) {
          printf("\n");
        }
        i++;
      }
      
      printf("\n");

      newseq = newseq->next;
      count++;
    }
}
