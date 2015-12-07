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
    Sequence *seq1;
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

    int j = 1;

    seq  = seqptr;
    seq1 = seqptr;

    i = 0;
      
    int curpos = -1;

    while (i < seq1->length) {
      char c1 = seq1->sequence[i];
      
      if (c1 != '-') {

	curpos++;

	j = 0;
	
	seq = seqptr;

	while (j < numseqs) {
	  seq->sequence[curpos] = seq->sequence[i];
	  j++;
	  seq = seq->next;
	}
      }
      i++;
    }

    j = 0;

    seq = seqptr;

    curpos++;

    while (j < numseqs) {
      seq->sequence[curpos] = '\0';
      j++;
      seq = seq->next;
    }

    count = 0;
    
    seq = seqptr;

    while (count < numseqs) {
      printf (">%s\n",seq->id);
      
      i = 0;
      
      int len = strlen(seq->sequence);
      
      while (i < len) {
	printf("%c",seq->sequence[i]);
	
	if (i > 0 && i%72 == 0) {
	  printf("\n");
	}
	i++;
      }
      
      printf("\n");
    
      seq = seq->next;
      count++;
    }
}

