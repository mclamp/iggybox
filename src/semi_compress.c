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
    int minnum    = atoi(argv[2]);
    int num = 0;

    char c;

    seqptr = read_fasta(argv[1]);

    seq = seqptr;

    while (seq != NULL) {
      numseqs++;
      seq = seq->next;
    }

    // Go through the alignment and pick out columns with > 3 species in
  

    while (i < seqptr->length) {

      num =  0;

      seq = seqptr;

      while (seq != NULL) {
	c = seq->sequence[i];

	if (c != '-') {
	  num++;
	}

	seq = seq->next;
      }

      if (num > minnum) {

	count = 0;
	
	if (newseqptr == NULL) {
	  newseqptr = (Sequence *)malloc(sizeof (Sequence));

	  newseq = newseqptr;
	  seq    = seqptr;

	  while (count < numseqs) {
	    newseq->name = (char *)malloc((1+strlen(seq->id))*sizeof(char));

	    strcpy(newseq->name,seq->id);
	    newseq->sequence = (char *)malloc(2 * sizeof(char));
	    
	    c = seq->sequence[i];

	    newseq->sequence[0] = c;
	    newseq->sequence[1] = '\0';
	    newseq->length      = 1;

	    count++;

	    if (count == numseqs) {
	      newseq->next = NULL;
	    } else {
	      newseq->next = (Sequence *)malloc(sizeof(Sequence));
	    }

	    newseq = newseq->next;
	    seq    = seq->next;

	  }
	  

	} else {

	  newseq = newseqptr;
	  seq    = seqptr;

	  while (count < numseqs) {
	    
	    newseq->sequence = (char *)realloc(newseq->sequence,(newseq->length+1)*sizeof(char));

	    c = seq->sequence[i];

	    newseq->sequence[newseq->length-1] = c;
	    newseq->sequence[newseq->length]   = '\0';
	    newseq->length++;

	    newseq = newseq->next;
	    seq    = seq->next;

	    count++;
	  }
	}
      }

      i++;
    }

    count = 0;
    
    newseq = newseqptr;

    while (count < numseqs) {

      printf (">%s\n",newseq->name);

      i = 0;

      while (i < newseq->length-1) {
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

