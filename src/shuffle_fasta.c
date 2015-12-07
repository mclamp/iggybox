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
    Sequence *newseqptr = NULL;
    Sequence *newseq;

    int r;
    int numseqs = 0;
    int i = 0;
    int count = 0;
    int numcols = 0;
    int num    = 0;
    int length;
    int realcol;

    char c;
    int  col;
    int  *done;

    seqptr = read_fasta(argv[1]);

    seq = seqptr;

    while (seq != NULL) {

      numseqs++;
      seq = seq->next;
    }

    length = seqptr->length;

    numcols = length;

    // Initialize the array to keep track of used columns

    done = (int *)malloc(length*sizeof(int));

    i = 0;
    
    while (i < numcols) {
      done[i] = 0;
      i++;
    }

    // Start the random number generator

    srand(8);
    
    newseqptr = (Sequence *)malloc(sizeof (Sequence));
	
    int j = 0;

    newseq = newseqptr;
    seq    = seqptr;

    while (j < numseqs) {

      newseq->name = (char *)malloc((1+strlen(seq->id))*sizeof(char));
      newseq->sequence = NULL;
      strcpy(newseq->name,seq->id);

      newseq->next = (Sequence *)malloc(sizeof (Sequence));
    
      newseq = newseq->next;
      seq    = seq->next;

      j++;
    }




    while (numcols > 0) {

      col = (int)((double)length*rand()*1.0/(1.0*RAND_MAX));


      while (done[col] == 1) {
	col = (int)((double)length*rand()*1.0/(1.0*RAND_MAX));
      }

      if (numcols % 10000 == 0) {
	//fprintf(stderr,"Found col %d %d\n",col,numcols);
      }

      done[col] = 1;

      // Now we take the chars at column col and add them onto the seuqences

      count = 0;


      newseq = newseqptr;
      seq    = seqptr;

      while (count < numseqs) {

	if (newseq->sequence == NULL) {

	  newseq->sequence = (char *)malloc(2 * sizeof(char));
	  
	  c = seq->sequence[col];

	  newseq->sequence[0] = c;
	  newseq->sequence[1] = '\0';
	  newseq->length      = 1;

	} else {

	  newseq->sequence = (char *)realloc(newseq->sequence,(newseq->length+1)*sizeof(char));
	  
	  c = seq->sequence[col];

	  newseq->sequence[newseq->length-1] = c;
	  newseq->sequence[newseq->length]   = '\0';
	  newseq->length++;
	}
	newseq = newseq->next;
	seq    = seq->next;

	count++;
      }

      numcols--;
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

