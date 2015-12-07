#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "datamodel.h"

extern Sequence *read_fasta(char *filename);

int main(int argc, char * argv[]) {

    Sequence *seq;
    Sequence *seqptr;
    Sequence *mouse;
    Sequence *dog;

    int i = 0;
    int numseqs = 0;

    seqptr = read_fasta(argv[1]);

    seq = seqptr;
    
    int cov[200];
    int md_cov[200];

    int len = seqptr->length;

    Sequence *prev;

    while (seq != NULL) {
      if (strcmp(seq->id,"mouse") == 0) {
	mouse = seq;
      } else if (strcmp(seq->id,"dog") == 0) {
	dog = seq;
      } else {
	prev->next = seq;
	numseqs++;
	prev = seq;
      }
      seq = seq->next;
    }

    while (i < 200) {
      cov[i] = 0;
      md_cov[i] = 0;
      i++;
    }

    i = 0;

    int count = 0;
    int md_count = 0;
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

      if (mouse->sequence[i] != '-' &&
	  dog->sequence[i]   != '-') {
	md_cov[nongap]++;
	md_count++;
      } else {
	cov[nongap]++;
	count++;
      }
      
      i++;
    }

    i = 0;

    while (i < numseqs) {
      printf("%d\t%d\t%d\t%d\t%d\n",i,cov[i],md_cov[i],count,md_count);
      i++;
    }
}

