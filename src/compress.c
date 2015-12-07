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

    // The format is
    // >human1_gapped
    // >species1_gapped
    // >human1_gapped
    // >species2_gapped


    int j = 0;

    seq = seqptr;

    while (j < numseqs) {

      Sequence *seq1;
      Sequence *seq2;

      // human seq
      seq1 = seq;
      seq = seq->next;

      // species seq
      seq2 = seq;

      i = 0;
      
      char *str1;
      char *str2;

      str1 = NULL;
      str2 = NULL;

      int len1 = 0;
      int len2 = 0;

      int strlen1 = 0;
      int strlen2 = 0;

      while (i < seq1->length) {
	char c1 = seq1->sequence[i];
	char c2 = seq2->sequence[i];

	if (c1 != '-') {

	  if (str1 == NULL) {
	    len1 = 1000;
	    strlen1 = 1;

	    str1 = (char *)malloc((len1+1)*sizeof(char));

	    str1[0] = c1;
	    str1[1] = '\0';

	  } else {
	    int len = strlen1;

	    if (strlen1 > len1) {
	      
	      len1 += 1000;
	      str1 = (char *)realloc(str1,(1 + len1)*sizeof(char));
	    }
	    str1[len] = c1;
	    str1[len+1] = '\0';
	    
	    strlen1++;

	  }

	  if (str2 == NULL) {
	    len2 = 1000;
	    strlen2 = 1;
	    str2 = (char *)malloc((len2 + 1)*sizeof(char));

	    str2[0] = c2;
	    str2[1] = '\0';

	  } else {
	    int len = strlen2;
	    if (strlen2 > len2) {

	      len2 += 1000;
	      str2 = (char *)realloc(str2,(1+ len2)*sizeof(char));
	    }
	    str2[len] = c2;
	    str2[len+1] = '\0';
	    strlen2++;
	  }
	}
	i++;
      }

      free(seq1->sequence);
      free(seq2->sequence);

      seq1->sequence = str1;
      seq2->sequence = str2;
     
      seq = seq->next;

      j += 2;
    }

    count = 0;
    
    newseq = seqptr;

    while (count < numseqs) {
      if (count%2 == 1 || count == 0){
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
     } 
      newseq = newseq->next;
      count++;
    }
    
}
