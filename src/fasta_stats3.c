#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "datamodel.h"

extern Sequence *read_fasta(char *filename);

int main(int argc, char * argv[]) {

    Sequence *seq;
    Sequence *seqptr;

    int i       = 0;
    int numseqs = 0;

    seqptr = read_fasta(argv[1]);
    seq    = seqptr;

    while (seq != NULL) {

      numseqs++;
      
      long len = seq->length;
      long i   = 0;
      int nongap = 0;
      int repbases = 0;

      while (i < seq->length) {
	if (seq->sequence[i] != '-' &&
	    seq->sequence[i] != 'N') {
	  nongap++;
	}
	if (seq->sequence[i] == 'a' ||
	    seq->sequence[i] == 't' ||
	    seq->sequence[i] == 'c' ||
	    seq->sequence[i] == 'g') {
	    repbases++;
	}
	i++;
      }
      
      printf("Sequence\t%s\t%s\tLength\t%d\tNongap\t%d\tPercent_nongap\t%d\tRepeat\t%d\tPercent_repeat\t%d\n",argv[1],seq->id,seq->length,nongap,(int)(100*nongap/seq->length),repbases,(int)(100*repbases/seq->length));

      seq = seq->next;
      
    }

    printf("Total seqs %d\n",numseqs);
}


Sequence *next_seq(FILE *file) {
  char c;

  Sequence *seq;
  char     *id;
  char     *str;

  seq = NULL;
  id  = NULL;
  str = NULL;

  long idchunk  = 80;
  long idlen    = 80;
  long seqchunk = 1000;
  long seqlen   = 1000;

  while ((c = fgetc(file)) != EOF) {
    if (c == '>') {

      if (id != NULL) {

	seq = (Sequence *)malloc(sizeof(Sequence));

	id  = (char *)realloc(id,(strlen(id)+1)*sizeof(char));
	str = (char *)realloc(str,(strlen(str)+1)*sizeof(char));

	seq->id       = id;
	seq->sequence = str;
	seq->length   = strlen(str);

	ungetc('>',file);
	
	return seq;
      } else {

	while ((c = fgetc(file)) != EOF && c != '\n') {
	  
	  if (id == NULL) {
	    id = (char *)malloc(idlen * sizeof(char));

	    id[0] = c;
	    id[1] = '\0';
	  } else {
	    long len = strlen(id);

	    if (len > idlen+1) {
	      idlen = idlen + idchunk;

	      id = (char *)realloc(id,(idlen)*sizeof(char));
	    }    
	    id[len] = c;
	    id[len+1] = '\0';
	    
	  }
	}
      }
    } else {
      if (c != '\n') {

	if (str == NULL) {
	  str = (char *)malloc(seqlen*sizeof(char));

	  str[0] = c;
	  str[1] = '\0';

	} else {
	  long len = strlen(str);
	  if (len > seqlen-1) {
	    seqlen = seqlen + seqchunk;
	    str = (char *)realloc(str,(seqlen)*sizeof(char));
	  } 
	  str[len] = c;
	  str[len+1] = '\0';
	}
      }
    }
  }
  if (id != NULL) {
    
    seq = (Sequence *)malloc(sizeof(Sequence));


    id  = (char *)realloc(id,(strlen(id)+1)*sizeof(char));
    str = (char *)realloc(str,(strlen(str)+1)*sizeof(char));
    
    seq->id = id;
    seq->sequence = str;
  }

  return seq;
}
