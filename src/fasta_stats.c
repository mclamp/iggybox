#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "datamodel.h"

//extern Sequence *read_fasta(char *filename);
Sequence *next_seq(FILE *file);

int main(int argc, char * argv[]) {

    Sequence *seq;

    int i       = 0;
    int numseqs = 0;

    FILE *file;

    file = fopen(argv[1],"r");

    while ((seq = next_seq(file)) != NULL) {

      numseqs++;
      
      int len = seq->length;
      int i   = 0;
      int nongap = 0;
      
      while (i < seq->length) {
	if (seq->sequence[i] != '-' &&
	    seq->sequence[i] != 'N') {
	  nongap++;
	}
	i++;
      }
      
      printf("Sequence\t%s\t%s\t%d\t%d\t%d\n",argv[1],seq->id,seq->length,nongap,(int)(100*nongap/seq->length));

      free(seq->sequence);
      free(seq->id);
      free(seq);
      
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

  int idchunk  = 80;
  int idlen    = 80;
  int seqchunk = 1000;
  int seqlen   = 1000;

  while ((c = fgetc(file)) != EOF) {
    if (c == '>') {

      if (id != NULL) {

	seq = (Sequence *)malloc(sizeof(Sequence));

	id  = (char *)realloc((strlen(id)+1)*sizeof(char));
	seq = (char *)realloc((strlen(seq)+1)*sizeof(char));

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
	    int len = strlen(id);

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
	  int len = strlen(str);
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


    id  = (char *)realloc((strlen(id)+1)*sizeof(char));
    seq = (char *)realloc((strlen(seq)+1)*sizeof(char));
    
    seq->id = id;
    seq->sequence = str;
    return seq;
  }

  return seq;
}
