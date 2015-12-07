#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>

#include "datamodel.h"
#include "pog_utils.h"
#include "seqio.h"

char * read_fastq_line(FILE *file);

Sequence *read_fastq_sequence(FILE *file) {
  char       c;
  
  int idlen  = 0;
  int seqlen = 0;
  
  int idchunk  = 20;
  int seqchunk = 10000;
  
  int curr_idchunk  = idchunk;
  int curr_seqchunk = seqchunk;
  
  
  Sequence *seq;

  seq = (Sequence *)malloc(sizeof(Sequence));
  
  if (seq == NULL) {
    fprintf(stderr,"Can't allocate memory for sequence\n");
    exit(0);
  }

  seq->id       = NULL;
  seq->sequence = NULL;
  seq->next     = NULL;
  seq->qual     = NULL;
  seq->comment  = NULL;

  while ((c = fgetc(file)) != EOF) {
    
    // Beginning of read sequence

    if (c == '@') {
      
      seq->id = read_fastq_line(file);
      
      if (seq->id == NULL) {
        fprintf(stderr,"Can't allocate memory for sequence id\n");
        exit(0);
      }
      
      seq->sequence = read_fastq_line(file);

      if (seq->sequence == NULL) {
        fprintf(stderr,"Can't allocate memory for sequence string\n");
        exit(0);
      }

      seq->comment = read_fastq_line(file);

      if (seq->comment == NULL) {
        fprintf(stderr,"Can't allocate memory for sequence comment\n");
        exit(0);
      }

      seq->qual    = read_fastq_line(file);

      if (seq->qual == NULL) {
        fprintf(stderr,"Can't allocate memory for quality string\n");
        exit(0);
      }
      
      seq->length = strlen(seq->sequence);
     
      return seq;
    }
  }
  return seq;
}

char * read_fastq_line(FILE *file) {

  char *str = NULL;
  char  c;
  int   n = 0;

  while ((c = getc(file)) != EOF && c != '\n') {

    if (str == NULL) {

      str = (char *)malloc(sizeof(char)*(2));
      n = 0;
      str[n] = c;
      str[n+1]   = '\0';

    } else {

      str = (char *)realloc(str,sizeof(char)*(strlen(str)+2));

      str[n] = c;
      str[n+1]   = '\0';

    }

    n++;

  }

  return str;
}


