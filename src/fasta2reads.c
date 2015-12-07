#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "datamodel.h"

extern Sequence *read_sequence(FILE *file);
extern void      free_seqptr(Sequence *seq);

int main(int argc, char *argv[]) {
  
  char *infilename = argv[1];


  FILE *infile = fopen(infilename,"r");

  if (infile == NULL) {
    printf("No file [%s] found\n",infilename);
    exit(0);
  }

  
  Sequence *seq;

  int readlen = 35;

  while ((seq = read_sequence(infile)) != NULL) {

    int len = seq->length;

    int i = 0;

    while (i < len-readlen) {
      char *newstr = substring(seq->sequence,i,i+readlen-1);

      printf(">%s.%d-%d\n%s\n",seq->id,i,i+readlen-1,newstr);

      free(newstr);

      i += readlen;
    }
    free_seqptr(seq);

  }

  return 1;
}


