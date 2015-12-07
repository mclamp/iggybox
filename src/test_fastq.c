#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "datamodel.h"

extern Sequence *read_fastq_sequence(FILE *file);

char *substring(char *seq,int start, int end);

int main(int argc, char *argv[]) {

  char *filename = argv[1];
  int   samplefreq = atoi(argv[2]);

  //printf("File and freq %s %d\n",filename,samplefreq);

  FILE     *file   = fopen(filename,"r");

  Sequence *seq = read_fastq_sequence(file);
  int       count = 0;

  while (seq != NULL) {
    if (count % samplefreq == 0) {
      printf("@%s\n",seq->id);
      printf("%s\n",seq->sequence);
      printf("%s\n",seq->comment);
      printf("%s\n",seq->qual);
    }
    seq = read_fastq_sequence(file);

    count++;
  }
  return 1;
}


