/*
 *  random_corr.c
 *  pogview
 *
 *  Created by mclamp on 10/3/08.
 *  Copyright 2008 The Broad Institute. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "datamodel.h"

extern Sequence *read_fasta(FILE *file);

char *substring(char *seq,int start, int end);

int main(int argc, char *argv[]) {

  char *filename = "/Users/mclamp/cvs/pogview/data/chr1_120000000-120999999.fa";
  //char *filename = "/ahg/scr3/cvs/pogview/data/chr1_120000000-120999999.fa";

  int j = 0;
  
  while (j < 1000) {

    FILE     *file   = fopen(filename,"r");

    Sequence *seqptr = read_fasta(file);
    Sequence *seq    = seqptr;
   
    while (seq != NULL) {
      printf("Sequence %s %d\n",seq->id,j);
      seq = seq->next;
    }

    seq = seqptr;
      
    while (seq != NULL) {
	
      Sequence *tmpseq = seq->next;

      pog_free(seq->sequence,"Sequence2",NO_DEBUG);
      pog_free(seq->id,      "Id",       NO_DEBUG);
      pog_free(seq,          "Seq",      NO_DEBUG);

      seq = tmpseq;
	
    }
    printf("Done\n");
    j++;
  }
  return 1;
}


