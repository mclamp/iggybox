/*
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

extern Pi    *read_pi(FILE *file);
extern void   print_pi(Pi *pi);

char *substring(char *seq,int start, int end);

int main(int argc, char *argv[]) {

  //char *filename = "/Users/mclamp/cvs/pogview/data/chr1_120000000-120999999.fa.pi";
  char *filename = "/ahg/scr3/cvs/pogview/data/chr1_120000000-120999999.fa.pi";

  int j = 0;

  while (j < 100000) {
    Pi     *pi = NULL;
    
    FILE     *file   = fopen(filename,"r");
    
    if (file == NULL) {
      printf("ERROR: Can't open file %s\n",filename);
      exit(0);
    }
    
    int count   = 0;
    int pichunk = 1000;

    Pi **pis = (Pi **)malloc(pichunk * sizeof(Pi *));
    
    int prev = -1;
    
    while ((pi = read_pi(file)) != NULL) {
      count++;
    
      // Increase pis if necessary
      int inc_chunk = 0;
      
      while (pi->start > pichunk) {
        pichunk += pichunk;
        inc_chunk = 1;
      }
      if (inc_chunk == 1) {
        pis = (Pi **)realloc(pis,pichunk*sizeof(Pi *));
      }
      
      // Fill in gaps with nulls
      while (pi->start - prev > 0) {
        prev++;
        pis[prev] = NULL;
      }
      
      // Add pi into array
      pis[pi->start] = pi;
    }

    // Fill in end of array
    while (prev < pichunk-1) {
      prev++;
      pis[prev] = NULL;
    }
    
    // Free memory
    int i = 0;
    
    while (i < pichunk) {

      if (pis[i] != NULL) {
        free(pis[i]->chr);
        free(pis[i]);
      } else {
        free(pis[i]);
      }
      i++;
    }
    
    //
    free(pis);
    
    printf("Count is %d\n",count);

    fclose(file);
    j++;
    
  }
  return 1;
}


