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

extern Pimer *read_pimer(FILE *file);

char *substring(char *seq,int start, int end);

int main(int argc, char *argv[]) {

  char *filename = "/Users/mclamp/cvs/pogview/data/chr1_120000000-120999999.fa.pi.12mer";

  int j = 0;


  
  while (j < 1000) {
    Pimer **pis = (Pimer **)malloc(10000*sizeof(Pimer *));
    
    Pimer *pi = NULL;
    
    FILE     *file   = fopen(filename,"r");
    
    while ((pi = read_pimer(file)) != NULL) {
      if (pi->start > 400000 && pi->start <= 410000) {
        pis[pi->start - 400000] = pi;
      } else {
        free(pi->chr);
        free(pi);
      }

    }

    int i = 0;
    
    while (i < 10000) {
      if (pis[i] != NULL) {
        if (i < 100) {
          printf("Pimer %d %s %d %f\n",j,pis[i]->chr,pis[i]->start,pis[i]->logodds);
        }
        free(pis[i]->chr);
        free(pis[i]);
      }
      i++;
    }
    free(pis);
    j++;
    
    fclose(file);
  }
  return 1;
}


