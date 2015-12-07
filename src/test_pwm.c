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

extern Pwm*read_pwm(FILE *file);

char *substring(char *seq,int start, int end);

int main(int argc, char *argv[]) {

  char *filename = "/Users/mclamp/cvs/pogvue/data/tf9.4.matrix.pwmline";
  //char *filename = "/ahg/scr3/cvs/pogview/data/chr1_120000000-120999999.fa";

  int k = 0;

  while (k < 1000) {
    
    FILE     *file   = fopen(filename,"r");
    Pwm *pwm;
    
    int j = 0;
    
    while ((pwm = read_pwm(file)) != NULL) {
      j++;

      printf("Pwm is %s %f\n",pwm->name, pwm->vals[pwm->len-1]);
      free_pwm(pwm);
    }
    printf("Found %d %d\n",j,k);
    
    fclose(file);
    k++;
  }
  return 1;
}


