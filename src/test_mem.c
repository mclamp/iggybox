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
#include <math.h>

#include "datamodel.h"

float find_correlation(float *vals, int start, int end, int offset);
int  *find_coverage(Sequence *seqptr);
int   find_numseqs(Sequence *seqptr);

int main(int argc, char *argv[]) {
  int totmem = 1000000000;
  int num = 10;
  int i = 0;

  while (i < num) {
  int j = 0;
  int chunk = 10000;

  int count = totmem/(sizeof(int)*chunk);

  int **arr = (int **)malloc(count*sizeof(int *));

  while (j < count) {
   
    arr[j] = (int *)malloc(chunk * sizeof(int));
    printf("Allocated %d\n",j);
    j++;
  }

  j = 0;

  while (j < 10000) {
     int arr = (int)(count*rand());
     printf("Found %d\n",j);
     j++;
  }

  j = 0;

  while (j < count) {
     free(arr[j]);
     j++;
  }
  free(arr);
  i++;
}

}
