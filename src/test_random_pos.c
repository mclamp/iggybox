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


int main(int argc, char *argv[]) {

  int j = 0;
  
  Genome *human = fetch_human_genome();
  printf("Size of int %lu\n",sizeof(unsigned long long int));


  
  while (j < 1000000) {

    GenomePosition *pos = get_random_genome_position(human,1000);
    
    printf("position %d %s %llu\n",j,pos->chr->name,pos->pos);
    
    free_random_genome_position(pos);

    j++;
  }
  
  return 0;
}


