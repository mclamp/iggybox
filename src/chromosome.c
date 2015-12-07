/*
 *  chromosome.c
 *  pogview
 *
 *  Created by mclamp on 10/3/08.
 *  Copyright 2008 The Broad Institute. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>


#include "datamodel.h"

char *human_chrnames[] = {"chr1",
"chr2",
"chr3",
"chr4",
"chr5",
"chr6",
"chr7",
"chr8",
"chr9",
"chr10",
"chr11",
"chr12",
"chr13",
"chr14",
"chr15",
"chr16",
"chr17",
"chr18",
"chr19",
"chr20",
"chr21",
"chr22",
"chrX"};

unsigned long long int hg17_chrlen[] = {
245522847,
243018229,
199505740,
191411218,
180857866,
170975699,
158628139,
146274826,
138429268,
135413628,
134452384,
132449811,
114142980,
106368585,
100338915,
88827254,
78774742,
76117153,
63811651,
62435964,
46944323,
49554710,
154824264};
//57701691


Genome *fetch_human_genome();
GenomePosition *get_random_genome_position(Genome *genome,int offset);

Genome *fetch_human_genome() {
  Genome *human = (Genome *)malloc(sizeof(Genome));
	
  human->chr = (Chromosome **)malloc(HUMAN_CHRNUM * sizeof(Chromosome *));
	
	int i  = 0;
  
	human->species = "Human";
	human->version = "hg17";
	human->chrnum  = HUMAN_CHRNUM;
	
	unsigned long long int len = 0;
	
	while (i < human->chrnum) {
	  Chromosome *chr = (Chromosome *)malloc(sizeof(Chromosome));
		
		chr->name   = human_chrnames[i];
		chr->length = hg17_chrlen[i];
    
		human->chr[i] = chr;
		len += chr->length;
		
		i++;
	}
  
	human->length = len;
  
	return human;
}

Chromosome *find_chromosome_by_name(Genome *genome, char *name) {
  int i = 0;

  while (i < genome->chrnum) {
    if (strcmp(genome->chr[i]->name,name) == 0) {
      return genome->chr[i];
    }
    i++;
  }
  return NULL;
}

void free_random_genome_position(GenomePosition *pos) {
  free(pos);
}

GenomePosition *get_random_genome_position(Genome *genome,int offset) {
	
	unsigned long long int tmppos = 0;
  
  int found = 0;
  
  while (found == 0) {
    // Find a random position
    double r = rand();
		
    unsigned long long int pos = (unsigned long long int)((2*r)*genome->length/(2.0*RAND_MAX));
		
    int i = 0;
		
    unsigned long long int cumul_len  = 0;
    unsigned long long int prev_len   = 0;
		
    Chromosome *tmpchr = NULL;
    
    while (i < genome->chrnum) {
      
      cumul_len += genome->chr[i]->length;
      
      if (cumul_len >= pos) {
        tmpchr = genome->chr[i];
        tmppos = prev_len;
        i = genome->chrnum;
      }	else {			
        prev_len = cumul_len;
      }
      i++;
    }
    
    //printf("\nChromosome %s len %llu pos %llu\n\n",tmpchr->name,tmpchr->length,pos-tmppos);
    
    if (tmpchr->length- pos + tmppos > offset) {
      found = 1;
      
      GenomePosition *gpos = (GenomePosition *)malloc(sizeof(GenomePosition));
      
      gpos->pos = (pos-tmppos);
    //
      
      gpos->chr = tmpchr;
      
      return gpos;
    }
  }
  return NULL;
}
