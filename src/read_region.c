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

char *orgs[] = {"Human",
                "Chimp",
                "RheMac",
                "MouseLemur",
                "Bushbaby",
                "TreeShrew",
                "Mouse7",
                "Rat",
                "Cavia",
                "Squirrel",
                "Rabbit",
                "Pika",
                "Cow",
                "Cat",
                "Dog2",
                "Hedgehog",
                "Shrew",
                "Elephant",
                "Tenrec",
                "Armadillo"};

int main(int argc, char *argv[]) {

  char *chrname = argv[1];
  int   start = atoi(argv[2]);
  int   end   = atoi(argv[3]);

  Genome         *human = fetch_human_genome();
  Chromosome     *chr   = find_chromosome_by_name(human,chrname);
  GenomePosition *pos   = (GenomePosition *)malloc(sizeof(GenomePosition));

  pos->chr = chr;
  pos->pos = start;
	
  Sequence *seqptr = read_fasta_test(pos,end-start+1);

  int i = 0;

  while (i < 20) {

    Sequence *tmpseq = seqptr;
    
    while (tmpseq != NULL) {

      if (strcmp(tmpseq->id,orgs[i]) == 0) {
	printf(">%s\n",tmpseq->id);
	printf("%s\n",tmpseq->sequence);
      }

      tmpseq = tmpseq->next;

    }
    i++;
  }
  return 1;
}


