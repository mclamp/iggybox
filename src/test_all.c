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
  int len = 100;

  if (argc > 1) {
     len = atoi(argv[1]);
  }
  Genome *human = fetch_human_genome();
    
  int j   = 0;
  
  while (j < 10000) {
    
    GenomePosition *pos     = get_random_genome_position(human,2*len);
    
    int offset1 = (int)(pos->pos/1000000);
    int offset2 = (int)((pos->pos+2*len)/1000000);
    
    while (offset1 != offset2) {
      free_random_genome_position(pos);
      pos     = get_random_genome_position(human,2*len);

      offset1 = (int)(pos->pos/1000000);
      offset2 = (int)((pos->pos+2*len)/1000000);
    }
    Pimer         **pis     = read_pimers_test(pos,2*len);    
    Sequence       *seqptr  = read_fasta_test(pos,2*len);
    int            *cov     = find_coverage(seqptr);
      
      
    printf("position %s %llu\n",pos->chr->name,pos->pos);
    
    // Pimers
    
    if (pis != NULL) {
      
      float *vals = (float *)calloc(2*len, sizeof(float));
      
      int i = 0;
      
      while (i < 2*len) {
        if (pis[i] != NULL) {
          vals[i] = pis[i]->logodds;
        }           
        i++;
      }
        
      int offset = 1;
        
      while (offset < 50) {
        
        float corr = find_correlation(vals,1,len,offset);
        
        printf("Corr %s\tcoord\t%llu\toffset\t%d\tcorr\t%f\tcov\t%d\n",pos->chr->name,pos->pos,offset,corr,cov[0]);
        offset++;
      }
      free(vals);
    }
    
    
    // Fasta 
    
    Sequence *seq    = seqptr;
      
    while (seq != NULL) {
      //printf("Sequence %20s %5d\n",seq->id,j);
      seq = seq->next;
    }


    // Coverage
    
    int i = 0;
    
    while (i < len) {
        //printf("Cov %d %d\n",i,cov[i]);
      i++;
    }

    // Now free
    
    free_random_genome_position(pos);
    free_pimers_test(pis,len);
    free_seqptr(seqptr);
    free(cov);
    
    printf("Done\n");

    
    j++;
    
  }
  return 1;
}

int find_numseqs(Sequence *seqptr) {
  Sequence *seq = seqptr;
  
  int numseqs = 0;
  
  while (seq != NULL) {
    numseqs++;
    seq = seq->next;
  }
  return numseqs;
}

float find_correlation(float *vals, int start, int end, int offset) {
	
	float sum_sq_x = 0.0;
	float sum_sq_y = 0.0;
	float sum_coproduct = 0.0;
	float mean_x = vals[start];
	float mean_y = vals[start+offset];
  
  
	int len = (end-start+1);
  
	int i = 1;
  
	while (i < len) {
    
		float sweep = (i - 1.0) / i;
		float delta_x = vals[i+start] - mean_x;
		float delta_y = vals[i+start+offset] - mean_y;
		sum_sq_x += delta_x * delta_x * sweep;
		sum_sq_y += delta_y * delta_y * sweep;
		sum_coproduct += delta_x * delta_y * sweep;
		mean_x += delta_x / i;
		mean_y += delta_y / i;
		i++;
	}
  
	float pop_sd_x = sqrt( sum_sq_x / len);
	float pop_sd_y = sqrt( sum_sq_y / len );
	float cov_x_y  = sum_coproduct / len;
  
	if (pop_sd_x > 0 && pop_sd_y > 0) {
		return cov_x_y / (pop_sd_x * pop_sd_y);
	} else {
		return 0.0;
	}
}

int *find_coverage(Sequence *seqptr) {
  int numcols = seqptr->length;
  
  int *cov;
  
  cov = (int *)calloc(numcols,sizeof(int));
  
  int i = 0;
  
  while (i < numcols) {
    Sequence *seq    = seqptr;
    
    int tmp   = 0;
    int count = 0;
      
    while (seq != NULL) {
      if (seq->sequence[i] != '-') {
        tmp++;
      }
      count++;
      seq = seq->next;
    }
    cov[i] = tmp;

    i++;
  }
  
  return cov;
}


