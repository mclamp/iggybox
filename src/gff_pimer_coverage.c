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
char **find_char_coverage(Sequence *seqptr,int numseqs);
int  *find_coverage(Sequence *seqptr);
int   find_numseqs(Sequence *seqptr);

int main(int argc, char *argv[]) {

  Genome *human = fetch_human_genome();
    
  int j   = 0;

  
  while (j < human->chrnum) {
    int chrlen = human->chr[j]->length;
    
    int i = 0;
    
    while (i < chrlen) {
      GenomePosition *pos   = (GenomePosition *)malloc(sizeof(GenomePosition));
      pos->chr = human->chr[j];
      pos->pos = i+1;
      int len = 1000000;
      
      if (i + len - 1 > chrlen) {
        len = chrlen - i;
      }
      
      Sequence       *seqptr  = read_fasta_test(pos,len);

      int k = 0;
      Sequence *seq = seqptr;
      while (seq != NULL) {
         k++;
	 seq = seq->next;
      }

      int numseqs = k;
      Pimer         **pis     = read_pimers_test(pos,len);
      char           **cov     = find_char_coverage(seqptr,numseqs);
      int            *fcov    = find_coverage(seqptr);

      len = seqptr->length;

      printf("Lenght %d\n",len);

      k = 0;
      while (k < len){
        if (pis[k] != NULL) {
         // printf("%s\tpimer\tpimer\t%llu\t%llu\t%f\t1\t%d\t%s\n",human->chr[j]->name,pos->pos+k,pos->pos+k,pis[k]->logodds,fcov[k],cov[k]);
        }
        k++;
      }
      free_pimers_test(pis,len);
      free_seqptr(seqptr);
      
      
      k = 0;
      
      while (k < len) {
        free(cov[k]);
        k++;
      }
      free(cov);
      free(fcov);
      free(pos);
      
      i += 1000000;
    }
    
    printf("Done\n");

    
    j++;
    
  }
  return 1;
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
char **find_char_coverage(Sequence *seqptr,int numseqs) {
  int numcols = seqptr->length;
  
  char **cov = (char **)calloc(numcols,sizeof(char *));
  
  int i = 0;
  
  while (i < numcols) {
    Sequence *seq    = seqptr;
    
    char *hmrd = (char *)malloc((numseqs+1) * sizeof(char));

    int j = 0;
    while (j < numseqs) {
      hmrd[j] = '-';
      j++;
    }
    hmrd[numseqs] = '\0';
    j = 0;
    while (seq != NULL) {
      if (seq->sequence[i] != '-') {
        hmrd[j] = seq->id[0];
      } 
      j++;
      seq = seq->next;
    }
    cov[i] = hmrd;
    
    i++;
  }
  
  return cov;
}


