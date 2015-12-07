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
char **find_char_coverage(Sequence *seqptr,int numseqs);
char *find_mercov(Sequence *seqptr,int numseqs,char **cov,int *fcov,Pimer *pi,int offset1,int offset2);


int main(int argc, char *argv[]) {
  
  Genome *human = fetch_human_genome();
  
  char *gfffilename = argv[1];//"/Users/mclamp/cvs/pogview/data/exons.gff";
  
  FILE *gfffile     = fopen(gfffilename,"r");
  
  if (gfffile == NULL) {
    printf("Can't open file %s\n",gfffilename);
    exit(1);
  }
  
  Gff    *gff   = read_gff(gfffile);
  
  Gff    *tmpgff = gff;
  
  while (gff != NULL) {
    GenomePosition *pos     = (GenomePosition *)malloc(sizeof(GenomePosition));
    Chromosome *chr         = find_chromosome_by_name(human, tmpgff->seqname);
    
    pos->pos = tmpgff->start;
    pos->chr = chr;
    
    int len = tmpgff->end - tmpgff->start + 1;
    
    if (chr == NULL) {
      printf("Can't find chromosome %s in genome %s\n",tmpgff->seqname,human->species);
    } else {
      
      int offset1 = (int)(pos->pos/1000000);
      int offset2 = (int)((pos->pos+len)/1000000);
      
      if (offset1 == offset2 && len > 50) {
        Pimer            **pis     = read_pimers_test(pos,len);
        Sequence       *seqptr     = read_fasta_test(pos,len);
        int            *fcov       = find_coverage(seqptr);

	if (pis != NULL && seqptr != NULL && fcov != NULL) {
	  int numseqs = 0;
	  
	  Sequence *seq = seqptr;
	  
	  while (seq != NULL) {
	    numseqs++;
	    printf("Id is %s\n",seq->id);
	    seq = seq->next;
	  }
	  
	  char          **cov     = find_char_coverage(seqptr,numseqs);
	  
	  print_gff(stdout,tmpgff);
	  
	  int i = 0;
	  
	  while (i < len - 50) {
	    if (pis[i] != NULL) {
	      char *mercov = find_mercov(seqptr,numseqs,cov,fcov,pis[i],tmpgff->start,offset1);
	      
	      printf("%s\tpimer\tpimer\t%d\t%d\t%f\t1\t.\t%d\t%s\t%s\n",tmpgff->seqname,tmpgff->start+i,tmpgff->end+i,pis[i]->logodds,fcov[i],cov[i],mercov);   
	      
	      free(mercov);
	    }
	    i++;
	  }
	  
        
	  int k = 0;
	  
	  while (k < len) {
	    free(cov[k]);
	    k++;
	  }
        free(cov);
	}
        free(fcov);
        free(seqptr);
        free_pimers_test(pis,len);

      }
      free(pos);
    }
    tmpgff = tmpgff->next;
  }
  
  free_gffptr(gff);
  
  return(1);
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

char *find_mercov(Sequence *seqptr,int numseqs,char **cov,int *fcov,Pimer *pi,int offset1,int offset2) {
  
  char *mercov = (char *)malloc((numseqs+1)*sizeof(char));
  
  int start = pi->start;
  int end   = pi->end;
  
  int i = 0;
  
  while (i < numseqs) {
    mercov[i] = '-';
    i++;
  }
  mercov[numseqs] = '\0';
  
  i = 0;
  
  Sequence *seq = seqptr;  
  int numcov = 0;
  
  while (i < numseqs) {
    int j = start;
    int tot = 0;
    
    while (j < end) {
      // Check if the char in cov at position j,i is "-"
      if (cov[j-offset1+1000000*offset2][i] != '-') {
        tot++;
      }
      j++;
    }
    if (tot > 0.8*50) {
      mercov[i] = seq->id[0];
      numcov++;
    }
    seq = seq->next;
    i++;
  }
  
  fcov[start-offset1+1000000*offset2] = numcov;
  return mercov;
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


