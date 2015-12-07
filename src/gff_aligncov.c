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


OrgHash *find_orgcov(Sequence *seqptr,int numseqs, int start, int end);
void free_orghash(OrgHash *hash);

Sequence *reorder_seqs(Sequence *seq);
int main(int argc, char *argv[]) {

  Genome *human = fetch_human_genome();
    
  Gff *gff;

  int bin;
  int prevbin=-1;
  char *prevchr = (char *)malloc(20*sizeof(char));
  char *gfffilename = argv[1];

  FILE *gfffile = fopen(gfffilename,"r");
  Sequence *seqptr = NULL;

  while ((gff = read_gff(gfffile)) != NULL) {  

    bin = (int)(gff->start/1000000);

    Chromosome     *chr   = find_chromosome_by_name(human,gff->seqname);
    GenomePosition *pos   = (GenomePosition *)malloc(sizeof(GenomePosition));
      
    if (chr != NULL && (int)(gff->end/1000000) == bin) {

      pos->chr = chr;
      pos->pos = 1000000*bin;
	
      if (seqptr  == NULL ||
	  prevbin != bin  ||
	  strcmp(gff->seqname,prevchr)!=0) {
	  
	if (seqptr != NULL) {
	  printf("New seq chunk %s %s %d %d %d\n",gff->seqname,prevchr,bin,prevbin,seqptr->length);
	  free_seqptr(seqptr);	
	} 
   
	if (pos->pos+1000000 > chr->length) {
	  seqptr  = read_fasta_test(pos,chr->length-pos->pos);	
	} else {
	  seqptr  = read_fasta_test(pos,1000000);	

	}
	seqptr = reorder_seqs(seqptr);

      }
      printf("read %d\n",seqptr->length);
      int k = 0;
      Sequence *seq = seqptr;
      while (seq != NULL) {
         k++;
        seq = seq->next;
      }

      int numseqs = k;

      OrgHash *cov = find_orgcov(seqptr,numseqs,gff->start-bin*1000000,gff->end-bin*1000000);

      k = 0;
      double avcov = 0;

      int len = gff->end-gff->start+1;
      while (k < numseqs) {
	avcov += cov->vals1[k];
	k++;
      }
      avcov /= numseqs;
      k = 0;

      while (k < numseqs) {
	
	printf("%20s\t%10s\t%10d\t%10d\t%5d\t%7d\t%7d\t%7d\t%7d\t%d\t%5.2f\n",cov->orgs[k],gff->seqname,gff->start,gff->end,len,cov->vals1[k],cov->vals2[k],k,cov->mdcount,cov->no_md,avcov);
	k++;

      }
      k = 0;
      while (k < len) {
	printf("COV\t%10d\t%10s\t%10d\t%7d\t%7d\t%s\n",k,gff->seqname,gff->start+k,cov->vals3[k],cov->xcov[k],cov->mercov[k]);
	k++;
      }
      free_orghash(cov);
      free(pos);
      
    }
    prevbin = bin;
    strcpy(prevchr,gff->seqname);
      
  }
}


void free_orghash(OrgHash *hash) {
  int i = 0;

  while (i < hash->count) {
    free(hash->orgs[i]);
    i++;
  }
  free(hash->vals1);
  free(hash->vals2);
  free(hash);
}
       
Sequence *reorder_seqs(Sequence *seq) {

  
  int i = 0;
  Sequence *tmpseq = seq;

  while (tmpseq != NULL) {
    tmpseq = tmpseq->next;
    i++;
  }

  
  Sequence *currseq = NULL;

  Sequence **newseq = (Sequence **)malloc(i * sizeof(Sequence *));
  

  int j = 0;

  int jj = 0;
  char *gapseq = (char *)malloc((seq->length+1)*sizeof(char));
  
  while (jj < seq->length) {
    gapseq[jj] = '-';
    jj++;
  }
  gapseq[jj] = '\0';
  
  while (j < ORGS) {
    printf("Looking for %d %s\n",j,orgs[j]);
    char *org = orgs[j];

    Sequence *found = NULL;
    tmpseq = seq;

    while (tmpseq != NULL) {
      if (strcmp(tmpseq->id,org) == 0) {
	found = tmpseq;
	tmpseq = NULL;
      } else {
	tmpseq = tmpseq->next;
      }
    }

    if (found != NULL) {
      printf("Found %s\n",found->id);
      newseq[j] = found;
    } else {
      printf("Not found = gapseq\n");
      newseq[j]->id = orgs[j];
      newseq[j]->length = seq->length;
      newseq[j]->sequence = gapseq;
    }
    j++;
  }

  seq = newseq[0];

  tmpseq = seq;

  j = 0;

  while (j < ORGS-1) {
    tmpseq->next = newseq[j+1];
    tmpseq = tmpseq->next;
    j++;
  }

  tmpseq->next = NULL;
    
  tmpseq = seq;

  while (tmpseq != NULL) {
    printf("Seq %s\n",tmpseq->id);
    tmpseq = tmpseq->next;
  }
    
  free(newseq);
  return seq;
}
      



OrgHash *find_orgcov(Sequence *seqptr, int numseqs, int start, int end) {
  OrgHash *cov = (OrgHash *)malloc(sizeof(OrgHash));

  int numcols = end-start+1;
  cov->count = numseqs;
  
  cov->orgs  = (char **)malloc(numseqs*sizeof(char *));
  cov->vals1  = (int *)calloc(numseqs,sizeof(int));
  cov->vals2  = (int *)calloc(numseqs,sizeof(int));
  cov->vals3  = (int *)calloc(numcols,sizeof(int));
  cov->mercov = (char **)malloc(numcols*sizeof(char *));
  cov->xcov  = (int *)calloc(numcols,sizeof(int));

  int i = 0;
  
  
  i=0;

  Sequence *seq = seqptr;

  printf("Fetching %d %d\n",start,end);

  while (i < numseqs) {
    cov->orgs[i] = (char *)malloc((strlen(seq->id)+1)*sizeof(char));
    strcpy(cov->orgs[i],seq->id);

    seq = seq->next;

    i++;
  }
  cov->mdcount = 0;
  cov->cov  = 0;
  cov->no_md = 0;
  int j = start;

  while (j <= end) {
    
    i = 0;
    seq    = seqptr;
    
    Sequence *topseq = seqptr;

    int mdcount = 0;
    int mandcount = 0;

    char *mercov = (char *)malloc((numseqs+1)*sizeof(char));

    cov->mercov[j-start] = mercov;
    
    while (i < numseqs) {
      mercov[i] = '-';
      i++;
    }
    
    mercov[numseqs] = '\0';

    i=0;

    while (i < numseqs) {
      
      if (seq->sequence[j] != '-') {
	cov->vals1[i]++;
	cov->vals3[j-start]++;

	mercov[i] = seq->id[0];

	if (seq->sequence[j] == topseq->sequence[j]) {
	  cov->vals2[i]++;
	}

	if (!(strcmp(seq->id,"Mouse7") == 0 ||
	      strcmp(seq->id,"Rat") ==0 ||
	      strcmp(seq->id,"RheMac") == 0 ||
	      strcmp(seq->id,"Chimp") == 0 ||
	      strcmp(seq->id,"Dog2") == 0 ||
	      strcmp(seq->id,"Human") == 0))
	  {
	    // printf("Mdcount %s\n",seq->id);
	    mdcount++;
	}
      }
      seq = seq->next;
      
      i++;

    }

      cov->xcov[j-start] = mdcount;
    if (mdcount == 2) {
      cov->mdcount++;
    } else if (mdcount == 0) {
      cov->no_md++;

    }
    j++;
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

char *find_mercov(Sequence *seqptr,int numseqs,char **cov,int *fcov,Gff *gff) {

  char *mercov = (char *)malloc((numseqs+1)*sizeof(char));
  
  int start = gff->start;
  int end   = gff->end;

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
      if (cov[j][i] != '-') {
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
  
  fcov[start] = numcov;
  return mercov;
}
    
      
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  







