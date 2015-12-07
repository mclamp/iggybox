#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "datamodel.h"

typedef struct _pimer {
  char  *chr;
  int    start;
  int    end;
  double logodds;
} Pimer;

extern char **read_tokens(FILE *file,char token, int *ntok);
extern Gff *read_gff(FILE *file);
extern void print_gff(FILE *file,Gff *gff);
Pimer *read_pimer(FILE *file) {

  int ntok;

  char **tokens;

  int count = 0;
  
  while ((tokens = read_tokens(file,'\t',&ntok))) {
    count++;
    if (count % 10000 == 0) {
      printf("Tokens : %d : %d : %s : %s : %s : %s\n",count,ntok,tokens[0],tokens[1],tokens[2],tokens[3],tokens[5]);
    }
    if (ntok == 6) {
      Pimer *pimer = (Pimer *)malloc(sizeof(Pimer));
      
      pimer->chr = tokens[0];
      pimer->start = atoi(tokens[1]);
      pimer->end   = atoi(tokens[2]);
      pimer->logodds = (double)atof(tokens[4]);
      
      return pimer;
    }
  }
  return NULL;
}

extern Sequence *read_fasta(char *filename);


int main(int argc, char * argv[]) {
    FILE *file;

    Sequence *seq;
    Sequence *seqptr;

    int numseqs = 0;
    int i       = 0;
    int numcols = 0;
    int length;

    int *cov;

    
    seqptr = read_fasta(argv[1]);

    seq = seqptr;

    while (seq != NULL) {
      numseqs++;
      seq = seq->next;
    }

    length = seqptr->length;

    numcols = length;

    cov  = (int *)malloc(numcols*sizeof(int));

    i = 0;

    while (i < numcols) {
      cov[i] = 0;
      i++;
    }

    i = 0;

    while (i < numcols) {
      seq    = seqptr;

      int tmp   = 0;
      int count = 0;

      while (count < numseqs) {

        if (seq->sequence[i] != '-') { 
          tmp++;
        }
        seq    = seq->next;
	
        count++;
      } 

      cov[i] = tmp;

      i++;
    }

    file = fopen(argv[2],"r");

    Pimer *pimer;
  
  int **bins = (int **)malloc(numseqs*sizeof(int *));
  
  i = 0;
  
  while (i < numseqs) {
    int j = 0;
    bins[i] = (int *)malloc(200*sizeof(int));
    while (j < 200) {
      bins[i][j] = 0;
      j++;
    }
    i++;
  }

  Pimer **pimers = (Pimer **)malloc(numcols*sizeof(Pimer *));
  
  i = 0;
  
  while (i < numcols) {
    pimers[i] = NULL;
    i++;
  }
  int count = 0;
  int bincount = -1;
  while ((pimer = read_pimer(file)) != NULL)  {
    pimers[count] = pimer;
    if (count % 10000 == 0) {
      fprintf(stderr,"Pimer %s %d %d %f %d\n",pimer->chr,pimer->start,pimer->end,pimer->logodds,cov[pimer->start]);
    }
    count++;
    int bin1 = (int)(pimer->logodds);
    int bin2 = cov[pimer->start];
  
    bins[bin2][bin1]++;
 //   printf("Bin\t%d\t%d\n",bin1,bin2);
    if (bin1 > bincount) {
      bincount = bin1;
    }
  }

  i = 0;
  while (i < bincount) {
    int j = 0;
    while (j < numseqs) {
      if (j > 0) {
       // bins[j][i] += bins[j-1][i];
      }
   //   printf("%d\t%d\n",i,bins[j][i]);
      j++;
    }

    i++;
  }
    
  FILE *gfffile = fopen(argv[3],"r");
  
  Gff *gff = read_gff(gfffile);
  
  while (gff != NULL) {

    
    int i = gff->start;
    
    if (i > numcols) {
      printf("Coordinate out of range %d [range is 0-%d]\n",i,numcols);
      exit(0);
    }
    float score = 0;
    
    while (i <= gff->end) {
      Pimer *pimer = pimers[i];
      
      if (pimer != NULL) {
        if (pimer->logodds > 80 && cov[i] > 6) {
          printf("chr1\tpimer\tpimer\t%d\t%d\t%f\t1\t.\n",i-6,i+6,pimer->logodds);
        }
        score += pimer->logodds;
        //printf("%d\t%d\t%f\t%s\n",i,cov[i],pimer->logodds,gff->hseqname);
      }
      i++;
    }
 //   if (score/(gff->end-gff->start+1) > 70) {
      gff->score = score/(gff->end-gff->start+1);
      print_gff(stdout,gff);    
 //   }
    gff = gff->next;
  }
  
    
    
}

