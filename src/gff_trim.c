#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "datamodel.h"

extern char **read_tokens(FILE *file,char token, int *ntok);
extern Gff *read_gff(FILE *file);
extern void print_gff(FILE *file,Gff *gff);
extern Sequence *read_fasta(char *filename);

extern Pimer *read_pimer(FILE *file);

int main(int argc, char * argv[]) {
    FILE *file;

    Sequence *seq;
    Sequence *seqptr;
    Sequence **seqarr;

    int numseqs = 0;
    int i       = 0;
    int numcols = 0;
    int length;

    int *cov;

    
    seqptr = read_fasta(argv[1]);

    seq = seqptr;


  
    while (seq != NULL) {
      numseqs++;
      //printf("Sequence %s %d %d\n",seq->id,seq->length,numseqs);
      seq = seq->next;

    }
    seqarr = (Sequence **)malloc(numseqs * sizeof(Sequence *));
  
  i = 0;
  seq = seqptr;
  
  while (i < numseqs) {
    seqarr[i] = seq;
    seq = seq->next;
    i++;
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
      if (i % 10000 == 0) {
	//printf("Coverage for %d is %d\n",i,cov[i]);
      }
      i++;
    }

  FILE *gfffile = fopen(argv[2],"r");

  char **trimstr = (char **)malloc(numseqs*sizeof(char *));
    
  int j = 0;
    
  while (j < numseqs) {
    trimstr[j] = (char *)malloc(numcols*sizeof(char *));
    j++;
  }
  
  Gff *gff = read_gff(gfffile);
  int len = 0;
  
  int tmp = 0;
  int count = 0;
  while (gff != NULL) {
    
    int i = gff->start;
    // print_gff(stdout,gff);
    tmp += gff->end - gff->start+1;
    count++;
    // if (count % 10000 == 0) {
    //printf("Count %d %d\n",count,tmp);
      //}
    if (i > numcols || gff->end < 1) {
      printf("Coordinate out of range %d [range is 0-%d]\n",i,numcols);
      exit(0);
    }

    if (gff->end > numcols) {
      gff->end = numcols-1;
    } 
    if (gff->start < 1) {
      gff->start = 1;
    }
    int j   = 0;
    
    int allpos = 0;

    while (j < numseqs) {
      Sequence *tmpseq = seqarr[j];
      
      int k = 0;
      int pos = 0;

      while (k <= (gff->end - gff->start + 1)) {
	  trimstr[j][len+pos] = tmpseq->sequence[gff->start+k];
	  pos++;
          k++;
      }
      if (allpos == 0) {
	allpos = pos;
      }
      trimstr[j][len+pos] = '\0';    
      j++;
  }
  len = len + allpos;
  
    gff = gff->next;
  }
  
  i = 0;
  
  while (i < numseqs) {
    printf(">%s\n",seqarr[i]->id);
    
    int j = 0;
    
    while (j < 300000 && j < len) {
      printf("%c",trimstr[i][j]);
      if (j % 72 == 0) {
        printf("\n");
      }
      j++;
    }
    printf("\n");
    i++;
  }
            
    
    
}

