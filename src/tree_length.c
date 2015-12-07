#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>

#include "datamodel.h"

extern Sequence *read_fasta(char *filename);

int main(int argc, char * argv[]) {
    FILE *file;

    Sequence  *seq;
    Sequence  *seqptr;
    Sequence **seqarr;

    int numseqs = 0;
    int numcols = 0;
    int length;
    
    long long *val;
    int   valnum = 0;
    int  *valcount;



    seqptr = read_fasta(argv[1]);

    seq = seqptr;

    while (seq != NULL) {
      numseqs++;
      seq = seq->next;
    }

    length  = seqptr->length;
    numcols = length;

    seqarr    = (Sequence **)malloc(numseqs*sizeof(Sequence *));

    int i = 0;

    seq = seqptr;

    while (i < numseqs) {
      seqarr[i] = seq;
      seq = seq->next;
      i++;
    }

    i = 0;

    while (i < numcols) {
      int count = 0;
      long long tmpval = 0;

      while (count < numseqs) {

	if (seqarr[count]->sequence[i] != '-') { 
	  tmpval += (long long)pow(2,count);
	  //	  printf("I count %d  %d %lld\n",i,count,tmpval);
	}
	
	count++;

      }

      // printf("Val %d\t%lld\t%d\n",i,tmpval,valnum);

      int j     = 0;
      int found = 0;
      
      while (found == 0  && j < valnum) {

	//	printf("Val2 %lld  %lld %d\n",tmpval,val[j],j);
	if (val[j] == tmpval) {
	  found = 1;
	  valcount[j]++;

	  //	  printf("Found existing for %lld  count %d  num %d\n",val[j],valcount[j],j);
	}
	j++;
      }

      if (found == 0) {
	if (valnum == 0) {
	  val      = (long long *)malloc(sizeof(long long));
	  valcount = (int *)malloc(sizeof(int));

	  valnum = 1;
	  val[0] = tmpval;
	  valcount[0] = 1;

	  //	  printf("Val3 for %d %lld %d\n",valnum,val[0],valcount[valnum-1]);
	  
	} else {
	  val      = (long long *)realloc(val,     (valnum+1)*sizeof(long long));
	  valcount = (int *)realloc(valcount,(valnum+1)*sizeof(int));

	  val[valnum] = tmpval;
	  valcount[valnum] = 1;

	  //  printf("Val4 for %d %lld %d\n",valnum,val[valnum],valcount[valnum]);

	  valnum++;
	}
      }

      if (i % 10000 == 0) {
	fprintf(stderr,"On %d %d %lld %d\n",i,valnum,val[i],valcount[i]);
      }
      i++;
    }

    i = 0;

    while (i < valnum) {
      printf("Val %d\t%lld\t%d\n",i,val[i],valcount[i]);
      i++;
    }
}
