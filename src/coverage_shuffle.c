#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "datamodel.h"

extern Sequence *read_fasta(char *filename);

int main(int argc, char * argv[]) {
    FILE *file;

    Sequence  *seq;
    Sequence  *seqptr;
    Sequence **seqarr;

    char     **newstr;

    int numseqs = 0;
    int i       = 0;
    int count   = 0;
    int numcols = 0;
    int length;


    int **covpos;
    int **done;

    int *covcount;
    int *cov;
    int *tmppos;

    char c;
    int  col;

    seqptr = read_fasta(argv[1]);

    seq = seqptr;

    while (seq != NULL) {
      numseqs++;
      seq = seq->next;
    }

    length  = seqptr->length;
    numcols = length;


    covpos    = (int **)malloc(numseqs*sizeof(int *));    //  list of cols with n sequences
    done      = (int **)malloc(numseqs*sizeof(int *));    //  stores the taken columns for each coverage no.

    seqarr    = (Sequence **)malloc(numseqs*sizeof(Sequence *));

    covcount  = (int *)malloc(numseqs*sizeof(int));       //  number of cols with n sequences
    cov       = (int *)malloc(numcols*sizeof(int));       //  number of seqs at each column 
    tmppos    = (int *)malloc(numseqs*sizeof(int));       //  End position of each coverage levels cov storarge

    newstr    = (char **)malloc(numseqs*sizeof(char *));   // Store the new strings

    i = 0;

    seq = seqptr;

    while (i < numseqs) {
      seqarr[i] = seq;
      covcount[i] = 0;
      tmppos[i] = 0;

      seq = seq->next;
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

      cov[i]= tmp;

      covcount[tmp]++;

      i++;
    }


    i = 0;

    while (i < numseqs) {
      covpos[i] = (int *)malloc(covcount[i]*sizeof(int));
      done[i]   = (int *)malloc(covcount[i]*sizeof(int));

      i++;
    }

    // Fill in the coverage arrays

    i = 0;

    while (i < numcols) {
      int seqs = cov[i];
      int pos  = tmppos[seqs];

      done[seqs][pos] = 0;

      covpos[seqs][pos] = i;

      tmppos[seqs]++;

      //printf("Position %d\tcoverage\t%d\tArr pos\t%d\n",i,seqs,pos);

      i++;
    }


    i = 0;

    while (i < numseqs) {
      newstr[i] = (char *)malloc((numcols+1)*sizeof(char));
      tmppos[i] = 0;

      i++;
    }

    i = 0;

    // Now loop over each column 

    while (i < numcols) {

      // Find coverage at column i
      int seqs = cov[i];                                 

      // Pick a rand column newcol from 0 to covcount[seqs]
      int newcol = (int)((double)covcount[seqs]*rand()*1.0/(1.0*RAND_MAX));

      while (done[seqs][newcol] == 1) {
	newcol = (int)((double)covcount[seqs]*rand()*1.0/(1.0*RAND_MAX));
      }


      int newcol2 = covpos[seqs][newcol];

      int j = 0;
      int pos =0;

      while (j < numseqs) {

	Sequence* tmpseq  = seqarr[j];

	char c = tmpseq->sequence[newcol2];

	if (c != 'A' &&
	    c != 'C' &&
	    c != 'T' &&
	    c != 'G' &&
	    c != 'a' &&
	    c != 'c' &&
	    c != 't' &&
	    c != 'g' &&
	    c != 'N' && 
	    c != 'n' && 
	    c != '-') {
	  fprintf(stderr,"AAAAAARGGGG! %d %d %c\n",newcol,newcol2,c);
	}

//	if (c != '-') {

//	  while (seqarr[pos]->sequence[i] == '-') {
//	    newstr[pos][i] = '-';


//	    pos++;
//	  }
	  //fprintf(stderr,"Old col %d new col %d - old char %c new char %c - old seq %d new seq %d cov %d\n",newcol2,i,seqarr[pos]->sequence[i],c,j,pos,seqs);
	  //newstr[pos][i] = c;
	  newstr[j][i] = c;
//	  pos++;
//	} 


	j++;

      }

     // while (pos < j) {
//	newstr[pos][i] = '-';
//	pos++;
//      }

      tmppos[seqs]++;
      done[seqs][newcol] = 1;

      if (i % 10000 == 0) {
	fprintf(stderr,"Found col %d %d\n",newcol,i);
      }

      // Newcov[i] = cov[newpos];

      i++;
    }

    i = 0;
    while (i < numseqs) {
      newstr[i][numcols] = '\0';
      i++;
    }

    i = 0;

    while (i < numseqs) {
      printf(">%s\n",seqarr[i]->id);

      int j = 0;
      while (j < numcols && j < 300000) {
	printf("%c", newstr[i][j]);
	if (j % 72 == 0) {
	  printf("\n");
	}
	j++;
      }
      printf("\n");
      i++;
    }

}

