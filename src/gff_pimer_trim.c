#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "datamodel.h"

extern char    **read_tokens(FILE *file,char token, int *ntok);
extern Gff      *read_gff   (FILE *file);
extern void      print_gff  (FILE *file,Gff *gff);
extern Sequence *read_fasta (char *filename);
extern Pimer    *read_pimer (FILE *file);

int *find_coverage(Sequence *seqptr, int numseqs, int numcols, int *covcount,Pimer **pimers);

Pimer **shuffle_pimers(Pimer **pimers, Sequence* seqptr, int numseqs, int numcols);


int main(int argc, char * argv[]) {
	
    Sequence *seq;
    Sequence *seqptr;
    Sequence **seqarr;
	
    int *covcount;
    int numseqs = 0;
    int i       = 0;
    int numcols = 0;
    int length;
	
    // A structure holding a Sequence array would be good here - would have numseqs, maybe coverage?
    
    seqptr = read_fasta(argv[1]);
	
    seq = seqptr;
	
    while (seq != NULL) {
		numseqs++;
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
	
    covcount  = (int *)malloc(numseqs*sizeof(int));       //  number of cols with n sequences   
    length  = seqptr->length;
    numcols = length;
	
    printf("Got seqs %d %d\n",numseqs,numcols);
	

    
    FILE *pimerfile = fopen(argv[2],"r");
    printf("Got pimers\n");    
    Pimer **pimers = (Pimer **)malloc(numcols*sizeof(Pimer *));
	
    i = 0;
    
    while (i < numcols) {
      pimers[i] = NULL;
      i++;
    }
    
    Pimer *pimer;
    int count = 0;
    int maxcol = 0;
    
    while ((pimer = read_pimer(pimerfile)) != (Pimer *)NULL) {
      pimers[pimer->start] = pimer;
      count++;
      if (pimer->end > maxcol) {
	maxcol = pimer->end;
      }
    }
	
    printf("Got pimers %d %d\n",maxcol,count);
    int *cov = find_coverage(seqptr,numseqs,numcols,covcount,NULL);	
    //Pimer **shuffled_pimers = shuffle_pimers(pimers, seqptr, numseqs, numcols);
	
    i = 0;
	

    if (argc == 4) {
      FILE *gfffile = fopen(argv[3],"r");
      
      Gff *gff = read_gff(gfffile);
      
      while (gff != NULL) {
	
	int i = gff->start;
	//  print_gff(stdout,gff);
	
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
	
	int k = gff->start;
	
	while (k <= gff->end && k <= numcols ) {
	  if (pimers[k] != NULL && pimers[k]->start > 0) {
	    Pimer *tmp = (Pimer *)pimers[k];
	    printf("chr1\t%s\tpimer\t%d\t%d\t%d\t1\t.\t%f\t%d\n", gff->hseqname,tmp->start,tmp->end,cov[k],tmp->logodds,k);
		    
	    // printf("%s\tpimer\tpimer\t%d\t%d\t%f\t1\t.\n",gff->seqname,k,k,pimers[k]->logodds);
	  }
	  k++;
	  
	}
		
	gff = gff->next;
      }
    } else {
      
      int k = 0;
      while (k < numcols) {
	if (pimers[k] != NULL && pimers[k]->start > 0) {
	  Pimer *tmp = (Pimer *)pimers[k];
	  printf("chr1\t%s\tpimer\t%d\t%d\t%d\t1\t.\t%f\t%d\n", "pimer",tmp->start,tmp->end,cov[k],tmp->logodds,k);
	}
	k++;
      }
    }

}

Pimer **shuffle_pimers(Pimer **pimers, Sequence* seqptr, int numseqs,int numcols) {
    int **covpos    = (int **)malloc(numseqs*sizeof(int *));    //  list of cols with n sequences
    int **done      = (int **)malloc(numseqs*sizeof(int *));    //  stores the taken columns for each coverage no.
	
    Sequence **seqarr = (Sequence **)malloc(numseqs*sizeof(Sequence *));
	
    int *covcount  = (int *)malloc(numseqs*sizeof(int));       //  number of cols with n sequences
    int *cov       = (int *)malloc(numcols*sizeof(int));       //  number of seqs at each column 
    int *tmppos    = (int *)malloc(numseqs*sizeof(int));       //  End position of each coverage levels cov storarge
	
    Sequence *seq = seqptr;
	
	Pimer **shuffled_pimers = (Pimer **)malloc(numcols*sizeof(Pimer *));
	
    int i = 0;
	
    while (i < numseqs) {
		seqarr[i]   = seq;
		covcount[i] = 0;
		tmppos[i]   = 0;
		
		seq = seq->next;
		
		i++;
		
    }
    
	i = 0;
	while (i < numcols) {
		shuffled_pimers[i] = (Pimer *)NULL;
		i++;
	}
    cov = find_coverage(seqptr,numseqs,numcols,covcount,pimers);       // Index is column number
	
    i = 0;
	
    while (i < numseqs) {
		covpos[i] = (int *)malloc(covcount[i]*sizeof(int));
		done[i]   = (int *)malloc(covcount[i]*sizeof(int));
		printf("Covcount for %d is %d\n",i,covcount[i]);
		i++;
    }
	
    // Fill in the coverage arrays
	
    i = 0;
	
    while (i < numcols) {
		int seqs = cov[i];
		//printf("Seqs %d %d\n",seqs,i);
		int pos  = tmppos[seqs];
		if (pimers[i] != NULL) {
			//printf("Pos %d %d\n",pos,covcount[seqs]);
			done[seqs][pos]   = 0;
			covpos[seqs][pos] = i;
			//printf("Covpos %d %d %d\n",seqs,pos,i);
			tmppos[seqs]++;
			//printf("Tmppos %d %d\n",seqs,tmppos[seqs]);
			//printf("Position %d\tcoverage\t%d\tArr pos\t%d\n",i,seqs,pos);
		}
		i++;
    }
	
    i = 0;
	
    printf("done\n");
    // Now loop over each column 
	
    while (i < numcols) {
		
		// Find coverage at column i
		int seqs = cov[i];                                 
		
		if (covcount[seqs] > 0) {
			//      printf("Seqs %d\n",seqs);
			// Pick a rand column newcol from 0 to covcount[seqs]
			int newcol = (int)((double)covcount[seqs]*rand()*1.0/(1.0*RAND_MAX));
			
			while (done[seqs][newcol] == 1) {
				newcol = (int)((double)covcount[seqs]*rand()*1.0/(1.0*RAND_MAX));
			}
			
			int newcol2 = covpos[seqs][newcol];
			
			double score = pimers[newcol2]->logodds;
			shuffled_pimers[i] = pimers[newcol2];

			done[seqs][newcol] = 1;
			
			if (i % 10000 == 0) {
				fprintf(stderr,"Score is %d %d %f %d\n",i,newcol2,score,seqs);
				fprintf(stderr,"Found col %d %d\n",newcol,i);
				fprintf(stderr,"Pimer %f\n",shuffled_pimers[i]->logodds);
			}
		}
		i++;
    }
	i = 0;
	while (i < 100) {
		if (shuffled_pimers[i] != NULL) {
			printf("Shuff %d %f\n",i,shuffled_pimers[i]->logodds);
		}
		i++;
	}
	return shuffled_pimers;
}

int *find_coverage(Sequence *seqptr, int numseqs, int numcols,int *covcount, Pimer **pimers) {
	
	int *cov = (int *)malloc(numcols*sizeof(int));
	
	int i = 0;
	
	Sequence *seq;
	
	while (i < numcols) {
		cov[i] = 0;
		covcount[i] = 0;
		i++;
	}
	
	i = 0;
	
	while (i < numcols) {
		seq    = seqptr;
		
		if (pimers == NULL || pimers[i] != NULL) {
			int tmp   = 0;
			int count = 0;
			
			while (count < numseqs) {
				if (seq->sequence[i] != '-') { 
					tmp++;
				}
				count++;
				seq    = seq->next;
			} 
			cov[i] = tmp;
			covcount[tmp]++;
		}
		i++;
	}
	
	return cov;
}
