#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "datamodel.h"

extern SeqArray *read_seq_array(char *filename);
extern PimerArray *read_pimer_array(char *filename, int numcols);
extern Pimer *read_pimer(FILE *file);
extern Sequence *read_fasta(char *filename);
int *find_coverage(Sequence *seqs, int start, int end, int numseqs);
int main(int argc, char **argv) {
float find_average_coverage(Sequence *seqs, int *cov,int start, int end, int numseqs);
float find_correlation(Sequence *seqs, float *vals, int start, int end, int numseqs, int offset);
	
  //SeqArray *seq_array = read_seq_array(argv[1]);

  //printf("Array %d\n",seq_array->len);
  //printf("Array %d\n",seq_array->seq[1]->length);

//	int numcols = seq_array->seq[1]->length;
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
    
	float *vals = (float *)malloc(sizeof(float)*numcols);
	i = 0;
	while (i < numcols) {
		vals[i] = 0.0;
		i++;
	}
    while ((pimer = read_pimer(pimerfile)) != (Pimer *)NULL) {
		pimers[pimer->start] = pimer;
		vals[pimer->start] = pimer->logodds;
		count++;
		if (pimer->end > maxcol) {
			maxcol = pimer->end;
		}
    }
	
    printf("Got pimers %d %d\n",maxcol,count);

	
	int *cov = find_coverage(seqptr,0,numcols,numseqs);
	
	printf("Cov %d\n",cov[0]);
	
	//PimerArray *pis = read_pimer_array(argv[2],seq_array->seq[1]->length);
	//printf("Pimers %d\n",pis->len);
	//printf("Pimers %f\n",pis->pimer[0]->logodds);
	
	
	i = 0;
	
	int len = 50;
	

	while (i < numcols - 2*len) {
		if (pimers[i] != NULL) {	
		// Get average coverage for i - i+50
		float icov = find_average_coverage(seqptr,cov,i,i+len,numseqs);
			//printf("Icov %d %f\n",i,icov);
		
		int l = 1;
			
		while (l < len) {
			float lcov = find_average_coverage(seqptr,cov,i+l,i+l+len,numseqs);	
			//printf("Lcov %d %f\n",i,lcov);
			float corr = find_correlation(seqptr,vals,i,i+len,numseqs,l);
			
			float odds = 0.0;
			if (pimers[i] != NULL) {
				odds = pimers[i]->logodds;
			}
			printf("I\t%d\t%f\t%f\t%d\t%f\t%f\n",i,icov,lcov,l,corr,odds);
			l++;
		}
		// print i, l, icov, lcov
		}
		i++;
	}
  return 0;
}
float find_average_coverage(Sequence *seqs, int *cov,int start, int end, int numseqs) {
	
	int tot = 0;
	int i   = 0;
	
	while (i < end-start+1) {
//		printf("Cov %d\n",cov[i]);
		tot += cov[i+start];
		i++;
	}
	float avcov = (tot*1.0/(1.0*(end-start+1)));
	
	return avcov;
}
int *find_coverage(Sequence *seqs, int start, int end, int numseqs) {
	
	//printf("Num %d %d %d\n",numseqs,start, end);
	
	int *cov = (int *)malloc((end-start+1)*sizeof(int));
	
	//printf("Cov done %d\n",cov);
	
	int i = 0;
	
	while (i < (end-start+1)) {
		cov[i] = 0;
		i++;
	}
	
	i = 0;
	
	while (i < (end-start+1)) {
		Sequence *seq = seqs;
		int num = 0;
		
		while (seq->next != NULL) {
			if (seq->sequence[i+start] != '-') {
				num++;
			}
			seq = seq->next;
		}
		cov[i] = num;
		i++;
	}
	
	return cov;
}

float find_correlation(Sequence *seqs, float *vals, int start, int end, int numseqs, int offset) {
	
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