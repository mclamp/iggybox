#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "datamodel.h"

char *int2char(int i);
Pwm *read_pwm(FILE *file);

double score_matrix(Sequence *seq, Pwm *pwm, int start, int offset, char *chr,FILE *file);

int main(int argc, char *argv[]) {

  Genome *human = fetch_human_genome();

  int j   = 0;

  int offset = atoi(argv[1]);
  int chunk  = atoi(argv[2]);
  char *stub = argv[3];

  
  //FILE *matfile = fopen("/Users/mclamp/cvs/pogvue/data/tf9.4.matrix.pwmline","r");

  FILE *matfile = fopen("/home/mclamp/iccpog/src/tf9.4.matrix.pwmline","r");
  
  Pwm **pwm = (Pwm **)malloc(800*sizeof(Pwm *));
  Pwm *tmppwm;
  int  pwmcount = 0;

  FILE **outfiles = (FILE **)malloc(800*sizeof(FILE *));
  char **outfilestr = (char **)malloc(800*sizeof(char *));

  while ((tmppwm = read_pwm(matfile))!= NULL) {


    if (tmppwm != NULL) {
      pwm[pwmcount] = tmppwm;
      char *countstr = int2char(pwmcount);
      
      outfilestr[pwmcount] = (char *)malloc((strlen(stub)  + strlen(".") + strlen(tmppwm->name) + strlen(".") + strlen(countstr) + strlen(".dat") + 1)*sizeof(char));
      
      strcpy(outfilestr[pwmcount],stub);  
      strcat(outfilestr[pwmcount],".");
      strcat(outfilestr[pwmcount],tmppwm->name);
      strcat(outfilestr[pwmcount],".");
      strcat(outfilestr[pwmcount],countstr);
      strcat(outfilestr[pwmcount],".dat");
      
      printf("File is %s\n",outfilestr[pwmcount]);
      
      FILE *tmpfile = fopen(outfilestr[pwmcount],"w");
      
      outfiles[pwmcount] = tmpfile;
      printf("File %d\n",pwmcount);
      pwmcount++;
    }

  }

  fclose(matfile);
  printf("Done pwms\n");

  while (j < human->chrnum) {

    if (j >= offset && j < offset + chunk) {
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

        int ii = 0;
	
        while (ii < pwmcount) {
	  Pwm *tmppwm = pwm[ii];
	  
	  len = seqptr->length;
	  //printf("pwm %s\n", tmppwm->name,ii);
	  int k = 0;
	  
	  while (k < len){
	    Sequence *seq = seqptr;
	    
	    double score = 0;//score_matrix(seq,tmppwm,k,i,pos->chr->name,outfiles[ii]);
	    
	    k += 50000;
	  }
	  ii++;
	}
	free_seqptr(seqptr);
	free(pos);
        int k = 0;	
        while (k < len){
          Sequence *seq = seqptr;

          double score = 0;//score_matrix(seq,pwm,k,i,pos->chr->name);

          k += 10000;
        }
        free_seqptr(seqptr);
        free(pos);
      
        i += 1000000;
      }
    }
    j++;
    
  }
  return 1;
}

double score_matrix(Sequence *seq, Pwm *pwm, int start, int offset, char *chr, FILE *file) {

  Sequence *tmpseq = seq;

  while (tmpseq != NULL) {
    double score = 0.0;    
    int    j     = 0;

    char *sub = substring(tmpseq->sequence,start,start+pwm->len-1);

    int k = 0;
    int count = 0;

    while (k < pwm->len) {
      if (sub[k] == '-' || sub[k] == 'N') {
        count++;
      }
      if (sub[k] > 95) {
        sub[k] -= 32;
      }
      k++;
    }

    double gcval = 0.333333;
    double atval = 0.2;

    double log2 = log(2);

    if (count < 3) {
      while (j < pwm->len) {
	
        if (sub[j] == 'A') {
          score += log(pwm->vals[j*4]/atval)/log2;
	  
        } else if (sub[j] == 'T') {
          score += log(pwm->vals[j*4+1]/atval)/log2;
	  
        } else if (sub[j] == 'C') {
          score += log(pwm->vals[j*4+2]/gcval)/log2;
	  
        } else if (sub[j] == 'G') {
          score += log(pwm->vals[j*4+3]/gcval)/log2;
        }
	
        j++;
      }
      
      fprintf(file,"S\t%s\t%s\t%d\t%s\t%7.2f\t%s\n",tmpseq->id,chr,start+offset,pwm->name,score,sub);

    }
    free(sub);
    tmpseq = NULL;
    //tmpseq = tmpseq->next;
  } 
      
    

  return 0.0;
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

char *find_mercov(Sequence *seqptr,int numseqs,char **cov,int *fcov,Pimer *pi) {

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
    
      
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  







