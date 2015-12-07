#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "datamodel.h"
#include "pog_utils.h"

//james change 4


#define LOG2 log(2)
#define LOG25 log(0.25)

static const char dna_complement[] =
"                                                                "
" TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      "
"                                                                "
"                                                                ";
/* ................................................................ */
/* @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~. */
/* ................................................................ */
/* ................................................................ */


Pwm *read_pwm(FILE *file);
Pwm **fetch_pwms(int *pwmcount,char *file);
void revcomp(char *s);

double score_matrix(Sequence *seq, Pwm *pwm, int start, int *strand);

int main(int argc, char *argv[]) {

  Genome *human   = fetch_human_genome();

  int pwmcount;

  char *matfile = argv[1];
  int   chroffset = atoi(argv[2]);

  Pwm **pwm = fetch_pwms(&pwmcount,matfile);

  int chunk = 1000000;


  int chrnum = chroffset;

  while (chrnum < HUMAN_CHRNUM) {
    int   len = human->chr[chrnum]->length;
    
    int startcoord = 0;

    while (startcoord < len) {
      int size = chunk;
      if (startcoord + chunk-1 > len) {
         size  = len - startcoord;
      }
      if (size > 5000) {

      GenomePosition *pos   = (GenomePosition *)malloc(sizeof(GenomePosition));

      pos->chr = human->chr[chrnum];
      pos->pos = startcoord;

      Sequence *seqptr  = read_fasta_test(pos,size);             printf("Read sequence\n");	

      /////////////////////////////////////////////// Loop over pwms //////////////////////////

      int pwmnum = 0;
	
      while (pwmnum < pwmcount) { 

	////////////////////////////////////////////  Find all scores ////////////////////////

	double *scores = (double *)malloc((size-20)*sizeof(double));
	int    *strand = (int *)malloc((size-20)*sizeof(int));

	int coord = 0;

	while (coord < size - 20) {
	  int snum;

	  double score = score_matrix(seqptr,pwm[pwmnum],coord,&snum);  

	  scores[coord] = score;
	  strand[coord] = snum;

	  coord++;
	}

	//////////////////////////////////////////// Find the max score in 5kb windows ///////

	int bincount = (int)(size/2500) - 1;

	double *maxes     = (double *)malloc(bincount*sizeof(double));
	int    *maxespos  = (int    *)malloc(bincount*sizeof(int));

	int i = 0;

	while (i < bincount) {

	  int j = 0;

	  double maxscore = -10000000;
	  int    maxpos   = 0;

	  while (j < 5000) {
	    if (scores[i*2500+j] > maxscore) {
	      maxscore = scores[i*2500+j];
	      maxpos   = i*2500+j;
	    }
	    j++;
	  }

	  //printf("Maxes %f\t%d\t%d\n",maxscore,maxpos,i);

	  maxes[i]    = maxscore;
	  maxespos[i] = maxpos;

	  i++;
	}

	//////////////////////////////////////////// Find the positions within score of 2 of the max ////


	coord = 0;

	while (coord < size-20) {

	  int bin = (int)(coord/2500);

	  if (bin >= bincount) {
	    bin = bincount-1;
	  }

	  double maxscore = maxes[bin];
	  int    maxpos   = maxespos[bin];

	  if (scores[coord]  >= maxscore && maxscore != 0) {
	    printf("%s\ttfscan\ttfscan\t%d\t%d\t%f\t%d\t.\t%s\t%d\t\t%f\t%d\t%d\n",human->chr[chrnum]->name,coord+startcoord,coord+startcoord + pwm[pwmnum]->len - 1,scores[coord],strand[coord],pwm[pwmnum]->name,pwmnum,maxscore,maxpos,bin);
	  }

	  coord++;
	}

	free(maxes);
	free(maxespos);
	free(scores);
	free(strand);

	//////////////////////////////////////////// End of pwm loop ////////////////////////////////

	pwmnum++;
      }

      free(pos);
      free_seqptr(seqptr);
      }
      startcoord += chunk;

      printf("Coord %d\n",startcoord);
    }
    chrnum++;
  }
  return 1;
}
  

double score_matrix(Sequence *seq, Pwm *pwm, int start, int *strand) {

  Sequence *tmpseq = seq;

  double totscore1=0;
  double totscore2=0;

  double numseqs=0;

  // start - position in seqstring to start
  char *sub = substring(tmpseq->sequence,start,start+pwm->len-1);
  
  int k     = 0;
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
  

  if (count < 3) {
    
    // First the forward strand
    int j = 0;
    double score = 0;

    double log2 = log(2);    


    while (j < pwm->len) {
    
      double tmpscore = 0;
      if (sub[j] == 'A') {
	tmpscore += log(pwm->vals[j*4]/0.25);
	
      } else if (sub[j] == 'T') {
	tmpscore += log(pwm->vals[j*4+1]/0.25);
	
      } else if (sub[j] == 'C') {
	tmpscore += log(pwm->vals[j*4+2]/0.25);
	
      } else if (sub[j] == 'G') {
	tmpscore += log(pwm->vals[j*4+3]/0.25);
      }
      tmpscore *= pwm->inf[j];      
      score += tmpscore;

      j++;
    }
    score /= LOG2;

    totscore1 += score;
    numseqs++;
    
    // Now the reverse strand
    revcomp(sub);
    j = 0;
    score = 0;

    
    while (j < pwm->len) {
      double tmpscore = 0;
      if (sub[j] == 'A') {
	tmpscore += log(pwm->vals[j*4]/0.25)/log2;
	
      } else if (sub[j] == 'T') {
	tmpscore += log(pwm->vals[j*4+1]/0.25)/log2;
	
      } else if (sub[j] == 'C') {
	tmpscore += log(pwm->vals[j*4+2]/0.25)/log2;
	
      } else if (sub[j] == 'G') {
	tmpscore += log(pwm->vals[j*4+3]/0.25)/log2;
      }
      
      tmpscore *= pwm->inf[j];
      score += tmpscore;
      j++;
    }
    totscore2 += score;
  }

  totscore2 /= (2*pwm->len);
  totscore1 /= (2*pwm->len);

  free(sub);
  
  double outscore = totscore1;
  *strand = 1;
  if (totscore2 > totscore1) {
    outscore = totscore2;
    *strand = -1;
  }
  
  return outscore;
}

int *find_coverage(Sequence *seqptr) {
  int numcols = seqptr->length;
  
  int *cov;
  
  cov = (int *)calloc((size_t) numcols,sizeof(int));
  
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
  
  char **cov = (char **)calloc((size_t) numcols,sizeof(char *));
  
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
    
      
void revcomp(char *s) {
  char *p = s + strlen(s) - 1;

  while (s <= p) {
    char c = dna_complement[(int)*s];
    *s = dna_complement[(int)*p];
    *p = c;
    ++s;
    --p;
  }
}


Pwm **fetch_pwms(int *pwmcount,char *mat) {

  FILE *matfile   = fopen(mat,"r"); //"/Users/mclamp/cvs/pogvue/data/tf9.4.matrix.pwmline","r");

  Pwm   **pwm     = (Pwm **)malloc(800 * sizeof(Pwm *));
  Pwm    *tmppwm;

  *pwmcount = 0;

  while ((tmppwm = read_pwm(matfile))!= NULL) {

    if (tmppwm != NULL) {

      pwm[*pwmcount] = tmppwm;
      (*pwmcount)++;

    }
  }

  fclose(matfile);
  
  return pwm;
}


