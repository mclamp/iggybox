#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "datamodel.h"

Pwm *read_pwm(FILE *file);
void revcomp(char *s);
//Thresh *read_threshfile();
char **read_tokens(FILE *file,char token, int *ntok);

static const char dna_complement[] =
"                                                                "
" TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      "
"                                                                "
"                                                                ";
/* ................................................................ */
/* @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~. */
/* ................................................................ */
/* ................................................................ */


double score_matrix(Sequence *seq, Pwm *pwm, int start);

int main(int argc, char *argv[]) {
  printf("Hello World!\n");
  Genome *human   = fetch_human_genome();

<<<<<<< scan_transfac_chunk.c
  FILE *matfile   = fopen("/ahg/scr3/cvs/pogvue/data/tf9.4.matrix.pwmline","r");
=======
  FILE *matfile   = fopen("/Users/mclamp/cvs/pogvue/data/tf9.4.matrix.pwmline","r");
>>>>>>> 1.3

  printf("Found matfile %s\n",matfile);
  Pwm   **pwm     = (Pwm **)malloc(800 * sizeof(Pwm *));
  Pwm    *tmppwm;

  int pwmcount = 0;

  while ((tmppwm = read_pwm(matfile))!= NULL) {

    if (tmppwm != NULL) {

      pwm[pwmcount] = tmppwm;
      pwmcount++;

    }
  }

  fclose(matfile);

  printf("Done pwms\n");

  int chunk = 100000;

  double **scores = (double **)malloc(pwmcount*sizeof(double *));
  double *maxes   = (double *)malloc(chunk*sizeof(double)/2500);
  int    *maxpos  = (int *)malloc(chunk*sizeof(int)/2500);

  int p = 0;

  while (p < pwmcount) {
    scores[p] = (double *)malloc(chunk * sizeof(double));

    int pp = 0;
    while (pp < chunk) {
      scores[p][pp] = -100;
      pp++;
    }
    p++;
  }



  // Now let's find the maxes

  int chrnum = 0;

  while (chrnum < HUMAN_CHRNUM) {
    int   len = human->chr[chrnum]->length;
    
    int startcoord = 0;

    while (startcoord < len) {

      int i = 0;

      while (i < chunk*sizeof(double)/2500) {
	maxes[i]  = -100000;
	maxpos[i] = -100000;
	i++;
      }

      p = 0;

      while (p < pwmcount) {
	int pp = 0;
	while (pp < chunk) {
	  scores[p][pp] = -100;
	  pp++;
	}
	p++;
      }

      GenomePosition *pos   = (GenomePosition *)malloc(sizeof(GenomePosition));

      pos->chr = human->chr[chrnum];
      pos->pos = startcoord;

      printf("Reading sequence %d %s\n",startcoord,human->chr[chrnum]->name);

      Sequence *seqptr  = read_fasta_test(pos,chunk);	

      printf("Done\n");

      // Find all the scores first
      int pwmnum = 0;
	
      while (pwmnum < 10) { //pwmcount) {
	int coord = 0;

	while (coord < chunk - 100) {
         // printf("JAMES: Seq len %d\n",seqptr->length);

	  scores[pwmnum][coord] = score_matrix(seqptr,pwm[pwmnum],coord);	  
	  
	  coord++;
	}
	
	coord = 0;
	
	printf("Maxes\n");

	while (coord < chunk - 5000) {
	  int count = (int)(coord/2500);

	  maxes[count] = -100000;
	  int windowcoord      = 0;

	  while (windowcoord < 5000) {

	    if (scores[pwmnum][coord+windowcoord] > maxes[count]) {

	      maxes[count]    = scores[pwmnum][coord+windowcoord];
	      maxpos[count]   = coord + windowcoord;

	      //printf("Setting max %d %d %f\n",count,coord+windowcoord,scores[pwmnum][coord+windowcoord]);
	    }
	    windowcoord++;
	  }

	  printf("Max_score\t%s\t%f\t%d\t%d\t%d\n",pwm[pwmnum]->name,maxes[count],coord,count,maxpos[count]);
	  coord +=2500;

	}

	// Now score for every base

	coord = 0;  

	while (coord < chunk - 5000) {

	  int count = (int)(coord/2500);
	  
	  double score = scores[pwmnum][coord];

	  double maxscore = maxes[count];
	  int    pos      = maxpos[count];

	  if (score +2 > maxscore && maxscore != 0) {
	    printf("Max\t%s\tpos\t%d\tscore\t%f\tmax\t%f\tmaxpos\t%d\t%d\n",pwm[pwmnum]->name,startcoord+coord,score,maxscore,pos,count);
	  }
	  coord++;
	}

	
	pwmnum++;
      }
      free(pos);
      free_seqptr(seqptr);
      startcoord += chunk;
      printf("Coord %d\n",startcoord);
    }
    chrnum++;
  }
return 1;
}
  

double score_matrix(Sequence *seq, Pwm *pwm, int start) {

  Sequence *tmpseq = seq;

  double totscore1=0;
  double totscore2=0;

  double numseqs=0;

  // start - position in seqstring to start
  char *sub = substring(tmpseq->sequence,start,start+pwm->len-1);
  
  //printf("Substring %s\n",sub);
  
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
  
  double log2 = log((float) 2);

  // Weight by information content at each position

  if (count < 3) {
    
    // First the forward strand
    int j = 0;
    double score = 0;
    
    while (j < pwm->len) {
      
      if (sub[j] == 'A') {
	score += log(pwm->vals[j*4]/0.25)/log2;
	
      } else if (sub[j] == 'T') {
	score += log(pwm->vals[j*4+1]/0.25)/log2;
	
      } else if (sub[j] == 'C') {
	score += log(pwm->vals[j*4+2]/0.25)/log2;
	
      } else if (sub[j] == 'G') {
	score += log(pwm->vals[j*4+3]/0.25)/log2;
      }
      
      j++;
    }
    
    totscore1 += score;
    numseqs++;
    
    // Now the reverse strand
    revcomp(sub);
    j = 0;
    score = 0;
    
    while (j < pwm->len) {
      
      if (sub[j] == 'A') {
	score += log(pwm->vals[j*4]/0.25)/log2;
	
      } else if (sub[j] == 'T') {
	score += log(pwm->vals[j*4+1]/0.25)/log2;
	
      } else if (sub[j] == 'C') {
	score += log(pwm->vals[j*4+2]/0.25)/log2;
	
      } else if (sub[j] == 'G') {
	score += log(pwm->vals[j*4+3]/0.25)/log2;
      }
      
      j++;
    }
    
    totscore2 += score;
  }

  free(sub);
  
  int outscore = totscore1;
  
  if (totscore2 > totscore1) {
    outscore = totscore2;
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
