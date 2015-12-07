#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "datamodel.h"

Pwm *read_pwm(FILE *file);
void revcomp(char *s);
Thresh *read_threshfile();
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


double score_matrix(Sequence *seq, Pwm *pwm, int start, int offset, char *chr,FILE *file, int gffstart, int threshold,int num);

int main(int argc, char *argv[]) {
  Genome *human   = fetch_human_genome();

  //FILE *matfile   = fopen("/home/mclamp/iccpog/src/tf9.4.matrix.pwmline","r");
  FILE *matfile   = fopen("/ahg/scr3/cvs/pogvue/data/tf9.4.matrix.pwmline","r");
  FILE *gfffile   = fopen(argv[1],"r");

  Pwm   **pwm     = (Pwm **)malloc(800 * sizeof(Pwm *));
  Pwm    *tmppwm;

  int pwmcount = 0;

  Thresh *thresh = read_threshfile();

  int offset = atoi(argv[2]);
  int chunk  = atoi(argv[3]);
  char *stub = argv[4];

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



  GenomePosition *pos   = (GenomePosition *)malloc(sizeof(GenomePosition));
  Sequence       *seqptr = NULL;
  int             count = 0;
  Gff            *gff;


  int bin;
  int prevbin;

  char *prevchr = (char *)malloc(20*sizeof(char));

  while ((gff = read_gff(gfffile)) != NULL) {

    if (count >= offset && count < offset+chunk) {
      printf("GFF\t%s\t%d\n",gff->hseqname,count);
      
      bin = (int)(gff->start/1000000);

      Chromosome     *chr   = find_chromosome_by_name(human,gff->seqname);
      
      if (chr != NULL && (int)(gff->end/1000000) == bin) {

	pos->chr = chr;
	pos->pos = 1000000*bin;
	
	if (seqptr  == NULL ||
	    prevbin != bin  ||
	    strcmp(gff->seqname,prevchr)!=0) {
	  
	  if (seqptr != NULL) {
	    printf("New seq chunk %s %s %d %d\n",gff->seqname,prevchr,bin,prevbin);
	    free_seqptr(seqptr);	
	  } 
          int len = 1000000;
	  if (pos->pos + len > chr->length) {
	     len = chr->length - pos->pos;
	  }

	  seqptr  = read_fasta_test(pos,len);	

	}

	int i = 0;
	
	while (i < pwmcount) {
	  int threshold = 0;
	  int k = 0;
	  while (k < thresh->count) {

	    if (i == thresh->num[k]) {
	      threshold = thresh->thresh[k];
	    }
	    k++;
	  }
	
	  
	  Pwm *tmppwm = pwm[i];

	  int j = gff->start - bin*1000000;
	  
	  int end = gff->start - bin*1000000 + (gff->end-gff->start+1);// - tmppwm->len;
	  

	  while (j < end) {
	    double score = score_matrix(seqptr,tmppwm,j,bin*1000000,gff->seqname,outfiles[i],gff->start-1000000*bin,threshold,i);
	    j++;
	    //j = end;
	  }

	  i++;
	  
	}
      }

      prevbin = bin;
      strcpy(prevchr,gff->seqname);

    }

    free(gff->seqname);
    free(gff->source);
    free(gff->feature);
    free(gff->hseqname);
    free(gff);
    count++;

  }
  free(pos);
}

double score_matrix(Sequence *seq, Pwm *pwm, int start, int offset, char *chr,FILE *file,int gffstart, int threshold, int num) {

  Sequence *tmpseq = seq;

  double totscore1=0;
  double totscore2=0;

  double numseqs=0;

  while (tmpseq != NULL) {
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

    double log2 = log(2);

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
//      fprintf(file,"F\t%s\t%s\t%d\t%s\t%7.2f\t%d\t%s\n",tmpseq->id,chr,start+offset,pwm->name,score,start-gffstart+1,sub);
      
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
 //     fprintf(file,"R\t%s\t%s\t%d\t%s\t%7.2f\t%d\t%s\n",tmpseq->id,chr,start+offset,pwm->name,score,start-gffstart+1,sub);

      totscore2 += score;


	
    
    }

    free(sub);

    if (numseqs >1  || (totscore1 >= threshold || totscore2 >= threshold)){
      tmpseq = tmpseq->next;
    } else {
      tmpseq = NULL;
    }

  } 

  if (totscore1/numseqs >= threshold) {
    fprintf(file,"AVF\t%s\t%d\t%d\t%s\t%7.2f\t%f\t%d\n",chr,start+offset,num,pwm->name,(totscore1/numseqs),numseqs,threshold);    
    printf("AVF\t%s\t%d\t%d\t%s\t%7.2f\t%f\t%d\n",chr,start+offset,num,pwm->name,(totscore1/numseqs),numseqs,threshold);    
  }
  if (totscore2/numseqs >= threshold) {
    fprintf(file,"AVR\t%s\t%d\t%d\t%s\t%7.2f\t%f\t%d\n",chr,start+offset,num,pwm->name,(totscore2/numseqs),numseqs,threshold);    
    printf("AVR\t%s\t%d\t%d\t%s\t%7.2f\t%f\t%d\n",chr,start+offset,num,pwm->name,(totscore2/numseqs),numseqs,threshold);    
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
    
      
void revcomp(char *s) {
  char *p = s + strlen(s) - 1;

  while (s <= p) {
    char c = dna_complement[*s];
    *s = dna_complement[*p];
    *p = c;
    ++s;
    --p;
  }
}

Thresh *read_threshfile() {

  char *threshfile = "/ahg/scr3/cvs/pogview/src/scan_transfac.all.thresh";

  FILE *file = fopen(threshfile,"r");

  if (file == NULL) {
    printf("File %s not found\n",threshfile);
    exit(0);
  }
  char **tokens;
  int   ntok;

  Thresh *thresh = (Thresh *)malloc(sizeof(Thresh));

  int count = 0;
  int chunk = 100;

  thresh->ids    = (char **)malloc(chunk*sizeof(char*));
  thresh->thresh = (int *)malloc(chunk*sizeof(int));
  thresh->num    = (int *)malloc(chunk*sizeof(int));
  
  while ((tokens = read_tokens(file,'\t',&ntok)) != NULL) {
    if (count > chunk) {
      chunk = 2*chunk;
      thresh->ids    = (char **)realloc(thresh->ids,chunk*sizeof(char *));
      thresh->thresh = (int *)realloc(thresh->thresh,chunk*sizeof(int));
      thresh->num    = (int *)realloc(thresh->num,chunk*sizeof(int));
    }

    thresh->ids[count] = tokens[0];
    thresh->thresh[count] = atoi(tokens[2]);
    thresh->num[count]    = atoi(tokens[1]);

    printf("Thresh %s\t%d\t%d\n",thresh->ids[count],thresh->thresh[count],thresh->num[count]);
    free(tokens[1]);
    free(tokens[2]);

    free(tokens);

    count++;
  }

  thresh->count = count;

  return thresh;

  fclose(file);
}
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  







