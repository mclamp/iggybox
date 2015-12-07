#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "datamodel.h"

char **read_tokens(FILE *file,char token, int *ntok) {
  char   c;
  
  int    toklen = 0;
  int    toknum;
  
  char **tokens   = NULL;
  char prev       = (char)NULL;
  
  int    toknumchunk = 200;
  int    toklenchunk = 200;
  
  int    curr_toknum = toknumchunk;
  int    curr_toklen = toklenchunk;
  
  tokens = (char **)malloc(toknumchunk*sizeof(char *));
  toknum = 0;
  
  while ((c = fgetc(file)) != EOF) {
    //printf("New char %c\n",c);
    if (c == '\n') {
      *ntok = toknum;
      return tokens;
    }
    
    if (c != token) {
      // Start of line or new token
      if (prev == (char)NULL || prev == token) {

        // New token
        if (toknum > curr_toknum) {
	  // This has a problem!!!
	  printf("New token %d\n",toknum);
          curr_toknum += toknum;
          tokens = (char **)realloc(tokens,(curr_toknum+1)*sizeof(char *));

        }

        char *newtok = (char *)malloc(toklenchunk*sizeof(char));
        
        newtok[0] = c;
        newtok[1] = '\0';
        
        toknum++;
        toklen = 1;
        
        tokens[toknum-1] = newtok;
        
      } else {        // Continue existing token
        toklen++;
	//printf("Token %s %d %d\n",tokens[toknum-1],toknum,curr_toklen);
        if (toklen > curr_toklen) {
          curr_toklen += toklen;
          tokens[toknum-1] = (char *)realloc(tokens[toknum-1],(curr_toklen + 1)*sizeof(char));
        }
        tokens[toknum-1][toklen-1] = c;
        tokens[toknum-1][toklen] = '\0';

      }
    }
    prev = c;
  }
  free(tokens);
  return NULL;
}



Gff *read_gff(FILE *file) {
  int ntok;
  char **tokens = read_tokens(file,'\t', &ntok);

  
  if (tokens != NULL) {
    Gff *tmp = (Gff *)calloc(1,sizeof(Gff));

    tmp->seqname = (char *)malloc((strlen(tokens[0])+1)*sizeof(char));
    tmp->source  = (char *)malloc((strlen(tokens[1])+1)*sizeof(char));
    tmp->feature = (char *)malloc((strlen(tokens[2])+1)*sizeof(char));

    strcpy(tmp->seqname,tokens[0]);
    strcpy(tmp->source, tokens[1]);
    strcpy(tmp->feature,tokens[2]);

    tmp->start   = atoi(tokens[3]);
    tmp->end     = atoi(tokens[4]);
    tmp->score   = atof(tokens[5]);
    
    if (tokens[6][0] == '+') {
      tmp->strand = 1;
    } else if (tokens[6][0] == '-') {
      tmp->strand = -1;
    } else if (tokens[6][0] != '.') {
      tmp->strand  = atoi(tokens[6]);
    }
    
    tmp->phase    = (int)NULL;
    tmp->hseqname = (char *)NULL;
    tmp->hstart   = (int)NULL;
    tmp->hend     = (int)NULL;
    tmp->hstrand  = (int)NULL;
    
    if (tokens[7][0] != '.') {
      //tmp->phase = atoi(tokens[7]);
    }
    
    if (ntok > 8) {
      tmp->hseqname = (char *)malloc((strlen(tokens[8])+1)*sizeof(char));
      strcpy(tmp->hseqname,tokens[8]);
    }
    
    if (ntok > 10) {
      tmp->hstart   = atoi(tokens[9]);
      tmp->hend     = atoi(tokens[10]);
    }
    if (ntok > 11) {
      tmp->hstrand  = atoi(tokens[11]);
    }
    tmp->next = NULL;

    int i = 0;

    while (i < ntok) {
	free(tokens[i]);
	i++;
    }
    free(tokens);
    return tmp;
  }
  
  return NULL;
}

Gff *read_genome_position_gff(GenomePosition *pos, int len) {
  
  char command[300];
  char outfile[50];
  
  int start = (int)pos->pos;
  int end   = start + len - 1;
  
  snprintf(outfile,50,"/tmp/bsearch.%s.%d-%d.out",pos->chr->name,start,end);
  
  //snprintf(command,300,"perl /ahg/scr3/cvs/twox/bsearch_genome_file -file  /ahg/scr3/mammals/local_data/human.31_35d.exon.gff.coordsort -chr %s -start %d -end %d -startcol 3 -chrcol 0 > %s\n",
  snprintf(command,300,"perl /ahg/scr3/cvs/twox/bsearch_genome_file -file /ahg/scr3/mammals/ucsc/rmsk.gff.coordsort -chr %s -start %d -end %d -startcol 3 -chrcol 0 > %s\n",
  //snprintf(command,300,"perl /Volumes/rg/cvs/twox/bsearch_genome_file -file /Volumes/rg/mammals/ucsc/rmsk.gff.coordsort -chr %s -start %d -end %d -startcol 3 -chrcol 0 > %s\n",
           pos->chr->name,start,end,
           outfile);
  
  printf("Command is %s\n",command);
  
  int status = system(command);
  
  if (status == 0) {
    FILE *file = fopen(outfile,"r");
    if (file != NULL) {
      Gff  *gff = read_gff(file);
      fclose(file);
      remove(outfile);
      return gff;
    } else {
      remove(outfile);
      return NULL;
    }
  } else {
    remove(outfile);
    return NULL;
  }
}

void free_gffptr(Gff *gff) {
  
  while (gff != NULL) {
    free(gff->seqname);
    free(gff->source);
    free(gff->feature);
    free(gff->hseqname);
    
    Gff *tmpgff = gff->next;
    
    free(gff);
    
    gff = tmpgff;
  }
}
void print_gff(FILE *file,Gff *gff) {
  
  fprintf(file,"%s\t%s\t%s\t%d\t%d\t%4.2f\t%d",
          gff->seqname,
          gff->source,
          gff->feature,
          gff->start,
          gff->end,
          gff->score,
          gff->strand);
  
  if (gff->phase != (int)NULL) {
    fprintf(file,"\t.\t");
  } else {
    fprintf(file,"\t%d",gff->phase);
  }
  
  if (gff->hseqname != NULL) {
    fprintf(file,"\t%s",gff->hseqname);
    
    if (gff->hstart != (int)NULL) {
      fprintf(file,"\t%d\t%d",gff->hstart,gff->hend);
      
      if (gff->hstrand != (int)NULL) {
        fprintf(file,"\t%d",gff->hstrand);
      }
    }
  }
  fprintf(file,"\n");
}

int gffmain(int argc,char *argv[]) {
  FILE *file;
  
  Gff *head;
  
  file = fopen(argv[1],"r");
  
  head = read_gff(file);
  
  while (head != NULL) {
    print_gff(stdout,head);
    head = head->next;
  }
  return 1;
}

Pwm *read_pwm(FILE *file) {
  int    ntok;

  char **tokens = read_tokens(file,'\t',&ntok);
  
  Pwm *pwm = NULL;

  if (tokens != NULL) {
    pwm = (Pwm *)malloc(sizeof(Pwm));

    pwm->name = tokens[0];

    int i = 3;

    pwm->vals = (double *)malloc((ntok-3)*sizeof(double));
    pwm->len = (ntok-3)/4;

    while (i < ntok) {

      pwm->vals[i-3] = atof(tokens[i]);

      if (pwm->vals[i-3] < 0) {
	printf("EEEEEK %d\t%f\n",i-3,pwm->vals[i-3]);
      }
      if (pwm->vals[i-3] == 0) {
        pwm->vals[i-3] = 1e-6;  
      }
      i++;
    }

    // Weight by information content at each position

    double *inf = (double *)malloc(pwm->len*sizeof(double));

    i = 0;
  
    while (i < pwm->len) {

      int j = 0;
      double infval = 0.0;

      while (j < 4) {
	//printf("Val is %f\t%f\n",pwm->vals[i*4+j],log(pwm->vals[i*4+j]));
	infval += pwm->vals[i*4+j]*log(pwm->vals[i*4+j])/log(2);
	j++;
      }
      
      infval += 2;
      
      inf[i] = infval;
      
      //printf("Inf is %f\n",infval);
      i++;
    }
    pwm->inf = inf;

    i = 1;

    while (i < ntok) {
      free(tokens[i]);
      i++;
    }
    free(tokens);

  }
  return pwm;

}


void free_pwm(Pwm *pwm) {
  free(pwm->name);
  free(pwm->vals);
  free(pwm);
}
  
