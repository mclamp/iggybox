#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "pog_utils.h"
#include "datamodel.h"

extern char **read_tokens(FILE *file,char token, int *ntok);

Pimer *read_pimer(FILE *file);

Pimer *read_pimer(FILE *file) {
  
  int ntok;
  
  char **tokens;
  
  int count = 0;
  
  
  while ((tokens = read_tokens(file,'\t',&ntok))) {
    Pimer *pimer = NULL;
    
    count++;
    
	  if (ntok == 6) {
	    pimer = (Pimer *)malloc(sizeof(Pimer));
      
	    pimer->chr = (char *)malloc(sizeof(char)*(strlen(tokens[0])+1));
      
		  strcpy(pimer->chr,tokens[0]);
      
		  pimer->start = atoi(tokens[1]);
		  pimer->end   = atoi(tokens[2]);
		  pimer->logodds = (double)atof(tokens[4]);
    }
    int i = 0;
    
    while (i < ntok) {
      free(tokens[i]);
      i++;
    }
    
    free(tokens);
    
    if (pimer != NULL) {
      return pimer;
    }
  }
  return NULL;
}

char *get_pi_filename(GenomePosition *pos) {
  char *pidir    = PI_RAW;
  //char *pidir    = "/Volumes/rg/genome/2x_aln/estimation/mammals21/runs/pi_ar_model/";
  //char *pidir    = "/ahg/scr3/genome/2x_aln/estimation/mammals21/runs/pi_ar_model/";
  
  int start = 1000000* (int)(pos->pos/1000000);
  int end   = start+999999;
  
  if (end > pos->chr->length) {
    end = pos->chr->length;
  }
  
  char startbuf[20];
  char endbuf[20];
  
  sprintf(startbuf, "%d", start);
  sprintf(endbuf,   "%d", end);  
  
  
  char *pifilename = (char *)malloc((strlen(pidir)     + strlen(pos->chr->name) + strlen("/") + strlen(pos->chr->name) + strlen("_") 
                                     + strlen(startbuf) + strlen("-")  + strlen(endbuf) + strlen(".pi") + 1)*sizeof(char));
  
  
  strcpy(pifilename,pidir);
  strcat(pifilename,pos->chr->name);
  strcat(pifilename,"/");
  strcat(pifilename,pos->chr->name);
  strcat(pifilename,"_");
  strcat(pifilename,startbuf);
  strcat(pifilename,"-");
  strcat(pifilename,endbuf);
  strcat(pifilename,".pi");
  
  return pifilename;
}

Pi **read_pis_test(GenomePosition *pos, int len) {
  
  Pi **pis = (Pi **)malloc((len+1)*sizeof(Pi *));
  Pi  *pi  = NULL;
  
  char *filename = get_pi_filename(pos);
  
  int offset = 1000000 * (int)(pos->pos / 1000000);
  
  printf("Filename %s %d\n",filename,offset);
  
  FILE     *file   = fopen(filename,"r");
  
  if (file == NULL) {
    printf("Filename for %s %s %llu is not found\n",filename,pos->chr->name,pos->pos);
    return NULL;
  }
  
  int start = pos->pos;
  int end   = pos->pos + len - 1;
  
  int prev = start - 1;
  
  while ((pi = read_pi(file)) != NULL  && prev < end) {
    
    pi->start = pi->start + offset;
    pi->end   = pi->end   + offset;
    
    
    if (pi->start >= start && pi->start <= end) {
      
      if (pi->start > prev) {
        int i = prev+1;
        
        while (i < pi->start) {
          pis[i - start] = NULL;
          
          //printf("Allocating1 for %d %d\n",i,prev);
          
          prev = i;
          i++;
        }
      } else {
        printf("Skipping %d\n",pi->start - start);
      }
      pis[pi->start - start] = pi;
      
      //printf("Allocating2 for %d %d\n",pi->start - start,prev);
      
      prev = pi->start;
      
      pi->start -= offset;
      pi->end   -= offset;
      
    } else {
      free(pi->chr);
      free(pi);
    }
  }
  
  int i = prev+1;
  
  while (i <= end) {
    //printf("Allocating3 for %d %d\n",i-start,i);
    pis[i - start] = NULL;
    prev = i;
    i++;
  }
  pis[end+1-start] = NULL;
  
  fclose(file);
  free(filename);
  
  return pis;
  
}

char *get_pimer_filename(GenomePosition *pos) {
  char   *pidir    = PI_MER;
  
  //char *pidir    = "/Volumes/rg 1/21mammals/pi/12mer/";
  //char *pidir = "/Volumes/rg/split_tree/21mammals.50mer/";
  //char *pidir    = "/ahg/scr3/split_tree/21mammals.50mer/";
  //char *pidir    = "/ahg/scr3/split_tree/HMRD_AR_aligns.50mer/";
  //char *pidir    = "/ahg/scr3/split_tree/HMRD.50mer/";
  //char *pidir    = "/seq/mgscratch/21mammals/pi/12mer/";
  
  int start = 1000000* (int)(pos->pos/1000000);
  int end   = start+999999;
  
  if (end > pos->chr->length) {
    end = pos->chr->length;
  }
  
  char startbuf[20];
  char endbuf[20];
  
  sprintf(startbuf, "%d", start);
  sprintf(endbuf,   "%d", end);  
  
  
  char *pifilename = (char *)malloc((strlen(pidir)     + strlen(pos->chr->name) + strlen("/") + strlen(pos->chr->name) + strlen("_") 
                                     + strlen(startbuf) + strlen("-")  + strlen(endbuf) + strlen(".pi.") + strlen(MER) + 1)*sizeof(char));
  
  
  strcpy(pifilename,pidir);
  strcat(pifilename,pos->chr->name);
  strcat(pifilename,"/");
  strcat(pifilename,pos->chr->name);
  strcat(pifilename,"_");
  strcat(pifilename,startbuf);
  strcat(pifilename,"-");
  strcat(pifilename,endbuf);
  strcat(pifilename,".pi.");
  strcat(pifilename,MER);
  
  return pifilename;
}

Pimer **read_pimers_test(GenomePosition *pos, int len) {
  
  Pimer **pis = (Pimer **)malloc((len+1)*sizeof(Pimer *));
  Pimer  *pi  = NULL;
  
  char *filename = get_pimer_filename(pos);
  
  int offset = 1000000 * (int)(pos->pos / 1000000);
  
  //char *filename = "/Users/mclamp/cvs/pogview/data/chr1_120000000-120999999.fa.pi.12mer";
  
  printf("Filename %s %d\n",filename,offset);
  
  FILE     *file   = fopen(filename,"r");

  if (file == NULL) {
    printf("Filename for %s %s %llu is not found\n",filename,pos->chr->name,pos->pos);
    return NULL;
  }
  
  

  int start = pos->pos;
  int end   = pos->pos + len - 1;

  int prev = start - 1;
  
  while ((pi = read_pimer(file)) != NULL  && prev < end) {
    
    pi->start = pi->start + offset;
    pi->end   = pi->end   + offset;
    

    if (pi->start >= start && pi->start <= end) {
      
      if (pi->start > prev) {
        int i = prev+1;
        
        while (i < pi->start) {
          pis[i - start] = NULL;
          
          //printf("Allocating1 for %d %d\n",i,prev);
          
          prev = i;
          i++;
        }
      } else {
        printf("Skipping %d\n",pi->start - start);
      }
      pis[pi->start - start] = pi;
      
      //printf("Allocating2 for %d %d\n",pi->start - start,prev);
      
      prev = pi->start;
      
      pi->start -= offset;
      pi->end   -= offset;
      
    } else {
      free_pimer(pi);
    }
  }
  
  int i = prev+1;
  
  while (i <= end) {
    //printf("Allocating3 for %d %d\n",i-start,i);
    pis[i - start] = NULL;
    prev = i;
    i++;
  }
  pis[end+1-start] = NULL;
  
  fclose(file);
  free(filename);
  
  return pis;
  
}

void free_pimers_test(Pimer **pis, int len) {
  int i = 0;
 
  if (pis == NULL) {
    return;
  }
  while (i < len) {
    if (pis[i] != NULL) {
      free_pimer(pis[i]);
    }
    i++;
  }
  free(pis);
}

void free_pis_test(Pi **pis, int len) {
  int i = 0;
  
  if (pis == NULL) {
    return;
  }
  while (i < len) {
    if (pis[i] != NULL) {
      free(pis[i]->chr);
      free(pis[i]);
    }
    i++;
  }
  free(pis);
}



void free_pimer(Pimer *pimer) {
  if (pimer->chr != NULL) {
    free(pimer->chr);
  }
  free(pimer);
}

void print_pi(Pi *pi) {
	printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t1\n",pi->start,pi->aval,pi->cval,pi->gval,pi->tval,pi->logodds,pi->tree_length);
}

Pi *read_pi(FILE *file) {
	
	int ntok;
	
	char **tokens;
	Pi    *pi = NULL;
  
	int count = 0;
	
	while ((tokens = read_tokens(file,'\t',&ntok))) {
		count++;
		if (ntok == 8) {
		  pi = (Pi *)malloc(sizeof(Pi));
			pi->chr   = '\0';
			pi->start = atoi(tokens[0]);
			pi->end   = pi->start;
			pi->aval  = atof(tokens[1]);
			pi->cval  = atof(tokens[2]);
			pi->gval  = atof(tokens[3]);
			pi->tval  = atof(tokens[4]);
			pi->logodds = (double)atof(tokens[5]);
			pi->tree_length = (double)atof(tokens[6]);
    }
    
    if (tokens != NULL) {
			int i = 0;
      while (i < ntok) {
        free(tokens[i]);
        i++;
      }
      free(tokens);
    }
    return pi;
	}
	return NULL;
}
