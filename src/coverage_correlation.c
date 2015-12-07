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
extern Pi       *read_pi    (FILE *file);
extern void      print_pi   (Pi *pi);

extern SeqArray   *read_seq_array(char *filename);
extern PimerArray *read_pimer_array(char *filename, int numcols);

int *find_coverage(SeqArray *seqs ,PimerArray *pis);
int *coverage_pis (SeqArray *seqs, PimerArray *pis);

int main(int argc, char * argv[]) {
  int i;


  PimerArray *pimer_array = read_pimer_array(argv[2],1000000);
    
  fprintf(stderr,"Got pis %d\n",pimer_array->len);

  SeqArray   *seq_array   = read_seq_array(argv[1]);


  int numseqs = seq_array->len;
  int numcols = seq_array->seq[0]->length;

  fprintf(stderr,"Got seqs %d %d\n",numseqs,numcols);
    
  int *cov = find_coverage(seq_array, pimer_array);
	
  i = 0;
   
  Pimer **pis = pimer_array->pimer;

  while (i < numcols) {
    Pimer *tmp = pis[i];
    tmp->start = i;
    tmp->end   = i;
    printf("chr1\tpi\tpi\t%d\t%d\t%d\t1\t.\t%f\n",tmp->start,tmp->end,cov[i],tmp->logodds);
    //print_pi(tmp);
    //printf("Original %f %d\n",tmp->chr,pis[i]->tree_length,cov[i]);
  }
  i++;
}

int *find_coverage(SeqArray *seqs ,PimerArray *pis) {
  int numseqs = seqs->len;
  int numcols = seqs->seq[0]->length;
  
  printf("Num %d %d\n",numseqs,numcols);

  int *cov;

  cov = (int *)malloc(numcols*sizeof(int));

  printf("Cov done %d\n",cov);

  int i = 0;
  
  while (i < numcols) {
    cov[i] = 0;
    i++;
  }
  
  i = 0;

  
  while (i < numcols) {
    Sequence *seq    = seqs->seq[i];
    
    printf("Seq2 %d\n",seq->length);
    printf("Seq2 %c\n",seq->sequence[0]);
    printf("Pimer %d - %d\n",i,pis->pimer[i]->start);
    if (pis->pimer[i] != (Pimer *)NULL) {
      int tmp   = 0;
      int count = 0;
      
      char *species = (char *)NULL;
      
      int mouse = 0;
      int dog   = 0;
      int cow   = 0;
      int hedge = 0;

      //printf("I %d\n",i);
      while (count < numseqs) {
	//	printf("Count %d\n",count);
	if (seq->sequence[i] != '-') { 
	  if (strcmp(seq->id ,"Mouse7") == 0) {
	    mouse = 1;
	  }
	  if (strcmp(seq->id, "Dog2") == 0) {
	    dog = 1;
	  }
	  if (strcmp(seq->id, "Cow") == 0) {
	    cow = 1;
	  }
	  if (strcmp(seq->id,"Tenrec") == 0) {
	    hedge = 1;
	  }
	  if (strcmp(seq->id,"Chimp") != 0 &&
	      strcmp(seq->id,"Cow")  != 0 &&
	      strcmp(seq->id,"RheMac") != 0 &&
	      strcmp(seq->id,"Rat") != 0 &&
	      strcmp(seq->id,"Mouse7") != 0 &&
	      strcmp(seq->id,"Dog2") &&
	      strcmp(seq->id,"Human") != 0)  {
	    
	    if (species == NULL) {
	      
	      species = (char *)malloc(20*(strlen(seq->id)+1)*sizeof(char));
	      strcpy(species,seq->id);
	      //		printf("Malloc %s\n",species);
	    } else {
	      //			printf("Realloc %s : %s\n",species,seq->id);
	      //		species = (char *)realloc(species,(strlen(species) + strlen(seq->id) + 1)*sizeof(char));
	      strcat(species,seq->id);
	    }
	  }
	    //		printf("Species %s\n",species);
	  tmp++;
	  
	}
	
	count++;	
      }
      if (mouse == 1 && dog == 1 && cow == 1 && hedge == 1) {
	//cov[tmp]++;	
	//printf("Species %d %s\n",tmp,species);
      }
      
      cov[i] = tmp;
    //  covcount[tmp]++;
    }
    i++;
  }

  return cov;
}
