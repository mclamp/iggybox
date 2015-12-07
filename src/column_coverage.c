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

    Sequence *seq;
    Sequence *seqptr;

    int numseqs = 0;
    int i       = 0;
    int count   = 0;
    int numcols = 0;
    int length;

    int *cov;

    char c;
    int  col;
    int  *done;

    seqptr = read_fasta(argv[1]);

    seq = seqptr;

    while (seq != NULL) {
      numseqs++;
      seq = seq->next;
    }

    length = seqptr->length;

    numcols = length;

    cov  = (int *)malloc(numseqs*sizeof(int));

    i = 0;

    while (i < numseqs) {
      cov[i] = 0;
      i++;
    }
	
    i = 0;
	
	while (i < numcols) {
		int mouse = 0;
		int dog   = 0;
		int cow   = 0;
		
		
		seq    = seqptr;
		
		int tmp   = 0;
		int count = 0;
		
		char *species = (char *)NULL;
		
		while (count < numseqs) {
			
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
				if (strcmp(seq->id,"Chimp") != 0 &&
					strcmp(seq->id,"Cow")  != 0 &&
					strcmp(seq->id,"RheMac") != 0 &&
					strcmp(seq->id,"Rat") != 0 &&
					strcmp(seq->id,"Mouse7") != 0 &&
					strcmp(seq->id,"Dog2") &&
					strcmp(seq->id,"Human") != 0)  {
					
				    if (species == NULL) {
						species = (char *)malloc((strlen(seq->id)+1)*sizeof(char));
						strcpy(species,seq->id);
					} else {
						species = (char *)realloc(species,(strlen(species) + strlen(seq->id) + 1)*sizeof(char));
						strcat(species,seq->id);
					}
				  tmp++;
				}
			}
			seq    = seq->next;
			
			count++;	
		}
		if (mouse == 1 && dog == 1 && cow == 1) {
			cov[tmp]++;	
			printf("Species %d %s\n",tmp,species,);
		}
		i++;
    }
    i = 0;

    while (i < numseqs) {
      printf("Cov\t%d\t%d\n",i,cov[i]);
      i++;
    }

}

