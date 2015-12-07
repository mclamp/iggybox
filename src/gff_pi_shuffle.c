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

int *find_coverage(Sequence *seqptr, int numseqs, int numcols, int *covcount,Pi **pis);

Pi **shuffle_pis(Pi **pis, Sequence* seqptr, int numseqs, int numcols);


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
	
    fprintf(stderr,"Got seqs %d %d\n",numseqs,numcols);
    
    FILE *pifile = fopen(argv[2],"r");
	
    Pi **pis = (Pi **)malloc(numcols*sizeof(Pi *));
    Pi **newpis = (Pi **)malloc(numcols*sizeof(Pi *));
	
    i = 0;
    
    while (i < numcols) {
		pis[i] = NULL;
		newpis[i] = NULL;
		i++;
    }
    
    Pi *pi;
    int count = 0;
    int maxcol = 0;
	
    while ((pi = read_pi(pifile)) != (Pi *)NULL) {
      pis[pi->start] = pi;
      count++;
      if (pi->end > maxcol) {
	maxcol = pi->end;
      }
    }
	
    fprintf(stderr,"Got pis %d %d\n",maxcol,count);
	
    
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
	if (pis[k] != NULL) {
	  newpis[k] = pis[k];
	  //printf("%s\tpimer\tpimer\t%d\t%d\t%f\t1\t.\n",gff->seqname,k,k,pis[k]->logodds);
	}
	k++;
	
      }
      
      gff = gff->next;
    }
    
	
    Pi **shuffled_pis = shuffle_pis(newpis, seqptr, numseqs, numcols);
	
    i = 0;
    
    while (i < numcols) {
      if (shuffled_pis[i] != NULL) {
	Pi *tmp = (Pi *)shuffled_pis[i];
	tmp->start = i;
	tmp->end   = i;
	//	fprintf(stderr,"%s\t%d\t%d\tw1\t%f\t+\n", tmp->chr,tmp->start,tmp->end,tmp->logodds);
	print_pi(tmp);
	//printf("Original %f\n",pis[i]->tree_length);
      }
      i++;
    }
}

Pi **shuffle_pis(Pi **pis, Sequence* seqptr, int numseqs,int numcols) {
  int **covpos    = (int **)malloc(numseqs*sizeof(int *));    //  list of cols with n sequences
  int **done      = (int **)malloc(numseqs*sizeof(int *));    //  stores the taken columns for each coverage no.
  
  Sequence **seqarr = (Sequence **)malloc(numseqs*sizeof(Sequence *));
  
  int *covcount  = (int *)malloc(numseqs*sizeof(int));       //  number of cols with n sequences
  int *cov       = (int *)malloc(numcols*sizeof(int));       //  number of seqs at each column 
  int *tmppos    = (int *)malloc(numseqs*sizeof(int));       //  End position of each coverage levels cov storarge
  
  Sequence *seq = seqptr;
  
  Pi **shuffled_pis = (Pi **)malloc(numcols*sizeof(Pi *));
  
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
    shuffled_pis[i] = (Pi *)NULL;
    i++;
  }

  cov = find_coverage(seqptr,numseqs,numcols,covcount,pis);       // Index is column number
  
  i = 0;
	
  while (i < numseqs) {
    covpos[i] = (int *)malloc(covcount[i]*sizeof(int));
    done[i]   = (int *)malloc(covcount[i]*sizeof(int));
    fprintf(stderr,"Covcount for %d is %d\n",i,covcount[i]);
    i++;
  }
	
  // Fill in the coverage arrays
  
  i = 0;
  
  while (i < numcols) {
    int seqs = cov[i];
    //printf("Seqs %d %d\n",seqs,i);
    int pos  = tmppos[seqs];
    if (pis[i] != NULL) {
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
  
  // Now loop over each column 
  
  while (i < numcols) {
    
    // Find coverage at column i
    int seqs = cov[i];                                 
    
    if (seqs < 23 && covcount[seqs] > 0) {
      //      printf("Seqs %d\n",seqs);
      // Pick a rand column newcol from 0 to covcount[seqs]
      int newcol = (int)((double)covcount[seqs]*rand()*1.0/(1.0*RAND_MAX));
      
      while (done[seqs][newcol] == 1) {
	newcol = (int)((double)covcount[seqs]*rand()*1.0/(1.0*RAND_MAX));
      }
      
      int newcol2 = covpos[seqs][newcol];
      
      double score = pis[newcol2]->logodds;
      shuffled_pis[i] = pis[newcol2];
      
      done[seqs][newcol] = 1;
      
      if (i % 10000 == 0) {
	fprintf(stderr,"Score is %d %d %f %d\n",i,newcol2,score,seqs);
	fprintf(stderr,"Found col %d %d\n",newcol,i);
	fprintf(stderr,"Pi %f\n",shuffled_pis[i]->logodds);
      }
    }
    i++;
  }
	
	
  i = 0;
  
  return shuffled_pis;
}

int *find_coverage(Sequence *seqptr, int numseqs, int numcols,int *covcount, Pi **pis) {
  
  int *cov = (int *)malloc(numcols*sizeof(int));
  
  int i = 0;
  
  Sequence *seq;
  
  while (i < numcols) {
    cov[i] = 0;
    i++;
  }
  
  i = 0;
  
  while (i < numcols) {
    seq    = seqptr;
    
    
    if (pis[i] != NULL) {
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
	seq    = seq->next;
	
	count++;	
      }
      if (mouse == 1 && dog == 1 && cow == 1 && hedge == 1) {
	//cov[tmp]++;	
	//printf("Species %d %s\n",tmp,species);
      }
      
      cov[i] = tmp;
      covcount[tmp]++;
    }
    i++;
  }

  return cov;
}
