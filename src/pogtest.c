#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "viewport.h"
#include "seqio.h"

#define _ (-1)

#define MAXCOMBSIZE 50

static void print_l (FILE *fp, unsigned long x,int n);

  static void print_ll(FILE *fp, unsigned long long x,int m) 
  {
  int i;
  for (i = (2*(m-1)); i >= 0; i -= 2) {
  //	fprintf(fp, "%d ",i);
  fprintf(fp, "%c", "ACGT"[(x>>i)&3]);
  //	fprintf(fp, "X0 = %u ",x>>0 & 3);
  }
  }


static void print_l (FILE *fp, unsigned long x,int n)
{
    int i;
    for (i = (2*(n-1)); i >= 0; i -= 2) {
      //	fprintf(fp, "%d ",i);
	fprintf(fp, "%c", "ACGT"[(x>>i)&3]);
    }
}

Hash **read_fasta_and_hash(char *filename, int n, int m,Comb *comb) {

  int i;

  int e;

  FILE *file;

  char *nmerstring;

  static unsigned long long ei = 0;
  static unsigned int ecode;

  struct rusage *res;
  int           *bits;

  char           c;

  Hash          **hash;
  Hash          *tmphash;

  hash    = (Hash **)malloc(nmer2pow(n) * sizeof (Hash *));

  printf("Size of int %ld\n"      ,sizeof(unsigned int));
  printf("Size of long %ld\n"     ,sizeof(unsigned long));
  printf("Size of long long %ld\n",sizeof(unsigned long long));

  res = (struct rusage *)malloc(sizeof(struct rusage));

  getrusage(0,res);

  printf("Time %d %d\n",res->ru_utime.tv_sec,res->ru_utime.tv_usec);

  bits = (int *)malloc(sizeof(int)*19);
  bits = generate_bitstring(n,m);

  memset(hash, 0, nmer2pow(n)*sizeof(Hash *));

  nmerstring = (char *)malloc((m+1) * sizeof(char));

  ei = 0;
  
  file = fopen(filename,"r");
  
  i = 0;

  while (i < m-1) {
    c = getc(file);


    if (c == EOF) {
      printf("Got end of file %s\n",filename);
      return hash;
    }

    if (c == '\n' || c == ' ') {
      while ((c == getc(file)) != EOF && (c != ' ' && c != '\n')) {
      }
    }
    
    if (c == '>') {
      while ((c = getc(file)) != EOF && c != '\n') {
	printf("%c",c);
      }
      printf("\n");
      c = getc(file);
    } 
    printf("Char %c\n",c);
    e = fasta_encoding[(int)c];
    ei = (ei << 2) | e;
    i++;
  }

  if (i < 25) {
    print_ll(stdout,ei,m);
    printf("\n");
  }
  getrusage(0,res);
  
  printf("Time before hash %d %d\n",res->ru_utime.tv_sec,res->ru_utime.tv_usec);
  
  i = m-1;

  while (( c = getc(file)) != EOF) {
    if (c != ' ' && c != '\n') {

      if (c == '>') {
	while ((c = getc(file)) != EOF && c != '\n') {
	  printf("%c",c);
	}
	printf("\n");
      }

      if (c == EOF) {
	printf("End of file %s\n",filename);
	return hash;
      }

      e = fasta_encoding[(int)c];
      
      if (e >= 0) {
	ei = (ei << 2) | e;
	
	ecode = map_comb(comb,ei,0);
	
	tmphash = (Hash *)malloc(sizeof(Hash));
	tmphash->pos = i;
	
	if (hash[ecode] != NULL) {
	  tmphash->count = hash[ecode]->count + 1;
	} else {
	  tmphash->count = 1;
	}
	
	tmphash->next = hash[ecode];
	hash[ecode] = tmphash;
      }
      i++;
    }
    
  }

  printf("Built hash\n"); 
  getrusage(0,res);
  
  printf("Time %d %d\n",res->ru_utime.tv_sec,res->ru_utime.tv_usec);
  fclose(file);
  return hash;
}  


int main(int argc, char *argv[]) {
  
  int i;
  int j;

  int e;
  static unsigned long long ei = 0;
  static unsigned int ecode;

  Sequence *query;

  Hash **hash;
  Comb *comb;

  Hash *tmphash;

  int n              = atoi(argv[3]);
  int m              = atoi(argv[4]);
  int ncombs         = atoi(argv[5]);
  //int do_transitions = atoi(argv[6]);

  int count    = 0;
  int combhits = 0;

  int *bits;
  int combcount;

  int found;
  int maxcount;
  int maxindex;
  int liscount;

  int l_qstart;
  int l_qend;
  int l_tstart;
  int l_tend;

  int tothits;

  Hash **hithash;
  Lis  **lis;

  Lis *tmplis;
  Lis *newlis;

  int do_blastz = 0;

  //unsigned int t;

  int falsecount = 0;
  int truecount  = 0;

  char *blastz_bitstring = "1110100110010101111";

  struct rusage *res;

  srand((unsigned)time(NULL));

  res = (struct rusage *)malloc(sizeof(struct rusage));

  getrusage(0,res);

  printf("Time %d %d\n",res->ru_utime.tv_sec,res->ru_utime.tv_usec);

  query    = (Sequence *)read_fasta(argv[2]);
  
  hash    = (Hash **)malloc(nmer2pow(n) * sizeof (Hash *));


  printf("Size of int %ld\n"      ,sizeof(unsigned int));
  printf("Size of long %ld\n"     ,sizeof(unsigned long));
  printf("Size of long long %ld\n",sizeof(unsigned long long));

  getrusage(0,res);

  printf("Time %d %d\n",res->ru_utime.tv_sec,res->ru_utime.tv_usec);
  printf("Hash size %ld\n",sizeof(Hash));

  hithash = (Hash **)malloc(query->length*sizeof(Hash *));

  printf("Hash size %ld\n",sizeof(Hash));

  getrusage(0,res);

  printf("Time %d %d\n",res->ru_utime.tv_sec,res->ru_utime.tv_usec);

  lis     = (Lis  **)malloc(query->length*20*sizeof(Lis *));

  getrusage(0,res);

  printf("Time %d %d\n",res->ru_utime.tv_sec,res->ru_utime.tv_usec);
  printf("Hash memory %ld %d\n",sizeof(*hash),nmer2pow(n));
  printf("Hit  memory %ld %d\n",sizeof(*hithash),query->length);
  printf("Lis  memory %ld %d\n",sizeof(*lis),query->length*20);

  memset(hithash, 0, query->length*sizeof(Hash *));

  getrusage(0,res);
  
  printf("Time %d %d\n",res->ru_utime.tv_sec,res->ru_utime.tv_usec);

  memset(lis    , 0, query->length*20*sizeof(Lis  *));
  
  getrusage(0,res);
  
  printf("Time %d %d\n",res->ru_utime.tv_sec,res->ru_utime.tv_usec);

  printf("Ncombs %d\n",ncombs);

  if (ncombs == 0) {
    ncombs = 1;
    do_blastz = 1;

    bits = (int *)malloc(sizeof(int)*19);
    i=0;
    while (i < 19) {
      bits[i] = (int)(blastz_bitstring[i]-48);
      i++;
    }
  }


  while (count < ncombs) {
    printf("\nComb number %d\n\n",count);
    
    if (count > 0) {
      /* blocks should be freed in the comb as well*/
      free(comb);
      free(bits);
    }

    if (do_blastz == 0) {
      printf("Generating bitstring %d %d\n",n,m);
      bits = generate_bitstring(n,m);
    } else {
      printf("Using blastz comb\n");
    }

    comb = bits2comb(bits,n,m);
    
    combhits = 0;
    
    memset(hash, 0, nmer2pow(n)*sizeof(Hash *));
    
    hash    = read_fasta_and_hash(argv[1],n,m,comb);

    ei = 0;
    combcount = 0;

    getrusage(0,res);

    printf("Time before lookup %d %d\n",res->ru_utime.tv_sec,res->ru_utime.tv_usec);
    tothits = 0;

    for (i=0; i < (m-1) ;i++) {
      e = fasta_encoding[(int)query->sequence[i]];
      ei = (ei << 2) | e;
    }

    for (i = (m-1); i < query->length; i++) {

      e = fasta_encoding[(int)query->sequence[i]];
      
      if (e >= 0) {
	ei = (ei << 2) | e;

	ecode = map_comb(comb,ei,0);
	
	//printf("position %d ecode %d string ",i,ecode);
	//print_l(stdout,ecode,n);
	//printf("\n");
	
	if (hash[ecode] != NULL) {
	  combhits++;
	  
	  tmphash = hash[ecode];

	  hithash[i] = tmphash;
	  
	  //printf("Hithash %d %d\n",i,hithash[i]->count);
	}
      }
    }
  

  getrusage(0,res);

  printf("Time after lookup %d %d\n",res->ru_utime.tv_sec,res->ru_utime.tv_usec);

    count++;
  }

  for (i = (m-1); i < query->length; i++) {
    tmphash = hithash[i];
    if (tmphash != NULL && tmphash->count < 5000) {
      while (tmphash != NULL) {
	if (tmphash->pos >= query->start &&
	    tmphash->pos <= query->end) {
	  printf("%s\ttrue\tsimilarity\t%d\t%d\t100\t+\t.\tsequence2\t%d\t%d\t%d\n",
		 query->id,i-(m/2),i-(m/2),tmphash->pos-m/2,tmphash->pos-m/2,tmphash->count);
	  truecount++;
	} else {
	  falsecount++;
	  //printf("%s\tfalse\tsimilarity\t%d\t%d\t100\t+\t.\tsequence2\t%d\t%d\t%d\n",
	//	 query->id,i-(m/2),i-(m/2),tmphash->pos-m/2,tmphash->pos-m/2,tmphash->count);
	}
	tmphash = tmphash->next;
      }
    }
  }
  
  printf ("Number of comb hits %d\n",combhits);

  /* Now do the lis
      - loop over all query positions
      - loop through all existing lis
      - are any of the query positions 
          - within a query distance d
	  - within a target distance t
      - if not start a new lis
  */

  
  // Loop over all query positions 
  for (i = (m-1); i < query->length; i++) {
    
    //printf("Lis %d\n\n",i);
    
    found = 0;

    tmphash = hithash[i];

    if (tmphash != NULL && tmphash->count < 5000) {    

    // Loop over all lis previously built 
      for (j = (m-1); j < query->length*20; j++) {
	tmphash = hithash[i];	

	tmplis  = lis[j];
      
	// Look at all hit positions for this query position 
	while (tmphash != NULL && tmplis != NULL) {
	  if (((i - tmplis->query_end) < 2000) &&
	      ((tmphash->pos - tmplis->target_end) < 2000) &&
	      ((tmphash->pos - tmplis->target_end) > 0)   &&
	      (abs(i - tmplis->query_end - tmphash->pos + tmplis->target_end) < 50)) {
	    
	    
	    // We've found an extension 
	    
	    found = 1;
	    printf("Extending lis %d %d %d %d %d\n",i,j,tmphash->pos,tmplis->target_end,(tmphash->pos-tmplis->target_end));
	    newlis = (Lis *)malloc(sizeof(Lis));
	    newlis->query_start  = i;
	    newlis->query_end    = i;
	    newlis->target_start = tmphash->pos;
	    newlis->target_end   = tmphash->pos;

	    newlis->next = tmplis;
	    
	    lis[j] = newlis;
	    lis[j]->count = tmplis->count++;
	  }
	
	
	  tmphash = tmphash->next;
	}
      }

      if (found == 0) {
	
	tmphash = hithash[i];
	
	while (tmphash != NULL) {
	  newlis = (Lis *)malloc(sizeof(Lis));
	  
	  newlis->query_start = i;
	  newlis->query_end   = i;
	  newlis->target_start = tmphash->pos;
	  newlis->target_end   = tmphash->pos;
	  newlis->next = NULL;
	  newlis->count = 1;	

	  //printf("Making new lis %d %d\n",i,tmphash->pos);
	  lis[i] = newlis;

	  tmphash = tmphash->next;
	}
      }
    }
  }


  maxcount = 0;

  for (j = (m-1); j < query->length*20; j++) {
    tmplis = lis[j];
    liscount = 0;
    
    //if (tmplis != NULL && tmplis->count > 4) {
      
      while (tmplis != NULL) {
	
	printf("%s_%s\ttmplis\tsimilarity\t%d\t%d\t100\t+\t.\tsequence2\t%d\t%d\n",
	       argv[2],query->id,tmplis->query_start,tmplis->query_end,tmplis->target_start,tmplis->target_end);
	
	liscount++;
	tmplis = tmplis->next;
      }
      
      if (liscount > 0) {
	printf("\n");
      }

      if (liscount >= maxcount) {
	maxcount = liscount;
	maxindex = j;
      }
      //}
  }

  printf("Max count %s %d %d ",query->id,maxcount,maxindex);
  
  tmplis = lis[maxindex];

  l_qstart = tmplis->query_start;
  l_qend   = tmplis->query_start;
  l_tstart = tmplis->target_start;
  l_tend   = tmplis->target_end;

  while (tmplis != NULL) {
    if (tmplis->query_start < l_qstart) {
      l_qstart = tmplis->query_start;
    }
    if (tmplis->query_end > l_qend) {
      l_qend = tmplis->query_end;
    }
    if (tmplis->target_start < l_tstart) {
      l_tstart = tmplis->target_start;
    }
    if (tmplis->query_start > l_tend) {
      l_tend = tmplis->query_start;
    }
    tmplis = tmplis->next;
  }

  printf("%d %d %d %d\n",l_qstart,l_qend,l_tstart,l_tend);

  tmplis = lis[maxindex];

    while (tmplis != NULL) {
    printf("%s_%s\tlis\tsimilarity\t%d\t%d\t100\t+\t.\tsequence2\t%d\t%d\n",
  	   argv[2],query->id,tmplis->query_start,tmplis->query_end,tmplis->target_start,tmplis->target_end);
    tmplis = tmplis->next;
  }
  
  getrusage(0,res);

  printf("Queryid/start/end/true/false/total/n/m/ncombs/query_len %s %d %d %d %d %d %5.5f %5.5f %d %d %d %d\n",query->id,query->start,query->end,truecount,falsecount,truecount+falsecount,truecount*100.0/(truecount+falsecount),falsecount*query->length/1000000.0,n,m,ncombs,query->length);
  printf("Time after lis %d %d\n",res->ru_utime.tv_sec,res->ru_utime.tv_usec);

  return 0;
}
  
    /*    ei = 0;
	  
    for (i=0; i < (m-1) ;i++) {
    e = fasta_encoding[(int)seq[i]];
    ei = (ei << 2) | e;
    }
    
    getrusage(0,res);

    printf("Time before hash %d %d\n",res->ru_utime.tv_sec,res->ru_utime.tv_usec);

    for (i = (m-1); i < database->length; i++) {
      
    e = fasta_encoding[(int)seq[i]];
      
      if (e >= 0) {
	ei = (ei << 2) | e;
	
	ecode = map_comb(comb,ei,0);

	tmphash = (Hash *)malloc(sizeof(Hash));
	tmphash->pos = i;
	  
	if (hash[ecode] != NULL) {
	  tmphash->count = hash[ecode]->count + 1;
	} else {
	  tmphash->count = 1;
	}

	tmphash->next = hash[ecode];
	hash[ecode] = tmphash;
	
      }
    }
    printf("Built hash\n"); 
    getrusage(0,res);

    printf("Time %d %d\n",res->ru_utime.tv_sec,res->ru_utime.tv_usec);
    
    */ 

	// transitions
    /*	if (do_transitions) {
	  for (j=0; j < n; ++j) {
	    t = ecode ^ (2 << (2 * j));  A-G; C-T 
	    
	    if (hash[t] != NULL) {
	      combhits++;
	      
	      tmphash = hash[t];
	      
	      tmphash->next = hithash[i];

	      if (hithash[i] != NULL) {
	      	tmphash->count += hithash[i]->count;
	      } 
	      
	      hithash[i] = tmphash;
	    }
	  }
	  for (j=0; j < n; ++j) {
	    t = ecode ^ (1 << (2 * j));  A-C; T-G 
	    
	    if (hash[t] != NULL) {
	      combhits++;
	      
	      tmphash = hash[t];
	      
	      tmphash->next = hithash[i];
	      

	      if (hithash[i] != NULL) {
	      	tmphash->count += hithash[i]->count;
	      } 
	      
	      hithash[i] = tmphash;
	    }
	  }
	  for (j=0; j < n; ++j) {
	    t = ecode ^ (3 << (2 * j));  A-T; C-G 
	    
	    if (hash[t] != NULL) {
	      combhits++;
	      
	      tmphash = hash[t];
	      
	      tmphash->next = hithash[i];

	      if (hithash[i] != NULL) {
	      	tmphash->count += hithash[i]->count;
	      } 
	      
	      hithash[i] = tmphash;
	    }
	  }

	  }
	}
	}*/




