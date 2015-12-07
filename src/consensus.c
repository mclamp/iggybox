#include <strings.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "datamodel.h"
#include "seqio.h"

void  get_partial_Sequence_hash(Sequence * head,int **consensus,int start, int end) {

  char c;

  int i = 0;
  int j = 0;

  int height = Sequence_count_Sequences(head);

  Sequence *tmpseq;

  int **tmp;
  int  *tmpint;

  tmp = consensus;

  for (i = start-1; i < end; i++) {

    for (j = 0; j < height; j++) {
      
      tmpseq = Sequence_get_Sequence_at((Sequence *)head,(int)j,height);
      
      c = Sequence_get_char_at(tmpseq,i);
      
      if (c-45 < 0) {
	c = 45;
      }
      if (c-45 > 45) {
	c = 45;
      }

      tmpint = (int *)(consensus + (i-start+1)*46 + c - 45);      
      
      *tmpint += 1;

    }

    if (i%1000000 == 0) {
      printf("column reached %d\n",i);
    }
  }
}

Sequence *get_partial_consensus_Sequence(Sequence* head, int start, int end) {
  char  *seq;
  Sequence *s;

  int i;
  int j;

  int maxindex = -1;
  int max;
  
  int *tmpint;
  int **consensus;

  consensus = (int **)malloc((end-start+1) * sizeof(int) * 100);

  if (consensus == NULL) {
    printf("Can't allocate memory for consensus\n");
    return NULL;
  }

  get_partial_Sequence_hash(head,consensus,start,end);

  seq = (char *)malloc((end-start+2) *sizeof(char));

  for (i=0 ; i < (end-start); i++) {
    max = -1;
    for (j=0; j <= 45; j++) {
      tmpint = (int*)(consensus + (i-start+1)*46 + j);

      if (*tmpint >= max) {
	max = *tmpint;
	maxindex = j;
      }
    }

    //printf("Max %d %d %c\n",max,maxindex,maxindex+45);

    seq[i] = maxindex+45;

    if (i%1000000 == 0) {
      printf("Calculated %d\n",i);
    }
  }

  seq[end - start+1] = '\0';

  //printf("Sequence %s\n",seq);

  free(consensus);

  s = (Sequence *)malloc(sizeof(Sequence));

  s->id = "consensus";
  s->sequence = seq;
  s->length = strlen(seq);

  return s;
}
void  get_Sequence_hash(Sequence * head,int **consensus) {
  int start    = 1;
  int end      = Sequence_max_sequence_length(head);
  
  get_partial_Sequence_hash(head,consensus,start,end);
}

/* Now I don't really understand this bit but it is
   getting a unique integer from a string */

int get_key(char *key_string) {
  int i;

  int key = -1;

  for (i = 0; i < strlen(key_string); i++) {

     key = (key << 4) + (key ^ (int)key_string[i]);
   }
  return key;
}

  
Sequence  *get_consensus_Sequence(Sequence *head) {
  Sequence *consensus;

  consensus = get_partial_consensus_Sequence(head,1,Sequence_max_sequence_length(head));

  return consensus;
}

Sequence  *pogget_consensus_Sequence(Sequence *head) {

  int i;
  int j;

  int **consensus;
  int width;
  int height;

  char *maxseq;

  int max;
  int maxindex = -1;

  Sequence *seq;

  int *tmpint;

  width  = Sequence_max_sequence_length(head);
  height = Sequence_count_Sequences(head);

  consensus = (int **)malloc(width * sizeof(int) * 100);

  if (consensus == NULL) {
    fprintf(stderr,"Can't allocate memory for consensus\n");
    exit(0);
  }

  get_Sequence_hash(head,consensus);

  maxseq = (char *)malloc((1+width)*sizeof(char));

  for (i=0 ; i <= width; i++) {
    max = 0;
    for (j=0; j < 46; j++) {
      tmpint = (int*)(consensus + i*46 + j);

      if (*tmpint >= max) {
	max = *tmpint;
	maxindex = j;

      }
    }
  
    if (max >= height) { 
      maxseq[i] = maxindex+45;
    } else {
      maxseq[i] = '-';
    } 

    if (i%1000000 == 0) {
      printf("Calculated %d\n",i);
    }
  }

  maxseq[width] = '\0';

  seq = make_Sequence();

  Sequence_set_seqstring(seq,maxseq);
  Sequence_set_id(seq,"consensus");

  seq->next = NULL;

  return seq;
}


  

