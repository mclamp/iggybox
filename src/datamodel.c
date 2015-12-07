#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
//#include <gtk/gtk.h>

#include "datamodel.h"


Sequence    *make_Sequence() {
  Sequence *seq;

  seq = (Sequence *)malloc(sizeof(Sequence));

  return seq;
}

void        *Sequence_set_name(Sequence *seq,char *name) {
  int length = strlen(name);

  seq->name = (char *)malloc(length*sizeof(char)+1);

  if (seq->name == NULL) {
    fprintf(stderr,"Couldn't allocate memory for sequence name\n");
    exit(0);
  }
  strcpy(seq->name,name);

  return NULL;

}

void        *Sequence_set_id(Sequence *seq, char   *id) {
  int length = strlen(id);

  seq->id = (char *)malloc(length*sizeof(char)+1);

  if (seq->id == NULL) {
    fprintf(stderr,"Couldn't allocate memory for sequence id\n");
  }

  strcpy(seq->id,id);
  
  return NULL;
}

void        *Sequence_set_seqstring(Sequence *seq,char *str) {
  int length = strlen(str);

  seq->sequence = (char *)malloc(length*sizeof(char)+1);
  seq->length   = length;

  if (seq->sequence == NULL) {
    fprintf(stderr,"Couldn't allocate memory for sequence string\n");
  }

  strcpy(seq->sequence,str);

  return NULL;

}

void Sequence_add_Sequence(Sequence *seq,Sequence *newseq) {

  while (seq->next != NULL) {
    seq = seq->next;
  }

  seq->next = newseq;

  newseq->next = NULL;
}

char Sequence_get_char_at(Sequence *seq,int position) {

  if (seq == NULL) {
    return ' ';
  } else {
    return seq->sequence[position];
  }
}

Sequence *Sequence_get_Sequence_at(Sequence *list,int pos,int max) {
  int count = 0;

  Sequence *tmp;

  tmp = list;

   while (count < pos && tmp != NULL) {
    tmp = tmp->next;
    count++;
  }
  return tmp;
}

int Sequence_count_Sequences(Sequence *list) {
  int count = 0;

  Sequence *tmp;

  tmp = list;

  while (tmp != NULL) {
    tmp = tmp->next;
    count++;
  }

  return count;
}

int Sequence_max_sequence_length(Sequence *list) {

  int max = 0;

  Sequence *tmp;

  tmp = list;

  while (tmp != NULL) {
    if (tmp->length > max) {
      max = tmp->length;
    }
    tmp = tmp->next;
  }
  return max;
}

Block **find_Blocks(Sequence *head) {
  Block **blockhead = NULL;
  Block *tmp       = NULL;
  
  Sequence *seq;

  int count = 0;
  int prev  = 0;
  int blen  = 0;

  int length = Sequence_max_sequence_length(head);
  
  char c;
  int  found = 0;
  int  seqnum;
  int  tmpcount;

  blockhead = (Block **)malloc(length*sizeof(Block *));

  while (count < length) {
    blockhead[count] = NULL;
    seq = head;
    
    seqnum = 0;
    found  = 0;

    while (seq != NULL) {

      if (seq->length > count && strcmp(seq->id,"consensus") != 0) {

	if (seqnum == 0) {

	  //if (seq->sequence[count] != '-') {
	    c = seq->sequence[count];
	    seqnum++;
	    found = 1;
	    //}

	} else {
	  //	  if (seq->sequence[count] != '-') {
	    if (seq->sequence[count] != c) {
	      found = 0;
	    } else {
	      seqnum++;
	    }
	    //}
	}
      }
      seq = seq->next;
    }

    //if (found == 1 && seqnum > 5) {
    //  printf("Position %d found %d seqnum %d char %c\n",count,found,seqnum,c);
    //}

    //Is it a start of a new block
    if (found == 1 && prev == 0) {
      
      prev   = 1;
      blen   = 1;
      
    } else if (found ==1 && prev == 1) {
      // Increase length
      blen++;
    } else if (found == 0 && prev == 1) {
      // Got to end of block


      // Make block and add to list
      
      tmp = (Block *)malloc(sizeof(Block));
     
      if (tmp == NULL) {
	printf("Can't allocate memory for blocks\n");
	return NULL;
      }

      tmp->position = count-blen;
      tmp->length   = blen;
      tmp->count    = seqnum;
      tmp->next = NULL;

      tmpcount = count-blen;

      while (tmpcount < count) {
	blockhead[tmpcount] = tmp;
	tmpcount++;
      }
      
      prev = 0;
      blen = 0;

    }
    count++;
  }
  return blockhead;
}




      

    

  




  
