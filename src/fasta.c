#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>

#include "datamodel.h"
#include "pog_utils.h"
#include "seqio.h"

Sequence *read_sequence(FILE *file) {
  char       c;
  
  int idlen  = 0;
  int seqlen = 0;
  
  int idchunk  = 20;
  int seqchunk = 10000;
  
  int curr_idchunk  = idchunk;
  int curr_seqchunk = seqchunk;
  
  Sequence *seq;
    

  seq = (Sequence *)malloc(sizeof(Sequence));
  
  if (seq == NULL) {
    fprintf(stderr,"Can't allocate memory for sequence\n");
    exit(0);
  }
  seq->id       = NULL;
  seq->sequence = NULL;
  seq->next     = NULL;


  
  while ((c = fgetc(file)) != EOF) {
    
    // Beginning of sequence - return old sequence if exists
    if (c == '>') {
      
      if (seq != NULL && seq->id != NULL) {
        ungetc(c,file);
        
        seq->sequence[seqlen] = '\0';
        seq->sequence = (char *)realloc(seq->sequence,(seqlen+1)*sizeof(char));
        seq->length = seqlen;
        
        return seq;
        
      }
      
      // Allocate space for id

      seq->id = (char *)malloc(idchunk * sizeof(char));
      
      if (seq->id == NULL) {
        fprintf(stderr,"Can't allocate memory for sequence id\n");
        exit(0);
      }
      
      idlen = 0;
      
      // Read id
      while ((c = fgetc(file)) != EOF && (c != '\n') && (c != ' ')) {
        
        if (idlen < idchunk) {
          seq->id[idlen] = c;
        } else {
          idchunk += curr_idchunk;
          seq->id = (char *)realloc(seq->id,idchunk*sizeof(char));
          if (seq->id == NULL) {
            fprintf(stderr,"Can't allocate memory for sequence id\n");
            exit(0);
          }
          seq->id[idlen] = c;
        }
        idlen++;
      }
      
      // Finish off id string and reallocate memory to the exact size for id
      seq->id[idlen] = '\0';
      seq->id = (char*)realloc(seq->id,(idlen+1)*sizeof(char));
      
      seq->sequence = (char *)malloc(seqchunk*sizeof(char));
      
      if (seq->sequence == NULL) {
        fprintf(stderr,"Can't allocate memory for sequence string\n");
        exit(0);
      }
      
      seqlen = 0;
      
    } else if (c != '\n') {
      if (seq->sequence != NULL) {
        
        if (seqlen < seqchunk) {
          seq->sequence[seqlen] = c;
        } else {
          seqchunk += curr_seqchunk;
          
          seq->sequence = (char *)realloc(seq->sequence,seqchunk*sizeof(char));
          seq->sequence[seqlen] = c;
        }
        seqlen++;
      }
    }
  }
  
  if (seq->id != NULL) {
    seq->sequence[seqlen] = '\0';
    seq->sequence = (char *)realloc(seq->sequence,(seqlen+1)*sizeof(char));
    seq->length = seqlen;
    return seq;
  } else {
    free(seq);
    return NULL;
  }
  
  return seq;
}

void free_seqptr(Sequence *seqptr) {
  Sequence *seq = seqptr;
  
  while (seq != NULL) {
    
    Sequence *tmpseq = seq->next;
    
    pog_free(seq->sequence,"Sequence2",NO_DEBUG);
    pog_free(seq->id,      "Id",       NO_DEBUG);
    pog_free(seq,          "Seq",      NO_DEBUG);
    
    seq = tmpseq;
    
  }
  
}

char *get_fasta_filename(GenomePosition *pos) {
  char *aligndir = FASTA_STUB;
  //char *aligndir = "/Volumes/rg/genome/2x_aln/filtered_alignments/";
  //char *aligndir    = "/ahg/scr3/genome/2x_aln/filtered_alignments/";
  //char *aligndir    = "/seq/mgscratch/21mammals/AR_aligns/";
  
  int start = 1000000* (int)(pos->pos/1000000);
  int end   = start+999999;
  
  if (end > pos->chr->length) {
    end = pos->chr->length;
  }
  
  char startbuf[20];
  char endbuf[20];
  
  sprintf(startbuf, "%d", start);
  sprintf(endbuf,   "%d", end);  
  
  
  char *filename = (char *)malloc((strlen(aligndir)     + strlen(pos->chr->name) + strlen("/") + strlen(pos->chr->name) + strlen("_") 
                                     + strlen(startbuf) + strlen("-")  + strlen(endbuf) + strlen(".fa") + 1)*sizeof(char));
  
  
  strcpy(filename,aligndir);
  strcat(filename,pos->chr->name);
  strcat(filename,"/");
  strcat(filename,pos->chr->name);
  strcat(filename,"_");
  strcat(filename,startbuf);
  strcat(filename,"-");
  strcat(filename,endbuf);
  strcat(filename,".fa");
  
  return filename;
}


Sequence *read_fasta_test(GenomePosition *pos, int len) {


  // This can't read things > 1000000
  char *filename = get_fasta_filename(pos);
  
 int offset = 1000000 * (int)(pos->pos/1000000);
  
  FILE *file = fopen(filename,"r");
 
  printf("File is %s\n",filename);

  if (file == NULL) {
    printf("Can't read fasta file %s\n",filename);
    return NULL;
  }
 
  printf("Position is %llu - offset = %d - length = %d\n",pos->pos,offset,len);

  Sequence *seq = read_fasta(file);
  
  Sequence *tmp = seq;
  
  while (tmp != NULL) {

    char *sub = substring(tmp->sequence,pos->pos - offset, pos->pos - offset + len - 1);

    //    printf("Getting substring %d\t%d\t%d\t%d\n",pos->pos - offset, pos->pos - offset + len - 1, strlen(sub),len);
    free(tmp->sequence);
    
    tmp->sequence = sub;
    tmp->length   = len;
    
    tmp = tmp->next;
  }
  
  free(filename);
  return seq;
}
  
Sequence *read_fasta(FILE *file) {
  
  Sequence    *head = NULL;
  Sequence    *tmpseq = NULL;
  Sequence    *curr_seq;
  
  Sequence    *p;
  
  int       count = 0;
  
  while ((curr_seq = read_sequence(file)) != NULL) {
    
    if (count == 0) {
      head = curr_seq;
      head->next = NULL;
      tmpseq = curr_seq;
      
    } else {
      tmpseq->next = curr_seq;
      curr_seq->next = NULL;
      
      p = head;
      
      tmpseq = curr_seq;
      
      while (p->next != NULL) {
        p = p->next;
      }
    }
    count++;
  }
  
  fclose(file);
  //fprintf(stderr,"Total seqs %d\n",count);
  
  return head;
}

void save_fasta(char *filename,Sequence * head) {
  
  FILE *file;
  
  Sequence *seq;
  char     *line;
  
  int       line_len = 72;
  int       pos      = 0;
  int       i;
  /* Should check whether the file exists */
  
  file = fopen(filename,"w");
  
  seq = head;
  
  line = (char *)malloc(sizeof(char) * (line_len+1));
  
  if (line == NULL) {
    fprintf(stderr,"Can't allocate memory for line\n");
    exit(0);
  }
  
  
  while (seq != NULL) {
    
    fprintf(file,">%s\n",seq->id);
    
    pos = 0;
    
    while (pos <= strlen(seq->sequence)) {
      
      for (i = 0; i < line_len; i++) {
        line[i] = seq->sequence[pos+i];
      }
      
      line[line_len] = '\0';
      
      fprintf(file,"%s\n",line);
      
      pos += line_len;
    }
    seq = seq->next;
  }
  
  fclose(file);
  
}

