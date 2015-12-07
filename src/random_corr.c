/*
 *  random_corr.c
 *  pogview
 *
 *  Created by mclamp on 10/3/08.
 *  Copyright 2008 The Broad Institute. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "datamodel.h"

extern Genome *fetch_human_genome();
extern GenomePosition *get_random_genome_position(Genome *genome,int offset);
extern Sequence *read_fasta(FILE *file);
extern Pimer    *read_pimer(FILE *file);

int main(int argc, char *argv[]) {
  Genome *human = fetch_human_genome();
	
  char *aligndir = "/Volumes/rg/genome/2x_aln/filtered_alignments/";
  char *pidir    = "/Volumes/rg 1/21mammals/pi/12mer/";
  //char *aligndir = "/ahg/scr3/genome/2x_aln/filtered_alignments/";
  //char *pidir    = "/seq/mgscratch/21mammals/pi/12mer/";
  
  int j = 0;
  
  while (j < 1000) {
    
    GenomePosition *gpos = get_random_genome_position(human,1000);
		
    printf("Position is %s %llu %llu %d\n",gpos->chr->name,gpos->pos,gpos->chr->length-gpos->pos,(int)(100*gpos->pos/gpos->chr->length));
    
    int start = 1000000* (int)(gpos->pos/1000000);
    int end   = start+999999;
    
    if (end > gpos->chr->length) {
      end = gpos->chr->length;
    }
    
    char startbuf[20];
    char endbuf[20];
    
    sprintf(startbuf, "%d", start);
    sprintf(endbuf,   "%d", end);
    
    char *filename   = (char *)malloc((strlen(aligndir)  + strlen(gpos->chr->name) + strlen("/") + strlen(gpos->chr->name) + strlen("_") 
                                       + strlen(startbuf) + strlen("-")  + strlen(endbuf) + strlen(".fa") + 1)*sizeof(char));
    char *pifilename = (char *)malloc((strlen(pidir)     + strlen(gpos->chr->name) + strlen("/") + strlen(gpos->chr->name) + strlen("_") 
                                       + strlen(startbuf) + strlen("-")  + strlen(endbuf) + strlen(".pi.12mer") + 1)*sizeof(char));
    
    strcpy(filename,aligndir);
    strcat(filename,gpos->chr->name);
    strcat(filename,"/");
    strcat(filename,gpos->chr->name);
    strcat(filename,"_");
    strcat(filename,startbuf);
    strcat(filename,"-");
    strcat(filename,endbuf);
    strcat(filename,".fa");
    
    strcpy(pifilename,pidir);
    strcat(pifilename,gpos->chr->name);
    strcat(pifilename,"/");
    strcat(pifilename,gpos->chr->name);
    strcat(pifilename,"_");
    strcat(pifilename,startbuf);
    strcat(pifilename,"-");
    strcat(pifilename,endbuf);
    strcat(pifilename,".pi.12mer");
    
    printf("Fasta file is %s\n",filename);
    
    FILE *file = fopen(filename,"r");
    
    if (! file) {
      printf("Filename %s can't be opened\n",filename);
      exit(0);
    }
    
    Sequence *seqptr = read_fasta(file);
    Sequence *seq    = seqptr;
    
    while (seq != NULL) {
      char *subseq = substring(seq->sequence,0,100);
      
      printf("Sequence %20s %s\n",seq->id,subseq);
      
      pog_free(subseq,"Subseq",DEBUG);
      
      //subseq = substring(seq->sequence,gpos->pos - start,gpos->pos+1000 - start);
      //pog_free(seq->sequence,"Sequence",DEBUG);
      //seq->sequence = subseq;
      seq = seq->next;
      
    }
    
    printf("Got here\n");	
    
    FILE *pifile = fopen(pifilename,"r");
    
    if (! pifile ) {
      printf("No file %s, Exit\n",pifilename);
    } else {
      Pimer **pis  = (Pimer **)malloc(sizeof(Pimer *)*(end-start+1));
      Pimer *pi;
      
      int prev     = 0;
      
      while ((pi = read_pimer(pifile)) != (Pimer *)NULL) {
        if (pi->start >= gpos->pos-start) {
          if (pi->start > prev) {
            prev++;
            while (prev < pi->start) {
              pis[prev] = (Pimer *)NULL;
              prev++;
            }
          }
          if (pi->start >= (gpos->pos - start) &&
              pi->start <= (gpos->pos - start + 1000)) {
            //printf("Pos %d %f\n",pi->start,pi->logodds);
            pis[pi->start] = pi;
            prev = pi->start;
            
          } else {
            free(pi->chr);
            free(pi);
          }
        
        } else {
          free(pi->chr);
          free(pi);
        }
      }
      
      free(gpos);
      
      pog_free(filename,"Filename",NO_DEBUG);
      pog_free(pifilename,"Pifilename",NO_DEBUG);
      
      seq = seqptr;
      
      while (seq != NULL) {
        
        Sequence *tmpseq = seq->next;
        
        
        pog_free(seq->sequence,"Sequence2",NO_DEBUG);
        pog_free(seq->id,      "Id",       NO_DEBUG);
        pog_free(seq,          "Seq",      NO_DEBUG);
        
        seq = tmpseq;
        
      }
      
      int i = 0;
      
      while (i < end-start+1) {
        
        if (pis[i] != (Pimer *)NULL) {      
          free(pis[i]->chr);
          free(pis[i]);
        }
        i++;
      }
    }
    
    j++;
  }
  return 1;
}

