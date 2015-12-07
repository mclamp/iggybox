#ifndef POG_UTILS

#include <stdio.h>
#include <dlfcn.h>
#include <stdlib.h>

#include "datamodel.h"
#include <math.h>

//#define MAC_OSX 1

#ifdef MAC_OSX

#define FASTA_STUB    "/Volumes/rg/mammals/data/filtered_alignments/"

#define PI_21_12MER   "/Volumes/rg 1/21mammals/pi/12mer/"
#define PI_21_50mer   "/Volumes/rg/genome/2x_align/estimation/mammals21/integration_runs/pi_ar_model_50mers/"

#define PI_4_12mer    "/Volumes/rg/genome/2x_align/estimation/mammals21/integration_runs/HMRD/12mer/"
#define PI_4_50mer    "/Volumes/rg/genome/2x_align/estimation/mammals21/integration_runs/HMRD/50mer/"

#define PI_4_RAW      "/Volumes/rg 1/21mammals/HMRD/pi_ar_model"
#define PI_21_RAW     "/Volumes/rg/genome/2x_aln/estimation/mammals21/runs/pi_ar_model"

# else

#define FASTA_STUB    "/home/mclamp/2x_aln/"

#define PI_21_12MER   "/seq/mgscratch/21mammals/pi/12mer/"
#define PI_21_50mer   "/ahg/scr3/genome/2x_align/estimation/mammals21/integration_runs/pi_ar_model_50mers/"

#define PI_4_12mer    "/ahg/scr3/genome/2x_align/estimation/mammals21/integration_runs/HMRD/12mer/"
#define PI_4_50mer    "/ahg/scr3/rg/genome/2x_align/estimation/mammals21/integration_runs/HMRD/50mer/"

#define PI_4_RAW      "/seq/mgscratch/21mammals/HMRD/pi_ar_model"
#define PI_21_RAW     "/ahg/scr3/genome/2x_aln/estimation/mammals21/runs/pi_ar_model"

#endif

char *int2char(int i);

void pog_free(void *ptr, char *desc, int debug) {

	if (debug == 1) {
		printf("free %s (%p)\n", desc,ptr);     
	}
    free(ptr);
}

char *substring(char *seq,int start, int end) {
  char *substr = (char *)malloc((end-start+2)*sizeof(char));
  
  int i = start;
  
  while (i <= end) {
    substr[i-start] = seq[i];
    i++;
  }
  substr[i-start] = '\0';
  return substr;
}


char *int2char(int i) {
  int len = 1;

  if (i > 0) {
    len = (int)(log10(i))+1;
  } else if (i < 0) {
    len = -1*(int)(log10(i))+3;
  }

  printf("Len is %d %d\n",len,i);
  char *buf = (char *)malloc((len+1)*sizeof(char));
  
  sprintf(buf, "%d", i);
  
  buf[len] = '\0';


  return buf;
}

static const char dna_complement[] =
"                                                                "
" TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      "
"                                                                "
"                                                                ";
/* ................................................................ */
/* @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~. */
/* ................................................................ */
/* ................................................................ */


#endif
