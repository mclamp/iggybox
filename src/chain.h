#include <stdio.h>

#include "EnsC.h"

#include "DBAdaptor.h"
#include "AnalysisAdaptor.h"
#include "BaseAdaptor.h"
#include "ChromosomeAdaptor.h"
#include "DNAAlignFeatureAdaptor.h"
#include "ExonAdaptor.h"
#include "GeneAdaptor.h"
#include "ProteinAlignFeatureAdaptor.h"
#include "PredictionTranscriptAdaptor.h"
#include "RawContigAdaptor.h"
#include "FeatureSet.h"
#include "StringHash.h"
#include "RepeatFeatureAdaptor.h"
#include "SequenceAdaptor.h"
#include "SimpleFeatureAdaptor.h"
#include "SliceAdaptor.h"
#include "TranscriptAdaptor.h"
#include "StrUtil.h"
#include "DBEntry.h"
#include "DBEntryAdaptor.h"
#include "DNAAlignFeature.h"
#include "DNAPepAlignFeature.h"
#include "PredictionTranscript.h"
#include "RawContig.h"
#include "RepeatFeature.h"
#include "SimpleFeature.h"

typedef struct _Seq {
  char    *id;
  char    *sequence;
  int     length;
  int     start;
  int     end;
  struct _Seq *next;
} Seq;


int Sequence_max_sequence_length(Seq *list);


typedef struct _chain {
  int     id;
  int     id2;
  char    *chr1;
  int     start;
  int     end;
  int     chr1len;
  char    *chr2;
  int     chr2len;
  int     strand;
  int     hstart;
  int     hend;
  int     chainid;
} Chain;

typedef struct _link {
  int     id;
  char     *chr;
  int     start;
  int     end;
  int     hstart;
  int     hend;
  int     chainid;
  struct _link *next;
} Link;


char **get_tokens(FILE *file,char token, int *ntok);
Chain *read_chain (FILE *file);
Link *read_links(FILE *file,int chainid,Link **newlink);
void print_feature(DNAAlignFeature *daf, int analysis);
int intlen(int num);
void pretty_print(char *seq1, char *seq2,char *id1,char *id2, int linelen);
void get_alignment(DNAAlignFeature *daf,DBAdaptor *dba1,DBAdaptor *dba2,  char *chr,int  chrstart,int  chrend,char *chr_filter);
void parse_chr(char *chr, int *start, int *end);
