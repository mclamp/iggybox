#include <string.h>
#include "chain.h"

typedef struct _PogSequence {
  char    *name;
  char    *id;
  char    *sequence;
  int     length;
  int     start;
  int     end;
  struct _PogSequence *next;
} PogSequence;

PogSequence *get_alignment2(DNAAlignFeature *daf,DBAdaptor *dba1,DBAdaptor *dba2,  char *chr,int  chrstart,int  chrend,char *chr_filter);
char *get_pad_str(int start,int end);

int main(int argc, char *argv[]) {
  DBAdaptor *dba1;
  DBAdaptor *dba2;
  DBAdaptor *dba3;
  DBAdaptor *dba4;

  FILE *outfile;
  char *chr     = argv[1];

  int  chrstart = atoi(argv[2]);
  int  chrend   = atoi(argv[3]);
  char *chr_filter = NULL;

  char *org = argv[5];

  PogSequence *seq;
  PogSequence *tmpseq;

  char *pad_str;
  char *bigstring2;

  int prev;

  if (argc == 7) {
     chr_filter = argv[6];
     printf("Setting chromosome filter to %s\n",chr_filter);
  } 

  outfile = fopen("pog.fa","w");

  
  initEnsC();
  
  printf("Opening connection...");

  
  dba1 = DBAdaptor_new("localhost","root",NULL,"homo_sapiens_core_18_34",3306,NULL);

  /*  if (strcmp("mouse",org) == 0) {
    printf("Connecting to mouse\n");
    dba2 = DBAdaptor_new("localhost","root",NULL,"mus_musculus_core_18_30",3306,NULL);
  } else if (strcmp("dog",org) == 0) {
    printf("Connecting to dog\n");
    dba2 = DBAdaptor_new("lead","ensrw","ensembl","dog_feb_1",3306,NULL);
    }*/
    dba2 = DBAdaptor_new("localhost","root",NULL,"mus_musculus_core_18_30",3306,NULL);
    dba3 = DBAdaptor_new("lead","ensrw","ensembl","dog_feb_1",3306,NULL);
  //dba2 = DBAdaptor_new("localhost","root",NULL,"homo_sapiens_core_18_34",3306,NULL);

  printf("Done\n");

  {
    SliceAdaptor *sa1 = DBAdaptor_getSliceAdaptor(dba1);

    Slice *slice1 = SliceAdaptor_fetchByChrStartEnd(sa1,chr,chrstart,chrend);
	
    Vector *features;

    int analysis;
    int i;
    int j;
    int k;
    char **keys;
    DNAAlignFeature *daf;

    int minquery = NULL;
    int maxquery = NULL;

    Vector *tmpv;
    StringHash *chrHash;

    chrHash  = StringHash_new(STRINGHASH_SMALL);

    features =  Slice_getAllDNAAlignFeatures(slice1,"",NULL);

    printf("Num features %d\n",Vector_getNumElement(features));
    
    Vector_sort(features, SeqFeature_startCompFunc);

    for (i=0;i<Vector_getNumElement(features);i++) {
      DNAAlignFeature *daf = Vector_getElementAt(features,i);

      analysis = Analysis_getDbID(DNAAlignFeature_getAnalysis(daf));
       
      if (analysis == 113 || analysis == 110 || analysis == 112) {

	if (minquery == NULL || FeaturePair_getStart(daf) < minquery) {
	  minquery = FeaturePair_getStart(daf);
	}
	if (maxquery == NULL || FeaturePair_getEnd(daf) > maxquery) {
	  maxquery = FeaturePair_getEnd(daf);
	}
	if (StringHash_contains(chrHash,FeaturePair_getHitSeqName(daf)) == NULL) {
	  tmpv = Vector_new();
	  Vector_addElement(tmpv,daf);
	  StringHash_add(chrHash,FeaturePair_getHitSeqName(daf),tmpv);
	} else {
	  tmpv = StringHash_getValue(chrHash,FeaturePair_getHitSeqName(daf));
	  Vector_addElement(tmpv,daf);
	}
      }
    }

    fprintf(outfile,">%s.%d-%d\n%s\n",chr,chrstart,chrend,Slice_getSubSeq(slice1,minquery,maxquery,1));
    {

      Vector *tmpv2;

      int qgap = 100000000;
      int hgap = 100000000;
      
      int l = 0;
      int found = 0;

      Vector *fvect;

      DNAAlignFeature *daf2;

      j = 0;
      
      keys = StringHash_getKeys(chrHash);
      
      while (j < StringHash_getNumValues(chrHash)) {
	fvect = Vector_new();
	tmpv = (Vector *)StringHash_getValue(chrHash,keys[j]);
	
	k = 0;
	
	Vector_sort(tmpv, SeqFeature_startCompFunc);
	
	while (k < Vector_getNumElement(tmpv)) {
	  daf = (DNAAlignFeature *)Vector_getElementAt(tmpv,k);
	
	  l = 0;
	  found = 0;

	    printf("FEAT %s\tblastz\tsimilarity\t%d\t%d\t%4.2f\t%d\t.\t%s\t%d\t%d\t%d\n",
		   FeaturePair_getSeqName(daf),
		   FeaturePair_getStart(daf),
		   FeaturePair_getEnd(daf),
		   FeaturePair_getScore(daf),
		   FeaturePair_getStrand(daf),
		   FeaturePair_getHitSeqName(daf),
		   FeaturePair_getHitStart(daf),
		   FeaturePair_getHitEnd(daf),
		   FeaturePair_getHitStrand(daf));

	  while (l < Vector_getNumElement(fvect) && found == 0) {
	    tmpv2 = (Vector *)Vector_getElementAt(fvect,l);

	    daf2 = (DNAAlignFeature *)Vector_getElementAt(tmpv2,Vector_getNumElement(tmpv2)-1);

	    if (FeaturePair_getHitStrand(daf) == 1) {
	      if ((FeaturePair_getHitStrand(daf) == FeaturePair_getHitStrand(daf2))   &&
		  (FeaturePair_getStart(daf) - FeaturePair_getEnd(daf2)) < qgap       && 
		  (FeaturePair_getStart(daf) - FeaturePair_getEnd(daf2)) > 0      
		) {
		printf("Adding element\n");
		Vector_addElement(tmpv2,daf);
		found = 1;
	      }
	    } else {
	      if ((FeaturePair_getHitStrand(daf) == FeaturePair_getHitStrand(daf2)) &&
		  (FeaturePair_getStart(daf) - FeaturePair_getEnd(daf2)) < qgap   &&
		  (FeaturePair_getStart(daf) - FeaturePair_getEnd(daf2)) > 0      
		  ) {
		printf("Adding element\n");
		Vector_addElement(tmpv2,daf);
		found = 1;
	      }
	    }
	    l++;
	  }

	  if (found == 0) {
	    printf("\nNew block\n");
	    tmpv2 = Vector_new();

	    Vector_addElement(tmpv2,daf);
	    Vector_addElement(fvect,tmpv2);
	  }

	  k++;
	}

	k = 0;
	
	while (k < Vector_getNumElement(fvect)) {
	  l = 0;
	  tmpv2 = (Vector *)Vector_getElementAt(fvect,k);
	  
	  printf("\nNew feature block\n");

	  // Pad the beginning
	  daf2 = (DNAAlignFeature *)Vector_getElementAt(tmpv2,0);

	  pad_str = get_pad_str(minquery,FeaturePair_getStart(daf2));

	  bigstring2 = (char *)calloc(1,(strlen((char*)pad_str)+1)*sizeof(char));

	  strcpy(bigstring2,pad_str);

	  while (l < Vector_getNumElement(tmpv2)) {
	  
	    daf2 = (DNAAlignFeature *)Vector_getElementAt(tmpv2,l);

	    printf("FEAT %s\tblastz\tsimilarity\t%d\t%d\t%4.2f\t%d\t.\t%s\t%d\t%d\t%d\n",
		   FeaturePair_getSeqName(daf2),
		   FeaturePair_getStart(daf2),
		   FeaturePair_getEnd(daf2),
		   FeaturePair_getScore(daf2),
		   FeaturePair_getStrand(daf2),
		   FeaturePair_getHitSeqName(daf2),
		   FeaturePair_getHitStart(daf2),
		   FeaturePair_getHitEnd(daf2),
		   FeaturePair_getHitStrand(daf2));

	    if (l > 0) {
	    //  if (FeaturePair_getStart(daf2) - prev > 1000) {
	    //		pad_str = get_pad_str(1,1000);
	    // } else {
		pad_str = get_pad_str(prev,FeaturePair_getStart(daf2));
		// }

	      bigstring2 = realloc(bigstring2,(strlen(bigstring2) + strlen(pad_str) + 1)*sizeof(char));
	      strcat(bigstring2,pad_str);
	    }

	    prev = FeaturePair_getEnd(daf2) + 1;

	    if (strncmp(FeaturePair_getHitSeqName(daf2),"UNK",3) == 0) {
	      seq = get_alignment2(daf2,dba1,dba3,chr,chrstart,chrend,NULL);
	    } else {
	      seq = get_alignment2(daf2,dba1,dba2,chr,chrstart,chrend,NULL);
	    }
	    bigstring2 = realloc(bigstring2,(strlen(bigstring2) + seq->length + 1)*sizeof(char));

	    tmpseq = seq;

	    tmpseq = tmpseq->next;	    
	    strcat(bigstring2,tmpseq->sequence);
	    l++;
	  }
	  fprintf(outfile,">%s\n%s\n",seq->next->id,bigstring2);	 
	  k++;
	}
	j++;
      }
    }
  }
  fclose(outfile);
  return 1;
}




PogSequence *get_alignment2(DNAAlignFeature *daf,DBAdaptor *dba1,DBAdaptor *dba2,  char *chr,int  chrstart,int  chrend,char *chr_filter) {

    int j;
    int k;

    Vector *ungapped;

    char   *alignstr1 = NULL;
    char   *alignstr2 = NULL;

    char   *tmpstr1;
    char   *tmpstr2;

    int start;
    int end;

    int hstart;
    int hend;

    int prev  = NULL;
    int hprev = NULL;

    char *pad;
    char *hpad;

    char *gappy;
    char *hgappy;

    char *id1;
    char *id2;

    char *hitid;

    int tmpstart = NULL;
    int tmpend   = NULL;

    SliceAdaptor *sa1 = DBAdaptor_getSliceAdaptor(dba1);
    SliceAdaptor *sa2 = DBAdaptor_getSliceAdaptor(dba2);

    Slice *slice1;
    Slice *slice2;

    PogSequence *seq1;
    PogSequence *seq2;

    prev = NULL;

    hitid = (char *)calloc(1,(strlen(DNAAlignFeature_getHitSeqName(daf))+1)*sizeof(char));


    if (hitid == NULL) {
      fprintf(stderr,"ERROR: Can't allocate memory for hitid\n");
    }
    seq1 = (PogSequence *)calloc(1,sizeof(PogSequence));
    seq2 = (PogSequence *)calloc(1,sizeof(PogSequence));

    if (seq1 == NULL) {
      fprintf(stderr,"ERROR: Can't allocate memory for sequence in get_alignment2\n");
      exit(0);
    }

    if (seq2 == NULL) {
      fprintf(stderr,"ERROR: Can't allocate memory for sequence in get_alignment2\n");
      exit(0);
    }

    strcpy(hitid,DNAAlignFeature_getHitSeqName(daf));

    parse_chr(hitid,&tmpstart,&tmpend);

    start = FeaturePair_getStart(daf);
    end   = FeaturePair_getEnd(daf);
    
    hstart = FeaturePair_getHitStart(daf);
    hend   = FeaturePair_getHitEnd(daf);

    //printf("DNA align feature: %d-%d id %d-%d " IDFMTSTR " %s\n", DNAAlignFeature_getStart(daf),
    //   DNAAlignFeature_getEnd(daf),
    //	   DNAAlignFeature_getHitStart(daf),
    //	   DNAAlignFeature_getHitEnd(daf),
    //	   DNAAlignFeature_getDbID(daf),
    //	   DNAAlignFeature_getHitSeqName(daf));

    if (tmpstart != NULL && tmpend != NULL) {
      hstart += tmpstart -1;
      hend   += tmpstart -1;
    }
    
    slice1 = SliceAdaptor_fetchByChrStartEnd(sa1,chr,
					     chrstart,
					     chrend);


    slice2 = SliceAdaptor_fetchByChrStartEnd(sa2,
					     hitid,
					     hstart,
					     hend);

    FeaturePair_setHitStart(daf,hstart);
    FeaturePair_setHitEnd(daf,hend);

    id1 = (char *)calloc(1,(strlen(chr) + intlen(chrstart+start-1) + intlen(end+chrstart-1) + 3) * sizeof(char));
    id2 = (char *)calloc(1,(strlen(FeaturePair_getHitSeqName(daf)) + 
			    intlen(hstart) + 
			    intlen(hend) + 3) * sizeof(char));

    if (id1 == NULL) {
      fprintf(stderr,"ERROR: Can't allocate memory for id1\n");
    }
    if (id2 == NULL) {
      fprintf(stderr,"ERROR: Can't allocate memory for id2\n");
    }

    sprintf(id1,"%s.%d-%d",chr,start+chrstart-1,end+chrstart-1);
    sprintf(id2,"%s.%d-%d",hitid,
	    hstart,
	    hend);

    ungapped = DNAAlignFeature_getUngappedFeatures((BaseAlignFeature *)daf);      
    //Vector_sort(ungapped, SeqFeature_startCompFunc);

    for (j = 0; j < Vector_getNumElement(ungapped); j++) {
      DNAAlignFeature *ugf = Vector_getElementAt(ungapped,j);


      tmpstr1 = Slice_getSubSeq(slice1,FeaturePair_getStart(ugf),FeaturePair_getEnd(ugf),1);
      if (slice2 == NULL) {
	printf("Slice2 is null\n");
      }

      tmpstr2 = Slice_getSubSeq(slice2,FeaturePair_getHitStart(ugf)-hstart+1,FeaturePair_getHitEnd(ugf)-hstart+1,FeaturePair_getHitStrand(ugf));

      if (prev != NULL) {
	if ((FeaturePair_getStart(ugf) - prev) == 1) {
	  
	  // mouse insert
	  
	  if (FeaturePair_getHitStrand(ugf) == 1) {
	    

	    hgappy = Slice_getSubSeq(slice2,hprev - hstart + 2,FeaturePair_getHitStart(ugf)  - hstart,FeaturePair_getHitStrand(ugf));  


	    k = 0;
	    
	    pad = (char *)calloc(1,(FeaturePair_getHitStart(ugf) - hprev )*sizeof(char));
	    
	    if (pad == NULL) {
	      fprintf(stderr,"Can't allocate memory for human pad string\n");
	      exit(0);
	    }
	    
	    
	    while (k < (FeaturePair_getHitStart(ugf) - hprev - 1)) {
	      pad[k] = '-';
	      k++;
	    }
	    
	    pad[k] = '\0';

	    //alignstr1 = (char *)StrUtil_appendString(alignstr1,pad);
	    alignstr1 = (char *)StrUtil_appendString(alignstr1,tmpstr1);
	    //alignstr2 = (char *)StrUtil_appendString(alignstr2,hgappy);
	    alignstr2 = (char *)StrUtil_appendString(alignstr2,tmpstr2);
	    
	  } else {
	    
	    hgappy = Slice_getSubSeq(slice2,FeaturePair_getHitEnd(ugf) - hstart +2,hprev  - hstart ,FeaturePair_getHitStrand(ugf));
	    
	    pad = (char *)calloc(1,(hprev - FeaturePair_getHitEnd(ugf))*sizeof(char));
	    
	    if (pad == NULL) {
	      fprintf(stderr,"Can't allocate memory for human pad string\n");
	      exit(0);
	    }
	    
	    k = 0;
	    
	    while (k < (hprev - FeaturePair_getHitEnd(ugf) - 1)) {
	      pad[k] = '-';
	      k++;
	    }
	    
	    pad[k] = '\0';

	    //printf("*************** %s\n",hgappy);
	    //printf("*************** %s\n",pad);
	    
	    //alignstr1 = (char *)StrUtil_appendString(alignstr1,pad);
	    alignstr1 = (char *)StrUtil_appendString(alignstr1,tmpstr1);
	    //alignstr2 = (char *)StrUtil_appendString(alignstr2,hgappy);
	    alignstr2 = (char *)StrUtil_appendString(alignstr2,tmpstr2);
	    
	  }
	  
	} else {
	  
	  // human insert
	  
	  hpad = (char *)calloc(1,(FeaturePair_getStart(ugf) - prev)*sizeof(char));
	  
	  if (hpad == NULL) {
	    fprintf(stderr,"Can't allocate memory for mouse pad string\n");
	    exit(0);
	  }
	  k = 0;
	  
	  while (k < (FeaturePair_getStart(ugf) - prev) - 1) {
	    hpad[k] = '-';
	    k++;
	  }
	  
	  hpad[k] = '\0';
	  
	  gappy = Slice_getSubSeq(slice1,prev + 1,FeaturePair_getStart(ugf)-1,1);
	  
	  //printf("Pad %s\n",hpad);
	  //printf("Gappy %d %d %d %s\n",prev,start,FeaturePair_getStart(ugf),gappy);
	  
	  alignstr1 = (char *)StrUtil_appendString(alignstr1,gappy);
	  alignstr1 = (char *)StrUtil_appendString(alignstr1,tmpstr1);
	  alignstr2 = (char *)StrUtil_appendString(alignstr2,hpad);
	  alignstr2 = (char *)StrUtil_appendString(alignstr2,tmpstr2);
	  
	} 
      } else {
	alignstr1 = (char*)calloc(1,(1+strlen(tmpstr1))*sizeof(char));
	alignstr2 = (char*)calloc(1,(1+strlen(tmpstr2))*sizeof(char));
	
	strcpy(alignstr1,tmpstr1);
	strcpy(alignstr2,tmpstr2);
	
      }
      
      
      prev = FeaturePair_getEnd(ugf);
      
      if (FeaturePair_getHitStrand(ugf) == 1) {
	hprev = FeaturePair_getHitEnd(ugf);
      } else {
	hprev = FeaturePair_getHitStart(ugf);
      }
    }

    seq1->id = id1;
    seq1->sequence = alignstr1;
    seq1->next = seq2;
    seq1->length = strlen(alignstr1);

    seq2->id = id2;
    seq2->sequence = alignstr2;
    seq2->next = NULL;
    seq2->length = strlen(alignstr2);

    return seq1;
    
    
    //pretty_print(alignstr1,alignstr2,id1,id2,80);

}

char *get_pad_str(int start,int end) {

  int i;

  char *pad_str;

  pad_str = (char *)calloc(1,sizeof(char)*(end-start+2));

  if (end < start) {
    fprintf(stderr,"Can't create pad string - end < start %d %d\n",end,start);
    exit(0);
  }

  i = start;

  while (i < end) {
    pad_str[i-start] = '-';
    i++;
  }

  pad_str[i-start] = '\0';

  return pad_str;
}
    
