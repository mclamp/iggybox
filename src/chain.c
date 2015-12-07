#include "chain.h"

int intlen(int num) {

    int i = 1;

    while ((int)num > 10) {
      i++;
      
      num = num/10;
    }
    return i;
  }

char **get_tokens(FILE *file,char token, int *ntok) {
  int intoken = 0;

  char   c;

  char  *prev     = NULL;

  char **tokens   = NULL;

  int    tokcount    = 0;
  int    toklen      = 0;
  int    tokchunk    = 12;
  int    toklenchunk = 24;

  char *   realtok = NULL;

  tokens = (char **)calloc(1,tokchunk * sizeof(char *));

  
  while ((c = fgetc(file)) != EOF) {
    //printf("Char is %c\n",c);
    if (c == token) {
      
      if (intoken == 0) {

	intoken = 1;
	
	if (prev != NULL) {

	  realtok = (char *)calloc(1,(toklen+1) * sizeof(char));
	  strncpy(realtok,prev,toklen);
	  realtok[toklen] = '\0';
	  free(prev);
	  tokens[tokcount] = realtok;

	  tokcount++;

	  prev = (char *)calloc(1,toklenchunk * sizeof(char));

	  toklen = 0;
	}
      }
    } else if ( c == '\n') {
      
      if (prev != NULL) {
	
	realtok = (char *)calloc(1,(toklen+1) * sizeof(char));
	strncpy(realtok,prev,toklen);
	realtok[toklen] = '\0';
	free(prev);
	tokens[tokcount] = realtok;
	tokcount++;
      }	

      *ntok = tokcount;

      return tokens;

    } else {
      
      if (intoken == 0) {
	if (prev == NULL) {
	  prev = (char *)calloc(1,toklenchunk * sizeof(char));
	  toklen = 0;
	}
	
	prev[toklen] = c;
	toklen++;
	//	printf("Setting char %c %d\n",c,toklen);
      } else {
	//prev = (char *)calloc(1,toklenchunk * sizeof(char));

	toklen = 0;
	prev[toklen] = c;
	toklen++;
	intoken = 0;
      }
    }
  }

  return NULL;
}

Chain *read_chain (FILE *file) {

  Chain *chain;

  char **tokens   = NULL;
  int ntok;
  int i;
  char *tmp;

  while ((tokens = get_tokens(file,'\t',&ntok)) != NULL) {
    chain = (Chain *)calloc(1,sizeof(Chain));
    
    if (chain == NULL) {
      fprintf(stderr,"Can't allocate memory for chain\n");
      exit(0);
    }
    
    chain->id      = atoi(tokens[0]);
    chain->id2     = atoi(tokens[1]);
    chain->chr1    = tokens[2];
    chain->chr1len = atoi(tokens[3]);
    chain->start   = atoi(tokens[4]);
    chain->end     = atoi(tokens[5]);
    chain->chr2    = tokens[6];
    chain->chr2len = atoi(tokens[7]);

					  
    if (strncmp("chr",chain->chr1,3) == 0) {
      chain->chr1 = (char *)calloc(1,(strlen(tokens[2])-2)*sizeof(char));
      
      if (chain->chr1 == NULL) {
	fprintf(stderr,"Can't allocate memory for chromosome name in chain\n");
	exit(0);
      }
      tmp  = tokens[2];
      tmp += 3;
      i    = 0;

      while (*tmp != '\0') {
	chain->chr1[i] = *tmp;
	tmp++;
	i++;
      }

      chain->chr1[i] = '\0';
    }

    if (strncmp("chr",chain->chr2,3) == 0) {
      chain->chr2 = (char *)calloc(1,(strlen(tokens[6])-2)*sizeof(char));

      if (chain->chr2 == NULL) {
	fprintf(stderr,"Can't allocate memory for chromosome2 name in chain\n");
	exit(0);
      }
      
      tmp  = tokens[6];
      tmp += 3;
      i    = 0;

      while (*tmp != '\0') {
	chain->chr2[i] = *tmp;
	tmp++;
	i++;
      }

      chain->chr2[i] = '\0';
    }
    
    if (strcmp(tokens[8],"-") == 0) {
      chain->strand =  -1;
    } else {
      chain->strand = 1;
    }
    chain->hstart  = atoi(tokens[9]);
    chain->hend    = atoi(tokens[10]);
    chain->chainid = atoi(tokens[11]);

    free(*tokens);
    free(tokens);

    return chain;
  }

  return NULL;
}
	  
Link *read_links(FILE *file,int chainid,Link **newlink) {

  char **tokens;
  int    count = 0;
  int    ntok;
  
  Link  *tmplink = NULL;
  Link  *link = NULL;
  Link  *prev = NULL;

  while ((tokens = get_tokens(file,'\t',&ntok)) != NULL) {

    if (ntok < 6) {
      fprintf(stderr,"Not enough tokens for a link line\n");
    } else {
      tmplink = (Link *)calloc(1,sizeof(Link));
      
      if (tmplink == NULL) {
	fprintf(stderr,"Can't allocate memory for link\n");
	exit(0);
      }

      tmplink->id      = atoi(tokens[0]);
      tmplink->chr     = tokens[1];
      tmplink->start   = atoi(tokens[2]);
      tmplink->end     = atoi(tokens[3]);
      tmplink->hstart  = atoi(tokens[4]);
      tmplink->chainid = atoi(tokens[5]);

      tmplink->next = NULL;

      if (tmplink->chainid != chainid) {
	//fprintf(stdout,"Chain id is wrong %d %d\n",chainid,tmplink->chainid);

	*newlink = tmplink;
	
	return link;
      }

      if (count == 0) {
	link = tmplink;
      } else {
	prev->next = tmplink;
      }
      prev = tmplink;
    }

    count++;
  }
  return NULL;
}

void print_feature(DNAAlignFeature *daf,int analysis) {
    Vector *newf;
    BaseAlignFeature *newfp;
    RawContig        *contig;

    printf("%d %d\n",FeaturePair_getStart(daf),FeaturePair_getEnd(daf));
    newf = DNAAlignFeature_transformToRawContig(daf);

    //printf("Num elements %d\n",Vector_getNumElement(newf));

    newfp = (BaseAlignFeature *)Vector_getElementAt(newf,0);

    contig = (RawContig *)DNAAlignFeature_getContig(newfp);

    printf("\\N\t%qd\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t0\tNULL\t0\t%s\n",
	   RawContig_getDbID(contig),
	   FeaturePair_getStart(newfp),
	   FeaturePair_getEnd(newfp),
	   FeaturePair_getStrand(newfp),
	   FeaturePair_getHitStart(newfp),
	   FeaturePair_getHitEnd(newfp),
	   FeaturePair_getHitStrand(newfp),
	   FeaturePair_getHitSeqName(newfp),
	   analysis,
	   BaseAlignFeature_getCigarString(newfp));

    Vector_free(newf);

  }



  void pretty_print(char *seq1, char *seq2,char *id1,char *id2, int linelen) {
    int i;
    int j;

    int tmpi;

    if (strlen(seq1) != strlen(seq2)) {
      fprintf(stderr,"Sequences not the same length - can't pretty print them\n%s\n%s\n",seq1,seq2);
      return;
    }

    i = 0;

    printf("\n");

    while (i < strlen(seq1)) {

      if (i % linelen == 0) {

	j = 0;

	printf("%-25s ",id1);

	tmpi = i;

	while (j < linelen && tmpi < strlen(seq1)) {
	  printf("%c",seq1[tmpi]);
	  j++;
	  tmpi++;
	}

	printf("\n");

	tmpi = i;

	j = 0;

	printf("%-25s "," ");

	while (j < linelen && tmpi < strlen(seq1)) {
	  if (seq1[tmpi] == seq2[tmpi]) {
	    printf("|");
	  } else {
	    printf(" ");
	  }
	  j++;
	  tmpi++;
	}

	
	printf("\n");

	j = 0;

	printf("%-25s ",id2);

	while (j < linelen && i < strlen(seq2)) {
	  printf("%c",seq2[i]);
	  j++;
	  i++;
	}

	printf("\n");
	
	printf("\n");
      }
    }
  }

void parse_chr(char *chr, int *start, int *end) {

  int i = 0;
  
  int s_start = NULL;
  int s_end   = NULL;
  int e_start = NULL;
  int e_end   = NULL;
  int c_end   = NULL;

  char *start_str;
  char *end_str;

  char c;

  while (i < strlen(chr)) {
    c = chr[i];
    if (c == '.') {
      c_end = i-1;
      s_start = i+1;
    }
    if (c == '-') {
      s_end = i-1;
      e_start = i+1;
    }
    i++;
  }

  e_end = strlen(chr) - 1;


  if (s_end != NULL) {
    start_str = (char *)calloc(1,(s_end - s_start + 2)*sizeof(char));
    end_str = (char *)calloc(1,(e_end - e_start + 2)*sizeof(char));
    
    i = 0;

    while (i < strlen(chr)) {
      if (i >= s_start && i <= s_end) {
	start_str[i-s_start] = chr[i];
      }
      if (i >= e_start && i <= e_end) {
	end_str[i-e_start] = chr[i];
      }
      i++;
    }

    start_str[s_end - s_start + 2] = '\0';
    end_str[e_end - e_start + 2]   = '\0';

    chr[s_start-1] = '\0';

    *start = atoi(start_str);
    *end   = atoi(end_str);

    printf("Final %s %d %d\n",chr,*start,*end);
  }


}

void get_alignment(DNAAlignFeature *daf,DBAdaptor *dba1,DBAdaptor *dba2,  char *chr,int  chrstart,int  chrend,char *chr_filter) {

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

    prev = NULL;

   


    hitid = (char *)calloc(1,(strlen(DNAAlignFeature_getHitSeqName(daf))+1)*sizeof(char));
    strcpy(hitid,DNAAlignFeature_getHitSeqName(daf));

    parse_chr(hitid,&tmpstart,&tmpend);

    start = FeaturePair_getStart(daf);
    end   = FeaturePair_getEnd(daf);
    
    hstart = FeaturePair_getHitStart(daf);
    hend   = FeaturePair_getHitEnd(daf);



    printf("DNA align feature: %d-%d id %d-%d " IDFMTSTR " %s\n", DNAAlignFeature_getStart(daf),
	   DNAAlignFeature_getEnd(daf),
	   DNAAlignFeature_getHitStart(daf),
	   DNAAlignFeature_getHitEnd(daf),
	   DNAAlignFeature_getDbID(daf),
	   DNAAlignFeature_getHitSeqName(daf));

    printf("Chr %s %d %d\n",hitid,tmpstart,tmpend);
    
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

    sprintf(id1,"%s.%d-%d",chr,start+chrstart-1,end+chrstart-1);
    sprintf(id2,"%s.%d-%d",hitid,
	    hstart,
	    hend);

    ungapped = DNAAlignFeature_getUngappedFeatures((BaseAlignFeature *)daf);      
    //Vector_sort(ungapped, SeqFeature_startCompFunc);

    for (j = 0; j < Vector_getNumElement(ungapped); j++) {
      DNAAlignFeature *ugf = Vector_getElementAt(ungapped,j);

      tmpstr1 = Slice_getSubSeq(slice1,FeaturePair_getStart(ugf),FeaturePair_getEnd(ugf),1);
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

	    //printf("*************** %s\n",hgappy);
	    //printf("*************** %s\n",pad);
	    
	    alignstr1 = (char *)StrUtil_appendString(alignstr1,pad);
	    alignstr1 = (char *)StrUtil_appendString(alignstr1,tmpstr1);
	    alignstr2 = (char *)StrUtil_appendString(alignstr2,hgappy);
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
	    
	    alignstr1 = (char *)StrUtil_appendString(alignstr1,pad);
	    alignstr1 = (char *)StrUtil_appendString(alignstr1,tmpstr1);
	    alignstr2 = (char *)StrUtil_appendString(alignstr2,hgappy);
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

    
    pretty_print(alignstr1,alignstr2,id1,id2,80);

}
