#include "chain.h"
#include "CigarStrUtil.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba1;
  DBAdaptor *dba2;

  char *chr     = argv[1];

  int  chrstart = atoi(argv[2]);
  int  chrend   = atoi(argv[3]);
  
  initEnsC();
  
  printf("Opening connection...");

  //dba1 = DBAdaptor_new("localhost","root",NULL,"homo_sapiens_core_18_34",3306,NULL);
  //dba2 = DBAdaptor_new("localhost","root",NULL,"mus_musculus_core_18_30",3306,NULL);

  dba1 = DBAdaptor_new("lead","ensrw","ensembl","dog_feb_1",3306,NULL);
  dba2 = DBAdaptor_new("localhost","root",NULL,"homo_sapiens_core_18_34",3306,NULL);

  printf("Done\n");

  {
    SliceAdaptor *sa1 = DBAdaptor_getSliceAdaptor(dba1);
    SliceAdaptor *sa2 = DBAdaptor_getSliceAdaptor(dba2);
    
    Slice *slice1 = SliceAdaptor_fetchByChrStartEnd(sa1,chr,chrstart,chrend);
    Slice *slice2;
	
    Vector *features;

    int analysis;
    int i;

    char *tmpname;
    char *newcigar;

    char *chr;
    char *hitid;

    int tmpstart = NULL;
    int tmpend   = NULL;
    
    char *seqname;
    char *hitseqname;

    int cigarlen;

    features =  Slice_getAllDNAAlignFeatures(slice1,"",NULL);

    printf("Num features %d\n",Vector_getNumElement(features));
    
    Vector_sort(features, SeqFeature_startCompFunc);

    for (i=0;i<Vector_getNumElement(features);i++) {
      DNAAlignFeature *daf = Vector_getElementAt(features,i);

      analysis = Analysis_getDbID(DNAAlignFeature_getAnalysis(daf));
      
      if (analysis == atoi(argv[4])) {
	
	printf("\nFeature %d\t%d\t%d\t%s\t%d\t%d\t%d\n",FeaturePair_getStart(daf),
	       FeaturePair_getEnd(daf),
	       FeaturePair_getStrand(daf),
	       FeaturePair_getHitSeqName(daf),
	       FeaturePair_getHitStart(daf),
	       FeaturePair_getHitEnd(daf),
	       FeaturePair_getHitStrand(daf));

	seqname = FeaturePair_getSeqName(daf);
	hitseqname = FeaturePair_getHitSeqName(daf);

	tmpstart  = FeaturePair_getStart(daf);
	tmpend    = FeaturePair_getEnd(daf);
	tmpname   = FeaturePair_getSeqName(daf);

	FeaturePair_setStart(daf,FeaturePair_getHitStart(daf));
	FeaturePair_setEnd  (daf,FeaturePair_getHitEnd(daf));

	FeaturePair_setHitStart(daf,tmpstart);
	FeaturePair_setHitEnd  (daf,tmpend);

	tmpstart = NULL;
	tmpend   = NULL;

	hitseqname = (char *)calloc(1,(strlen(DNAAlignFeature_getHitSeqName(daf))+1)*sizeof(char));
	strcpy(hitseqname,FeaturePair_getHitSeqName(daf));

	parse_chr(hitseqname,&tmpstart,&tmpend);

	if (tmpstart != NULL && 
	    tmpend   != NULL) {
	  
	  FeaturePair_setStart(daf,FeaturePair_getStart(daf) + tmpstart - 1);
	  FeaturePair_setEnd  (daf,FeaturePair_getEnd(daf) + tmpstart - 1);
	}


	printf("    Cigar %s\n",BaseAlignFeature_getCigarString(daf));

	cigarlen = strlen(BaseAlignFeature_getCigarString(daf));

	newcigar = CigarStrUtil_flip(BaseAlignFeature_getCigarString(daf),cigarlen);
	//newcigar = BaseAlignFeature_getCigarString(daf);
	
	if (FeaturePair_getHitStrand(daf) == -1) {
	  
	  newcigar = CigarStrUtil_reverse(newcigar,strlen(newcigar));


	}

	newcigar[cigarlen] = '\0';

	DNAAlignFeature_setCigarString(daf,newcigar);

	printf("New Cigar %s\n",newcigar);

	printf("Slice %s %d %d\n",hitseqname,FeaturePair_getStart(daf),
	       FeaturePair_getEnd(daf));

	slice2 = SliceAdaptor_fetchByChrStartEnd(sa2,hitseqname,
						 FeaturePair_getStart(daf),
						 FeaturePair_getEnd(daf));

	printf("\nFeature %d\t%d\t%d\t%s\t%d\t%d\t%d\n",FeaturePair_getStart(daf),
	       FeaturePair_getEnd(daf),
	       FeaturePair_getStrand(daf),
	       FeaturePair_getHitSeqName(daf),
	       FeaturePair_getHitStart(daf),
	       FeaturePair_getHitEnd(daf),
	       FeaturePair_getHitStrand(daf));



	FeaturePair_setStart(daf,1);
	FeaturePair_setEnd  (daf,FeaturePair_getEnd(daf) - FeaturePair_getStart(daf) + 1);
	FeaturePair_setHitSeqName(daf,seqname);

	FeaturePair_setContig(daf,slice2);

	print_feature(daf,analysis);
      }

    }

    return 1;
  }

}
