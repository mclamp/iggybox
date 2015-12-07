
#include "chain.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba1;
  DBAdaptor *dba2;

  char *chr     = argv[1];

  int  chrstart = atoi(argv[2]);
  int  chrend   = atoi(argv[3]);
  char *chr_filter = NULL;

  char *org = argv[5];

  if (argc == 7) {
     chr_filter = argv[6];
     printf("Setting chromosome filter to %s\n",chr_filter);
  } 
  initEnsC();
  
  printf("Opening connection...");

  
  dba1 = DBAdaptor_new("localhost","root",NULL,"homo_sapiens_core_18_34",3306,NULL);

  if (strcmp("mouse",org) == 0) {
    printf("Connecting to mouse\n");
    dba2 = DBAdaptor_new("localhost","root",NULL,"mus_musculus_core_18_30",3306,NULL);
  } else if (strcmp("dog",org) == 0) {
    printf("Connecting to dog\n");
    dba2 = DBAdaptor_new("lead","ensrw","ensembl","dog_feb_1",3306,NULL);
  }
  //dba2 = DBAdaptor_new("localhost","root",NULL,"homo_sapiens_core_18_34",3306,NULL);

  printf("Done\n");

  {
    SliceAdaptor *sa1 = DBAdaptor_getSliceAdaptor(dba1);

    Slice *slice1 = SliceAdaptor_fetchByChrStartEnd(sa1,chr,chrstart,chrend);
	
    Vector *features;

    int analysis;
    int i;

    features =  Slice_getAllDNAAlignFeatures(slice1,"",NULL);

    printf("Num features %d\n",Vector_getNumElement(features));
    
    Vector_sort(features, SeqFeature_startCompFunc);

    for (i=0;i<Vector_getNumElement(features);i++) {
      DNAAlignFeature *daf = Vector_getElementAt(features,i);

      analysis = Analysis_getDbID(DNAAlignFeature_getAnalysis(daf));
       
      if (analysis == atoi(argv[4])) {
        if (chr_filter == NULL || (chr_filter != NULL && strncmp(FeaturePair_getHitSeqName(daf),chr_filter,1)  == 0)) {
        	
	printf("\nFEAT %s\tblastz\tsimilarity\t%d\t%d\t%4.2f\t%d\t.\t%s\t%d\t%d\t%d\n",
	       FeaturePair_getSeqName(daf),
	       FeaturePair_getStart(daf),
	       FeaturePair_getEnd(daf),
	       FeaturePair_getScore(daf),
	       FeaturePair_getStrand(daf),
	       FeaturePair_getHitSeqName(daf),
	       FeaturePair_getHitStart(daf),
	       FeaturePair_getHitEnd(daf),
	       FeaturePair_getHitStrand(daf));
	
	get_alignment(daf,dba1,dba2,chr,chrstart,chrend,chr_filter);
        }
      }

    }

    return 1;
  }

}
