#include <string.h>
#include "chain.h"

/* TODO
 * 
 * Display sequences as proper multiple alignment
 *   - piece features together with padding sequence
 *   - tie hit features together with common chromosome name
 *   - deal with gaps in query sequence
 *
 * Display exons as sequence features in pogview
 *   - input features as gff
 *   - display them under a specified sequence - dealing with gaps
 *   - exclude the features from any calculations
 * 
 * Not display consensus
 * Exclude some sequences from calculations
 * Display coords when base is selected
 * Select a base
 * When zooming keep selected base centred
 * Show scale
 * Analyse msa for nmers or patterns 
 *    Find all nmers - nmer is variable
 *    Find common nmers within n aligned base pairs
 *    Find all non conserved base pairs.
 */
int main(int argc, char *argv[]) {
  DBAdaptor *dba1;
  DBAdaptor *dba2;

  FILE *file;

  char *filestub = argv[1];
  char *filename = NULL;
  int  filenum = 1;

  char *chrname;
  int   chrstart;
  int   chrend;

  char *name    = argv[1];
  
  DBEntry *dbentry;
  FeatureSet *tmpset;

  StringHash *hash = StringHash_new(STRINGHASH_SMALL);

  initEnsC();
  
  printf("Opening connection...");

  dba1 = DBAdaptor_new("18.103.3.97","root",NULL,"homo_sapiens_core_18_34",3306,NULL);
  dba2 = DBAdaptor_new("18.103.3.97","root",NULL,"mus_musculus_core_18_30",3306,NULL);
  //dba1 = DBAdaptor_new("localhost","root",NULL,"homo_sapiens_core_18_34",3306,NULL);
  //dba2 = DBAdaptor_new("localhost","root",NULL,"mus_musculus_core_18_30",3306,NULL);

  printf("Done\n");

  {


    DBEntryAdaptor *dbe = DBAdaptor_getDBEntryAdaptor(dba1);
    SliceAdaptor   *sa1 = DBAdaptor_getSliceAdaptor(dba1);

    Slice *slice1;
    Seq   *head;

    Seq   *all = NULL;
    Seq   *tmp;

    Vector *features;

    IDType *idtype;

    int analysis;
    int i;
    int j;

    FeaturePair *tmpfeat;
    char **keys;
    dbentry = DBEntry_new();
    
    DBEntry_setDisplayId(dbentry,name);
    DBEntry_setDbName(dbentry,argv[2]);

    if ((idtype = DBEntryAdaptor_exists(dbe,dbentry)) != NULL) {

      printf("Id type is %d\n",(int)idtype);

      dbentry = DBEntryAdaptor_fetchByDbID(dbe, (IDType)idtype);

      if (dbentry != NULL) {
	printf("Entry is %d %s %d %d %s %s %s\n", DBEntry_getPrimaryId(dbentry),
	       DBEntry_getDisplayId(dbentry),
	       DBEntry_getVersion(dbentry),
	       DBEntry_getRelease(dbentry),
	       DBEntry_getDbName(dbentry),
	       DBEntry_getEnsemblId(dbentry),
	       DBEntry_getObjectType(dbentry));


	if (DBEntry_getEnsemblId(dbentry) != NULL &&
	    DBEntry_getObjectType(dbentry) != NULL) {

	  
	  if (strcmp(DBEntry_getObjectType(dbentry),"Translation") == 0) {
	    int size;

	    size = 2000;

	    slice1 = SliceAdaptor_fetchByTranslationId(sa1,
						       atoi(DBEntry_getEnsemblId(dbentry)),
						       &size);
	    
	    slice1 = SliceAdaptor_fetchByChrStartEnd(sa1,
						     Slice_getChrName(slice1),
						     Slice_getChrStart(slice1) - 2000,
						     Slice_getChrStart(slice1));

	    chrname  = Slice_getChrName(slice1);
	    chrstart = Slice_getChrStart(slice1);
	    chrend   = Slice_getChrEnd(slice1);

	    printf("Slice %s %d %d\n",Slice_getChrName(slice1),
		   Slice_getChrStart(slice1),
		   Slice_getChrEnd(slice1));
	    
	    
	    features =  Slice_getAllDNAAlignFeatures(slice1,"",NULL);

	    printf("Num features %d\n",Vector_getNumElement(features));
    
	    Vector_sort(features, SeqFeature_startCompFunc);
	    
	    for (i=0;i<Vector_getNumElement(features);i++) {


	      DNAAlignFeature *daf = Vector_getElementAt(features,i);
	      
	      analysis = Analysis_getDbID(DNAAlignFeature_getAnalysis(daf));

	      if (analysis == 107) {
	      printf("\nFeature %d\t%d\t%d\t%s\t%d\t%d\t%d\n",FeaturePair_getStart(daf),
		     FeaturePair_getEnd(daf),
		     FeaturePair_getStrand(daf),
		     FeaturePair_getHitSeqName(daf),
		     FeaturePair_getHitStart(daf),
		     FeaturePair_getHitEnd(daf),
		     FeaturePair_getHitStrand(daf));
	      

	      if (StringHash_contains(hash,FeaturePair_getHitSeqName(daf)) == NULL) {
		tmpset = (FeatureSet *)calloc(1,sizeof(FeatureSet));
		StringHash_add(hash,FeaturePair_getHitSeqName(daf),tmpset);
	      } else {
		tmpset = (FeatureSet *)StringHash_getValue(hash,FeaturePair_getHitSeqName(daf));
	      }

	      FeatureSet_addFeature(tmpset,daf);
	      
	    }

	    i = 0;

	    keys = StringHash_getKeys(hash);

	    while (i < StringHash_getNumValues(hash)) {
	      tmpset = (FeatureSet *)StringHash_getValue(hash,keys[i]);

	      j = 0;

	      while (j < FeatureSet_getNumFeature(tmpset)) {
		tmpfeat = FeatureSet_getFeatureAt(tmpset,j);

		head = get_alignment(tmpfeat,dba1,dba2,Slice_getChrName(slice1),
				     Slice_getChrStart(slice1),
				     Slice_getChrEnd(slice1),NULL);
		
		printf("\nFeature %d\t%d\t%d\t%s\t%d\t%d\t%d\n",FeaturePair_getStart(daf),
		       FeaturePair_getEnd(daf),
		       FeaturePair_getStrand(daf),
		       FeaturePair_getHitSeqName(daf),
		       FeaturePair_getHitStart(daf),
		       FeaturePair_getHitEnd(daf),
		       FeaturePair_getHitStrand(daf));
		
		get_alignment(daf,dba1,dba2,Slice_getChrName(slice1),
			      Slice_getChrStart(slice1),
			      Slice_getChrEnd(slice1),NULL);
	      }
	      i++;
	    }
	    }
	      
	    /*
	    filename = (char*)calloc(1,(strlen(filestub) + intlen(filenum) + 5)*sizeof(char));
	    sprintf(filename,"%s.%d.fa",filestub,filenum);
	    
	    file = fopen(filename,"w");
	    
	    print_fasta(file,all,80);
	    
	    fclose(file);

	    */
	  }
	}
      }
    }
  }

  return 1;
}


