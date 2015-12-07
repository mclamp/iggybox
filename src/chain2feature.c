#include "chain.h"
    
int main(int argc, char *argv[]) {
  DBAdaptor *dba1;
  DBAdaptor *dba2;

  FILE *chainfile;
  FILE *linkfile;

  Link  *links;
  Link  *newlink;

  Link   *prev;

  int rand_flag = 0;
  char *randstr;

  DNAAlignFeature *daf = NULL;

  printf("Opening %s\n",argv[1]);
  chainfile = fopen(argv[1],"r");
  linkfile  = fopen(argv[2],"r");

  initEnsC();

  printf("Opening connection...");

  dba1 = DBAdaptor_new("localhost","root",NULL,"homo_sapiens_core_18_34",3306,NULL);
  dba2 = DBAdaptor_new("localhost","root",NULL,"mus_musculus_core_18_30",3306,NULL);

  printf("Done\n");

  {

    Chain * tmp = NULL;
    Link  * tmplink;
    Link  * prevlink = NULL;

    SliceAdaptor *sa1 = DBAdaptor_getSliceAdaptor(dba1);

    while ((tmp = read_chain(chainfile)) != NULL) {
      rand_flag = 0;
      links  = read_links(linkfile,tmp->chainid,&prev);      
      tmplink = links;

      
      if (prevlink != NULL) {
	prevlink->next = links;
	links = prevlink;
      }

      prevlink = prev;
      
      {
	int tmpcoord;
	
	tmp->start++;
	tmp->end++;

	if (tmp->strand == -1) {
	  tmpcoord = tmp->hstart;
	  tmp->hstart = tmp->chr2len - tmp->hend;
	  tmp->hend   = tmp->chr2len - tmpcoord;
	}
      }

      if (strlen(tmp->chr1) > 5) {
	randstr = StrUtil_copyNString(&randstr,tmp->chr1,strlen(tmp->chr1)-6,6);
	
	printf("Rand str %s\n",randstr);
	
	if (strcmp("random",randstr) == 0) {
	  rand_flag = 1;
	}
      }
      if (strlen(tmp->chr2) > 5) {
	randstr = StrUtil_copyNString(&randstr,tmp->chr2,strlen(tmp->chr2)-6,6);
	
	printf("Rand str %s\n",randstr);
	
	if (strcmp("random",randstr) == 0) {
	  rand_flag = 1;
	}
      }
      
      if (rand_flag == 0) {
      
	Slice *slice1 = SliceAdaptor_fetchByChrStartEnd(sa1,tmp->chr1,tmp->start,tmp->end);
	
	DNAAlignFeature *tmpfp = NULL;

	Link *prev       = NULL;

	Vector *features = NULL;

	printf("Chain is %d %d %d %d %d\n",tmp->start,tmp->end,tmp->strand,tmp->hstart,tmp->hend);

	while (links != NULL) {

	  links->start++;
	  links->hend = links->hstart + (links->end - links->start);
	  
	  if (tmp->strand == -1) {
	    links->hend   = tmp->chr2len - links->hstart;
	    links->hstart = links->hend - (links->end - links->start);
	  } else {
	    links->hstart++;
	    links->hend++;
	  }

	  tmpfp = DNAAlignFeature_new();

	  DNAAlignFeature_setContig(tmpfp,slice1); 
	  DNAAlignFeature_setStart     (tmpfp,links->start);
	  DNAAlignFeature_setEnd       (tmpfp,links->end);
	  DNAAlignFeature_setStrand    (tmpfp,1);
	  DNAAlignFeature_setHitStart  (tmpfp,links->hstart);
	  DNAAlignFeature_setHitEnd    (tmpfp,links->hend);
	  DNAAlignFeature_setHitSeqName(tmpfp,tmp->chr2);
	  DNAAlignFeature_setHitStrand (tmpfp,tmp->strand);

	  /*
	  printf("Feature\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t107\n",
		 FeaturePair_getStart(tmpfp),
		 FeaturePair_getEnd(tmpfp),
		 FeaturePair_getStrand(tmpfp),
		 FeaturePair_getHitStart(tmpfp),
		 FeaturePair_getHitEnd(tmpfp),
		 FeaturePair_getHitStrand(tmpfp),
		 FeaturePair_getHitSeqName(tmpfp));
	  */
	  if (prev != NULL) {
	    
	    if (tmp->strand == -1) {
	      if ((prev->end - links->start) < -1 && (prev->hstart - links->hend) > 1) {
		
		daf = DNAAlignFeature_new();

		BaseAlignFeature_parseFeatures((BaseAlignFeature *)daf,features);

		DNAAlignFeature_setContig(daf,slice1); 
		//DNAAlignFeature_setAnalysis(daf,analysis);
		DNAAlignFeature_setStart(daf,daf->start - slice1->start + 1);
		DNAAlignFeature_setEnd  (daf,daf->end   - slice1->start + 1);

		//get_alignment(daf,dba1,dba2,tmp->chr1,tmp->start,tmp->end);
		print_feature(daf);

		DNAAlignFeature_free(daf);
		daf = NULL;
		Vector_free(features);
		features = NULL;

		// reset the feature array.
	      }
	    } else if ((links->start - prev->end) > 1 && (links->hstart - prev->hend) > 1) {


	      daf = DNAAlignFeature_new();

	      BaseAlignFeature_parseFeatures((BaseAlignFeature *)daf,features);
	      DNAAlignFeature_setContig(daf,slice1); 
	      
	      DNAAlignFeature_setStart(daf,daf->start - slice1->start + 1);
	      DNAAlignFeature_setEnd  (daf,daf->end   - slice1->start + 1);

	      //get_alignment(daf,dba1,dba2,tmp->chr1,tmp->start,tmp->end);	      
	      print_feature(daf);

	      DNAAlignFeature_free(daf);
	      daf = NULL;
	      Vector_free(features);
	      features = NULL;
	      
	    }
	  }
	  
	  // Add the feature into the array
	  // set the prev

	  if (features == NULL) {
	    features = Vector_new();
	  }

	  Vector_addElement(features,tmpfp);

	  prev  = links;
	  links = links->next;
	}
	
        
	daf = DNAAlignFeature_new();

	BaseAlignFeature_parseFeatures((BaseAlignFeature *)daf,features);
	DNAAlignFeature_setContig(daf,slice1); 

	//DNAAlignFeature_setAnalysis(daf,analysis);
	DNAAlignFeature_setStart(daf,daf->start - slice1->start + 1);
	DNAAlignFeature_setEnd  (daf,daf->end   - slice1->start + 1);

	//get_alignment(daf,dba1,dba2,tmp->chr1,tmp->start,tmp->end);	
	print_feature(daf);



	DNAAlignFeature_free(daf);
	daf = NULL;
	Vector_free(features);

	free(slice1);
      }

      free(tmp->chr1);
      free(tmp->chr2);
      free(tmp);

      while (tmplink != NULL) {
	newlink = tmplink->next;
	free(tmplink->chr);
	free(tmplink);
	tmplink = newlink;
      }

    }
   
    return 1;
  }

    
}
