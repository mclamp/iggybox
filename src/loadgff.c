#include "chain.h"

void  *read_gff (FILE *file,SliceAdaptor *sa) {

  char **tokens   = NULL;
  int ntok;

  DNAAlignFeature *tmpfp;
  Slice *slice;

  while ((tokens = get_tokens(file,'\t',&ntok)) != NULL) {
    tmpfp = DNAAlignFeature_new();
  
    DNAAlignFeature_setStart     (tmpfp,atoi(tokens[3]));
    DNAAlignFeature_setEnd       (tmpfp,atoi(tokens[4]));
    DNAAlignFeature_setStrand    (tmpfp,atoi(tokens[6]));
    DNAAlignFeature_setHitStart  (tmpfp,atoi(tokens[9]));
    DNAAlignFeature_setHitEnd    (tmpfp,atoi(tokens[10]));
    DNAAlignFeature_setHitSeqName(tmpfp,tokens[8]);
    DNAAlignFeature_setHitStrand (tmpfp,atoi(tokens[11]));

    slice = SliceAdaptor_fetchByChrStartEnd(sa,tokens[0],
					    DNAAlignFeature_getStart(tmpfp),
					    DNAAlignFeature_getEnd(tmpfp));

    DNAAlignFeature_setContig(tmpfp,slice);

    print_feature(tmpfp,1160);

    
  }

}
	  
int main(int argc, char *argv[]) {
  FILE *file;
  DBAdaptor *dba;
  SliceAdaptor *sa;

  file = fopen(argv[1],"r");

  dba = DBAdaptor_new("18.103.3.97","root",NULL,"homo_sapiens_core_18_34",3306,NULL);
  sa  = DBAdaptor_getSliceAdaptor(dba);

  read_gff(file,sa);

  return 0;
}
