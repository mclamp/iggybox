#include "chain.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba1;

  char *chr     = argv[1];

  initEnsC();
  
  printf("Opening connection...");

  dba1 = DBAdaptor_new("18.103.3.97","root",NULL,"homo_sapiens_core_18_34",3306,NULL);

  printf("Done\n");

  {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba1);
    
    Slice *slice = SliceAdaptor_fetchByChrName(sa,chr);
    
    char * tmpstr;
    
    tmpstr = Slice_getSeq(slice);

    char prev = tmpstr[0];
    char c;
    int count = 1;
    int len = 0;
    int slen;
    tmpstr = Slice_getSeq(slice);

    slen = strlen(tmpstr);

    printf("Got seq\n");

    while (count < slen) {
      c = *tmpstr;
      tmpstr++;
      //  printf("Seq %c %d\n",c,count);

      if (c == 'N') {
	if (prev == 'N') {
	  len++;
	} else {
	  len = 1;
	}
      } else {
	if (prev == 'N') {
	  printf("Gap %d %d\n",count,len);
	  len = 0;
	}
      }

      prev = c;
      count++;
      
      if (count % 1000000 == 0) {
	printf("At base %d %c\n",count,c);
      }
    }
    if (len > 0) {
      printf("Gap %d %d\n",count,len);
    }
  }
}
