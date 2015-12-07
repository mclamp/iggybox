#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

char * read_line(FILE *file);
char * next_field(char *str, int *pos);
char **read_words(int verbose,char *str);

typedef struct _Gerp {
  int    coord;
  double expected;
  double observed;
  struct _Gerp *next;
} Gerp;


int main(int argc, char * argv[]) {
    FILE *file;

    file = fopen(argv[1],"r");


    char c;

    char  *line;
    char **words;

    Gerp *masterGerp;
    Gerp *tmpGerp;

    masterGerp = NULL;

    line = NULL;

    int count = 0;

    while ((line = read_line(file)) && line != NULL) {

      words = read_words(0,line);

      if (masterGerp == NULL) {
	masterGerp = (Gerp *)malloc(sizeof(Gerp));
	masterGerp->next = NULL;

	tmpGerp = masterGerp;
      } else {
	tmpGerp->next = (Gerp *)malloc(sizeof(Gerp));
	tmpGerp = tmpGerp->next;
	tmpGerp->next = NULL;
      }

      tmpGerp->coord = atoi(words[1]);
      tmpGerp->expected = atof(words[2]);
      tmpGerp->observed = atof(words[3]);

      tmpGerp->next = NULL;

      count++;

    }

    int numcols = count;

    srand(8);

    int *done;

    done = (int *)malloc(count*sizeof(int));

    int i = 0;

    while (i < count) {
      done[i] = 0;
      i++;
    }

    Gerp *gerps[count];

    i = 0;

    tmpGerp = masterGerp;

    while (i < count){ 
      gerps[i] = tmpGerp;
      tmpGerp = tmpGerp->next;
      i++;
    }
    int j = 1;

    while (numcols > 0) {
      
      int col = (int)((double)count*rand()*1.0/(1.0*RAND_MAX));

      while (done[col] == 1) {
	col = (int)((double)count*rand()*1.0/(1.0*RAND_MAX));
      }

      if (numcols % 1000 == 0) {
	fprintf(stderr,"Numcols %d\n",numcols);
	//printf("Numcols %d\n",numcols);
      }

      tmpGerp = gerps[col];

      printf("Lengths\t%d\t%f\t%f\n",j,tmpGerp->expected,tmpGerp->observed);

      done[col] = 1;

      numcols--;
      j++;
    }
}

char **read_words(int verbose,char *str) {

  char **words = NULL;
  char  *word  = NULL;

  int    count = 0;
  int    pos   = 0;

  if (verbose == 1) {
    printf("String %s\n",str);
  }

  words = (char **)malloc(sizeof(char *));
  
  while ((word = next_field(str,&pos)) != NULL) {
    if (verbose == 1) {
      printf("Field %s\n",word);
    }
    
    count++;
    
    words = (char **)realloc(words,sizeof(char *)*count);
    
    words[count-1] = word;
    
  }
  
  return words;
}
 
char * next_field(char *str, int *pos) {

  char *outstr = NULL;
  int   n      = 0;
  int   i      = *pos;

  while (i < strlen(str) && str[i]  != '\t') {
    if (outstr == NULL) {
      n = 0;
      outstr = (char *)malloc(sizeof(char)*(2));

      if (outstr == NULL) {
         printf("Can't allocate memory for outstr\n");
	 exit(0);
      }
      outstr[n] = str[i];
      n++;
      outstr[n] = '\0';

    } else {

      outstr = (char *)realloc(outstr,sizeof(char)*(strlen(outstr)+2));
      outstr[n] = str[i];
      n++;
      outstr[n] = '\0';
    }
    i++;
  }

  i++;
   
  *pos = i;

  return outstr;
}

char * read_line(FILE *file) {

  char *str = NULL;
  char  c;
  int   n = 0;

  while ((c = getc(file)) != EOF && c != '\n') {

    if (str == NULL) {

      str = (char *)malloc(sizeof(char)*(2));
      n = 0;
      str[n] = c;
      str[n+1]   = '\0';

    } else {

      str = (char *)realloc(str,sizeof(char)*(strlen(str)+2));

      str[n] = c;
      str[n+1]   = '\0';

    }

    n++;

  }

  return str;
}
