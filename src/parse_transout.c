#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "datamodel.h"

char **read_tokens(FILE *file,char token, int *ntok);

int main(int argc, char *argv[]) {

  FILE *file = fopen(argv[1],"r");

  int ntok;

  
  char **tokens;

  while (tokens = read_tokens(file,'\t',&ntok)) {
    
    if (ntok > 4) {
      if (strcmp(tokens[1],"Human") == 0 &&
	  strcmp(tokens[4],argv[2]) == 0) {
	
	int i = 0;
	
	while (i < ntok) {
	  printf("%s\t",tokens[i]);
	  i++;
	}
	printf("\n");
      } 
      int i = 0;
      
      while (i < ntok) {
	free(tokens[i]);
	i++;
      }
      free(tokens);
    }
  }
}
      
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  







