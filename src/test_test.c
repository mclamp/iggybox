/* test_test.c  - test if the test macros work */

#include <stdio.h>
#include "test.h"

int tests_run = 0;

int foo = 7;
int bar = 4;
 
 static char * test_foo() {
   test_assert("error, foo != 7", foo == 7);
   return 0;
 }
 
static char * test_bar() {
  test_assert("error, bar != 5", bar == 5);
  return 0;
}
 
 static char * all_tests() {
   test_run(test_foo);
   test_run(test_bar);
   return 0;
 }
 
int main(int argc, char **argv) {
  char *result = all_tests();
  if (result != 0) {
    printf("%s\n", result);
  }
  else {
    printf("ALL TESTS PASSED\n");
  }
  printf("Tests run: %d\n", tests_run);
  
  return result != 0;
}
