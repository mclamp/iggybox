/* test.h */

#define test_assert(message, test) do { if (!(test)) return message; } while (0)
#define test_run   (test)          do ( char *message = test(); tests_run++; if (message) return message; } while (0)

extern int tests_run;
