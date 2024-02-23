#include "Errors.h"

#include <iostream>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>

void fatalError(const char *function_name, ...)
//-----------------------------------------------------------------------------
{
  va_list ap;
  char *fmt;
  va_start(ap, function_name);

  // Displays the fatal error message
  fprintf(stderr, "\nFATAL ERROR in function: %s\n", function_name);

  // Displays a line of comments
  fmt = va_arg(ap, char *);
  vfprintf(stderr, fmt, ap);
  fprintf(stderr, "\n");
  va_end(ap);

  // Stop the execution of the code
  fprintf(stderr, "program aborted\n");
  exit(-1);
}