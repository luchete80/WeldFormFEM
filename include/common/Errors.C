/**********************************************************************************
 *                                                                                *
 *  DynELA Finite Element Code v.4.0                                              *
 *  by Olivier PANTALE                                                            *
 *  Olivier.Pantale@enit.fr                                                       *
 *                                                                                *
 *********************************************************************************/
//@!CODEFILE = DynELA-C-file
//@!BEGIN = PRIVATE

/*
  \file Errors.C
  Definition of the errors for the DynELA FEM code

  This file defines the errors and warning that are used in the DynELA FEM code.
  The functions defined in this files serves to display messages and warnings.
  \ingroup dnlKernel
*/

#include <Errors.h>
#include <iostream>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>

unsigned int totalWarningMessagesDisplayed = 1; // Total number of warnings actually displayed. This later is initialized to 1.

/*
  Displays a warning message

  This function displays a warning message on the stderr flux. \n
  Each call to this function will increment the internal counter \ref totalWarningMessagesDisplayed. \n
  When the maximal allowed value of \ref maxWarnings is reached, the execution of program is stopped.

  Example :
  \code
  warning("duplicate::operation(int i)", "the value of i is negative\ni=%d\n", i);
  \endcode

  - function_name Warning message to be displayed given as const char *
*/
//-----------------------------------------------------------------------------
void warning(const char *function_name, ...)
//-----------------------------------------------------------------------------
{
  va_list ap;
  char *fmt;
  va_start(ap, function_name);

  // Displays the warning message
  fprintf(stderr, "WARNING [%d] in function: %s\n",
          totalWarningMessagesDisplayed, function_name);

  // Displays the comment line
  fmt = va_arg(ap, char *);
  vfprintf(stderr, fmt, ap);
  fprintf(stderr, "\n");
  va_end(ap);

  // Force display now
  std::cerr.flush();

  // If total warning messages exceed the maximum allowed, stops the execution of the program
  if (totalWarningMessagesDisplayed++ >= maxWarnings)
    fatalError("Maximum number of warnings reached",
               "You may modify the source code of your application "
               "or increase the total number of allowed warning messages "
               "fixed here to %d\n",
               maxWarnings);
}

/*
  Displays a warning message and ask the user if he wants to continue the execution of the program

  This function displays a warning message on the stderr flux and ask the user if he wants to continue the execution of the program. \n
  Each call to this function will increment the internal counter \ref totalWarningMessagesDisplayed. \n
  When the maximal allowed value of \ref maxWarnings is reached, the execution of program is stopped.

  Exemple :
  \code
  warningQuestion("duplicate::operation(int i)", "the value of i is negative\ni=%d\n", i);
  \endcode
  - function_name Warning message to be displayed given as const char *
*/
//-----------------------------------------------------------------------------
void warningQuestion(char *function_name, ...)
//-----------------------------------------------------------------------------
{
  va_list ap;
  char *fmt;
  va_start(ap, function_name);

  // Displays the warning message
  fprintf(stderr, "WARNING [%d] in function: %s\n",
          totalWarningMessagesDisplayed, function_name);

  // Displays the comment line
  fmt = va_arg(ap, char *);
  vfprintf(stderr, fmt, ap);
  fprintf(stderr, "\n");
  va_end(ap);

  // Force display now
  std::cerr.flush();

  // If total warning messages exceed the maximum allowed, stops the execution of the program
  if (totalWarningMessagesDisplayed++ >= maxWarnings)
    fatalError("Maximum number of warnings reached",
               "You may modify the source code of your application "
               "or increase the total number of allowed warning messages "
               "fixed here to %d\n",
               maxWarnings);

  // Displays the question and ask for an answer
  std::cout << "\ndo you want to continue the parsing of the file (y/n)\n";
  std::cout.flush();
  std::fflush(stdin);

  // Exits if user ask for it
  if (getchar() != 'y')
  {
    exit(0);
  }
}

/*
  Displays an error message and stop the execution of the program.

  This function displays an error message on the stderr flux and the execution of program is stopped.

  Exemple :
  \code
  fatalError("duplicate::operation(int i)", "the value of i is negative\ni=%d\n", i);
  \endcode

  - function_name Error message to be displayed given as const char *
*/
//-----------------------------------------------------------------------------
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

/*
  Displays an error message and stop the execution of the program.

  This function displays an internal error message on the stderr flux and the execution of program is stopped.
  This function also displays a message asking to contact the main developper of the FEM code.

  Exemple :
  \code
  internalFatalError("duplicate::operation(int i)", "the value of i is negative\ni=%d\n", i);
  \endcode
  - function_name Error message to be displayed given as const char *
*/
//-----------------------------------------------------------------------------
void internalFatalError(const char *function_name, ...)
//-----------------------------------------------------------------------------
{
  va_list ap;
  char *fmt;
  va_start(ap, function_name);

  // Displays the fatal error message
  fprintf(stderr, "\nINTERNAL FATAL ERROR in function: %s\n", function_name);

  // Displays a line of comments
  fmt = va_arg(ap, char *);
  vfprintf(stderr, fmt, ap);
  fprintf(stderr, "\n");
  va_end(ap);

  // Displays the contact message
  fprintf(stderr, "You might contact main developper of this module\n");
  fprintf(stderr, "Or DynELA main developper\n");
  fprintf(stderr, "email this message to : pantale@enit.fr\n");

  // Stop the execution of the code
  fprintf(stderr, "\nprogram aborted\n");
  std::exit(-1);
}
