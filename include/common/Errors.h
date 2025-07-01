/**********************************************************************************
 *                                                                                *
 *  DynELA Finite Element Code v.4.0                                              *
 *  by Olivier PANTALE                                                            *
 *  Olivier.Pantale@enit.fr                                                       *
 *                                                                                *
 *********************************************************************************/
//@!CODEFILE = DynELA-H-file
//@!BEGIN = PRIVATE

/*
  \file Errors.h
  Declaration of the errors for the DynELA Finite Element code

  This file defines the errors and warning that are used in the DynELA Finite Element code.
  The functions defined in this files serves to display messages and warnings.
  \ingroup dnlKernel
*/

#ifndef __dnlKernel_Error_h__
#define __dnlKernel_Error_h__

#define maxWarnings 20 // Maximum number of warnings allowed before stopping the program

void fatalError(const char *st, ...);
void internalFatalError(const char *st, ...);
void warning(const char *st, ...);
void warningQuestion(char *st, ...);

extern unsigned int totalWarningMessagesDisplayed;
/*
  \def lineERROR
  Display an error message.

  This inline function display an error message including the name of the source file and the line number where the error has occured. \n
  This is a very usefull method for debugging. \n
  The execution of the program is not stopped, this is just a display of the error.
*/
#define lineERROR "source file %s line %d", __FILE__, __LINE__

/*
  \def fatalErrorLine
  Display an error message and stop the execution of the program.

  This inline function display an error message including the name of the source file and the line number where the error has occured. \n
  This is a very usefull method for debugging. \n
  The execution of the program is stopped.
*/
#define fatalErrorLine                                                          \
  std::cerr << "source file " << __FILE__ << " line " << __LINE__ << std::endl; \
  fatalError

/*
  \def internalFatalErrorLine
  Display an error message and stop the execution of the program.

  This inline function display an error message including the name of the source file and the line number where the error has occured. \n
  This function has to be used for internal fatal errors. \n
  This is a very usefull method for debugging. \n
  The execution of the program is stopped.
*/
#define internalFatalErrorLine                                                  \
  std::cerr << "source file " << __FILE__ << " line " << __LINE__ << std::endl; \
  internalFatalError

/*
  \def warningLine
  Display a warning message.

  This inline function display a warning message including the name of the source file and the line number where the error has occured. \n
  This is a very usefull method for debugging. \n
  The execution of the program is not stopped, this is just a display of the warning.
*/
#define warningLine                                                             \
  std::cerr << "source file " << __FILE__ << " line " << __LINE__ << std::endl; \
  warning

#endif
