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

/////#### DO NOT CALL THIS ERROR.H --> INTERFERES WITH MINGW ERROR.H FILE

#ifndef __dnlKernel_Error_h__
#define __dnlKernel_Error_h__

#define maxWarnings 20 // Maximum number of warnings allowed before stopping the program

void fatalError(const char *st, ...);

#endif