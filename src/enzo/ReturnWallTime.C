/***********************************************************************
/
/  RETURN THE CURRENT ELAPSED CPU TIME
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/

#ifdef USE_MPI
#include "communicators.h"
#endif
#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"

double ReturnWallTime()
{
#ifdef USE_MPI

  return MPI_Wtime();

#else /* USE_MPI */

  return (double)(clock()) / ((double)CLOCKS_PER_SEC);
  
#endif /* USE_MPI */
}
