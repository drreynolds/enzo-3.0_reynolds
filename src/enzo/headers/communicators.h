/***********************************************************************
/
/  COMMUNICATORS 
/
/  written by: Samuel Skillman 
/  date:       June, 2012
/  modified1:
/
/  PURPOSE:
/
/    This file is dual-purposed:
/        1) read with    DEFINE_STORAGE defined for the (single) definition
/        2) read without DEFINE_STORAGE defined for external linkage
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif

#ifndef COMMUNICATORS__
#define COMMUNICATORS__

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

#ifdef USE_MPI
EXTERN MPI_Comm EnzoTopComm;
#endif

#endif /* COMMUNICATORS__ */
