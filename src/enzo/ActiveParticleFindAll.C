/***********************************************************************
/
/  COMMUNICATION ROUTINE: SYNCHRONIZE GLOBAL ACTIVE PARTICLE LIST ACROSS
/                         PROCESSORS
/
/  written by: John Wise
/  date:       March 2009
/  modified1:  Nathan Goldbaum, March 2012
/              Porting to active particles
/
/
/  PURPOSE: Generates a list of all active particles in the simulation
/           with id = ActiveParticleIDToFind.
/
*************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <map>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "ActiveParticle.h"
#include "SortCompareFunctions.h"

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

int ActiveParticleFindAll(LevelHierarchyEntry *LevelArray[], ActiveParticleType** GlobalList, 
			  int &GlobalNumberOfActiveParticles, int ActiveParticleIDToFind)
{
  int i, level, type, ap_id, GridNum, LocalNumberOfActiveParticles, proc, buffer_size, 
    LocalNumberOfActiveParticlesOnThisLevel, header_size, element_size, count, offset;
  ActiveParticleType **LocalActiveParticlesOfThisType, **LocalActiveParticlesOnThisLevel;
  ActiveParticleType **ParticlesOnThisProc;
  
  HierarchyEntry **Grids;
  int NumberOfGrids, *NumberOfActiveParticlesInGrids;
  ActiveParticleType_info *ap_info;

  GlobalNumberOfActiveParticles = 0;
  LocalNumberOfActiveParticles = 0;

  for (type = 0; type < EnabledActiveParticlesCount; type++) {
      
    ap_info = EnabledActiveParticles[type];
    ap_id = ap_info->GetEnabledParticleID();
    
    if (ap_id == ActiveParticleIDToFind) {
      
      /* Traverse the hierarchy and generate a buffer of all of the active 
	 particles of this type on this processor */
      
      for (level = 0; level <= MaximumRefinementLevel; level++) {
	NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
	NumberOfActiveParticlesInGrids = new int[NumberOfGrids];
	LocalNumberOfActiveParticlesOnThisLevel = 0;

	/* In a first pass, find the number of active particles on each grid */
	for(GridNum = 0; GridNum < NumberOfGrids; GridNum++) {
	  
	  NumberOfActiveParticlesInGrids[GridNum] = Grids[GridNum]->GridData->
	    ReturnNumberOfActiveParticlesOfThisType(ActiveParticleIDToFind);
	  LocalNumberOfActiveParticlesOnThisLevel += NumberOfActiveParticlesInGrids[GridNum];
	} /* ENDFOR grids */
	
	offset = 0;
	
	if (LocalNumberOfActiveParticlesOnThisLevel > 0) {

	  /* If we've already found active particles, save the list */

	  if (LocalNumberOfActiveParticles != 0) 
	    ParticlesOnThisProc = LocalActiveParticlesOfThisType;
	    
	  LocalActiveParticlesOnThisLevel = new ActiveParticleType*[LocalNumberOfActiveParticlesOnThisLevel];

	  /* In a second pass, fill up the active particle list for this level*/
	  for(GridNum = 0; GridNum < NumberOfGrids; GridNum++) {
	    Grids[GridNum]->GridData->
	      AppendActiveParticlesToList(LocalActiveParticlesOnThisLevel,offset,ActiveParticleIDToFind);
	    offset += NumberOfActiveParticlesInGrids[GridNum];
	  } 

	  LocalActiveParticlesOfThisType = new ActiveParticleType*[offset];
	  
	  /* If we've already found active particles, copy the cached
	     list to the new one and delete the old list */
	  if (LocalNumberOfActiveParticles != 0) {
	    for (i = 0; i < LocalNumberOfActiveParticles; i++)
	      LocalActiveParticlesOfThisType[i] = ParticlesOnThisProc[i];
	    delete [] ParticlesOnThisProc;
	  }

	  /* Finally, append the new active particles to the list */
	  for(i = LocalNumberOfActiveParticles; i < offset; i++)
	    LocalActiveParticlesOfThisType[i] = LocalActiveParticlesOnThisLevel[i - LocalNumberOfActiveParticles];

	  LocalNumberOfActiveParticles += LocalNumberOfActiveParticlesOnThisLevel;
	    
	  /* Delete the list for this level */
	  delete [] LocalActiveParticlesOnThisLevel;
	
	} 
	
	delete [] Grids;
	delete [] NumberOfActiveParticlesInGrids;
	
      } /* ENDFOR level */
      
    } /* ENDIF id == id to search for

  } /* ENDFOR active particle type */
  
    /**************************************************/
    /*                                                */
    /* Share active particle counts on all processors */
    /*                                                */
    /**************************************************/

    Eint32 *nCount = NULL;
        
    if (NumberOfProcessors > 1) {
#ifdef USE_MPI
      
      nCount = new Eint32[NumberOfProcessors];
      
      MPI_Allgather(&LocalNumberOfActiveParticles, 1, MPI_INT, 
		    nCount, 1, MPI_INT,MPI_COMM_WORLD);
      
      for (i = 0; i < NumberOfProcessors; i++) {
	GlobalNumberOfActiveParticles += nCount[i];
      }
      
#endif /* USE_MPI */
    } /* ENDIF Number of processors > 1 */
    else {
      GlobalNumberOfActiveParticles = LocalNumberOfActiveParticles;
    }
    /**************************************************/
    /*                                                */
    /* Gather the active particles on all processors  */
    /*                                                */
    /**************************************************/
    
    if (GlobalNumberOfActiveParticles > 0) {

      Eint32 *displace = new Eint32[NumberOfProcessors];
  
      GlobalList = new ActiveParticleType*[GlobalNumberOfActiveParticles];
      
      if (NumberOfProcessors > 1) {
	
#ifdef USE_MPI
	/* Construct the MPI packed  buffer from the list of local particles*/
	Eint32 total_buffer_size, local_buffer_size, position = 0;
	int mpi_buffer_size;
	char *send_buffer, *recv_buffer;
	header_size = ap_info->return_header_size();
	element_size = ap_info->return_element_size();
	
	displace[0] = position;
	for (i = 1; i < NumberOfProcessors; i++) {
	  if (nCount[i-1] > 0)
	    position += header_size + nCount[i-1]*element_size;
	  displace[i] = position;
	}

	position = 0;

	local_buffer_size = header_size+LocalNumberOfActiveParticles*element_size;
	total_buffer_size = NumberOfProcessors*header_size+GlobalNumberOfActiveParticles*element_size;
	send_buffer = new char[local_buffer_size];
	recv_buffer = new char[total_buffer_size];
	
	ap_info->allocate_buffer(LocalActiveParticlesOfThisType, LocalNumberOfActiveParticles, send_buffer,
				 total_buffer_size, mpi_buffer_size, position, ap_id, -1);

	/* Share all data with all processors */

	MPI_Allgatherv(send_buffer, LocalNumberOfActiveParticles, MPI_PACKED,
		       recv_buffer, nCount, displace, MPI_PACKED, MPI_COMM_WORLD);

	/* Unpack MPI buffers, generate global active particles list */

	count = 0;
	for (proc = 0; proc < NumberOfProcessors; proc++) {
	  if (nCount[proc] > 0) {
	    ParticlesOnThisProc = new ActiveParticleType*[nCount[proc]];
	    buffer_size = header_size+(nCount[proc]*element_size);
	    ap_info->unpack_buffer(recv_buffer+displace[proc],buffer_size,nCount[proc],
				   ParticlesOnThisProc, count);
	    for (i = count; i < nCount[proc]; i++) 
	      GlobalList[count+i] = ParticlesOnThisProc[i];
	    count += nCount[proc];
	    delete [] ParticlesOnThisProc;
	  }
	}

	delete [] nCount;
	delete [] displace;
	delete [] send_buffer;
	delete [] recv_buffer;

#endif /* USE_MPI */
       
      } /* ENDIF multi-processor */
      else {
	GlobalList = LocalActiveParticlesOfThisType;
      } // ENDIF serial
      
    }  /* ENDIF number of active particles > 0 */
    else { /* GlobalNumberOfActiveParticles = 0 */
      // need to clean up still

      delete [] nCount;
    }
  } /* ENFOR Active particle types */

  return SUCCESS;
}
  
