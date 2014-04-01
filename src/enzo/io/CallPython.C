/***********************************************************************
/
/  INITIALIZE PYTHON INTERFACE AND START INTERPRETER
/
/  written by: Matthew Turk
/  date:       September, 2008
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#ifdef USE_PYTHON
#ifndef NEW_PROBLEM_TYPES
#error "Sorry, you need to have the new problem types enabled."
#endif
#endif

#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "TopGridData.h"
#include "ProblemType.h"
#ifdef USE_PYTHON
#include "message_passing.h"
#endif

int ExposeDataHierarchy(TopGridData *MetaData, HierarchyEntry *Grid, 
		       int &GridID, FLOAT WriteTime, int reset, int ParentID, int level);
void ExposeGridHierarchy(int NumberOfGrids);
void ExportParameterFile(TopGridData *MetaData, FLOAT CurrentTime, FLOAT OldTime, float dtFixed);
void CommunicationBarrier();

int CallPython(LevelHierarchyEntry *LevelArray[], 
               HierarchyEntry *Grids[],
               TopGridData *MetaData,
               int level, int from_topgrid)
{
#ifndef USE_PYTHON
    return SUCCESS;
#else
    if(access("user_script.py", F_OK) == -1) return SUCCESS;

    NumberOfPythonCalls++;
    if (from_topgrid) {
      NumberOfPythonTopGridCalls++;
      if (!(PythonTopGridSkip) ||
	  (NumberOfPythonTopGridCalls % PythonTopGridSkip) != 0) return SUCCESS;
    }
    else {
      if (LevelArray[level+1] != NULL) return SUCCESS;
      NumberOfPythonSubcycleCalls++;
      if (!(PythonSubcycleSkip) ||
	  (NumberOfPythonSubcycleCalls % PythonSubcycleSkip) != 0) return SUCCESS;
    }

	/* Initialize our message passer */
    static int message_passing_initialized = FALSE;
    if (message_passing_initialized == FALSE) {
        initmessage_passing();
        message_passing_initialized = TRUE;
    }
	create_coordinator(Grids, *MetaData);

    FLOAT CurrentTime, OldTime;
    float dtFixed;
    int num_grids, start_index;
    num_grids = 0; start_index = 1;

    PyDict_Clear(grid_dictionary);
    PyDict_Clear(old_grid_dictionary);

    LevelHierarchyEntry *Temp2 = LevelArray[0];
    /* Count the grids */
    /* I think there is a better idiom for this somewhere
       but I couldn't find it, and I think this works     */
    for (int lc = 0; LevelArray[lc] != NULL; lc++){
        Temp2 = LevelArray[lc];
        while (Temp2 != NULL) {
            num_grids++; Temp2 = Temp2->NextGridThisLevel;
        }
    }
    ExposeGridHierarchy(num_grids);
    Temp2 = LevelArray[0];
    while (Temp2->NextGridThisLevel != NULL)
        Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */
    CurrentTime = LevelArray[level]->GridData->ReturnTime();
    OldTime = LevelArray[level]->GridData->ReturnOldTime();
    dtFixed = LevelArray[level]->GridData->ReturnTimeStep();
    if (ExposeDataHierarchy(MetaData, Temp2->GridHierarchyEntry, start_index,
                CurrentTime, 1, 0, 0) == FAIL) {
        fprintf(stderr, "Error in ExposeDataHierarchy\n");
        return FAIL;
    }

    ExportParameterFile(MetaData, CurrentTime, OldTime, dtFixed);

    CommunicationBarrier();
    if(PythonReloadScript == TRUE) PyRun_SimpleString("reload(user_script)\n");
    PyRun_SimpleString("user_script.main()\n");

    PyDict_Clear(grid_dictionary);
    PyDict_Clear(old_grid_dictionary);
    PyDict_Clear(hierarchy_information);
    PyDict_Clear(yt_parameter_file);
    PyDict_Clear(conversion_factors);
    PyRun_SimpleString("gc.collect()\n");
    return SUCCESS;
#endif
}

