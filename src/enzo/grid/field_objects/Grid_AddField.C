/***********************************************************************
/
/  GRID CLASS (RECONSTRUCT MAPPING FOR FIELDS)
/
/  written by: Matthew Turk, Greg Bryan
/  date:       March, 2013
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
//
//  Assign basic values to a grid (allocate fields)
//

#include "preincludes.h"
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
FieldDescriptor* Grid::AddField(char *Name, char *UnitsName, int NumberOfFieldGhostZones,
                          InterpolationType InterpolationMethod,
                          CenteringType ValueCentering) {
  // We own this field.
  int Dimensions[MAX_DIMENSIONS];
  long long LeftEdge[MAX_DIMENSIONS];
  int dim;
  for (dim = 0; dim < MAX_DIMENSIONS; dim++) {
    Dimensions[dim] = this->Dimensions[dim] + 2*NumberOfFieldGhostZones;
    LeftEdge[dim] = this->LeftEdge[dim] - NumberOfFieldGhostZones;
  }
  if (this->GridLevel < 0)  ENZO_FAIL("Illegal Grid Level");
  FieldDescriptor *NewField = new FieldDescriptor(
        ValueCentering, 3, Dimensions, LeftEdge,
        InterpolationMethod, Name, UnitsName,
        this->BaryonFields + this->NumberOfFields++, this->GridLevel);
  this->Fields[Name] = NewField;
  return NewField;
}

