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
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
void grid::ReconstructFieldMapping(int ForceReconstruction)
{
 
  int field;
  FieldDescriptor *fd_base;
  std::string name;
  if (ForceReconstruction == FALSE && this->Fields.size() > 0) return;
  if (this->GridLevel < 0)  ENZO_FAIL("Illegal Grid Level");

  // determine the global "LeftEdge" of this field's active data
  long_int LeftEdge[3] = {0, 0, 0};
  FLOAT dx;
  for (int dim = 0; dim < this->GridRank; dim++) {
    dx = (this->GridRightEdge[dim] - this->GridLeftEdge[dim])
       / (this->GridEndIndex[dim]-this->GridStartIndex[dim]+1);
    LeftEdge[dim] = nlongint( (this->GridLeftEdge[dim] - DomainLeftEdge[dim]) / dx);
  }

  for (int field = 0; field < NumberOfBaryonFields; field++) {
    // We now do a double map lookup
    name = BaseFieldIDs[FieldType[field]];
    fd_base = BaseFieldTypes[name];
    this->Fields[name] = new FieldDescriptor(
	     fd_base, this->GridDimension, LeftEdge, this->BaryonField + field, 1, this->GridLevel);
#ifdef FIELD_DEBUG
    fprintf(stderr, "Locked %s to %"ISYM" (%"ISYM")\n", name.c_str(), field, FieldType[field]);
#endif
  }
 
}
