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
  static long_int Zero[3] = {0, 0, 0}; // Will add more later
  if (ForceReconstruction == FALSE && this->Fields.size() > 0) return;
  for (int field = 0; field < NumberOfBaryonFields; field++) {
    // We now do a double map lookup
    name = BaseFieldIDs[FieldType[field]];
    fd_base = BaseFieldTypes[name];
    delete this->Fields[name];
    this->Fields[name] = new FieldDescriptor(
            fd_base, this->GridDimension, Zero, this->BaryonField + field, 1);
#ifdef FIELD_DEBUG
    fprintf(stderr, "Locked %s to %"ISYM" (%"ISYM")\n", name.c_str(), field, FieldType[field]);
#endif
  }
 
}
