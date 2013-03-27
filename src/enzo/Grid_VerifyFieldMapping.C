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
 
void grid::VerifyFieldMapping()
{
#ifdef FIELD_DEBUG
  int field;
  int FoundError = FALSE;
  std::string name;
  FieldDescriptor *fd;
  for (int field = 0; field < NumberOfBaryonFields; field++) {
    // We now do a double map lookup
    name = BaseFieldIDs[FieldType[field]];
    fd = this->Fields[name];
    if ((fd == NULL) || (fd->GetValues() != this->BaryonField[field])) {
      fprintf(stderr, "PROBLEM: %s has become unlocked (%"ISYM" / %"ISYM")\n",
              name.c_str(), field, FieldType[field]);
      FoundError = TRUE;
    }
  }
  if (FoundError == FALSE) {
    fprintf(stderr, "ERROR-FREE\n");
  }
#endif
 
}
