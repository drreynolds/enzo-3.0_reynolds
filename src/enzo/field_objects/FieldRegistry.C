//
// Field Descriptor Object
//
// Authors: Matthew Turk
//          Greg Bryan

#include "FieldObjects.h"
#include "malloc.h"
#include "assert.h"
#include <map>
#include <string>

#define NEWFIELD(N, U, C, I) (*Fields)[N] = \
     new FieldDescriptor(C, Rank, ZeroDims, ZeroLeftEdge, I, N, U, NULL);

void FillFieldRegistry(int Rank, FieldRegistry* Fields) {
  // We provide a routine that will return a map of base field descriptors.
  // Part of the reasoning behind this is that we want to be able to
  // initialize from basic fields with their own centering, names, units, and
  // so on.  However, we don't really need this to be run in a .h file.

  int ZeroDims[MAX_DIMENSIONS] = {1, 1, 1};
  long long ZeroLeftEdge[MAX_DIMENSIONS] = {0,0,0};

  // First the main fields we know about
  
  NEWFIELD("Density",           "g / cm^3",     CellCentered, InterpolateDirectly);
  NEWFIELD("TotalEnergy",       "erg / g",      CellCentered, MultiplyByDensity);
  NEWFIELD("GasEnegy",          "erg / g",      CellCentered, MultiplyByDensity);
  NEWFIELD("x-velocity",        "cm / s",       CellCentered, MultiplyByDensity);
  NEWFIELD("y-velocity",        "cm / s",       CellCentered, MultiplyByDensity);
  NEWFIELD("z-velocity",        "cm / s",       CellCentered, MultiplyByDensity);

  // Species Fields
 
  NEWFIELD("Electron_Density",  "g / cm^3",     CellCentered, InterpolateDirectly);
  NEWFIELD("HI_Density",        "g / cm^3",     CellCentered, InterpolateDirectly);
  NEWFIELD("HII_Density",       "g / cm^3",     CellCentered, InterpolateDirectly);
  NEWFIELD("HeI_Density",       "g / cm^3",     CellCentered, InterpolateDirectly);
  NEWFIELD("HeII_Density",      "g / cm^3",     CellCentered, InterpolateDirectly);
  NEWFIELD("HeIII_Density",     "g / cm^3",     CellCentered, InterpolateDirectly);
  NEWFIELD("HM_Density",        "g / cm^3",     CellCentered, InterpolateDirectly);
  NEWFIELD("H2I_Density",       "g / cm^3",     CellCentered, InterpolateDirectly);
  NEWFIELD("H2II_Density",      "g / cm^3",     CellCentered, InterpolateDirectly);
  NEWFIELD("DI_Density",        "g / cm^3",     CellCentered, InterpolateDirectly);
  NEWFIELD("DII_Density",       "g / cm^3",     CellCentered, InterpolateDirectly);
  NEWFIELD("HDI_Density",       "g / cm^3",     CellCentered, InterpolateDirectly);
  NEWFIELD("Metal_Density",     "g / cm^3",     CellCentered, InterpolateDirectly);

}
