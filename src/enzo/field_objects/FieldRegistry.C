//
// Field Descriptor Object
//
// Authors: Matthew Turk
//          Greg Bryan

#include "FieldObjects.h"

#define NEWFIELD(N, U, C, I, E) (*Fields)[N] = \
     new FieldDescriptor(C, Rank, ZeroDims, ZeroLeftEdge, I, N, U, NULL); \
     (*FieldIDs)[E] = N;

void FillFieldRegistry(int Rank, FieldRegistry* Fields, FieldNumbers *FieldIDs) {
  // We provide a routine that will return a map of base field descriptors.
  // Part of the reasoning behind this is that we want to be able to
  // initialize from basic fields with their own centering, names, units, and
  // so on.  However, we don't really need this to be run in a .h file.

  int ZeroDims[MAX_DIMENSIONS] = {1, 1, 1};
  long_int ZeroLeftEdge[MAX_DIMENSIONS] = {0,0,0};

  // First the main fields we know about
  
  NEWFIELD("Density",           "g / cm^3", CellCentered, InterpolateDirectly, Density);
  NEWFIELD("TotalEnergy",       "erg / g",  CellCentered, MultiplyByDensity, TotalEnergy);
  NEWFIELD("GasEnegy",          "erg / g",  CellCentered, MultiplyByDensity, InternalEnergy);
  NEWFIELD("x-velocity",        "cm / s",   CellCentered, MultiplyByDensity, Velocity1);
  NEWFIELD("y-velocity",        "cm / s",   CellCentered, MultiplyByDensity, Velocity2);
  NEWFIELD("z-velocity",        "cm / s",   CellCentered, MultiplyByDensity, Velocity3);

  // Species Fields
 
  NEWFIELD("Electron_Density",  "g / cm^3", CellCentered, InterpolateDirectly, ElectronDensity);
  NEWFIELD("HI_Density",        "g / cm^3", CellCentered, InterpolateDirectly, HIDensity);
  NEWFIELD("HII_Density",       "g / cm^3", CellCentered, InterpolateDirectly, HIIDensity);
  NEWFIELD("HeI_Density",       "g / cm^3", CellCentered, InterpolateDirectly, HeIDensity);
  NEWFIELD("HeII_Density",      "g / cm^3", CellCentered, InterpolateDirectly, HeIIDensity);
  NEWFIELD("HeIII_Density",     "g / cm^3", CellCentered, InterpolateDirectly, HeIIIDensity);
  NEWFIELD("HM_Density",        "g / cm^3", CellCentered, InterpolateDirectly, HMDensity);
  NEWFIELD("H2I_Density",       "g / cm^3", CellCentered, InterpolateDirectly, H2IDensity);
  NEWFIELD("H2II_Density",      "g / cm^3", CellCentered, InterpolateDirectly, H2IIDensity);
  NEWFIELD("DI_Density",        "g / cm^3", CellCentered, InterpolateDirectly, DIDensity);
  NEWFIELD("DII_Density",       "g / cm^3", CellCentered, InterpolateDirectly, DIIDensity);
  NEWFIELD("HDI_Density",       "g / cm^3", CellCentered, InterpolateDirectly, HDIDensity);
  NEWFIELD("Metal_Density",     "g / cm^3", CellCentered, InterpolateDirectly, Metallicity);

#ifdef FIELD_DEBUG
  FieldRegistry::iterator iter;
  for (iter = Fields->begin(); iter != Fields->end(); ++iter) {
    fprintf(stderr, "Registered %s\n", iter->second->GetName());
  }
  FieldNumbers::iterator iter2;
  for (iter2 = FieldIDs->begin(); iter2 != FieldIDs->end(); ++iter2) {
    fprintf(stderr, "Locked %d to %s\n", iter2->first, iter2->second.c_str());
  }
#endif
}
