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
  NEWFIELD("colour",            "unk",      CellCentered, InterpolateDirectly, SNColour);
  NEWFIELD("Metal_Density",     "g / cm^3", CellCentered, InterpolateDirectly, Metallicity);
  NEWFIELD("MetalSNIaDensity",  "unk",      CellCentered, InterpolateDirectly, MetalSNIaDensity);

  NEWFIELD("CI_Density",        "g / cm^3", CellCentered, InterpolateDirectly, CIDensity);
  NEWFIELD("CII_Density",       "g / cm^3", CellCentered, InterpolateDirectly, CIIDensity);
  NEWFIELD("OI_Density",        "g / cm^3", CellCentered, InterpolateDirectly, OIDensity);
  NEWFIELD("OII_Density",       "g / cm^3", CellCentered, InterpolateDirectly, OIIDensity);
  NEWFIELD("SiI_Density",       "g / cm^3", CellCentered, InterpolateDirectly, SiIDensity);
  NEWFIELD("SiII_Density",      "g / cm^3", CellCentered, InterpolateDirectly, SiIIDensity);
  NEWFIELD("SiIII_Density",     "g / cm^3", CellCentered, InterpolateDirectly, SiIIIDensity);
  NEWFIELD("CHI_Density",       "g / cm^3", CellCentered, InterpolateDirectly, CHIDensity);
  NEWFIELD("CH2I_Density",      "g / cm^3", CellCentered, InterpolateDirectly, CH2IDensity);
  NEWFIELD("CH3II_Density",     "g / cm^3", CellCentered, InterpolateDirectly, CH3IIDensity);
  NEWFIELD("C2I_Density",       "g / cm^3", CellCentered, InterpolateDirectly, C2IDensity);
  NEWFIELD("COI_Density",       "g / cm^3", CellCentered, InterpolateDirectly, COIDensity);
  NEWFIELD("HCOII_Density",     "g / cm^3", CellCentered, InterpolateDirectly, HCOIIDensity);
  NEWFIELD("OHI_Density",       "g / cm^3", CellCentered, InterpolateDirectly, OHIDensity);
  NEWFIELD("H2OI_Density",      "g / cm^3", CellCentered, InterpolateDirectly, H2OIDensity);
  NEWFIELD("O2I_Density",       "g / cm^3", CellCentered, InterpolateDirectly, O2IDensity);

  NEWFIELD("ExtraType0",        "unk", CellCentered, InterpolateDirectly, ExtraType0);
  NEWFIELD("ExtraType1",        "unk", CellCentered, InterpolateDirectly, ExtraType1);
  NEWFIELD("kphHI",             "unk", CellCentered, MultiplyByDensity, kphHI);
  NEWFIELD("PhotoGamma",        "unk", CellCentered, MultiplyByDensity, PhotoGamma);
  NEWFIELD("kphHeI",            "unk", CellCentered, MultiplyByDensity, kphHeI);
  NEWFIELD("gammaHeI",          "unk", CellCentered, MultiplyByDensity, gammaHeI);
  NEWFIELD("kphHeII",           "unk", CellCentered, MultiplyByDensity, kphHeII);
  NEWFIELD("gammaHeII",         "unk", CellCentered, MultiplyByDensity, gammaHeII);
  NEWFIELD("kdissH2I",          "unk", CellCentered, MultiplyByDensity, kdissH2I);
  NEWFIELD("GravPotential",     "unk", CellCentered, InterpolateDirectly, GravPotential);
  NEWFIELD("Acceleration0",     "unk", CellCentered, InterpolateDirectly, Acceleration0);
  NEWFIELD("Acceleration1",     "unk", CellCentered, InterpolateDirectly, Acceleration1);
  NEWFIELD("Acceleration2",     "unk", CellCentered, InterpolateDirectly, Acceleration2);
  NEWFIELD("RadPressure0",      "unk", CellCentered, InterpolateDirectly, RadPressure0);
  NEWFIELD("RadPressure1",      "unk", CellCentered, InterpolateDirectly, RadPressure1);
  NEWFIELD("RadPressure2",      "unk", CellCentered, InterpolateDirectly, RadPressure2);
  NEWFIELD("Emissivity0",       "unk", CellCentered, InterpolateDirectly, Emissivity0);

  NEWFIELD("Bfield1",           "unk", CellCentered, InterpolateDirectly, Bfield1);
  NEWFIELD("Bfield2",           "unk", CellCentered, InterpolateDirectly, Bfield2);
  NEWFIELD("Bfield3",           "unk", CellCentered, InterpolateDirectly, Bfield3);
  NEWFIELD("PhiField",          "unk", CellCentered, InterpolateDirectly, PhiField);
  NEWFIELD("Phi_pField",        "unk", CellCentered, MultiplyByDensity, Phi_pField);
  NEWFIELD("DebugField",        "unk", CellCentered, InterpolateDirectly, DebugField); 

  NEWFIELD("DrivingField1", "unk", CellCentered, InterpolateDirectly, DrivingField1); 
  NEWFIELD("DrivingField2", "unk", CellCentered, InterpolateDirectly, DrivingField2); 
  NEWFIELD("DrivingField3", "unk", CellCentered, InterpolateDirectly, DrivingField3);

  NEWFIELD("AccelerationField1", "unk", CellCentered, MultiplyByDensity, AccelerationField1); 
  NEWFIELD("AccelerationField2", "unk", CellCentered, MultiplyByDensity, AccelerationField2); 
  NEWFIELD("AccelerationField3", "unk", CellCentered, MultiplyByDensity, AccelerationField3);

  NEWFIELD("Galaxy1Colour",       "unk", CellCentered, InterpolateDirectly, Galaxy1Colour);
  NEWFIELD("Galaxy2Colour",       "unk", CellCentered, InterpolateDirectly, Galaxy2Colour);
  NEWFIELD("Mach",                "unk", CellCentered, InterpolateDirectly, Mach);
  NEWFIELD("PreShockTemperature", "unk", CellCentered, InterpolateDirectly, PreShockTemperature);
  NEWFIELD("PreShockDensity",     "unk", CellCentered, InterpolateDirectly, PreShockDensity);  

  NEWFIELD("MBHColour",           "unk", CellCentered, InterpolateDirectly, MBHColour);
  NEWFIELD("ForbiddenRefinement", "unk", CellCentered, InterpolateDirectly, ForbiddenRefinement);

  NEWFIELD("RadiationFreq0", "unk", CellCentered, MultiplyByDensity, RadiationFreq0);
  NEWFIELD("RadiationFreq1", "unk", CellCentered, MultiplyByDensity, RadiationFreq1);
  NEWFIELD("RadiationFreq2", "unk", CellCentered, MultiplyByDensity, RadiationFreq2);
  NEWFIELD("RadiationFreq3", "unk", CellCentered, MultiplyByDensity, RadiationFreq3);
  NEWFIELD("RadiationFreq4", "unk", CellCentered, MultiplyByDensity, RadiationFreq4);
  NEWFIELD("RadiationFreq5", "unk", CellCentered, MultiplyByDensity, RadiationFreq5);
  NEWFIELD("RadiationFreq6", "unk", CellCentered, MultiplyByDensity, RadiationFreq6);
  NEWFIELD("RadiationFreq7", "unk", CellCentered, MultiplyByDensity, RadiationFreq7);
  NEWFIELD("RadiationFreq8", "unk", CellCentered, MultiplyByDensity, RadiationFreq8);
  NEWFIELD("RadiationFreq9", "unk", CellCentered, MultiplyByDensity, RadiationFreq9);

  NEWFIELD("RaySegments",    "unk", CellCentered, MultiplyByDensity, RaySegments);
  NEWFIELD("FieldUndefined", "unk", CellCentered, InterpolateDirectly, FieldUndefined);
   
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
