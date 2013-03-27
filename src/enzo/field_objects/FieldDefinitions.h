//
// Field Definitions
//
// Authors: Matthew Turk
//          Greg Bryan

#ifndef __FIELD_DEFINITIONS_H__
#define __FIELD_DEFINITIONS_H__

// These enums may not be 64-bit safe

#define MAX_DIMENSIONS 3
#define MAX_FIELDS 40

#define TEST_huge_number 1e300
#define TEST_tiny_number 1e-300

enum CenteringType {
  CellCentered,
  CornerCentered,
  FaceCenteredX,
  FaceCenteredY,
  FaceCenteredZ,
  EdgeCenteredX,
  EdgeCenteredY,
  EdgeCenteredZ
};

enum InterpolationType {
  MultiplyByDensity,   // This quantity should be multiplied by Density before interpolating
  InterpolateDirectly, // This quantity should not be multiplied by anything before interpolating
  Undefined            // Interpolation behavior is undefined
};

#endif
