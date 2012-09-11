/***********************************************************************
/
/  DETERMINE IF TWO BOXES OVERLAP
/
/  written by: Stephen Skory
/  date:       September, 2012
/  modified1:
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/

#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"

int one_overlap(FLOAT a[MAX_DIMENSION], FLOAT b[MAX_DIMENSION],
                FLOAT s[MAX_DIMENSION], FLOAT t[MAX_DIMENSION],
                int dims) {
  int dim;
  int overlap = SUCCESS;
  for (dim = 0; dim < dims; dim++) {
    if ((a[dim] >= t[dim]) || (b[dim] <= s[dim])) {
      overlap = FAIL;
      break;
    }
  }
  if (overlap) {return SUCCESS;}
  else {return FAIL;}
}

int check_overlap(FLOAT a[MAX_DIMENSION], FLOAT b[MAX_DIMENSION],
                  FLOAT s[MAX_DIMENSION], FLOAT t[MAX_DIMENSION],
                  int dims, FLOAT period[MAX_DIMENSION],
                  int *shift) {
  // a (left), b (right) - corners of box 1
  // s, t - corners of box 2
  FLOAT a_temp[MAX_DIMENSION], b_temp[MAX_DIMENSION];
  FLOAT max_shift[MAX_DIMENSION];
  int max_sum = -1, this_sum;
  int shift0, shift1, shift2, dim, overlap, any_overlap, max1, max2;
  any_overlap = 0;
  // Test the simplest case, where they already overlap.
  overlap = one_overlap(a, b, s, t, dims);
  if (overlap) return SUCCESS;
  // If we're here, we need to test all cases.
  // Keep checking until max_sum = dim or we've exhausted all cases.
  // This is to ensure that the box is shifted in all the dimensions it needs
  // to for reproducibility, meaning if we shift box 1 by an amount, we would
  // (swapping box 1 and box 2 in the algorithm) shift box 2 by the exact
  // opposite amount.
  // Shift box 1 around, keeping grid 2 static.
  overlap = SUCCESS;
  // Decide how many dimensions we move in.
  max1 = (dims > 1) ? 2 : 0;
  max2 = (dims > 2) ? 2 : 0;
  for (shift0 = -1; shift0 < 2; shift0++) {
    shift[0] = shift0;
    for (shift1 = -1; shift1 < max1; shift1++) {
      if (max1 > 0) shift[1] = shift1;
      for (shift2 = -1; shift2 < max2; shift2++) {
        if (max2 > 0) shift[2] = shift2;
        // We can skip [0,0,0]...
        if ((shift0 == 0) && (shift1 == 0) && (shift2 == 0)) continue;
        this_sum = 0;
        for (dim = 0; dim < dims; dim++) {
          a_temp[dim] = a[dim] + shift[dim] * period[dim];
          b_temp[dim] = b[dim] + shift[dim] * period[dim];
          this_sum += abs(shift[dim]);
        }
        overlap = one_overlap(a_temp, b_temp, s, t, dims);
        if (overlap) {
          if (this_sum == dim) return SUCCESS;
          any_overlap = 1;
          max_sum = max(max_sum, this_sum);
          if (max_sum == this_sum) {
            for (dim = 0; dim < dims; dim++) {
              max_shift[dim] = shift[dim];
            }
          }
        } // if overlap
      } // shift2
    } // shift1
  } // shift0
  // If we've had any overlap at all, copy max_shift back to shift and return
  // positively.
  if (any_overlap) {
    for (dim = 0; dim < dims; dim++) {
      shift[dim] = max_shift[dim];
    }
    return SUCCESS;
  }
  // If we get here, they don't overlap.
  return FAIL;
}

