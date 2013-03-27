#include "FieldObjects.h"
#include <iostream>

int main(int argc, char **argv) {
  int Dimensions[MAX_DIMENSIONS] = {16, 16, 16};
  long long LeftEdge[MAX_DIMENSIONS] = {0, 0, 0};
  FieldDescriptor *fd = new FieldDescriptor(
    CellCentered, 3, Dimensions, LeftEdge,
    InterpolateDirectly, "Density", "g/cc", NULL);
  double *v = fd->GetValues();
  delete fd;
}
